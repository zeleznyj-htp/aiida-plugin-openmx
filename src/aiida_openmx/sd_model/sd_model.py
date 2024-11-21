"""
This is a program for constructing an sd tight-binding models.

The model is stored in class Model, which also provides function
for writing the Hamiltonian to a tb_hr.dat file in a format that
can be read by the fortran linear response code.

The model is defined in input file model.in, which has an .ini
syntax.

This program when it is run reads the model.in file and writes the
tb_hr.dat file, but it can also be used as a library.
"""

import sys
import copy
from datetime import datetime

import configparser
import argparse
import itertools
import sympy as sp
from sympy import sympify as spf
from sympy import I, N
from sympy.physics.matrices import msigma
from sympy.physics.quantum import TensorProduct
from sympy.functions import re, im, exp, Abs, cos, sin

import numpy as np
from numpy.linalg import norm
from scipy.constants import hbar
import matplotlib.pyplot as plt

def sym2num(mat):
    return np.array(mat).astype(complex)

def symmetrize_nns(nns):
    nns_out = copy.deepcopy(nns)
    for i,nn in enumerate(nns):
        for n in nn:
            j = n[0]
            cell = n[1]
            dist = n[2]
            cell_op = tuple([-c for c in cell])
            is_in = False
            for n2 in nns[j]:
                if n2[0] == i and n2[1] == cell_op:
                    if abs( n2[2] - dist ) > 1e-5:
                        raise Exception('Something is wrong with symmetrize_nns')
                    is_in = True
            if not is_in:
                new_n = [i,cell_op,dist,-1]
                nns_out[j].append(new_n)
                
    return nns_out

def input_from_pymatgen(struct,output=None,nns=1,t=-1,J=1.7,nn_scaling=2):
    
    config = configparser.ConfigParser()
    
    config['structure'] = {}
    config['structure']['dim'] = '3'
    config['structure']['lattice_scale'] = '1'
    config['structure']['lattice'] = '\n' + np.array2string(struct.lattice.matrix).replace('[','').replace(']','')
    config['structure']['n_atoms'] = str(len(struct.sites))
    
    atoms = '\n'
    for site in struct.sites:
        atoms += np.array2string(site.frac_coords).replace('[','').replace(']','') + '\n'
    config['structure']['atoms'] = atoms
    
    config['tb'] = {}
    config['tb']['hopping'] = str(t)
    config['tb']['nns'] = str(nns)
    config['tb']['nn_scaling'] = str(nn_scaling)
    config['tb']['exchange'] = str(J)

    moms = '\n'
    for mom in struct.site_properties['magmom']:
        moms += np.array2string(mom).replace('[','').replace(']','') + '\n'
    config['tb']['mag_moms'] = moms
    
    if output is not None:
        with open(output,'w') as f:
            config.write(f)

    return config

class Model:
    """
    This class stores the model.

    Parameters of the model are read from an input file.
    """

    def __init__(self, inp_file=None, config=None):
        """
        The parameters of the model are read from inp_file.

        The input file has an .ini like syntax.
        """
        if inp_file is None and config is None:
            raise Exception('inp_file or config must be provided')
        if inp_file is not None:
            config = configparser.ConfigParser()
            config.read(inp_file)

        dim = int(config['structure']['dim'])

        self.dim = dim

        latt = config['structure']['lattice'].split('\n')
        B = np.zeros((dim, dim))
        for i in range(dim):
            for j in range(dim):
                B[i, j] = spf(latt[i + 1].split()[j])

        B = B * spf(config['structure']['lattice_scale'])

        B = np.array(B).astype(np.float64)

        self.B = B

        if config.has_option('structure', 'periodicity'):
            self.periodicity = [int(x) for x in config['structure']['periodicity'].split()]
            self.fully_periodic = False
        else:
            self.periodicity = [1, 1, 1]
            self.fully_periodic = True

        n_atoms = int(config['structure']['n_atoms'])
        atoms = sp.zeros(n_atoms, dim)
        at = config['structure']['atoms'].split('\n')
        for i in range(n_atoms):
            for j in range(dim):
                atoms[i, j] = spf(at[i + 1].split()[j])

        atoms = np.array(atoms).astype(np.float64)

        self.n_atoms = atoms.shape[0]
        self.n_orbs = self.n_atoms * 2
        self.atoms = atoms

        # implement different onsite energies
        if config.has_option('tb', 'onsite_energies'):
            onsite = config['tb']['onsite_energies'].split('\n')[1:]
            onsite = [np.float(on) for on in onsite]
        else:    # ignore if onsite energies are not defined in the input file
            onsite = np.nan
        
        self.onsite = onsite

        t = np.float64(spf(config['tb']['hopping']))
        J = np.float64(spf(config['tb']['exchange']))
        self.t = t
        self.J = J
        mags = sp.zeros(n_atoms, 3)
        ms = config['tb']['mag_moms'].split('\n')
        for i in range(n_atoms):
            for j in range(3):
                mags[i, j] = spf(ms[i + 1].split()[j])

        # we normalize the magnetic moments
        for at in range(n_atoms):
            if mags[at, :].norm() > 1e-6:
                mags[at, :] = mags[at, :] / mags[at, :].norm()

        mags = np.array(mags).astype(np.float64)

        self.mags = mags

        if config.has_option('tb', 'nns'):
            self.nns = spf(config['tb']['nns'])
            if self.nns == 0 or self.nns < -1:
                raise Exception('nns must be > 1 or -1')
        else:
            self.nns = 1

        if self.nns > 1 or self.nns == -1:
            self.nn_scaling = config['tb']['nn_scaling']

            if self.nn_scaling == '3':
                self.nn_scaling_factor = float(config['tb']['nn_scaling_factor'])

        if config.has_option('tb', 'nn_prec'):
            self.nn_prec = float(config['tb']['nn_prec'])
        else:
            self.nn_prec = 1e-5

        if config.has_option('tb', 'nn_n_cells'):
            self.nn_n_cells = float(config['tb']['nn_n_cells'])
        else:
            self.nn_n_cells = 2

        self.debug = True

        print('{} starting to create Hr'.format(datetime.now().strftime("%H:%M:%S")))
        self.Hr = self.create_Hr()

    def t_scale(self, dist, dist_nn):
        """
        If we consider also neighbors beyond the nearest neighbor than we scale the hopping for these.

        """

        if self.nns == 1:
            return self.t
        else:
            if self.nn_scaling == '2':
                return self.t * dist_nn ** 2 / dist ** 2
            elif self.nn_scaling == '1':
                return self.t
            elif self.nn_scaling == '3':
                return self.t * (dist_nn ** (self.nn_scaling_factor) ) / ( dist ** (self.nn_scaling_factor) )

    def find_nn(self):
        """
        Finds nearest neighbours.

        Only finds the nearest neighbourse, but could be easily modified for
        next to nearest etc.

        Returns:
            nn: A list, which for each atom contains a list of all nearest neighbours
                up to nth-nearest neighbors, where n=self.nns.
                Nearest neighbours have format:
                [atom_index,unit cell coordinates of the neighbor,distance,n],
                where n means order of the neighbor, i.e., it's nth nearest neighbor.
        """

        dim = self.dim

        ncells = 2  # ncells determines how many cells we use to look for nearest neighbors

        B = self.B
        atoms = self.atoms

        def find_position(n, atom_n):
            """
            Finds position of atom with index atom_n in unit cell with coordinates given by n.
            n is given in fractional coordinates.
            """

            position = np.matmul(B.transpose(), atoms[atom_n, :])
            # now move to correct unit cell:
            for i in range(dim):
                position += B[i, :] * n[i]
            return position

        # this builds a list, which runs over nearest cells
        iter1 = range(-ncells, ncells + 1, 1)
        iterlist = []
        for i in range(dim):
            iterlist.append(iter1)
        cells_iter = list(itertools.product(*iterlist))

        nn = []
        # loop over all atoms
        for a in range(self.n_atoms):
            neighbs = []  # neighbs contains all neighbors in the cells we consider for each atom
            for n in cells_iter:  # loop over nearest cells
                for a2 in range(self.n_atoms):  # loop over all atoms in given cell
                    dist = find_position(n, a2) - find_position([0] * dim, a)
                    neighbs.append([a2, n, norm(dist)])

            # we sort the neighbors by distance
            neighbs = sorted(neighbs, key=lambda tup: tup[2])
            neighbs.pop(0)  # remove the first one since it is the atom itself
            # count the order of each neighbor
            for i, ne in enumerate(neighbs):
                if i == 0:
                    n_nn = 1
                else:
                    if norm(ne[2] - neighbs[i - 1][2]) > self.nn_prec:
                        n_nn += 1
                ne.append(n_nn)

            # now we select only the nearest neighbors up to order self.nns
            nn_a = []
            for ne in neighbs:
                # if self.debug:
                # if ne[3] < 10:
                #    print('printing ne')
                #    print ne
                if ne[3] <= self.nns:
                    nn_a.append(ne)
            nn.append(nn_a)

        return nn

    def find_nn_new(self):
        """
        Finds nearest neighbours.

        Only finds the nearest neighbourse, but could be easily modified for
        next to nearest etc.

        Returns:
            nn: A list, which for each atom contains a list of all nearest neighbours
                up to nth-nearest neighbors, where n=self.nns.
                Nearest neighbours have format:
                [atom_index,unit cell coordinates of the neighbor,distance,n],
                where n means order of the neighbor, i.e., it's nth nearest neighbor.
        """

        dim = self.dim

        ncells = self.nn_n_cells  # ncells determines how many cells we use to look for nearest neighbors

        B = self.B
        atoms = self.atoms

        def find_position(n, atom_n):
            """
            Finds position of atom with index atom_n in unit cell with coordinates given by n.
            n is given in fractional coordinates.
            """

            position = np.matmul(B.transpose(), atoms[atom_n, :])
            # now move to correct unit cell:
            for i in range(dim):
                position += B[i, :] * n[i]
            return position

        # this builds a list, which runs over nearest cells
        iter1 = range(-ncells, ncells + 1, 1)
        iterlist = []
        for i in range(dim):
            iterlist.append(iter1)
        cells_iter = list(itertools.product(*iterlist))

        # loop over all atoms
        neighbs_all = []
        for a in range(self.n_atoms):
            neighbs = []  # neighbs contains all neighbors in the cells we consider for each atom
            for n in cells_iter:  # loop over nearest cells
                for a2 in range(self.n_atoms):  # loop over all atoms in given cell
                    dist = find_position(n, a2) - find_position([0] * dim, a)
                    neighbs.append([a2, n, norm(dist)])

            # we sort the neighbors by distance
            neighbs = sorted(neighbs, key=lambda tup: tup[2])
            neighbs.pop(0)  # remove the first one since it is the atom itself
            neighbs_all.append(neighbs)

        dists_all = []
        for neighbs in neighbs_all:
            for n in neighbs:
                dists_all.append(n[-1])
        dists_all.sort()

        unique_dists = []
        for i,d in enumerate(dists_all):
            if i == 0:
                unique_dists.append(d)
            else:
                if abs(d-dists_all[i-1]) > self.nn_prec:
                    unique_dists.append(d)
                if self.nns > 0:
                    if len(unique_dists) == self.nns:
                        max_dist = unique_dists[-1] + self.nn_prec
                        break

        def get_order(dist):
            difs = []
            for i,ud in enumerate(unique_dists):
                difs.append(abs(ud-dist))
            return difs.index(min(difs)) + 1

        nn = []

        if self.nns == -1:

            orders = []
            for neighbs in neighbs_all:
                n0_order = get_order(neighbs[0][-1])
                orders.append(n0_order)
            min_order = max(orders)
            print('min_order',min_order)
            max_dist = unique_dists[min_order-1] + self.nn_prec


        for neighbs in neighbs_all:
            nn_a = []
            for n in neighbs:
                d = n[-1]
                if d > max_dist:
                    break
                else:
                    n.append(get_order(d))
                    nn_a.append(n)

            nn.append(nn_a)

        #for ns in nn:
        #    print(ns)
        return nn

    def create_Hr(self):
        """
        This creates a tight-binding Hamiltonian in a similar format to that of wannier90.

        Returns: Hr: format is a dictionary, where each key has a format:
            (n,i,j), where are cell fractional coordinates and i,j are orbital numbers
            in this case, there are only two orbitals per atom, spin-up orbital on atom i
            has orbital number i+1, spin-down i+1+n_atoms (we start the counting from 1)
            Hr[n0,n1,i,j] is <0,0,i|H|n0,n1,j>
        """

        print('{} finding nearest neighbours'.format(datetime.now().strftime("%H:%M:%S")))
        nn = self.find_nn_new()

        #This shouldn't be necessary now, but just to be safe we keep it here.
        nn = symmetrize_nns(nn)

        print('{} nearest neighbors finished'.format(datetime.now().strftime("%H:%M:%S")))
        Hr = {}
        # first we add the hoppings
        for i in range(self.n_atoms):
            # finds the nearest neighbor distance
            dist_nn = nn[i][0][2]
            for neg in nn[i]:
                # construct the key
                cell = [0] * self.dim
                for j in range(self.dim):
                    cell[j] = neg[1][j]
                cell = tuple(cell)
                key = (cell, i + 1, neg[0] + 1)
                if self.fully_periodic:
                    skip_key = False
                else:
                    skip_key = False
                    for d in range(self.dim):
                        if self.periodicity[d] == 0:
                            if cell[d] != 0:
                                skip_key = True
                if skip_key:
                    continue
                # add the hopping
                # if we consider more than nearest neighbors than we scale the hopping according
                # to the distance
                if key in Hr:
                    Hr[key] += -self.t_scale(neg[2], dist_nn)
                else:
                    Hr[key] = -self.t_scale(neg[2], dist_nn)
                # the Hamiltonian has to be Hermitian
                key2 = (cell, i + 1 + self.n_atoms, neg[0] + 1 + self.n_atoms)
                if key2 in Hr:
                    Hr[key2] += -self.t_scale(neg[2], dist_nn)
                else:
                    Hr[key2] = -self.t_scale(neg[2], dist_nn)

            for site in np.arange(self.n_orbs):
                if self.onsite is np.nan:
                    break
                new_onsite = np.concatenate((self.onsite, self.onsite))
                cell = tuple([0] *self.dim)
                key3 = (cell, (site+1), (site+1))
                Hr[key3] = new_onsite[site]
                
        print('{} hoppings finished'.format(datetime.now().strftime("%H:%M:%S")))
        # this constructs the spin-operators projected on individual atoms

        msigma_np = []
        for i in range(3):
            msigma_np.append(sym2num(msigma(i + 1)))
        ms = np.zeros((self.n_atoms * 2, self.n_atoms * 2), dtype=complex)
        for at in range(self.n_atoms):
            proj = np.zeros((self.n_atoms, self.n_atoms))
            proj[at, at] = 1
            for i in range(3):
                ms += self.mags[at, i] * np.kron(msigma_np[i], proj)

        # add the exchange part
        for i in range(self.n_atoms * 2):
            for j in range(self.n_atoms * 2):
                cell = (0,) * self.dim
                key = (cell, i + 1, j + 1)
                if key in Hr:
                    Hr[key] += self.J * ms[i, j]
                else:
                    Hr[key] = self.J * ms[i, j]

        # print('Hr-----')
        # for key in Hr:
        #    if key[1] == key[2]:
        #        print(key,Hr[key])

        return Hr

    def write_wann_Hr(self):
        """
        This prints the Hr Hamiltonian in a file wannier90_hr.dat in the same format
        as wannier90 uses.

        !!!This is not the correct way how to do it and is left here only as a reference.!!!
        """

        Hr = self.create_Hr()
        rs = []
        for key in Hr:
            r = (key[0])
            if r not in rs:
                rs.append(r)
        n_rs = len(rs)

        with open('wannier90_hr.dat', 'w') as f:
            f.write('\n')
            f.write(str(self.n_atoms * 2) + '\n')
            f.write(str(n_rs) + '\n')
            nrlines = n_rs // 15 + 1
            for l in range(nrlines):
                for i in range(15):
                    if l * 15 + i < n_rs:
                        f.write('1  ')
                f.write('\n')
            for key in sorted(Hr):
                f.write(str(key[0][0]) + '  ')
                try:
                    f.write(str(key[0][1]) + '  ')
                except:
                    f.write(str(0) + '  ')
                try:
                    f.write(str(key[0][2]) + '  ')
                except:
                    f.write(str(0) + '  ')

                f.write(str(key[1]) + '  ')
                f.write(str(key[2]) + '  ')

                f.write(str(float(re(Hr[key]))) + '  ')
                f.write(str(float(im(Hr[key]))))

                f.write('\n')

        with open('POSCAR', 'w') as f:
            f.write('\n')
            f.write('1\n')
            Bout = sp.eye(3)
            Bout[0:self.dim, 0:self.dim] = self.B[:, :]
            for i in range(3):
                for j in range(3):
                    f.write(str(float(Bout[i, j])) + '  ')
                f.write('\n')

    def write_Hr(self,path='./'):
        """
        This prints the Hamiltonian in the tight-binding basis in a format
        that can be easily used for transforming to a k-dependent Hamiltonian.
        It writes a tb_hr.dat file.
        """

        Hr = self.Hr

        with open(path+'tb_hr.dat', 'w') as f:
            f.write('\n')
            f.write(str(self.n_atoms * 2) + '\n')
            f.write(str(len(Hr)) + '\n')
            for key in sorted(Hr):
                R = sp.zeros(1, 3)
                at_a = (key[2] - 1) % self.n_atoms
                at_b = (key[1] - 1) % self.n_atoms
                for i in range(self.dim):
                    R[i] += key[0][i] - self.atoms[at_b, i] + self.atoms[at_a, i]
                    f.write(str(float(R[i])) + '  ')
                for i in range(self.dim, 3):
                    f.write('0  ')

                f.write(str(key[1]) + '  ')
                f.write(str(key[2]) + '  ')

                f.write(str(float(re(Hr[key]))) + '  ')
                f.write(str(float(im(Hr[key]))))

                f.write('\n')

        with open(path+'input_latt', 'w') as f:
            f.write('\n')
            f.write('1\n')
            Bout = sp.eye(3)
            Bout[0:self.dim, 0:self.dim] = self.B[:, :]
            for i in range(3):
                for j in range(3):
                    f.write(str(float(Bout[i, j])) + '  ')
                f.write('\n')

    def create_Hk(self, k):

        #When the peridicioty is lower than the dimension we set the k along
        #the non-periodic dimensions to 0. Setting this to a different value would result
        #in a different Hamiltonian due to the presence of the (pa-pb) factor, however,
        #this Hamiltonian is related to the 0 one by a unitary transformation.
        kp = np.zeros(self.dim)
        for i in range(self.dim):
            if self.periodicity[i] == 1:
                kp[i] = k[i]

        Hr = self.Hr
        Hk = np.zeros((self.n_atoms * 2, self.n_atoms * 2), dtype=complex)

        for key in Hr:

            R = np.zeros(self.dim)
            for i in range(self.dim):
                R += self.B[i, :] * key[0][i]

            pa = np.zeros(self.dim)
            pb = np.zeros(self.dim)
            for i in range(self.dim):
                at_a = (key[2] - 1) % self.n_atoms
                at_b = (key[1] - 1) % self.n_atoms
                pa += self.B[i, :] * self.atoms[at_a, i]
                pb += self.B[i, :] * self.atoms[at_b, i]

            rho = R + pa - pb
            Hk[key[1] - 1, key[2] - 1] += exp(1j * kp.dot(rho)) * Hr[key]

        return Hk

    def plot_bands(self,kpts,nk,plot=True):
        Eks = np.zeros( ( (len(kpts)-1)*nk + 1,self.n_orbs))
        all_ks = []
        ki = 0
        for i in range(len(kpts)-1):
            for j in range(nk):
                k = kpts[i] + (kpts[i+1] - kpts[i]) * j/nk
                all_ks.append(k)
                Hk = self.create_Hk(k)
                w,v = np.linalg.eigh(Hk)
                Eks[ki,:] = w
                ki += 1
        k = kpts[-1]
        all_ks.append(k)
        Hk = self.create_Hk(k)
        w,v = np.linalg.eigh(Hk)
        Eks[ki,:] = w
        if plot:
            for i in range(Eks.shape[1]):
                plt.plot(Eks[:,i])
        return all_ks,Eks

    def get_spin_operators(self,atom=None):
        sigma_x = np.array([[0, 1],[1, 0]])
        sigma_y = np.array([[0, -1j],[1j, 0]])
        sigma_z = np.array([[1, 0],[0, -1]])
        nP = int(self.n_orbs/2)
        if atom is None:
            P = np.diag([1]*nP)
        else:
            P = np.zeros((nP,nP))
            P[atom,atom] = 1
        s_x = np.kron(sigma_x,P)
        s_y = np.kron(sigma_y,P)
        s_z = np.kron(sigma_z,P)    
        return s_x,s_y,s_z
        

    def get_spins(self,k,atom=None):
        Hk = self.create_Hk(k)
        s_x,s_y,s_z = self.get_spin_operators(atom)
        E,w = np.linalg.eigh(Hk)
        spins = []
        for i in range(len(E)):
            Sx = np.real(np.dot(np.conjugate(w[:,i]),np.dot(s_x,w[:,i])))
            Sy = np.real(np.dot(np.conjugate(w[:,i]),np.dot(s_y,w[:,i])))
            Sz = np.real(np.dot(np.conjugate(w[:,i]),np.dot(s_z,w[:,i])))
            spins.append(np.array([Sx,Sy,Sz]))
        return E,spins

    def spin_rot(self, atom, n, phi):

        n = n / n.norm()

        # creates the projected pauli matrices
        proj = sp.zeros(self.n_atoms)
        proj[atom, atom] = 1
        sigma_n = sp.zeros(self.n_atoms * 2, self.n_atoms * 2)
        for i in range(3):
            sigma_n += n[i] * TensorProduct(msigma(i + 1), proj)

        s_rot = cos(phi / 2) * TensorProduct(sp.eye(2), proj) + I * sigma_n * sin(phi / 2)

        return s_rot

    def lat2cart(self, vec):
        """converts a vector in lattice coordinates to the cartesian coordinates
        """
        vec_c = sp.zeros(1, self.dim)

        for i in range(self.dim):
            vec_c += self.B[i, :] * vec[i]

        return vec_c

    def get_reciprocal_vectors(self):
        R = np.zeros((3,3))
        V = np.dot(self.B[0,:],np.cross(self.B[1,:],self.B[2,:]))
        R[0,:] = 2 * np.pi * np.cross(self.B[1,:],self.B[2,:]) / V
        R[1,:] = 2 * np.pi * np.cross(self.B[2,:],self.B[0,:]) / V
        R[2,:] = 2 * np.pi * np.cross(self.B[1,:],self.B[2,:]) / V
        return R


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    args = parser.parse_args()

    print('{} starting'.format(datetime.now().strftime("%H:%M:%S")))
    model = Model('model.in')

    print('{} Model Initialized'.format(datetime.now().strftime("%H:%M:%S")))
    model.write_Hr()

    print('{} Hamiltonain written'.format(datetime.now().strftime("%H:%M:%S")))
