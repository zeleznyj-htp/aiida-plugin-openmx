def replace_dict(input_dict, to_be_replaced_string, new_string):
    return {key.replace(to_be_replaced_string, new_string): value for key, value in input_dict.items()}

def replace_dots(a):
    return replace_dict(a,'.','\\\\')

def replace_backslash(a):
    return replace_dict(a, '\\\\', '.')