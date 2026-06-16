from itertools import product
from  copy import deepcopy

def generate_arg_combinations(args):
    r"""Generator for argument dictionaries based on `args`.
        Every value which is a list is taken to be a list of options for its key.
        The generator goes through all combinations of these options, and yields
        dictionaries which contain them, and the other ones from `args`.
        A value of `None` is a proxy for removing that key.
        Example:
        ---------
        for a in generate_arg_combinations({'a':[0,1],'b':[None,''],'c':4}):
            print a
        {'a': 0, 'c': 4}
        {'a': 0, 'c': 4, 'b': ''}
        {'a': 1, 'c': 4}
        {'a': 1, 'c': 4, 'b': ''}
            """
    variable = {k:v for k,v in args.items() if type(v)==list}
    varkeys = variable.keys()
    for vartup in product(*tuple(variable[k] for k in varkeys)):  
        args_i = deepcopy(args)
        args_i.update(dict(zip(varkeys, list(vartup))))
        remove_keys = [k for k in args_i.keys() if args_i[k] is None]
        for k in remove_keys:
            args_i.pop(k)
        yield args_i

def build_arg_string(args):
    r"""Turn a dictionary into a string with `--key value`.
    """
    arg_string = " ".join([f'--{k} {str(v)}' for k,v in args.items()])
    return arg_string