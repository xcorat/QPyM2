""" Contains utulity classes and routines. """

class TempConfig:
    """ A temporary configuration class that will have all the keyword arguments 
        implemented as data members.

        Example:
            >>> cfg = TempConfig(name='hmu', threshold=threshold)
            >>> cfg.name
            'hmu'  

        TODO: find an existing implementation for this. 
    """
    def __init__(self, **kwds):
        self.__dict__.update(kwds)
        self._kwd_names = kwds.keys()

    def __repr__(self):
        return f'TempConfig({self._kwd_names})'

def get_ncores(print_state=False):
    """ Check the number of cores available on the node. Copied from NERSC doc.

    Returns:
        int: number of cores available on the node.
    """
    import multiprocessing as mp
    import os
    import psutil

    nthreads = psutil.cpu_count(logical=True)
    ncores = psutil.cpu_count(logical=False)
    nthreads_per_core = nthreads // ncores
    nthreads_available = len(os.sched_getaffinity(0))
    ncores_available = nthreads_available // nthreads_per_core

    assert nthreads == os.cpu_count()
    assert nthreads == mp.cpu_count()

    if print_state:
        print(f'{nthreads=}')
        print(f'{ncores=}')
        print(f'{nthreads_per_core=}')
        print(f'{nthreads_available=}')
        print(f'{ncores_available=}')

    return ncores_available  