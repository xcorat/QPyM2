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
  