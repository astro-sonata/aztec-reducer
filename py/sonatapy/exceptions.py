'''
Some custom exceptions
'''

class MissingDataException(Exception):
    '''
    General error to be thrown if the code can't find some of the data that should
    be there.
    '''
    pass

class FitFailed(Exception):
    '''
    General error to be thrown when a fit failed in the code
    '''
    pass
