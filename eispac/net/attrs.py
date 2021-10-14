from sunpy.net.attr import SimpleAttr

__all__ = ['FileType']


class FileType(SimpleAttr):
    """
    Specifies the type of EIS level 1 file
    
    Parameters
    ----------
    value: `str`
        Either "data" or "header"
    """
    
    def __init__(self, value):
        if not isinstance(value, str):
            raise ValueError('File type must be a string')
        value = value.lower()
        if value not in ['data', 'header']:
            raise ValueError(f'File type {value} must be either "data" or "header".')
        # The actual files are labeled "head" but "header" is more readable
        value = 'head' if value == 'header' else value
        super().__init__(value)
