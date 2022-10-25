from sunpy.net.attr import SimpleAttr

__all__ = ['FileType']


class FileType(SimpleAttr):
    """
    Specifies the type of EIS level 1 file
    
    Parameters
    ----------
    value: `str`
        Possible values are "HDF5 data" or "HDF5 header" to retrieve the 
        data and header files, respectively, in HDF5 format, or "FITS" to
        retrieve the FITS files. Inputs are not case sensitive.
    """
    
    def __init__(self, value):
        if not isinstance(value, str):
            raise ValueError('File type must be a string')
        value = value.lower()
        if 'hdf5' in value:
            value = '.'.join([value[5:], 'h5'])
        if value == 'header.h5':
            value = 'head.h5'
        if value not in ['data.h5', 'head.h5', 'fits']:
            raise ValueError(f'File type {value} must be either "HDF5 data", "HDF5 header", or "FITS".')
        super().__init__(value)
