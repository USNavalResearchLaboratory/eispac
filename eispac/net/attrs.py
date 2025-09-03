from sunpy.net.attr import SimpleAttr

__all__ = ['FileType']


class FileType(SimpleAttr):
    """
    Specifies the type of EIS level 1 file
    
    Parameters
    ----------
    value: `str`
        Select from "HDF5" or "HDF5 data" to retrieve both the data and header 
        files in HDF5 format, "HDF5 header" for ONLY the headers, "FITS" for
        for both the "l1" and "er" FITS files, or "ANY" to retrieve all
        file types available. Inputs are not case sensitive.
    """
    
    def __init__(self, value):
        if not isinstance(value, str):
            raise ValueError('File type must be a string')
        value = value.lower().strip()
        if 'hdf5' in value:
            if len(value) > 4:
                value = '.'.join([value[5:], 'h5'])
            else:
                # Selecting just 'HDF5' will search/download 'data.h5' files
                value = 'data.h5'
        if value == 'header.h5':
            value = 'head.h5'
        if value == 'all':
            value = 'any'
        if value not in ['data.h5', 'head.h5', 'fits', 'any']:
            raise ValueError(f'File type {value} must be either "HDF5",'
                            +f' "HDF5 data", "HDF5 header", "FITS", or "ANY".')
        super().__init__(value)
