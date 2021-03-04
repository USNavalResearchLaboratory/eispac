
__all__ = ['ccd_offset']

"""

This routine is a (nearly) verbatim translation of the IDL version of eis_ccd_offset written by
Peter Young. The IDL documentation is included here.

;+
; NAME
;
;     EIS_CCD_OFFSET()
;
; PROJECT
;
;     Hinode/EIS
;
; EXPLANATION
;
;     Evaluates the spatial offset of the specified wavelength relative
;     to He II 256. If you see a feature in the He II 256 image at Y
;     coordinate 100 arcsec (for example), then the correct Y coordinate
;     for any other wavelength is:
;
;     Y = 100 - eis_ccd_offset(wvl)
;
;     Note that the spatial coordinate system for EIS is evalulated for
;     He II 256 by cross-correlating with SOT.
;
;     The value of eis_ccd_offset(wvl) represents the number of pixels
;     that a feature seen in a wavelength wvl sits above the He II 256
;     image on the CCD.
;
;     The following goes into the calculation of eis_ccd_offset:
;
;     - the tilt of the grating is assumed to be linear with wavelength
;
;     - the tilt for SW was derived by Young et al. (2008, A&A, submitted)
;
;     - the tilt for LW was *assumed* to be the same as for SW
;
;     - the offset between SW and LW was measured by co-aligning images
;       from Fe VIII 185.21 and Si VII 275.35
; 
; INPUTS
;
;     WVL    Wavelength specified in angstroms. Can be an array.
;
; OUTPUTS
;
;     The spatial offset between the specified wavelength and He II 256.32.
;     The value represents how many pixels above He II 256.32 the specified
;     wavelength sits on the EIS detectors.
;
; EXAMPLES
;
;     IDL> offset=eis_ccd_offset(195.12)
;
;     IDL> l=findgen(120)+170
;     IDL> plot,l,eis_ccd_offset(l)
;
; HISTORY
;
;     Ver.1, 15-Aug-2008, Peter Young
;-

"""

import numpy as np

def ccd_offset(wvl):

    grating_tilt = -0.0792

    wvl = np.atleast_1d(wvl)
    offset = np.zeros(wvl.shape)

    for n in range(wvl.shape[0]):

        if wvl[n] > 230:
            offset[n] = grating_tilt*(wvl[n]-256.32)
        else:
            offset[n] = grating_tilt*(wvl[n]-185.21) + 18.5 + grating_tilt*(275.35-256.32)

    return offset

def test_ccd_offset():

    waves = [185.21, 195.119, 256.317, 258.375, 284.160]
    idl = [16.992825, 16.208033,  0.000239, -0.162755, -2.204928]
    py = ccd_offset(waves)

    print(f' {"wave":>7s} {"idl":>12s} {"python":>12s}')
    for n, w in enumerate(waves):
        print(f' {w:7.3f} {idl[n]:12.6f} {py[n]:12.6f}')

if __name__ == '__main__':

    test_ccd_offset()
