Pro Load_abs_coeffs

;+
; NAME:
;	LOAD_ABS_COEFFS, filename
; PURPOSE:
;	Loads mass absorption coefficients data from an XDR binary data file 
;	into the common block SXR_STUFF.
; CATEGORY:
;	X-ray calculations
; CALLING SEQUENCE:
;	LOAD_ABS_COEFFS
; INPUTS:
;	None
; OPTIONAL INPUT PARAMETERS:
;       filename
;           If this parameter is specified then it must be the name of a binary
;           XDR file containing the absorption coefficients.
;           If this parameter is not specified, then the filename is determined
;           from the environment variable IDL_ABS_COEFFS
; OUTPUTS:
;	None explicit, data put in common block
; OPTIONAL OUTPUT PARAMETERS:
;	None.
; COMMON BLOCKS:
;	SXR_STUFF.  Contains the following:
;	    SORLIST:	Structure.  List of X-ray sources and their parameters.
;	    ABCTAB:	Structure.  Mass absorption coefficient table.
;	    ENUN:	String.  Current energy units.
;	    DENUN:	String.  Default energy units.
;	    CURBREAM:	Structure.  Beam parameters.
; SIDE EFFECTS:
;	None.
; RESTRICTIONS:
;	The data file to be read from must be defined with the environment
; PROCEDURE:
;       Checks to see if it has been previously called. If so, it just returns
;       immediately.  This makes it convenient to call this routine in
;       every procedure which needs the absorption coefficient data.
;
;	Reads the data and creates an array of structures of type ABCO, one 
;	structure for each element.  Each ABCO structure containes, in order:
;	NAME  -	Element's name.
;	CSYM  -	Chemical symbol.
;	Z     -	Element's Z.
;	A     -	Element's A.
;	RO    -	Element's density in gr/cm^3.
;	COLEN -	Length of the following table.  Max. value is 63.
;	COTAB -	An (COLEN+1,3) array containing the spline coefficients that are
;		needed to calcualte the mass absorption coefficient for any
;		energy.  For information about the spline array, see routines
;		SPLIN_COEFFS and SPLIN_EVAL in MIDL.
;	EDLEN -	Length of the following table.  Max. value is 15.
;	EDTAB -	An (EDLEN+1,3) array containing, in order, the energies of the
;		edges (including first energy in table), the jump values at 
;		edges and -2 minus the values of the partial power coefficients
;		around the edge.
;
; MODIFICATION HISTORY:
;	Created 1-MARCH-1993 by Mati Meron.
;	Modified 1-MAY-1994 by Mati Meron.  Default read data directory set to
;	IDL_DIR:[USER]
;       Complete rewrite 6-SEP-1995 by Mark Rivers to make the routine portable
;       to non-VMS platforms.  Switched to use of XDR, environment variables.
;-

    common sxr_stuff, sorlist, abctab, enun, denun, curbeam

    if (n_elements(abctab) ne 0) then return
    if (n_elements(filename) eq 0) then filename = getenv('IDL_ABS_COEFFS')
    if (n_elements(filename) eq 0) then message, $
        'Must either define environment variable IDL_ABS_COEFFS' + $
        'or specify filename'
    restore, filename
end
