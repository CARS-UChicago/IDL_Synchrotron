Pro save_abs_coeffs, infile, outfile

;+
; NAME:
;	SAVE_ABS_COEFFS
; PURPOSE:
;	Reads an input ASCII file of absorption coefficients and writes it out
;       as a binary XDR file.
; CATEGORY:
;	X-ray calculations
; CALLING SEQUENCE:
;	SAVE_ABS_COEFFS INFILE [, OUTFILE]
; INPUTS:
;	INFILE:
;           The name of the ASCII file containing the absorption coefficients
; OPTIONAL INPUT PARAMETERS:
;    OUTFILE
;	Name of output XDR binary file.  Default is the environment variable
;       IDL_ABS_COEFFS
; OUTPUTS:
;	None.
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
;	None.
; PROCEDURE:
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
;	Uses calls to DEFAULT, STRPARSE, STRMATCH, SPLIN_COEFFS and SPLIN_EVAL 
;	in MIDL.
; MODIFICATION HISTORY:
;	Created 1-MARCH-1993 by Mati Meron.
;	Modified 1-MAY-1994 by Mati Meron.  Default read data directory set to
;	IDL_DIR:[USER]
;       Modified 1-SEP-1995 by Mark Rivers to make this routine only write
;       the new file. LOAD_ABS_COEFFS.PRO loads the binary file.
;-

    common sxr_stuff, sorlist, abctab, enun, denun, curbeam

    if (n_elements(outfile) eq 0) then outfile = getenv('IDL_ABS_COEFFS')
    if (n_elements(outfile) eq 0) then message, $
        'Must either define environment variable IDL_ABS_COEFFS' + $
        'or specify filename'

    cmx = 64
    emx = 16
    denun = Default(denun,'keV')
    enun  = Default(enun,denun)
    blank = {abco, name: '', csym: '', z: 0, a: 0., ro: 0., $
	colen: 0, cotab: fltarr(cmx,3), edlen: 0, edtab: fltarr(emx,3)}
    kwords = ['begin', 'end']

    on_ioerror, file_no_good
    openr, datun, infile, /get_lun
    nli = -1
    datfl = 0
    line = ''
    while not eof(datun) do begin
        readf, datun, line
        if line ne '' then begin
        mlast = Strparse(line,'	 ',list)
        mnum = Strmatch_mm(list(0),kwords,3)
        case (datfl + mnum + 1) of
            1:begin
                nli = nli + 1
                if nli eq 0 then abctab=blank else abctab=[abctab,blank]
                abctab(nli).name	= list(1)
                abctab(nli).csym	= list(2)
                abctab(nli).z	        = list(3)
                abctab(nli).a	        = list(4)
                abctab(nli).ro	        = list(5)
                hold = fltarr(cmx,2)
                datfl = 3
                npo = -1
            end
            2:message, 'END without BEGIN, aborting!'
            3:begin
                npo = npo + 1
                if npo eq cmx then message, 'Data overflow!'
                hold(npo,0:1) = float(list(0:1))
            end
            4:message, 'BEGIN without END, aborting!'
            5:begin
                abctab(nli).colen = npo
                use = alog(hold(0:npo,*))
                abctab(nli).cotab(0:npo,*) = $
                Splin_coeffs(use(*,0),use(*,1), /segment)
                ed = where(hold(0:npo-1,0) eq hold(1:npo,0), ned)
                abctab(nli).edlen = ned
                if ned eq 0 then ed = 0 else ed = [0,ed]
                abctab(nli).edtab(0:ned,0) = hold(ed,0)
                abctab(nli).edtab(0:ned,1) = hold(ed+1,1) - hold(ed,1)
                abctab(nli).edtab(0,1) = hold(0,1)
                led = use(ed,0)
                eps = 1e-5
                abctab(nli).edtab(0:ned,2) = -2. - 0.5* $
                (Splin_eval(led+eps,abctab(nli).cotab(0:npo,*),der=1)+ $
                Splin_eval(led-eps,abctab(nli).cotab(0:npo,*),der=1))
                datfl = 0
            end
            else:
            endcase
        endif
    endwhile
    free_lun, datun
    abctab.name = strcompress(abctab.name,/remove_all)
    abctab.csym = strcompress(abctab.csym,/remove_all)
    save, file=outfile, /XDR, abctab, enun, denun
    return

    file_no_good:
    free_lun, datun
    return
end
