Function Econv, en, from = ceun, to = deun

;+
; NAME:
;	ECONV
; PURPOSE:
;	Unit convertion of energy (or wavelength).
; CATEGORY:
;	X-ray calculations.
; CALLING SEQUENCE:
;	Result = ECONV (EN, FROM = CEUN, TO = DEUN)
; INPUTS:
;    EN
;	Numeric, scalar or array.  Input energy (or wavelength).
; OPTIONAL INPUT PARAMETERS:
;	None.
; KEYWORD PARAMETERS:
;    FROM
;	Character value, specifies input units.  Acceptable units are:
;	KEV, EV, NANOMETERS, ANGSTREM.  Only first 2 letters matter.
;    TO
;	Character value, specifies output units.  Acceptable units are:
;	KEV, EV, NANOMETERS, ANGSTREM.  Only first 2 letters matter.
; OUTPUTS:
;	Returns energy (or wavelength) in the specified output units.
; OPTIONAL OUTPUT PARAMETERS:
;	None.
; COMMON BLOCKS:
;	None.
; SIDE EFFECTS:
;	None.
; RESTRICTIONS:
;	None.
; PROCEDURE:
;	Straightforward.  Uses CAST and STRMATCH from MIDL.
; MODIFICATION HISTORY:
;	Created 1-MARCH-1993 by Mati Meron.
;	Modified 1-JUNE-1994 by Mati Meron.  Added nanometer option.
;-

    on_error, 1
    conv = 12.398424
    posib = ['keV','eV','Nano','Angs']

    fin = StrMatch_mm(ceun,posib,2)
    tin = StrMatch_mm(deun,posib,2)
    if fin lt 0 or tin lt 0 then message, 'Unknown units!'

    if fin ne tin then begin
	case fin of
	    1	:    begin
			fmul = 0.001
			fin = 0
		     end
	    2	:    begin
			fmul = 10.
			fin = 3
		     end
	    else:    fmul = 1.
	endcase
	case tin of
	    1	:    begin
			tmul = 1000.
			tin = 0
		     end
	    2	:    begin
			tmul = 0.1
			tin = 3
		     end
	    else:    tmul = 1.
	endcase
	if fin eq tin then res = (fmul*tmul)*en else res = (conv*tmul/fmul)/en
    endif else res = Cast(en,4,5)

    return, res
end
