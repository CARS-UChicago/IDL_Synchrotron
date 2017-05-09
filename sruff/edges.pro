Pro Edges, elem

;+
; NAME:
;		EDGES
; VERSION:
;		7.12
; PURPOSE:
;		Displays absorption edges for the selected element
; CATEGORY:
;		X-ray utility
; CALLING SEQUENCE:
;		EDGES, ELEM
; INPUTS:
;	ELEM
;		Scalar input, two possibilities:
;			1)	Numeric input - rounded to integer and taken as a Z-value.
;			2)	Character input - if the length is <= 2, taken as chemical
;				symbol, else taken as element name (full or partial).
; OPTIONAL INPUT PARAMETERS:
;		None.
; KEYWORD PARAMETERS:
;		None.
; OUTPUTS:
;		None other than screen output.
; OPTIONAL OUTPUT PARAMETERS:
;		None.
; COMMON BLOCKS:
;		SXR_STUFF.  See LOAD_ABS_COEFFS for more information.
; SIDE EFFECTS:
;		None.
; RESTRICTIONS:
;		None.
; PROCEDURE:
;		Straightforward search in the table ABCTAB.  Calls LOAD_ABS_COEFFS. 
;		Calls STRMATCH_MM, TABULATE and TYPE, from MIDL.
; MODIFICATION HISTORY:
;		Created 25-NOV-2009 by Mati Meron.
;-

	common sxr_stuff, sorlist, abctab, enun, denun, curbeam

	on_error, 1
	Load_abs_coeffs
	edlis = ['','K ','L1','L2','L3','M1','M2','M3','M4','M5',$
			 'N1','N2','N3','N4','N5','N6','N7']

	case n_elements(elem) of
		0	:	message, 'Missing input!'
		1	:	begin
					if Type(elem) eq 7 then begin
						slen = strlen(elem)
						if slen le 2 then comp = abctab.csym $
						else comp = abctab.name
						ind = Strmatch_mm(elem,comp,slen>2,/all,num=num)
						case num of
							0	:	message, 'Not found!'
							1	:	ind = ind[0]
							else:	message, 'Not unique!
						endcase
					endif else begin
						ind = where(abctab.z eq round(elem))
						if ind lt 0 then message, 'Not in table!'
					endelse
				end
		else:	message, 'Only scalar inputs accepted!'
	endcase

	print
	edl = abctab[ind].edlen
	if edl gt 0 then begin
		enam = reverse(edlis[1:edl])
		eval = (abctab[ind].edtab)[1:edl]
		tabulate, enam, eval, form = ['a8','f7.3'], $
		head = ['Edge', 'Energy (keV)'], tit = abctab[ind].name + ' edges'
	endif else print, '		No edges in table for ' + abctab[ind].name
	print

	return
end