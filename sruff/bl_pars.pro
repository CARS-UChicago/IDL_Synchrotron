Pro BL_pars, synch = syn, show = sho

;+
; NAME:
;		BL_PARS
; VERSION:
;		5.2
; PURPOSE:
;		Creates or updates a system variable named !BLPAR.  !BLPAR is a
;		structure the fields of which contain values of physical constants used
;		in source independent machine specific synchrotron calculations, as
;		follows:
;		!BLPAR
;			SYNC  -	Synchrotron's name
;			RENE  -	Ring energy in GeV.
;			RGAM  - The relativistic GAMMA of the machine.
;			CURR  -	The current in the machine in Amperes.
;			BMMF  -	Bending magnet magnetic field, in Tesla.
;			DLEN  -	The standard insertion device length in meters.
;			MGAP  -	Minimal insertion device gap, in mm.
;			RSIG  -	Source (ID) spatial dimensions in mm.  2-element vector
;					in a [sigma_x,sigma_y] format.
;			ASIG  -	Source (ID) angular dimensions in mrad.  2-element vector
;					in a [sigma_x',sigma_y'] format.
;			APER  -	The standard front end aperture size in mm.  2-element
;					vector in a [xsize,ysize] format.
;			DIST  -	The standard distance from source to the aperture in meters.
;		All the numerical values are given in a floating format.
; CATEGORY:
;		Utility.
; CALLING SEQUENCE:
;		BL_PARS
; INPUTS:
;		None.
; OPTIONAL INPUT PARAMETERS:
;		None.
; KEYWORD PARAMETERS
;	/SHOW
;		Switch.  If set, basic information about the structure !BLPAR is shown.
; OUTPUTS:
;		None.
; OPTIONAL OUTPUT PARAMETERS:
;		None.
; COMMON BLOCKS:
;		None.
; SIDE EFFECTS:
;		None.
; RESTRICTIONS:
;		None.
; PROCEDURE:
;		Either creates or, if already exists, updates the system variable
;		!BLPAR.  Calls STRMATCH_MM from MIDL.
; MODIFICATION HISTORY:
;		Created 25-JUN-2005 by Mati Meron.
;		Modified 5-FEB-2006 by Mati Meron.  Added fields RSIG and ASIG.
;		Modified 30-DEC-2006 by Mati Meron.  Added field MGAP.
;-

	on_error, 1

	posib = ['APS']
	energy = [7.]
	current = [.1]
	bm_field = [0.599]
	dev_len = [2.4]
	min_gap = [10.5]
	rsig = [0.2710,0.0086]
	asig = [0.0114,0.0029]
;	rsig = [0.2753,0.0085]
;	asig = [0.0113,0.0030]
	aperture = [3.,2.]
	distance = [25.]

	sfl = keyword_set(syn)
	if sfl then begin
		ind = Strmatch_mm(syn,posib)
		if ind eq -1 then message, 'Unknown Synchrotron!'
	endif else ind = 0

	defsysv, '!blpar', exists = exs
	if exs then if sfl then make = (!blpar.synch ne syn) else make = 0 $
	else make = 1

	if make then begin
		tem = {sync: posib[ind], rene: energy[ind], $
		rgam: float(1e6*energy[ind]/(!srcon.ee*!srcon.scal)), $
		curr: current[ind], bmmf: bm_field[ind], $
		dlen: dev_len[ind], mgap: min_gap[ind], $
		rsig: rsig[*,ind], asig: asig[*,ind], $
		aper: aperture[*,ind], dist: distance[ind]}
		defsysv, '!blpar', tem
	endif
	if keyword_set(sho) then help, /st, !blpar

	return
end