Function Grav_sag, sig, tetg, range = ran, consig = csig, precision = prec, $
	show = sho, full = ful, nosag = nos, _extra = _e

;+
; NAME:
;		GRAV_SAG
; VERSION:
;		6.2
; PURPOSE:
;		Calculates reflection broadening due to mirror gravitational sag.
; CATEGORY:
;		Mathematical, optics specific.
; CALLING SEQUENCE:
;		Result = GRAV_SAG( SIG, TETG, RANGE = RAN [, optional keywords])
; INPUTS:
;	SIG
;		The "sigma" of the angular distribution (assumed gaussian) of the
;	incoming radiation.
;	TETG
;		Gravitational sag angle, a parameter characterizing the magnitude of
;		the sag.
; OPTIONAL INPUT PARAMETERS:
;		None.
; KEYWORD PARAMETERS:
;	RANGE
;		Scalar, the angular range to be used in evaluation.  More precisely,
;		the evaluation is over the angular range [-RANGE,RANGE].  Mandatory.
;	CONSIG
;		The "sigma" of additional angular distribution (also assumed gaussian)
;		to be convoluted with the result.  Typically represents the angular
;		distribution resulting from the source's finite angular size.
;	PRECISION
;		Character scalar specifying the evaluation precision.  The three
;		possible values are "low", "med", "hig".  Default is "med".
;	/SHO
;		Switch.  If set, graphical display of calculation results is given.
;	/FULL
;		Switch.  Output modifier, see OUTPUTS below.
;	NOSAG
;		Optional output, see below.
; OUTPUTS:
;		In standard use returns the integral of the sag modified distribution
;		over the range [-RANGE,RANGE].  If the keyword /FULL is set, the whole
;		distribution is returned in the form of a [2,N] array, with the first
;		column being angle values and the second distribution values.
; OPTIONAL OUTPUT PARAMETERS:
;	NOSAG
;		Returns the results for the no-sag distribution i.e. the gaussian
;		distribution defined by CONSIG (see above).  This output is of the
;		same format as the main output and is modified by the keyword /FULL
;		in the same way.
; COMMON BLOCKS:
;		None.
; SIDE EFFECTS:
;		None.
; RESTRICTIONS:
;		None.
; PROCEDURE:
;		Generates the distribution as described in "Gravity sagging of mirrors"
;		and convolutes it as needed.  Calls CALCTYPE, CAST, CONVOL_MM, DEFAULT,
;		DIF, ERRORF_MM, FPU_FIX, MAKE_GRID, PREP_FT, SIGN and STRMATCH_MM, from
;		MIDL
; MODIFICATION HISTORY:
;		Created 10-SEP-2005 by Mati Meron.
;		Modified 10-JUL-2007 by Mati Meron.  Minor streamlining.
;-

	on_error, 1

	posib = ['low','med','hig']
	ppnt = [100l,1000l,10000l]
	whi = abs(Strmatch_mm(prec,posib,3))
	pnt = ppnt(whi)

	typ = Calctype(sig,tetg,0.)
	wsig = Cast(sig,5)
	wtetg = Cast(tetg,5)
	cfl = keyword_set(csig)
	wcsig = Default(csig,0d,/dtyp)
	psig = wsig^3/wtetg^2
	ssig = sqrt(psig^2 + wcsig^2)

	tet = Make_grid(ran*[-1d,1d],Prep_ft(2*round(ran*pnt/ssig)+1),dim=npo)
	tem = Sign(tet)*(abs(tet)/psig)^(1./3)

	dtet = Dif(tet,/cen,/lin)
	dist = Dif(Errorf_mm(tem/sqrt(2))/2,/cen,/lin)/dtet

	if cfl then begin
		cdist = Dif(Errorf_mm(tet/(wcsig*sqrt(2)))/2,/cen,/lin)/dtet
		if psig gt wcsig then dist = Convol_mm(dist,dtet*cdist,/edge_truncate) $
		else dist = FPU_fix(Convol_mm(dtet*cdist,dist,/edge_truncate))
	endif else begin
		cdist = 0*dist
		dum = min(abs(tet),cind)
		cdist[cind] = 1/dtet[cind]
	endelse

	if keyword_set(sho) then begin
		yran = [min(dist),max(cdist)]
		plot, tet, dist, yran = yran, _extra = _e
		oplot, tet, cdist, line = 2
	endif

	nosfl = arg_present(nos)
	if keyword_set(ful) then begin
		res = Cast(transpose([[tet],[dist]]),typ,typ)
		if nosfl then nos = Cast(transpose([[tet],[cdist]]),typ,typ)
	endif else begin
		res = Cast(total(dist*dtet),typ,typ)
		if nosfl then nos = Cast(total(cdist*dtet),typ,typ)
	endelse

	return, res
end