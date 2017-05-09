Pro Multi_gen, mult, formula = form, $
	toplay = tlay, tden = tden, twei = twei, tthi = tthi, trough = trof,$
	flayer = flay, fden = fden, fwei = fwei, fthi = fthi, frough = frof, $
	slayer = slay, sden = sden, swei = swei, sthi = sthi, srough = srof, $
	botlay = blay, bden = bden, bwei = bwei, bthi = bthi, brough = brof, $
	elements = elar, num_elems = numels, dens = rarr, weights = warr, $
	thicks = tarr, roughness = roff

;+
; NAME:
;		MULTIGEN
; VERSION:
;		5.5
; PURPOSE:
;		Generates multilayer parameters for use with MIRROR or REFLECT
; CATEGORY:
;		X-ray calculations.
; CALLING SEQUENCE:
;		MULTIGEN, MULT, FORMULA = FORM, $
;		[TOPLAY = TLAY, TDEN = TDEN, TWEI = TWEI, TTHI = TTHI, TROUGH = TROF,] $
;		FLAYER = FLAY, FDEN = FDEN, FWEI = FWEI, FTHI = FTHI, [FROUGH = FROF,] $
;		SLAYER = SLAY, SDEN = SDEN, SWEI = SWEI, STHI = STHI, [SROUGH = SROF,] $
;		[BOTLAY = BLAY, BDEN = BDEN, BWEI = BWEI, BTHI = BTHI, BROUGH = BROF,] $
;		ELEMENTS = ELAR, NUM_ELEMS = NUMELS, DENS = RARR, WEIGHTS = WARR, $
;		THICKS = TARR, ROUGHNESS = ROFF
; INPUTS:
;	MULT
;		Multiplicity, the number of repetitions in the multilayer.  Positive
;		integer, mandatory.
; OPTIONAL INPUT PARAMETERS:
;		None.
; KEYWORD PARAMETERS:
;	/FORMULA
;		Switch.  If set, specifies that all the weights to follow are "formula
;		weights" (see MIRROR or REFLECT for more details).  Irrelevant to the
;		operation of MULTIGEN, other that it sets the FORM flag in the common
;		block, for the benefit of any routine which may use the results of
;		MULTIGEN.
;	TOPLAY
;		Accepts a list of one or more elements which comprise the optional top
;		layer (coating) of the assembly.  This layer, if provided, is distinct
;		from the multilayer itself.  Provided as character array (each element
;		represented by its chemical symbol).
;	TDEN
;		The density of the top layer (in gr/cm^3).  Not needed (though allowed)
;		when TOPLAY consists of a single element, mandatory otherwise.
;	TWEI
;		The relative weights of the elements in TOPLAY.  Not needed (and ignored
;		if provided) when TOPLAY consists of a single element, mandatory
;		otherwise.  Provided as numerical array of same length as TOPLAY.
;
;		Note:	The routines MIRROR and REFLECT distinguish between weights by
;				mass and weights by chemical formula (see keyword /FORMULA
;				there).  MULTIGEN isn't concerned with this but the user should
;				rememeber which weights ewere intended when using the output of
;				MULTIGEN.
;	TTHI
;		The thickness of the top layer, in Angstrem.  Numeric scalar,mandatory.
;	TROUGH
;		The roughness of the top layer, in Angstrem.  Numeric scalar, optional.
;	FLAYER
;		Same as TOPLAY, for the first layer of the multilayer.  Mandatory.
;	FDEN
;		Same as TDEN, for the first layer of the multilayer.
;	FWEI
;		Same as TDEN, for the first layer of the multilayer.
;	FTHI
;		Same as TTHI, for the first layer of the multilayer.  Mandatory.
;	FROUGH
;		Same as TROUGH, for the first layer of the multilayer.
;	SLAYER
;		Same as TOPLAY, for the second layer of the multilayer.  Mandatory.
;	SDEN
;		Same as TDEN, for the second layer of the multilayer.
;	SWEI
;		Same as TDEN, for the second layer of the multilayer.
;	STHI
;		Same as TTHI, for the second layer of the multilayer.  Mandatory.
;	SROUGH
;		Same as TROUGH, for the second layer of the multilayer.
;	BOTLAY
;		Same as TOPLAY, for the bottom layer (substrate) of the multilayer.
;		Optional.
;	BDEN
;		Same as TDEN, for the bottom layer of the multilayer.
;	BWEI
;		Same as TDEN, for the bottom layer of the multilayer.
;	BTHI
;		Same as TTHI, for the bottom layer of the multilayer.  Optional,
;		defaults to 0.  The value doesn't matter except if the multilayer is
;		used for transmission.
;	BROUGH
;		Same as TROUGH, for the second layer of the multilayer.
;
;	Note:	In the case of the top and bottom layers, their presence is detected
;			through the presence of TOPLAY and BOTLAY, respectively.  If any of
;			these is absent, the subsequent parameters for its layer are ignored
;	ELEMENTS
;		Optional output, see below.
;	NUM_ELEMS
;		Optional output, see below.
;	DENS
;		Optional output, see below.
;	WEIGHTS
;		Optional output, see below.
;	THICKS
;		Optional output, see below.
;	ROUGHNESS
;		Optional output, see below.
; OUTPUTS:
;		No direct output other then through optional output parameters (see
;		below).  Fills the common block MULTILAYER_PARS, to be used by other
;		routines.
; OPTIONAL OUTPUT PARAMETERS:
;	ELEMENTS
;		Returns an ELEMENTS vector, representing the multilayer (with the
;		optional coating and substrate), in a form acceptable by MIRROR or
;		REFLECT.  See there for details.
;	NUM_ELEMS
;		Ditto, for NUM_ELEMS.
;	DENS
;		Ditto, for DENS.
;	WEIGHTS
;		Ditto, for WEIGHTS.
;	THICKS
;		Ditto, for THICKS.
;	ROUGHNESS
;		Ditto, for ROUGHNESS.
; COMMON BLOCKS:
;		SXR_STUFF.  See LOAD_ABS_COEFFS for more information.
;		MULTILAYER_PARS.  On return, contains the values of the six optional
;		output parameters listed above, as well as:
;			M_EXS	:	Flag, set to 1 when the common block is defined.
;			M_DEF	:	Flag array, 4 entries corresponding to tob, first layer
;						second layer and bottom, respectively, each set to 1
;						when the corresponding layer has been defined.
;			M_FORM	:	Flag, value of 1 indicating the the weights are to be
;						taken as formula weights.
; SIDE EFFECTS:
;		None.
; RESTRICTIONS:
;		Only elements defined in ABCTAB may be used.
; PROCEDURE:
;		Straightforward.  Uses density values from ABCTAB.  Calls ELECOMP and
;		LOAD_ABS_COEFFS from SRUFF.  Also calls ARREQ and DEFAULT from MIDL.
; MODIFICATION HISTORY:
;		Created 1-JUN-2006 by Mati Meron.
;-

	common sxr_stuff, sorlist, abctab, enun, denun, curbeam
	common multilayer_pars, m_exs, m_def, m_form, $
	m_elar, m_numels, m_rarr, m_warr, m_tarr, m_roff

	on_error, 1
	Load_abs_coeffs
	m_exs = 0
	m_def = intarr(4)

	fnum = n_elements(flay)
	frof = Default(frof,0.,/dtyp)
	case fnum of
		0	:	message, "missing first layer composition data!'
		1	:	begin
					fden = Default(fden,abctab[Elecomp(flay)].ro)
					fwei = 1.
				end
		else:
	endcase
	check= [n_elements(fwei),n_elements(fden),n_elements(fthi),n_elements(frof)]
	if Arreq(check,[fnum,1,1,1]) then m_def[1] = 1 $
	else message, 'Bad first layer data!'

	snum = n_elements(slay)
	srof = Default(srof,0.,/dtyp)
	case snum of
		0	:	message, "missing first layer composition data!'
		1	:	begin
					sden = Default(sden,abctab[Elecomp(slay)].ro)
					swei = 1.
				end
		else:
	endcase
	check= [n_elements(swei),n_elements(sden),n_elements(sthi),n_elements(srof)]
	if Arreq(check,[snum,1,1,1]) then m_def[2] = 1 $
	else message, 'Bad first layer data!'

	m_elar = strarr((fnum + snum)*mult)
	m_warr = fltarr((fnum + snum)*mult)
	m_numels = lonarr(2l*mult)
	m_rarr = (m_tarr = (m_roff = fltarr(2l*mult)))

	flo = (fnum + snum)*lindgen(mult)
	fhi = flo + fnum - 1
	slo = fhi + 1
	shi = slo + snum - 1
	for i = 0l, mult-1 do begin
		m_elar[flo[i]:fhi[i]] = flay
		m_elar[slo[i]:shi[i]] = slay
		m_warr[flo[i]:fhi[i]] = fwei
		m_warr[slo[i]:shi[i]] = swei

		m_numels[2*i] = fnum
		m_numels[2*i+1] = snum
		m_rarr[2*i] = fden
		m_rarr[2*i+1] = sden
		m_tarr[2*i] = fthi
		m_tarr[2*i+1] = sthi
		m_roff[2*i] = frof
		m_roff[2*i+1] = srof
	endfor

	tnum = n_elements(tlay)
	if tnum gt 0 then begin
		trof = Default(trof,0.,/dtyp)
		if tnum eq 1 then begin
			tden = Default(tden,abctab[Elecomp(tlay)].ro)
			twei = 1.
		endif
		check = $
		[n_elements(twei),n_elements(tden),n_elements(tthi),n_elements(trof)]
		if Arreq(check,[tnum,1,1,1]) then m_def[0] = 1 $
		else message, 'Bad top layer data!'
		m_elar = [tlay,m_elar]
		m_warr = [twei,m_warr]
		m_numels = [tnum,m_numels]
		m_rarr = [tden,m_rarr]
		m_tarr = [tthi,m_tarr]
		m_roff = [trof,m_roff]
	endif

	bnum = n_elements(blay)
	if bnum gt 0 then begin
		brof = Default(brof,0.,/dtyp)
		bthi = Default(bthi,0.,/dtyp)
		if bnum eq 1 then begin
			bden = Default(bden,abctab[Elecomp(blay)].ro)
			bwei = 1.
		endif
		check = $
		[n_elements(bwei),n_elements(bden),n_elements(bthi),n_elements(brof)]
		if Arreq(check,[bnum,1,1,1]) then m_def[3] = 1 $
		else message, 'Bad top layer data!'
		m_elar = [m_elar,blay]
		m_warr = [m_warr,bwei]
		m_numels = [m_numels,bnum]
		m_rarr = [m_rarr,bden]
		m_tarr = [m_tarr,bthi]
		m_roff = [m_roff,brof]
	endif

	m_exs = 1
	m_form = keyword_set(form)
	elar = m_elar
	numels = m_numels
	rarr = m_rarr
	warr = m_warr
	tarr = m_tarr
	roff = m_roff

	return
end






