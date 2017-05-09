Pro ID_conv, field= bfl, gap= gap, period= per, k= kvl, emax= emx, ene= ene, $
	mask= msk, ring_energy= ren, rgam= rgm, enhanced= enh, show= sho, num= ndat

;+
; NAME:
;		ID_CONV
; VERSION:
;		5.7
; PURPOSE:
;		Conversion between various insertion device parameters.
; CATEGORY:
;		SR specific
; CALLING SEQUENCE:
;		Result = ID_CONV, keywords
; INPUTS:
;		None
; OPTIONAL INPUT PARAMETERS:
;		None.
; KEYWORD PARAMETERS:
;	FIELD										|	These 6 parameters serve
;		The magnetic, in Tesla.					|	as both input and output
;	GAP											|	parameters.  Only two can
;		ID gap, in mm.							|	be provided as inputs, at
;	PERIOD										|	any given time (unless
;		ID period, in mm.						|	MASK is used, see below),
;	K											|	but all six will be outputs.
;		The K value of the ID.					|
;	EMAX										|
;		The maximal (saturation) energy, in keV.|
;	ENE											|
;		Current energy, in keV.					|
;
;	Note:	One (and no more than one) of the inputs may be an array.  In such
;			case all outputs will be arrays of same size.
;	MASK
;		An optional mask which can be used to force ID_CONV to ignore some
;		of the inputs.  Provided as numerical array, any 0 value translates to
;		"ignore" (and any nonzero value to "accept").  See the routine HOW_MANY
;		in MIDL_LIB for details.
;
;		Alternatively, the mask can also be provided as a character array,
;		specifying explicitly the two inputs to be used, out of the list
;		['FIELD','GAP','PERIOD','K','EMAX','ENE'].  Only the first 2 characters
;		are needed and case or order do not matter.
;
;		Note:	When MASK is used, one can have more than 2 defined inputs.
;				The number of inputs which are both "defined" and "accepted"
;				(i.e. corresponding to nonzero MASK entries) must still be 2.
;	RING_ENERGY													|
;		Storage ring energy, in GeV.  Optional, if not provided,| One and only
;		 the APS value (7 GeV) is used.							| one of these
;	RGAM														| may be used.
;		The relativistic Gamma of the stored electron.  Optional|
;	/ENHANCED
;		Switch.  If set, the calculation is for an "enhanced field", higher by
;		about 2% than the usual.  Details elsewhere.
;	/SHOW
;		Switch.  If set, the calculated results are displayed to the screen,
;		in tabular form.
;	NUM
;		Optional output, see below.
; OUTPUTS:
;		See Optional Output Parameters, below.
; OPTIONAL OUTPUT PARAMETERS:
;	FIELD, GAP, PERIOD, K, EMAX, ENE
;		All six input parameters serve as output parameters as well.  Given
;		any two inputs, the cooresponding values for the remaining four
;		parameters are calculated and returned.  If one of the inputs is an
;		array, all the calculated parameters are returned in this format and
;		the second (scaler) input is reformated as an array as well.
;	NUM
;		Returns the length (number of elements) of the data.  As mentioned
;		above, all six return variables have same length.
; COMMON BLOCKS:
;		None.
; SIDE EFFECTS:
;		As mentioned above, when one of the inputs is a scalar and the other an
;		arrray, the scalar is reformated as an array.  Also, on return all
;		parameters are of type no lower than 4 (float).
; RESTRICTIONS:
;		1)	The calculations of field values have a restricted range (see
;			ID_FIELD for details).
;		2)	Combinations of inputs which are equivalent to ENE > EMAX are
;			unphysical and do not allow for a calculation.
;		3)	The specific input combination containing PERIOD and EMAX is
;			inadequate as these two parameters are not independent.
;
;		In the first 2 cases no error message will be issued but NanN will be
;		returned for incalculable values.  In the last case, the routine will
;		exit with an error message.
; PROCEDURE:
;		Follows the basic definitions OF SR quantities.  Uses BL_DEFAULTS and
;		ID_FIELD from SRUFF_LIB, as well as ID_CONV_FUN (from same location)
;		for internal evaluations.  Calls CAST, HOW_MANY, ONE_OF, ROOT and
;		TABULATE from MIDL.
; MODIFICATION HISTORY:
;		Created 20-SEP-2005 by Mati Meron.
;		Modified 15-NOV-2005 by Mati Meron.  Added option of string valued mask.
;		Modified 20-NOV-2006 by Mati Meron.  Internal changes.
;-

	on_error, 1
	mult = 1e3
	posib = ['field','gap','period','k','emax','ene']

	case One_of(ren,rgm,val=val) of
		-1	:	BL_defaults, gam = wrgm
		0	:	wrgm = 1e6*ren/(!srcon.scal*!srcon.ee)
		1	:	wrgm = rgm
	endcase
	a = float(2*mult*!pi*!srcon.bere)
	b = float(4*mult*!pi*wrgm^2*!srcon.scal*!srcon.re*!srcon.ee/!srcon.alp)

	if Type(msk) eq 7 then begin
		wmsk = intarr(n_elements(posib))
		for i = 0, n_elements(msk) - 1 do begin
			j = Strmatch_mm(msk[i],posib,2)
			if j ge 0 then wmsk[j] = 1
		endfor
	endif else if Isnum(msk) then wmsk = msk

	check = [n_elements(bfl),n_elements(gap),n_elements(per),$
			 n_elements(kvl),n_elements(emx),n_elements(ene)]
	if Isnum(wmsk) then check = check*wmsk
	dum = where(check gt 1, ndum)
	if ndum le 1 then ndat = max(check) $
	else message, 'At most one input can be an array!'
	if ndat gt 1 then begin
		cc = check eq 1
		if cc[0] then bfl = replicate(bfl,ndat)
		if cc[1] then gap = replicate(gap,ndat)
		if cc[2] then per = replicate(per,ndat)
		if cc[3] then kvl = replicate(kvl,ndat)
		if cc[4] then emx = replicate(emx,ndat)
		if cc[5] then ene = replicate(ene,ndat)
	endif

	wha = $
	How_many(fir=bfl,sec=gap,thi=per,fou=kvl,fif=emx,six=ene,mask=wmsk,whi=whi)
	case whi of
		3	:	begin
					bfl = Cast(bfl,4)
					gap = Cast(gap,4)
					per = gap/ID_field(bfl,/inverse,enhanced=enh,/qui)
					kvl = bfl*per/a
					emx = b/per
					ene = emx/(1 + kvl^2/2)
				end
		5	:	begin
					bfl = Cast(bfl,4)
					per = Cast(per,4)
					gap = per*ID_field(bfl,/inverse,enhanced=enh,/qui)
					kvl = bfl*per/a
					emx = b/per
					ene = emx/(1 + kvl^2/2)
				end
		6	:	begin
					gap = Cast(gap,4)
					per = Cast(per,4)
					bfl = ID_field(gap/per,enhanced=enh,/qui)
					kvl = bfl*per/a
					emx = b/per
					ene = emx/(1 + kvl^2/2)
				end
		9	:	begin
					bfl = Cast(bfl,4)
					kvl = Cast(kvl,4)
					per = a*kvl/bfl
					gap = per*ID_field(bfl,/inverse,enhanced=enh,/qui)
					emx = b/per
					ene = emx/(1 + kvl^2/2)
				end
		10	:	begin
					gap = Cast(gap,4)
					kvl = Cast(kvl,4)
					par1 = gap/(a*kvl)
					bfl = 0*par1
					dum = ID_field(1,/inv,enh=enh,ran=ran)
					par = [whi,0.,0.,keyword_set(enh)]
					for i = 0l, ndat - 1 do begin
						par[1] = par1[i]
						bfl[i]= Root('id_conv_fun',ran,par=par,stat=stat)
						if not stat then bfl[i] = !values.f_nan
					endfor
					per = a*kvl/bfl
					emx = b/per
					ene = emx/(1 + kvl^2/2)
				end
		12	:	begin
					per = Cast(per,4)
					kvl = Cast(kvl,4)
					bfl = a*kvl/per
					gap = per*ID_field(bfl,/inverse,enhanced=enh,/qui)
					emx = b/per
					ene = emx/(1 + kvl^2/2)
				end
		17	:	begin
					bfl = Cast(bfl,4)
					emx = Cast(emx,4)
					per = b/emx
					gap = per*ID_field(bfl,/inverse,enhanced=enh,/qui)
					kvl = bfl*per/a
					ene = emx/(1 + kvl^2/2)
				end
		18	:	begin
					gap = Cast(gap,4)
					emx = Cast(emx,4)
					per = b/emx
					bfl = ID_field(gap/per,enhanced=enh,/qui)
					kvl = bfl*per/a
					ene = emx/(1 + kvl^2/2)
				end
		20	:	message, 'No solution for this case!'
		24	:	begin
					kvl = Cast(kvl,4)
					emx = Cast(emx,4)
					per = b/emx
					bfl = a*kvl/per
					gap = per*ID_field(bfl,/inverse,enhanced=enh,/qui)
					ene = emx/(1 + kvl^2/2)
				end
		33	:	begin
					bfl = Cast(bfl,4)
					ene = Cast(ene,4)
					kvl = 2*sqrt(2./3)*sinh($
					imaginary(asin(complex(0,(3./2)^(3./2)*b*bfl/(a*ene))))/3)
					per = a*kvl/bfl
					gap = per*ID_field(bfl,/inverse,enhanced=enh,/qui)
					emx = b/per
				end
		34	:	begin
					gap = Cast(gap,4)
					ene = Cast(ene,4)
					par1 = gap*ene/b
					par2 = gap/a
					bfl = 0*par1
					dum = ID_field(1,/inv,enh=enh,ran=ran)
					par = [whi,0.,0.,keyword_set(enh)]
					for i = 0l, ndat - 1 do begin
						par[1:2] = [par1[i],par2[i]]
						bfl[i] = Root('id_conv_fun',ran,par=par,stat=stat)
						if not stat then bfl[i] = !values.f_nan
					endfor
					per = gap/ID_field(bfl,/inverse,enhanced=enh,/qui)
					kvl = bfl*per/a
					emx = b/per
				end
		36	:	begin
					per = Cast(per,4)
					ene = Cast(ene,4)
					kvl = (bfl = (gap = (emx = fltarr(ndat))))
					arg = 2*(b/(per*ene) - 1)
					good = where(arg ge 0, ngood,comp=bad,ncomp=nbad)
					if ngood gt 0 then begin
						kvl[good] = sqrt(arg[good])
						bfl[good] = a*kvl[good]/per[good]
						gap[good] = per[good]*$
							ID_field(bfl[good],/inverse,enhanced=enh,/qui)
						emx[good] = b/per[good]
					endif
					if nbad gt 0 then $
					kvl[bad]=(bfl[bad]=(gap[bad]=(emx[bad]=!values.f_nan)))
				end
		40	:	begin
					kvl = Cast(kvl,4)
					ene = Cast(ene,4)
					per = b/(ene*(1 + kvl^2/2.))
					bfl = a*kvl/per
					gap = per*ID_field(bfl,/inverse,enhanced=enh,/qui)
					emx = b/per
				end
		48	:	begin
					emx = Cast(emx,4)
					ene = Cast(ene,4)
					bfl = (gap = (per = (kvl = fltarr(ndat))))
					arg = 2*(emx/ene - 1)
					good = where(arg ge 0, ngood,comp=bad,ncomp=nbad)
					if ngood gt 0 then begin
						kvl[good] = sqrt(arg[good])
						per[good] = b/emx[good]
						bfl[good] = a*kvl[good]/per[good]
						gap[good] = per[good]*$
							ID_field(bfl[good],/inverse,enhanced=enh,/qui)
					endif
					if nbad gt 0 then $
					bfl[bad]=(gap[bad]=(per[bad]=(kvl[bad]=!values.f_nan)))
				end
		else:	message, 'Unacceptable input, 2 and only 2 params needed!'
	endcase

	if ndat eq 1 then begin
		bfl = bfl[0]
		gap = gap[0]
		per = per[0]
		kvl = kvl[0]
		emx = emx[0]
		ene = ene[0]
	endif

	if keyword_set(sho) then begin
		print
		Tabulate, bfl, gap, per, kvl, emx, ene, form= replicate('f8.3',6),head=$
		['Field (T)','Gap (mm)','Period (mm)','K  ','Emax (keV)','E (keV)']
	endif

	return
end