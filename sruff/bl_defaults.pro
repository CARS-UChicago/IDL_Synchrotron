Pro BL_Defaults, synch = sync, energy = rene, gamma = rgam, current = curr,$
	bm_field = bmmf, dev_length = dlen, min_gap = mgap, $
	radial_size= rsig, angular_size= asig, aperture = aper, distance = dist

;+
; NAME:
;		BL_DEFAULTS
; VERSION:
;		5.0
; PURPOSE:
;		Allows for setting (or, alternatively, using defaults) values of
;		machine specific parameters used in synchrotron calculations.
; CATEGORY:
;		Utility.
; CALLING SEQUENCE:
;		BL_DEFAULTS [, keywords]
; INPUTS:
;		None.
; OPTIONAL INPUT PARAMETERS:
;		None.
; KEYWORD PARAMETERS
;	SYNCH
;		Allows to enter character value, a synchrotron name.  Currently only
;		'APS' is accepted.
;	ENERGY
;		An input or output parameter.  Ring energy, in GeV.	| At input, only
;	GAMMA													| one of these may
;		An input or output parameter.  Ring /Gamma.  		| be given.
;	CURRENT
;		An input or output parameter.  Ring current in Amp.
;	BM_FIELD
;		An input or output parameter.  Bending magnet field in Tesla.
;	DEV_LENGTH
;		An input or output parameter.  Standard length of an insertion device,
;		in meters.
;	MIN_GAP
;		An input or output parameter.  Minimal insertion device gap, in mm.
;	RADIAL_SIZE
;		An input or output parameter.  Two element vector containing the (sigma)
;		radial size of the source, in [x,y] order, in milimeters.
;	ANGULAR_SIZE
;		An input or output parameter.  Two element vector containing the (sigma)
;		angular size of the source, in [x,y] order, in miliradians.
;
;		Note:  RADIAL_SIZE and ANGULAR_SIZE shoud be 2-element vectors.
;	APERTURE
;		An input or output parameter.  Standard aperture, given as a 2-element
;		vector, in milimeters.
;	DISTANCE
;		An input or output parameter.  Distance of standard aperture from
;		source, in meters.
; OUTPUTS:
;		None.
; OPTIONAL OUTPUT PARAMETERS:
;		Any of the keyword parameters, other than SYNCH, when provided with
;		named (and undefined) variable, returns the approriate current default
;		value, as defined in the structure !BLPAR.  See the routine BL_PARS for
;		details.  If the provided variable is defined, it is returned as is.
; COMMON BLOCKS:
;		None.
; SIDE EFFECTS:
;		If a synchrotron name is given, through SYNCH, the structure !BLPAR is
;		reinitiated with the default values of the appropriate machine.  Else,
;		there are no side effects.
; RESTRICTIONS:
;		None.
; PROCEDURE:
;		Checks the inputs and, if any of them is a named (and undefined)
;		variable, it is replaced by the proper default.  Defined variables
;		and constants are not changed.  Calls BL_PARS.  Calls DEFAULT, ONE_OF
;		and TYPE from MIDL.
;		Note:	Since the parameters ENERGY and GAMMA are not independent,
;		providing one will change the other one and providing both will result
;		in error.
; MODIFICATION HISTORY:
;		Created 25-JUN-2005 by Mati Meron.
;		Modified 30-DEC-2006 by Mati Meron.  Some internal changes.  Added
;		keywords DEV_LENGTH, MIN_GAP, RADIAL_SIZE and ANGULAR_SIZE.
;-

	on_error, 1

	if Type(sync) eq 7 then BL_pars, sync = sync

	case One_of(rene,rgam) of
		-1	:	begin
					rene = !blpar.rene
					rgam = !blpar.rgam
				end
		0	:	rgam = float(1e6*rene/(!srcon.ee*!srcon.scal))
		1	:	rene = float(1e-6*rgam*(!srcon.ee*!srcon.scal))
	endcase

	curr = Default(curr,!blpar.curr,/dtyp)
	bmmf = Default(bmmf,!blpar.bmmf,/dtyp)
	aper = Default(aper,!blpar.aper,/dtyp)
	dist = Default(dist,!blpar.dist,/dtyp)
	dlen = Default(dlen,!blpar.dlen,/dtyp)
	mgap = Default(mgap,!blpar.mgap,/dtyp)
	rsig = Default(rsig,!blpar.rsig,/dtyp)
	asig = Default(asig,!blpar.asig,/dtyp)

	return
end
