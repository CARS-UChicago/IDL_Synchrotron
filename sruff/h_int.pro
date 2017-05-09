Function H_int, e_ec, gplim, pow = pow, prec = prec

;+
; NAME:
;		H_INT
; VERSION:
;		6.0
; PURPOSE:
;		Calculates the (vertical) angle integral of the bending magnet power
;		distribution function H.
; CATEGORY:
;		Mathematical, SR specific
; CALLING SEQUENCE:
;		Result = H_INT (E_EC,GPLIM [, keywords])
; INPUTS:
;	E_EC
;		Numeric scalar or vector.  Photon energy in units of the critical
;		energy E_c.
;	GPLIM
;		Scalar, integration limit in units of /gamma*/psi.  The integral is
;		performed over [-GPLIM,GPLIM].  If not given, an infinite range is used.
; OPTIONAL INPUT PARAMETERS:
;		None.
; KEYWORD PARAMETERS:
;	POW
;		Scalar.  If provided, an integral of H_FUN^POW is calculated.
;	PREC
;		Scalar, optional, specifies (relative) integration precision required.
;		If not given, internal defaults of the integration routine apply.
; OUTPUTS:
;		Returns integration result, type FLOAT or the highest of the types of
;		[E_EC, GPLIM].
; OPTIONAL OUTPUT PARAMETERS:
;		None
; COMMON BLOCKS:
;		None.
; SIDE EFFECTS:
;		None.
; RESTRICTIONS:
;		None.
; PROCEDURE:
;		Straighforward Romberg integration.  Uses CALCTYPE, CAST, DEFAULT,
;		ROMBERG and TOLER from MIDL.  Also uses (indirectly) RH_FUN from
;		SRUFF_LIB.
; MODIFICATION HISTORY:
;		Created 10-JAN-2005 by Mati Meron.
;		Modified 20-JAN-2007 by Mati Meron.  Internal changes.
;-

	on_error, 1

	typ = Calctype(e_ec,gplim,def=4)
	ulim = 4d/sqrt(e_ec > Toler())
	wlim = Default(gplim,ulim,/dtyp) < ulim
	pow = Default(pow,1)

	res = 0.*e_ec
	for i = 0l, n_elements(res) - 1 do res[i] = $
	2*Romberg('rh_fun',[0,wlim[i]],prec,/rel,par=[e_ec[i],pow])

	return, Cast(res,typ,typ)
end