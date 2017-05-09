Function Jsd_fun, x, k, l, negative = neg, positive = pos

;+
; NAME:
;		JSD_FUN
; VERSION:
;		5.6
; PURPOSE:
;		Calculates a linear combination of Bessel functions used in undulator
;		computations.
; CATEGORY:
;		Mathematical, SR specific
; CALLING SEQUENCE:
;		Result = JSD_FUN (X, K [,L] [NEGATIVE = NEG, POSITIVE = POS])
; INPUTS:
;	X
;		Numerical input, arbitrary form.  Mandatory.
;	K
;		Integer scalar, mandatory.
; OPTIONAL INPUT PARAMETERS:
;	L
;		Integer scalar.  If given, must be of same parity as K.  If not given,
;		defaults to (K mod 2).
; KEYWORD PARAMETERS:
;	/NEGATIVE												|	One and only one
;		Switch.  If set, a combination with negative sign	|	of these may be
;		is returned (see PROCEDURE below).					|	given. If non is
;	/POSITIVE												|	specified, the
;		Switch.  If set, a combination with positive sign	|	default is
;		is returned.										|	POSITIVE.
; OUTPUTS:
;		Returns calculation result, type FLOAT or higher (if X is higher).
; OPTIONAL OUTPUT PARAMETERS:
;		None
; COMMON BLOCKS:
;		None.
; SIDE EFFECTS:
;		None.
; RESTRICTIONS:
;		None other than the restrictions on K, L specified above.
; PROCEDURE:
;		Straighforward, evaluates the result as
;
;			Result = (-1)^((K-L)/2)*(J_((K-L)/2)(KX) +- J_((K+L)/2)(KX))
;
;		where the +- sign follows the POSITIVE/NEGATIVE above.
;		Calls DEFAULT, ISNUM and ONE_OF from MIDL.
; MODIFICATION HISTORY:
;		Created 25-DEC-2006 by Mati Meron.
;-

	on_error, 1

	sgn = abs((2*One_of(neg,pos))) - 1
	l = Default(l,k mod 2,/dtyp)
	if Isnum(k,/int) and Isnum(l,/int) and ((k + l) mod 2) eq 0 then $
	res = (-1)^((k-l)/2)*(beselj(k*x,(k-l)/2) + sgn*beselj(k*x,(k+l)/2)) $
	else message, 'Bad or missing input!'

	return, res
end