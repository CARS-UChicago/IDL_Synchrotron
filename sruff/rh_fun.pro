Function RH_fun, gpsi, pars

;+
; NAME:
;		RH_FUN
; VERSION:
;		5.0
; PURPOSE:
;		An interface to the function H_FUN, whid calculates bending magnet power
;		distribution (see details there).  For internal use only.
; CATEGORY:
;		Mathematical, SR specific
; CALLING SEQUENCE:
;		Result = RH_FUN (GPSI, PARS)
; INPUTS:
;	GPSI
;		Scalar or vector.  /Gamma*/Psi
;	PARS
;		Scalar or two-element vector.  The first element (mandatory) corresponds
;		to the E_EC input of H_FUN.  The second, optional element corresponds
;		to the POW input of H_FUN.
; OPTIONAL INPUT PARAMETERS:
;		None.
; KEYWORD PARAMETERS:
;		None.
; OUTPUTS:
;		Returns calculation result, type FLOAT or the highest of the types of
;		[GPSI, PAR].
; OPTIONAL OUTPUT PARAMETERS:
;		None.
; COMMON BLOCKS:
;		None.
; SIDE EFFECTS:
;		None.
; RESTRICTIONS:
;		None.
; PROCEDURE:
;		Calls H_FUN (see there).  RH_FUN serves simply as an interface to H_FUN
;		for numerical integration purposes.
; MODIFICATION HISTORY:
;		Created 10-JAN-2005 by Mati Meron.
;-

	on_error, 1

	case n_elements(pars) of
		1	:	res = H_fun(gpsi,pars)
		2	:	res = H_fun(gpsi,pars[0],pars[1])
		else:	message, 'Wrong number of variables!'
	endcase

	return, res
end