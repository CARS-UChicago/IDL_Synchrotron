Function ID_conv_fun, x, par

;+
; NAME:
;		ID_CONV_FUN
; VERSION:
;		5.1
; PURPOSE:
;		An internal function to be used by ROOT for ID_CONV purposes.
; CATEGORY:
;		Mathematical, SR specific
; CALLING SEQUENCE:
;		Result = ID_CONV_FUN( X, PAR)
; INPUTS:
;	X
;		Scalar or vector.  Represents magnetic field value.
;	PAR
;		A 4 element vector, including:
;			0	:	Calculation option number.
;			1-2	:	Calculation parameters.
;			3	:	Flag to be transferred to ID_FIELD.
; OPTIONAL INPUT PARAMETERS:
;		None.
; KEYWORD PARAMETERS:
;		None.
; OUTPUTS:
;		Returns calculation result, details depend on case.
; OPTIONAL OUTPUT PARAMETERS:
;		None
; COMMON BLOCKS:
;		None.
; SIDE EFFECTS:
;		None.
; RESTRICTIONS:
;		None.
; PROCEDURE:
;		Straighforward.  Calls ID_FIELD.

; MODIFICATION HISTORY:
;		Created 20-SEP-2005 by Mati Meron.
;-

	on_error, 1

	tem = ID_field(x,/inv,/qui,enh=par[3])
	case par[0] of
		10	:	res = tem - par[1]*x
		34	:	res = tem^3 - par[1]*(tem^2 + (par[2]*x)^2/2)
		else:	message, 'Option code ' + string(fix(par[0],form='(i3)')) + $
				'not acceptable!'
	endcase

	return, res
end