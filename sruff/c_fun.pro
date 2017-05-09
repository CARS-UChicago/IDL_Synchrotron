Function C_fun, y

;+
; NAME:
;		C_FUN
; VERSION:
;		8.33
; PURPOSE:
;		Calculates the C(y) function form X-ray Data Booklet, Sec. 2.1, Eq. 7.
;		However, see comment below in PROCEDURE.
; CATEGORY:
;		Mathematical, SR specific
; CALLING SEQUENCE:
;		Result = C_FUN (Y)
; INPUTS:
; 	Y
; 		Numeric, otherwise arbitrary.
; OPTIONAL INPUT PARAMETERS:
;		None.
; KEYWORD PARAMETERS:
; 		None.
; OUTPUTS:
;		Returns calculation result for the function C(y), type FLOAT or higher
;		(if x is higher).
; OPTIONAL OUTPUT PARAMETERS:
;		None
; COMMON BLOCKS:
;		None.
; SIDE EFFECTS:
;		None.
; RESTRICTIONS:
;		None.
; PROCEDURE:
;		Straighforward, follows the formula (mostly, see comment below).  Uses
;		BESELK_MM, CAST and TYPE, from MIDL.
;		
;		Comment:	The value of C(y) is defined up to a multiplicative constant
;					depending on convention used.  In the Orange Booklet, the
;					covention is set by defining:
;					
;					Sigma = sqrt(2/pi)*C(y)/gamma
;					
;					With this convention, C(1) = 0.812 (note: y = 1 means 
;					E = Ec) and, multiplying by sqrt(2/pi) yields
;					
;					Sigma = 0.648/gamma
;					
;					The convention used here is one where the sqrt(2/pi)
;					coefficient is folded into C(y), so that
;					
;					Sigma = C(y)/gamma
;					
;					With this convention, C(1) = 0.648, so that again, for y=1
;					
;					Sigma = 0.648/gamma
;					
;					As before.  A side effect of the current convention is that
;					the RHS of Eq. 7 in the booklet changes, we've
;					
;					0.408*C(y)[mr]/E[Gev] --> 0.511*C(y)[mr]/E[GeV]
; MODIFICATION HISTORY:
;		Created 17-OCT-2014 by Mati Meron.
;-

	on_error, 1

	typ = Type(y)
	res = sqrt(2d*!dpi/3d)*Beselk_mm(y,5d/3d,/int)/(y*Beselk_mm(y/2d,2d/3d)^2)

	return, Cast(res,4,typ,/fix)
end