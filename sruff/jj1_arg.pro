Function JJ1_arg, x, par

;+
; NAME:
;	JJ1_ARG
; PURPOSE:
;	Used as kernel in calculation of the US1 function (see there).
; CATEGORY:
;	Mathematical, Integration Kernel.
; CALLING SEQUENCE:
;	Result = JJ1_ARG (X, PAR)
; INPUTS:
;    X
;	Numeric, otherwise arbitrary.
;    PAR
;	3 element vector containing [k,u,v].
; OPTIONAL INPUT PARAMETERS:
;	None.
; KEYWORD PARAMETERS:
;	None.
; OUTPUTS:
;	Returns 1/4*[cos (k*(v*sin(x) - u*sin(2*x) - x)) - cos(k*x)] in a 
;	numerically convenient form.
; OPTIONAL OUTPUT PARAMETERS:
;	None.
; COMMON BLOCKS:
;	None.
; SIDE EFFECTS:
;	None.
; RESTRICTIONS:
;	None.
; PROCEDURE:
;	Returns the result as a product of sines.
; MODIFICATION HISTORY:
;	Created 30-MARCH-1994 by Mati Meron.
;-

    tem = (par(2)*sin(x) - par(1)*sin(2*x))/2
    return, sin(par(0)*tem)*sin(par(0)*(x-tem))
end
