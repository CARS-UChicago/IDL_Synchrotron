Pro EB_conv, energy = enr, field = mab, radius = rad, critical_energy = ecr, $
    gamma = gam

;+
; NAME:
;	EB_CONV
; PURPOSE:
;	Conversions between synchrotron radiation parameters.
; CATEGORY:
;	X-ray calculations.
; CALLING SEQUENCE:
;	EB_CONV, keywords
; INPUTS:
;	None.
; OPTIONAL INPUT PARAMETERS:
;	None.
; KEYWORD PARAMETERS:
;    ENERGY
;	Value of electron/positron beam energy.		|
;    FIELD						|    Two and ONLY two
;	Value of magnetic field.			|    out of the five
;    RADIUS						|    need to be 
;	Value of magnetic bend radius.			|    defined.
;    CRITICAL_ENERGY					|
;	Value of (synchrotron radiation) critical 	|    The remaining 
;	energy.						|    three (whichever
;    GAMMA						|    they are) serve as
;	Value of the relativistic gamma of the 		|    output parameters.
;	electron/positron beam.				|
; OUTPUTS:
;	None
; OPTIONAL OUTPUT PARAMETERS:
;	The three nondefined keywords serve as optional output parameters,
; COMMON BLOCKS:
;	None.
; SIDE EFFECTS:
;	The nondefined variables assigned to keywords become defined.
; RESTRICTIONS:
;	As mentioned above, two and only two of the parameters should be 
;	defined.  Morever, the combination of ENERGY and GAMMA is not allowed.
;	Also, all input parameters should be nonzero.  Zero inputs are 
;	considered non-defined.
; PROCEDURE:
;	Straightforward.  Uses DEFAULT and HOW_MANY from MIDL
; MODIFICATION HISTORY:
;	Created 15-AUG-1994 by Mati Meron.
;-

    on_error, 1
    cf = [333.56410, 0.66502546, 1956.9507]
    tot = How_many(fir=enr,sec=mab,thir=rad,fort=ecr,fif=gam,/nozero,which=cod)
    if tot lt 2 then message, 'Underspecified!' else $
    if tot gt 2 then message, 'Overspecified!'

    case cod of
	3    :
	5    :	mab = cf(0)*enr/rad
	6    :	enr = mab*rad/cf(0)
	9    :	mab = ecr/(cf(1)*enr^2)
	10   :	enr = sqrt(ecr/(cf(1)*mab))
	12   :	begin
		    enr = (ecr*rad/(cf(0)*cf(1)))^(1./3.)
		    mab = cf(0)*enr/rad
		end
	18   :	enr = gam/cf(2)
	20   :	begin
		    enr = gam/cf(2)
		    mab = cf(0)*enr/rad
		end
	24   :	begin
		    enr = gam/cf(2)
		    mab = ecr/(cf(1)*enr^2)
		end
	17   :	message, "Can't specify both energy and gamma!"
    endcase

    rad = Default(rad,cf(0)*enr/mab,/dtype)
    ecr = Default(ecr,cf(1)*enr^2*mab,/dtype)
    gam = Default(gam,cf(2)*enr,/dtype)

    return
end
