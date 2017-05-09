Pro SR_constants, show = sho

;+
; NAME:
;		SR_CONSTANTS
; VERSION:
;		5.0
; PURPOSE:
;		Creates or updates a system variable named !SRCON.  !SRCON is a
;		structure the fields of which contain values of physical constants used
;		in synchrotron calculations, as follows:
;		!SRCON
;			RE	  -	The classical electron radius r_c, in m.
;			TE	  -	The time corresponding to the classical electron radius,
;					i.e. r_c/c, in s.
;			EE	  -	Self energy of electron, m_e*c^2, in J.
;			BE    -	The field corresponding to classical electron radius, in T.
;			ALP	  -	The fine structure constant.
;			BERE  -	The product BE*RE, rigidity (in Tesla*meter).
;			HC	  - The product H*C (in Joule*meter).
;			SCAL  -	The scaling factor for converting Joule to keV or ma to
;					electron number.
;			CONV  -	The conversion factor from inverse wavelength (in Angstrem)
;					to energy (in keV).
;
;		All the values are given in a double precision format.  Of course, the
;		actual accuracy depends on how precisely they can be measured.
; CATEGORY:
;		Utility.
; CALLING SEQUENCE:
;		SR_CONSTANTS
; INPUTS:
;		None.
; OPTIONAL INPUT PARAMETERS:
;		None.
; KEYWORD PARAMETERS
;	/SHOW
;		Switch.  If set, basic information about the structure !SRCON is shown.
; OUTPUTS:
;		None.
; OPTIONAL OUTPUT PARAMETERS:
;		None.
; COMMON BLOCKS:
;		None.
; SIDE EFFECTS:
;		None.
; RESTRICTIONS:
;		None.
; PROCEDURE:
;		Either creates or, if already exists, updates the system variable
;		!SRCON.  Calls CONSTANTS and STREQ from MIDL.
; MODIFICATION HISTORY:
;		Created 5-JAN-2005 by Mati Meron.
;-

	defsysv, '!srcon', exists = exs
	if not exs then begin
		mksfl = Streq(!pcon.usys,'MKS')
		if mksfl then Constants, /cgs

		e = !pcon.e
		m = !pcon.me
		c = !pcon.c
		h = !pcon.hbar*2*!dpi
		emks = 10*e/c

		tem = {re: 1d-2*e^2/(m*c^2), te: e^2/(m*c^3), ee: 1d-7*m*c^2, $
				be: 1d-4*m^2*c^4/e^3, bere: 1d-6*m*c^2/e, hc: 1d-9*h*c, $
				alp: !pcon.alpha, scal: 1d-4*c/e, conv: 1d-3*h*c^2/e}

		if mksfl then Constants, /mks
		defsysv, '!srcon', tem
	endif
	if keyword_set(sho) then help, /st, !srcon

	return
end