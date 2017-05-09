Pro Multi_info

;+
; NAME:
;		MULTI_INFO
; VERSION:
;		5.5
; PURPOSE:
;		Displays info about currently defined multilayer.
; CATEGORY:
;		X-ray calculations.
; CALLING SEQUENCE:
;		MULTI_INFO
; INPUTS:
;		None
; OPTIONAL INPUT PARAMETERS:
;		None.
; KEYWORD PARAMETERS:
;		None.
; OUTPUTS:
;		None other then screen display.
; OPTIONAL OUTPUT PARAMETERS:
;		None.
; COMMON BLOCKS:
;		MULTILAYER_PARS.  See MULTI_GEN for more information.
; SIDE EFFECTS:
;		None.
; RESTRICTIONS:
;		None
; PROCEDURE:
;		Straightforward.  Calls DEFAULT and TABULATE from MIDL.
; MODIFICATION HISTORY:
;		Created 5-JUN-2006 by Mati Meron.
;-

	common multilayer_pars, m_exs, m_def, m_form, $
	m_elar, m_numels, m_rarr, m_warr, m_tarr, m_roff

	on_error, 1

	ftit = ['Top Film', 'First Layer', 'Second Layer', 'Substrate']
	head = ['Elements', 'Weights']

	print
	if Default(m_exs,0) then begin
		len = n_elements(m_numels)
		smult = string((len - m_def[0] - m_def[3])/2,form='(i6)')
		glob = strcompress(smult,/rem) + ' repetitions'
		print, '			' + glob
		print, '			' + strjoin(replicate('_',strlen(glob)))
		print
		print
		if m_form then head[1] = head[1] + ' (formula)'
		j = 0
		fir = 0
		for i = 0, 3 do begin
			if m_def[i] then begin
				if i eq 3 then j = len - 1
				stit = 'Thickness = ' + string(m_tarr[j],form='(f8.1)') + $
				', Density = ' + string(m_rarr[j],form = '(g9.4)') + $
				', Roughness = ' + string(m_roff[j], form = '(f6.2)')
				tit = ['			' + ftit[i],strcompress(stit)]
				if i eq 3 then begin
					las = round(total(m_numels)) - 1
					fir = las - m_numels[j] + 1
				endif else las = fir + m_numels[j]-1
				Tabulate, m_elar[fir:las], m_warr[fir:las], $
				tit = tit, head = head
				print
				print
				j = j + 1
				fir = las + 1
			endif
		endfor
	endif else print, 'Multilayer not defined!'

	return
end
