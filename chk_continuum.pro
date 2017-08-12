pro chk_continuum,obj_id,jmag=jmag;


	restore,obj_id+'_1d.idl'
	s1d = s1d*(1d-18)	;in erg s-1 cm-2 A-1
	fnu = s1d * lambda^2.d0/3.0d18
	nu = double(3.d18/lambda)

	if ~keyword_set(jmag) then jmag = 16.873200 ;star_14984
;	if ~keyword_set(jmag) then jmag = 23.4899998 ;20497 (i mag)

;	filter = '~/SED/filters/HST/WFC3/WFC3_F125W.txt'
	filter = '~/SED/filters/HST/ACS/ACS_F850lp.txt'
;	filter = '~/SED/filters/HST/ACS/ACS_F775w.txt'
	readcol,filter,lambda_filter,trans_filter,format='(d,d)',/silent
	interp_trans = double(interpol(trans_filter,lambda_filter,lambda))
	avg_flux = tsum(nu,interp_trans*fnu/nu)/tsum(nu,interp_trans/nu)
	print,avg_flux,fnu(0)
	mag = -48.6 - 2.5*alog10(avg_flux)	
	print,-48.6 - 2.5*alog10(fnu(0))	

;	print,'> Jmag     =',jmag
	print,'> Measured =',mag

end
