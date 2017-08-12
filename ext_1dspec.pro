;	construct 1D spectra
; 
;	1) Take science position from slitmask
;   2) Extract 1D spectrum within 9 pixels (0".1185/pixel * 9pixel = 1.0665 ~ Seeing)

pro ext_1dspec,obj_id,star=star


	chk_gz = file_test('spSlit*fits.gz')
	sltmak = strsplit(file_search('*.plan'),'.',/extract)
	readcol,'slit_objs_'+sltmak(0)+'.txt',slitnum,objectID,chipnum,dum,n_exp,$
		RA_objs,dec_objs,RA_slits,dec_slits,$
		format='(I,A,I,A,I,d,d,d,d)'
	tmp = where(objectID eq obj_id)
	num = slitnum(tmp)
	chipno = chipnum(tmp)
	nexp = n_exp(tmp)
	B_add = chipno(0)-1
	R_add = chipno(0)-1

	mskfile = sltmak(0)+'_sltmsk.fits'
	readcol,'obj_pos_in_slit_'+sltmak(0)+'.txt',dum,obj_list,flag_list,format='(i,a,i)',/silent
	ra_obj = RA_objs(tmp)
	dec_obj = dec_objs(tmp)
	ra_slit = RA_slits(tmp)
	dec_slit = dec_slits(tmp)
	gcirc,2,RA_obj,dec_obj,RA_slit,dec_slit,dis ;dis in arcsec
	yoffset=dis/0.1185 ; yoffset in pixels (y-axis)
	flag = flag_list(where(obj_list eq obj_id))
	if flag lt 0 then yoffset = yoffset*(-1)
	

;	[0] READ wavelength 
	print,'>>> BUILD 1D SPECTRA of '+obj_id,num

	s2d_B = mrdfits(obj_id+'/B_2dspec.fits',0,/silent)	
	n2d_B = mrdfits(obj_id+'/BN_2dspec.fits',0,/silent)	
	sky_B = mrdfits(obj_id+'/BS_01-out.fits',0,/silent)	
	if chk_gz then begin
		calib_B = mrdfits('calibSlit.'+sltmak(0)+'.'+string(num,format='(I03)')+'B.fits.gz',1,/silent)
	endif else begin
		calib_B = mrdfits('calibSlit.'+sltmak(0)+'.'+string(num,format='(I03)')+'B.fits',1,/silent)
	endelse
	lambda2d_B = lambda_eval(calib_B)

	s2d_R = mrdfits(obj_id+'/R_2dspec.fits',0,/silent)	
	n2d_R = mrdfits(obj_id+'/RN_2dspec.fits',0,/silent)	
	sky_R = mrdfits(obj_id+'/RS_01-out.fits',0,/silent)	
	if chk_gz then begin
		calib_R = mrdfits('calibSlit.'+sltmak(0)+'.'+string(num,format='(I03)')+'R.fits.gz',1,/silent)
	endif else begin
		calib_R= mrdfits('calibSlit.'+sltmak(0)+'.'+string(num,format='(I03)')+'R.fits',1,/silent)
	endelse
	lambda2d_R = lambda_eval(calib_R)

	restore,'response.idl'

	ny = n_elements(s2d_B(0,*))
	ycen = round((ny-1)/2.+yoffset)
;	ycen = round((n_elements(s2d_B(0,*))-1)/2.+yoffset)
	print,'>> ycen =',ycen
	print,'RA_obj,dec_obj,RA_slit,dec_slit'
	print,RA_obj,dec_obj,RA_slit,dec_slit,dis,yoffset

	nxB = n_elements(s2d_B(*,0))
	nxR = n_elements(s2d_R(*,0))
	nx = nxB+nxR
	npx = 9 ;

	exc1 = indgen(ycen-fix(npx/2.))
	exc2 = indgen(n_elements(s2d_B(0,*))-1-ycen+fix(npx/2.))+ycen+fix(npx/2.)+1
	exc_B = [exc1,exc2]

	exc2 = indgen(n_elements(s2d_R(0,*))-1-ycen+fix(npx/2.))+ycen+fix(npx/2.)+1
	exc_R = [exc1,exc2]

	std_n1d_exc_B = dblarr(nxB)
	std_n1d_exc_R = dblarr(nxR)
	for nn = 0,nxB-1 do std_n1d_exc_B(nn) = stddev(s2d_B(nn,exc_B))
	for nn = 0,nxR-1 do std_n1d_exc_R(nn) = stddev(s2d_R(nn,exc_R))
	std_n1d_exc = [std_n1d_exc_B,std_n1d_exc_R]

	s2d_B = reform(s2d_B(*,ycen-fix(npx/2.):ycen+fix(npx/2.)))
	lambda2d_B = reform(lambda2d_B(*,ycen-fix(npx/2.):ycen+fix(npx/2.)))
	n2d_B = reform(n2d_B(*,ycen-fix(npx/2.):ycen+fix(npx/2.)))
	s2d_R = reform(s2d_R(*,ycen-fix(npx/2.):ycen+fix(npx/2.)))
	lambda2d_R = reform(lambda2d_R(*,ycen-fix(npx/2.):ycen+fix(npx/2.)))
	n2d_R = reform(n2d_R(*,ycen-fix(npx/2.):ycen+fix(npx/2.)))

	sky_2d_B = reform(sky_B(*,ycen-fix(npx/2.):ycen+fix(npx/2.)))
	sky_2d_R = reform(sky_R(*,ycen-fix(npx/2.):ycen+fix(npx/2.)))

	s2d_ex = [s2d_B,s2d_R]
	n2d_ex = [n2d_B,n2d_R]
	sky_2d = [sky_2d_B,sky_2d_R]
	lambda2d_ex = [lambda2d_B,lambda2d_R]

	ycen_ex = fix(npx/2.)
;	lambda0 = reform(lambda2d_B(*,ycen))  ; lambda at the center of a science object
	lambda0_ex = reform(lambda2d_ex(*,ycen_ex)) ; lambda at the center of a science object
	nrows = n_elements(lambda2d_ex(0,*))
	npix = n_elements(lambda0_ex)
	dldx = (lambda0_ex[npix-1]-lambda0_ex[0])/npix
	tiltx = (lambda2d_ex[npix/2.,nrows-1]-lambda2d_ex[npix/2.,0])/nrows
	dxdp = tiltx/dldx
	
	s1d = dblarr(nx)
	s1d_num = dblarr(nx)
	s1d_denom = dblarr(nx)
	n1d_num = dblarr(nx)
	n1d_denom = dblarr(nx)

	n1d = dblarr(nx)
	n1d_std = dblarr(nx)
	sky_1d = dblarr(nx)
	for i=0l,nx-1 do sky_1d(i) = total(sky_2d(i,*))
	sky_1d = sky_1d - median(sky_1d)

	sigma = 4.31392d0 ; sigma in pixels for 2014Feb20 ........check!!!
;	sigma = 2.5d0 ; sigma in pixels for 2014Feb20 ........check!!!

	yy = indgen(npx) ; for 9pixel-extraction for 1D spectrum
	gprofile = (1.d0/(sigma*sqrt(2.d0*!dpi)))*exp((-1.d0)*((yy-fix(npx/2.))^2.d0)/(2.d0*sigma^2.d0))
	gprofile = gprofile/total(gprofile)
	gprofile = gprofile ## (fltarr(npix)+1.)
	weight = 1.d0/(n2d_ex^2.d0)
	ycen_ex = fix(npx/2.)
	for yi = 0l,nrows-1 do begin
		cshift = round((yi-ycen)*dxdp)
		s1d_num = s1d_num + shift(s2d_ex(*,yi)*gprofile(*,yi)*weight(*,yi),cshift)
		s1d_denom = s1d_denom + shift(gprofile(*,yi)^2.d0*weight(*,yi),cshift)
		n1d_num = n1d_num + shift(gprofile(*,yi),cshift)
		n1d_denom = n1d_denom + shift(gprofile(*,yi)^2.d0*weight(*,yi),cshift)
	endfor
	s1d = s1d_num/s1d_denom
	n1d = sqrt(n1d_num/n1d_denom)

	calib = interpol(calib_factor,calib_lam,lambda0_ex)
	s1d_calib = s1d * calib
	n1d_calib = n1d * calib
	n1d_exc_calib = std_n1d_exc * calib
	n1d_std_calib = n1d_std * calib

	s1d = s1d_calib * 1d18	;1e-18 erg ...
	n1d = n1d_calib * 1d18
	n1d_exc = n1d_exc_calib * 1d18
	lambda = lambda0_ex

	fin = where(finite(n1d))  & nan = where(~finite(n1d))
	n1d(nan) = interpol(n1d(fin),lambda(fin),lambda(nan))
	fin = where(finite(s1d))  & nan = where(~finite(s1d))
	s1d(nan) = interpol(s1d(fin),lambda(fin),lambda(nan))

	if ~keyword_set(star) then begin
;   calculate continuum flux
	restore,'~/SED/resolved/'+obj_id+'/mcmc_int_ref3.idl'
	pz = u_int.redshift
    pmass = median(u_int.mass)
    pmet = median(u_int.met)
    psfh = median(u_int.sfh)+1
	pmet = 10.d0^median(u_int.met)
	pdust = median(u_int.ebv)
    page = median(u_int.age)
    pfesc = 0.
    bc03 = makemodel(z=pz,mass=pmass,met=pmet,sfh=psfh,dust=pdust,page=page,fesc=pfesc,H0=67.8d0,omega_m=0.308d0,lambda0=0.692d0,/noplot,/silent)
    fnu = bc03.flux;*pmass 
    flam_rs = fnu*3.d18/(bc03.lambda^2.d0) * 1.d18 ; 10^-18 erg s-1 cm-2 A-1
    flam_rest = flam_rs*(1.d0+pz)
    lam_rest = bc03.lambda/(1.d0+pz) 
	continuum = mean(flam_rest(where(lam_rest gt 1216.d0 and lam_rest lt (1216.d0+100))))
    print,'>>> continuum =',continuum
	obj_continuum = continuum
	spawn,'rm -rf '+obj_id+'_1d.idl'
	save,filename=obj_id+'_1d.idl',s1d,n1d,lambda,sky_1d,obj_continuum
	endif
	spawn,'rm -rf '+obj_id+'_1d.idl'
	save,filename=obj_id+'_1d.idl',s1d,n1d,lambda,sky_1d
	


end
