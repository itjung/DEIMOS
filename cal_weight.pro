pro cal_weight

	sltmak = strsplit(file_search('*.plan'),'.',/extract)
	readcol,'slit_objs_'+sltmak(0)+'.txt',slitnum,objectID,chipnum,dum,n_exp,format='(I,A,I,A,I)'
	nslit = n_elements(slitnum)

	obj_list = ['']
	for i=0,nslit-1 do begin
		wrd = strsplit(objectid(i),'_',/extract)
		if wrd(0) eq 'star' then obj_list = [[obj_list],[objectID(i)]]
	endfor
	obj_list = obj_list(1:-1)
	obj_list = obj_list(where(obj_list eq 'star_14984' or obj_list eq 'star_15723'))
;	chk = where(obj_list eq 'star_8773')
;	if chk(0) ne -1 then obj_list = obj_list(where(obj_list ne 'star_8773'))
;	chk = where(obj_list eq 'star_15723')
;	if chk(0) ne -1 then obj_list = obj_list(where(obj_list ne 'star_15723'))
;	chk = where(obj_list eq 'star_44065')
;	if chk(0) ne -1 then obj_list = obj_list(where(obj_list ne 'star_44065'))
	n_stars = n_elements(obj_list)
	nexp = n_exp(0)
	print,obj_list

;	xcut = indgen(13)*300+150
	xcut = [254,530,1900,2150,2800,3150,3310,3665,3975]
	ncut = n_elements(xcut)
	maxf_arr = dblarr(n_stars,nexp,ncut)
	sigma_arr = dblarr(n_stars,nexp,ncut)
	ycen_arr = dblarr(n_stars,nexp,ncut)
	chisq_arr = dblarr(n_stars,nexp,ncut)

	tmp_ycen = dblarr(nexp)
	for ni=0,nexp-1 do begin
		s2d = mrdfits(obj_list(0)+'/Bf'+string(ni+1,format='(I02)')+'_sky.fits',0,/silent)
		tmp_add = dblarr(ncut)
		for ci=0,ncut-1 do begin
			maxf = max(s2d(xcut(ci),*),add)
			tmp_add(ci) = add
		endfor
		tmp_ycen(ni) = round(mean(tmp_add))
	endfor

	print,tmp_ycen
	for oi = 0,n_stars-1 do begin
		obj = obj_list(oi)
		print,obj
		for ni=0,nexp-1 do begin
			s2d = mrdfits(obj+'/Bf'+string(ni+1,format='(I02)')+'_sky.fits',0,/silent)	
			s2d_for_seeing = mrdfits(obj+'/BS_'+string(ni+1,format='(I02)')+'-out.fits',0,/silent)	
;			ny = n_elements(s2d(0,*))
			for ci=0,ncut-1 do begin
				yarr = reform(s2d(xcut(ci),tmp_ycen(ni)-5:tmp_ycen(ni)+5))	
				ny = n_elements(yarr)
				x = indgen(ny)
;				xx = dindgen(ny*10)/10
;				yyarr = frebin(yarr,ny*10)
				yfit = gaussfit(x,yarr,coeff,chisq=chisq,nterms=4)
				yy = coeff(0)*exp(-((x-coeff(1))/coeff(2))^2.d0/2.)+coeff(3);+coeff(4)*xx+coeff(5)*xx^2.
				maxf_arr(oi,ni,ci) = max(yy,add)

				yarr = reform(s2d_for_seeing(xcut(ci),*))	
				ny = n_elements(yarr)
				x = indgen(ny)
;				xx = dindgen(ny*10)/10
;				yyarr = frebin(yarr,ny*10)
				yfit = gaussfit(x,yarr,coeff,chisq=chisq,nterms=4)
				yy = coeff(0)*exp(-((x-coeff(1))/coeff(2))^2.d0/2.)+coeff(3);+coeff(4)*xx+coeff(5)*xx^2.
				tmp_maxf_arr = max(yy,add)

				sigma_arr(oi,ni,ci) = coeff(2)
				ycen_arr(oi,ni,ci) = round(x(add));+tmp_ycen(ni)-4
				chisq_arr(oi,ni,ci) = chisq
			endfor
		endfor
	endfor

	warr = dblarr(n_stars,nexp,ncut)
	for oi=0,n_stars-1 do for i=0,nexp-1 do for j=0,ncut-1 do warr(oi,i,j) = maxf_arr(oi,i,j)/maxf_arr(oi,0,j)
	avg_flux = dblarr(n_stars)
	for oi=0,n_stars-1 do avg_flux(oi) = mean(maxf_arr(oi,*,*))
;	bright = max(avg_flux,brightest)
	brightest = 0
	w = dblarr(nexp)
	sigma = dblarr(nexp)
	fwhm = dblarr(nexp)
	ycen = dblarr(nexp)
	for i=0,nexp-1 do begin
		meanclip,warr(brightest,i,*),sigclip_mean,clipsig=3.
		w(i) = sigclip_mean
		meanclip,sigma_arr(brightest,i,*),sigclip_sigma,clipsig=3.
		sigma(i) = sigclip_sigma
		fwhm(i) = sigclip_sigma * 2.35482d0
		ycen(i) = round(mean(ycen_arr(*,i,*)))
	endfor
	print,'>>> '+obj_list(brightest)
	print,'>> w'
	print,w
	print,'>> Seeing (FWHM) in pixels'
	print,fwhm
	print,'>> Seeing (FWHM) in arcsec'
	print,fwhm*0.1185
	print,'>> ycen'
	print,ycen

	save,filename='weight.idl',obj_list,maxf_arr,warr,fwhm,sigma,sigma_arr,ycen_arr,brightest,w,ycen
end
