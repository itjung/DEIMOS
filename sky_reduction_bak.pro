;	A' = A2 - ((A1+A3)/2)
;	consider weights & calculate noise (Kriek et al. 2015)
function noise1,SF1,SF2,G,R
	noise = sqrt(G*SF1+G*SF2+2.d0*R^2.d0)/G
	return,noise
end

function noise2,SF0,SF1,SF2,G,R
	noise = sqrt(G*SF1+G*(SF0+SF2)/4.d0+(3.d0*R^2.d0/2))/G
	return,noise
end

pro sky_reduction_bak,obj_id

	print,'>> reduce sky background ',obj_id
	
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

	restore,'weight_tmp.idl'
	dither = sort(ycen(0:2))   
	dither = [dither,dither,dither,dither,dither,dither,dither,dither,dither]

	bin = mrdfits(sltmak(0)+'.bintabs.fits',1)
	
	GA = [1.224d0,1.181d0,1.266d0,1.209d0,1.181d0,1.171d0,1.207d0,1.218d0]
	RO = [2.566d0,2.476d0,2.654d0,2.534d0,2.477d0,2.455d0,2.531d0,2.554d0]
	B_Gain = GA(B_add)
	B_Readout = RO(B_add)
	R_Gain = GA(R_add)
	R_Readout = RO(R_add)
;   from https://www2.keck.hawaii.edu/inst/deimos/deimos_detector_data.html	

	path  = './'+obj_id+'/'
	print,'>>>>>>'+path+'<<<<<<'
	print,'>>> Read sky-subtracted science images'

;	loop for total exposures
	for iexp=1,nexp(0) do begin
;		read CR-cleaned images
		Bim = path + 'BS_'+string(iexp,format='(I02)')+'-out.fits'
		Rim = path + 'RS_'+string(iexp,format='(I02)')+'-out.fits'
		pre_Bim = path + 'BS_'+string(iexp-1,format='(I02)')+'-out.fits'
		next_Bim = path + 'BS_'+string(iexp+1,format='(I02)')+'-out.fits'
		pre_Rim = path + 'RS_'+string(iexp-1,format='(I02)')+'-out.fits'
		next_Rim = path + 'RS_'+string(iexp+1,format='(I02)')+'-out.fits'

;		calculate noise from non-sky-subtracted science frame: eq(2) in Kriek et al. (2015)
;		sky_subtraction
		BF = mrdfits(Bim,0,/silent)
		RF = mrdfits(Rim,0,/silent)
		if file_test(pre_Bim) then begin
			pre_BF = mrdfits(pre_Bim,0,/silent)
			pre_RF = mrdfits(pre_Rim,0,/silent)
			if file_test(next_Bim) then begin	; A1' = A2 - ((A1+A3)/2)
				next_BF = mrdfits(next_Bim,0,/silent)
				next_RF = mrdfits(next_Rim,0,/silent)
				BN = noise2(pre_BF,BF,next_BF,B_Gain,B_Readout)
				RN = noise2(pre_RF,RF,next_RF,R_Gain,R_Readout)
				BF_sky = BF - ((pre_BF+next_BF)/2.)
				RF_sky = RF - ((pre_RF+next_RF)/2.)
			endif else begin
				BN = noise1(BF,pre_BF,B_Gain,B_Readout)
				RN = noise1(RF,pre_RF,R_Gain,R_Readout)
				BF_sky = BF - pre_BF
				RF_sky = RF - pre_RF
			endelse
		endif else begin
			next_BF = mrdfits(next_Bim,0,/silent)
			next_RF = mrdfits(next_Rim,0,/silent)
			BN = noise1(BF,next_BF,B_Gain,B_Readout)
			RN = noise1(RF,next_RF,R_Gain,R_Readout)
			BF_sky = BF - next_BF
			RF_sky = RF - next_RF
		endelse
		spawn,'rm -rf '+path+'Bf'+string(iexp,format='(I02)')+'_sky.fits'
		spawn,'rm -rf '+path+'Rf'+string(iexp,format='(I02)')+'_sky.fits'
		spawn,'rm -rf '+path+'BN'+string(iexp,format='(I02)')+'.fits'
		spawn,'rm -rf '+path+'RN'+string(iexp,format='(I02)')+'.fits'
		mwrfits,BF_sky,path+'Bf'+string(iexp,format='(I02)')+'_sky.fits',/silent
		mwrfits,RF_sky,path+'Rf'+string(iexp,format='(I02)')+'_sky.fits',/silent
		mwrfits,BN,path+'BN'+string(iexp,format='(I02)')+'.fits',/silent
		mwrfits,RN,path+'RN'+string(iexp,format='(I02)')+'.fits',/silent

;		residual sky-subtraction		
		nx = n_elements(BF_sky(*,0))
		nyB = n_elements(BF_sky(0,*))
		nyR = n_elements(RF_sky(0,*))
	

		ycen = round((nyB-1)/2.+yoffset(0))
		case dither(iexp-1) of 
		0: ycen = round(ycen-(1.d0/0.1185))
		1: ycen = ycen
		2: ycen = round(ycen+(1.d0/0.1185))
		endcase
		obj_ypx = indgen(13)-6 + ycen
		all_ypx = indgen(nyB)
		y1 = 0 & y2 = obj_ypx(0)-1 & y3 = obj_ypx(-1)+1 & y4 = nyB-1
		if y2 lt y1 then y2 = y1
		if y3 gt y4 then y3 = y4
		sky_ypx_B = [all_ypx(y1:y2),all_ypx(y3:y4)]
		all_ypx = indgen(nyR)
		y4 = nyR-1
		if y3 gt y4 then y3 = y4
		sky_ypx_R = [all_ypx(y1:y2),all_ypx(y3:y4)]

		for ix=0L,nx-1 do begin
			meanclip,BF_sky(ix,sky_ypx_B),mean_clipped,sig,clipsig=3.,SUBS=SUBS & BF_sky(ix,*) = BF_sky(ix,*) - mean_clipped
			meanclip,RF_sky(ix,sky_ypx_R),mean_clipped,sig,clipsig=3.,SUBS=SUBS & RF_sky(ix,*) = RF_sky(ix,*) - mean_clipped 
		endfor
		spawn,'rm -rf '+path+'Bf'+string(iexp,format='(I02)')+'_sky_res.fits'
		spawn,'rm -rf '+path+'Rf'+string(iexp,format='(I02)')+'_sky_res.fits'
		mwrfits,BF_sky,path+'Bf'+string(iexp,format='(I02)')+'_sky_res.fits',/silent
		mwrfits,RF_sky,path+'Rf'+string(iexp,format='(I02)')+'_sky_res.fits',/silent

	endfor

end
