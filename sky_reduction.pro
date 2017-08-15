function noise1,SF1,SF2,G,R
	noise = sqrt(G*SF1+G*SF2+2.d0*R^2.d0)/G
	return,noise
end

function noise2,SF0,SF1,SF2,G,R
	noise = sqrt(G*SF1+G*(SF0+SF2)/4.d0+(3.d0*R^2.d0/2))/G
	return,noise
end

pro sky_reduction,obj_id

;+NAME:
;	sky_reduction
;
; PURPOSE:
;   Reduce sky emission with adjacent science frames as A2' = A2 - (A1+A3)/2
;	Create noise maps 
;
; INPUTS:
;	./OBJECTNAME/B(R)S_##-out.fits
;   MASKNAME.bintabs.fits
;	MASKNAME.plan
;	slit_objs_MASKNAME.txt
;	
; OUTPUTS:
;  	./OBJECTNAME/B(R)f##_sky.fits
;  	./OBJECTNAME/B(R)N##.fits
;
; KEYWORD PARAMETERS:
;   No keyword.
;
; EXAMPLE:
;	IDL> sky_reduction,'object ID'
;
; MEMO:
;   READOUT & GAIN parameters are from 
;		https://www2.keck.hawaii.edu/inst/deimos/deimos_detector_data.html	
;
; MODIFICATION HISTORY:
;	Written by Intae Jung @ Aug 2017
;-

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

	bin = mrdfits(sltmak(0)+'.bintabs.fits',1)
	
;	https://www2.keck.hawaii.edu/inst/deimos/deimos_detector_data.html	
	GA = [1.224d0,1.181d0,1.266d0,1.209d0,1.181d0,1.171d0,1.207d0,1.218d0]
	RO = [2.566d0,2.476d0,2.654d0,2.534d0,2.477d0,2.455d0,2.531d0,2.554d0]
	B_Gain = GA(B_add)
	B_Readout = RO(B_add)
	R_Gain = GA(R_add)
	R_Readout = RO(R_add)

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
		mwrfits,BF_sky,path+'Bf'+string(iexp,format='(I02)')+'_sky.fits',/silent,/create
		mwrfits,RF_sky,path+'Rf'+string(iexp,format='(I02)')+'_sky.fits',/silent,/create
		mwrfits,BN,path+'BN'+string(iexp,format='(I02)')+'.fits',/silent,/create
		mwrfits,RN,path+'RN'+string(iexp,format='(I02)')+'.fits',/silent,/create

	endfor

end
