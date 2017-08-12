;	consider weights (Kriek et al. 2015)
;	weight = the maximum flux of the best-fit Gaussian from guide stars
pro weighted_combine,obj_id

	restore,'weight.idl'
;	w = weights
;	ycen = dither pattern
	dither = sort(ycen(0:2))	; 0 - bottom, 1 - middle, 2 - top
	dither = [dither,dither,dither,dither,dither,dither,dither,dither,dither]

	path  = './'+obj_id+'/'
	print,'>>>>>>'+path+'<<<<<<'

	chk_gz = file_test('spSlit*fits.gz')
	sltmak = strsplit(file_search('*.plan'),'.',/extract)
	readcol,'slit_objs_'+sltmak(0)+'.txt',slitnum,objectID,chipnum,dum,n_exp,format='(I,A,I,A,I)'
	
	subbin = 50
	gap = 1.d0/0.1185 ;in pixels
	gap = gap * subbin

	BF_in = path+'Bf01_sky.fits'
	RF_in = path+'Rf01_sky.fits'
	BF = mrdfits(BF_in,0,/silent)
	RF = mrdfits(RF_in,0,/silent)
	Bny = n_elements(BF(0,*))
	Bnx = n_elements(BF(*,0))
	Rny = n_elements(RF(0,*))
	Rnx = n_elements(RF(*,0))
	B_com = dblarr(Bnx,Bny*subbin+gap*2)
	BN_com = dblarr(Bnx,Bny*subbin+gap*2)
	R_com = dblarr(Rnx,Rny*subbin+gap*2)
	RN_com = dblarr(Rnx,Rny*subbin+gap*2)

	for ni=0,n_exp(0)-1 do begin
		BF_in = path+'Bf'+string(ni+1,format='(I02)')+'_sky.fits'
		BN_in = path+'BN'+string(ni+1,format='(I02)')+'.fits'
		RF_in = path+'Rf'+string(ni+1,format='(I02)')+'_sky.fits'
		RN_in = path+'RN'+string(ni+1,format='(I02)')+'.fits'
		BF = mrdfits(BF_in,0,/silent)
		BN = mrdfits(BN_in,0,/silent)
		RF = mrdfits(RF_in,0,/silent)
		RN = mrdfits(RN_in,0,/silent)
		print,'combine '+BF_in
		
		BF = frebin(BF,Bnx,Bny*subbin,/total)
		BN = frebin(BN,Bnx,Bny*subbin,/total)
		RF = frebin(RF,Rnx,Rny*subbin,/total)
		RN = frebin(RN,Rnx,Rny*subbin,/total)

		case dither(ni) of
		0:	begin
			B_com(*,gap*2:Bny*subbin+gap*2-1) = B_com(*,gap*2:Bny*subbin+gap*2-1) + w(ni)*BF
			R_com(*,gap*2:Rny*subbin+gap*2-1) = R_com(*,gap*2:Rny*subbin+gap*2-1) + w(ni)*RF
			BN_com(*,gap*2:Bny*subbin+gap*2-1) = BN_com(*,gap*2:Bny*subbin+gap*2-1) + (w(ni)*BN)^2. 
			RN_com(*,gap*2:Rny*subbin+gap*2-1) = RN_com(*,gap*2:Rny*subbin+gap*2-1) + (w(ni)*RN)^2. 
			end;bottom
		1:	begin
			B_com(*,gap:Bny*subbin+gap-1) = B_com(*,gap:Bny*subbin+gap-1) + w(ni)*BF 
			R_com(*,gap:Rny*subbin+gap-1) = R_com(*,gap:Rny*subbin+gap-1) + w(ni)*RF 
			BN_com(*,gap:Bny*subbin+gap-1) = BN_com(*,gap:Bny*subbin+gap-1) + (w(ni)*BN)^2. 
			RN_com(*,gap:Rny*subbin+gap-1) = RN_com(*,gap:Rny*subbin+gap-1) + (w(ni)*RN)^2. 
			end;middle
		2:	begin
			B_com(*,0:Bny*subbin-1) = B_com(*,0:Bny*subbin-1) + w(ni)*BF
			R_com(*,0:Rny*subbin-1) = R_com(*,0:Rny*subbin-1) + w(ni)*RF
			BN_com(*,0:Bny*subbin-1) = BN_com(*,0:Bny*subbin-1) +  (w(ni)*BN)^2. 
			RN_com(*,0:Rny*subbin-1) = RN_com(*,0:Rny*subbin-1) +  (w(ni)*RN)^2. 
			end;top
		endcase
	endfor

	gap = round(1.d0/0.1185)
	B_com = B_com/total(w)
	R_com = R_com/total(w)
	B_com = frebin(B_com,Bnx,Bny+gap*2,/total)
	R_com = frebin(R_com,Rnx,Rny+gap*2,/total)

	B_com = reform(B_com(*,gap:gap+Bny-1))
	R_com = reform(R_com(*,gap:gap+Bny-1))

	spawn,'rm -rf '+path+'B_2dspec.fits'
	spawn,'rm -rf '+path+'R_2dspec.fits'

	mwrfits,B_com,path+'B_2dspec.fits'
	mwrfits,R_com,path+'R_2dspec.fits'

	BN_com = sqrt(BN_com)/sqrt(total(w^2.))
	RN_com = sqrt(RN_com)/sqrt(total(w^2.))
	BN_com = frebin(BN_com,Bnx,Bny+gap*2,/total)
	RN_com = frebin(RN_com,Rnx,Rny+gap*2,/total)

	BN_com = reform(BN_com(*,gap:gap+Bny-1))
	RN_com = reform(RN_com(*,gap:gap+Bny-1))

	spawn,'rm -rf '+path+'BN_2dspec.fits'
	spawn,'rm -rf '+path+'RN_2dspec.fits'

	mwrfits,BN_com,path+'BN_2dspec.fits'
	mwrfits,RN_com,path+'RN_2dspec.fits'

end
