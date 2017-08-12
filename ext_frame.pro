;	extract individual science frames
pro ext_frame,obj_id

	chk_gz = file_test('spSlit*fits.gz')
	GA = [1.224d0,1.181d0,1.266d0,1.209d0,1.181d0,1.171d0,1.207d0,1.218d0]
    RO = [2.566d0,2.476d0,2.654d0,2.534d0,2.477d0,2.455d0,2.531d0,2.554d0]

	sltmak = strsplit(file_search('*.plan'),'.',/extract)
	readcol,'slit_objs_'+sltmak(0)+'.txt',slitnum,objectID,chipnum,format='(I,A,I)'
	num = slitnum(where(objectID eq obj_id))
    chipno = chipnum(where(objectID eq obj_id))
	print,'>> Extract individual science frames of object ',obj_id,' ',num,' ',chipno
   	B_add = chipno-1
  	R_add = chipno-1
	B_Gain = reform(GA(B_add))
	B_Readout = reform(RO(B_add))
	R_Gain = reform(GA(R_add))
	R_Readout = reform(RO(R_add))

	Bfilename='spSlit.'+sltmak(0)+'.'+string(num,format='(I03)')+'B.fits'
	Rfilename='spSlit.'+sltmak(0)+'.'+string(num,format='(I03)')+'R.fits'
	if chk_gz then Bfilename = Bfilename+'.gz'
	if chk_gz then Rfilename = Rfilename+'.gz'

	fits_info,Bfilename,N_ext = n_ext,/silent
	
	spawn,'mkdir '+obj_id
	path  = './'+obj_id+'/'
	spawn,'rm -rf '+path+'*'
	print,'>>>>>>'+path+'<<<<<<'
	k=1
	for i=1,n_ext,2 do begin
		SB = mrdfits(Bfilename,i,/silent)
		SR = mrdfits(Rfilename,i,/silent)
		Bflux = SB.flux
		Rflux = SR.flux
		Bout = path+'BS_'+string(k,format='(I02)')+'.fits'
		Rout = path+'RS_'+string(k,format='(I02)')+'.fits'
		mwrfits,Bflux,Bout,/silent
		mwrfits,Rflux,Rout,/silent
		K = K+1
	endfor

end
