pro cr_clean,obj_id,sigclip=sigclip

	if ~keyword_set(sigclip) then sigclip = 4.5
	chk_gz = file_test('spSlit*fits.gz')
	sltmak = strsplit(file_search('*.plan'),'.',/extract)
	readcol,'slit_objs_'+sltmak(0)+'.txt',slitnum,objectID,chipnum,format='(I,A,I)',/silent
	num = slitnum(where(objectID eq obj_id))
	Bfilename='spSlit.'+sltmak(0)+'.'+string(num,format='(I03)')+'B.fits'
	if chk_gz then Bfilename = Bfilename+'.gz'
	fits_info,Bfilename,N_ext = n_ext,/silent
	
	path  = './'+obj_id+'/'
	spawn,'rm -rf '+path+'*-out.fits'
	spawn,'rm -rf '+path+'*-mask.fits'
	k=1 
	for i=1,n_ext,2 do begin
		Braw = path+'BS_'+string(k,format='(I02)')+'.fits'
		Rraw = path+'RS_'+string(k,format='(I02)')+'.fits'
		la_cosmic,Braw,/isbig,sigclip=sigclip,gain=B_Gain,readn=B_Readout,sigfrac=0.5
		la_cosmic,Rraw,/isbig,sigclip=sigclip,gain=R_Gain,readn=R_Readout,sigfrac=0.5
		k = k+1
	endfor
end
