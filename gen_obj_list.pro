pro gen_obj_list

	chk_gz = file_test('spSlit*fits.gz') 

	binfile = file_search('*bintabs.fits')
	bin_obj = mrdfits(binfile,1) 
	bin_slit = mrdfits(binfile,3) 

	sltmak = strsplit(file_search('*.plan'),'.',/extract)
	openw,11,'slit_objs_'+sltmak(0)+'.txt'
	printf,11,'Slit#       ObjectID   Chip#  Obj_type  N_exp        obj_RA       obj_dec       slit_RA      slit_dec'
	for i=0,n_elements(bin_obj)-1 do begin
		Bfile = 'spSlit.'+sltmak(0)+'.'+string(i,format='(I03)')+'B.fits'
		Rfile = 'spSlit.'+sltmak(0)+'.'+string(i,format='(I03)')+'R.fits'
		if chk_gz then Bfile = Bfile+'.gz'
		if chk_gz then Rfile = Rfile+'.gz'
		if file_test(Bfile) and file_test(Rfile) then begin
			s = mrdfits(Bfile,1,hd,/silent)
			slitnum = SXPAR(hd,'SLITNO')
			chipnum = SXPAR(hd,'CHIPNO')
			obstype = SXPAR(hd,'obstype')
			fits_info,Bfile,N_ext = n_ext,/silent
			printf,11,i,bin_obj(i).object,chipnum,obstype,n_ext/2,bin_obj(i).ra_obj,bin_obj(i).dec_obj,bin_slit(i).SLITRA,bin_slit(i).slitdec,format='(I5,3x,A12,I8,3x,A9,I5,4F14.8)'
			print,i,bin_obj(i).object,chipnum,obstype,n_ext/2,bin_obj(i).ra_obj,bin_obj(i).dec_obj,bin_slit(i).SLITRA,bin_slit(i).slitdec,format='(I5,3x,A12,I8,3x,A9,I5,4F14.8)'
		endif
	endfor
	close,11

end
