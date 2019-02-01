;.r muser2fits.pro
;.compile package_muser.pro 
; 
; Name :
; muser2fits
;
; Purpose :
; transform the file .sav of clean image to fits file.
; attentions, no rotate.
;
; Inputs :
; inputfile = filename_sav
; 
; Outputs :
; outputfile = filename_fits
;
; History :
; Xingyao Chen, July 2017.
; Xingyao Chen, 05 Dec 2017.
;-------------------------begins-------------------------
data_name=['0.4-0.8_L',$
  '0.4-0.8_R',$
  '0.8-1.2_L',$
  '0.8-1.2_R',$
  '1.2-1.6_L',$
  '1.2-1.6_R',$
  '1.6-2.0_L',$
  '1.6-2.0_R']
temp=replicate('num2400',8)
data_name=data_name+temp

data_namet=data_name[6]

;sys=0 ;for windows system
sys=1 ;for linux   system

ktkt=['_test_20_100_1500','_test_20_100_2000', $
  '_test_20_180_1500','_test_20_180_2000', $
  '_test_30_100_1500','_test_30_100_2000', $
  '_test_30_180_1500','_test_30_180_2000', $
  '_test_30','_test_40']


nm_ind='17'+ktkt[9]

;ch=0
for ch=4,4 do begin
  ;nm_ind='sf'+strcompress(ch,/remove_all)


;;;;-------------------------file name for windows or linux-------------------------
;;;;============;;for multi-frequencies imaging data
if sys eq 0 then begin
  filename_sav=file_search('.\sav\'+data_namet+'\it_clean_img'+nm_ind+'\*.sav',count=count)
  data_fits='.\sav\'+data_namet+'\it_clean_img'+nm_ind+'_fits\'
  data_raw_fits='.\sav\'+data_namet+'\it_clean_img'+nm_ind+'_fits_raw\'
endif
if sys eq 1 then begin
  filename_sav=file_search('./sav/'+data_namet+'/it_clean_img'+nm_ind+'/*.sav',count=count)
  data_fits='./sav/'+data_namet+'/it_clean_img'+nm_ind+'_fits/'
  data_raw_fits='./sav/'+data_namet+'/it_clean_img'+nm_ind+'_fits_raw/'
endif

;;;;============;;for quiet_clean_img data
;if sys eq 0 then begin
;  filename_sav=file_search('.\sav\'+data_namet+'\quiet_clean_img\clean_'+nm_ind+'*.sav',count=count)
;  data_fits='.\sav\'+data_namet+'\quiet_clean_img_fits\'
;  data_raw_fits='.\sav\'+data_namet+'\quiet_clean_img_fits\'
;endif
;if sys eq 1 then begin
;  filename_sav=file_search('./sav/'+data_namet+'/quiet_clean_img/clean_'+nm_ind+'*.sav',count=count)
;  data_fits='./sav/'+data_namet+'/quiet_clean_img_fits/'
;  data_raw_fits='./sav/'+data_namet+'/quiet_clean_img_fits/'
;endif

;;;;============;;for residual image from quiet_clean_img data
;if sys eq 0 then begin
;  filename_sav=file_search('.\sav\'+data_namet+'\quiet_clean_img\clean_'+nm_ind+'*.sav',count=count)
;  data_fits='.\sav\'+data_namet+'\quiet_clean_img_fits_residual\'
;  data_raw_fits='.\sav\'+data_namet+'\quiet_clean_img_fits_residual\'
;endif
;if sys eq 1 then begin
;  filename_sav=file_search('./sav/'+data_namet+'/quiet_clean_img/clean_'+nm_ind+'*.sav',count=count)
;  data_fits='./sav/'+data_namet+'/quiet_clean_img_fits_residual/'
;  data_raw_fits='./sav/'+data_namet+'/quiet_clean_img_fits_residual/'
;endif

;;;============
;file_mkdir,data_fits
file_mkdir,data_raw_fits

;count=1
;-------------------------begin-------------------------
for cc=0,count-1 do begin;count-1
  ;muser2fits, filename_sav[cc], data_fits,1024,/fovaia
  muser2fits, filename_sav[cc], data_raw_fits,512,/left, /band4;,/fovaia
  ;muser2fits, filename_sav[cc], data_raw_fits,512,/right, /band4;,/fovaia
  ;muser2fits, filename_sav[cc], data_raw_fits,512,/right, /band3;,/fovaia
  ;muser2fits, filename_sav[cc], data_raw_fits,512,/left, /band3, /residual;,/fovaia ;;;for residual data

endfor

;-------------------------end-------------------------
endfor
print,'fits has saved to file :',data_fits
print,'it is ok.............'
end