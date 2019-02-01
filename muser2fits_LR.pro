;.r muser2fits_LR.pro
;.compile package_muser.pro 
; 
; Purpose :
; L+R
; transform the file .sav of clean image to fits file.
; attentions, no rotate.
;
; Keywords :
; /pixaia : transform pixels size to AIA image, [4096,4096]
;
; History :
; Xingyao Chen, July 2017.
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

nm_ind='05'
nm_ind='17'

;ch=0
for ch=4,4 do begin
  ;nm_ind='sf'+strcompress(ch,/remove_all)

;;;-------------------------step 1: generate sav file-------------------------
;if sys eq 0 then begin
;  data_namet=data_name[7]
;  filename_savR=file_search('.\sav\'+data_namet+'\it_clean_img'+nm_ind+'\*.sav',count=count)
;  data_namet=data_name[6]
;  filename_savL=file_search('.\sav\'+data_namet+'\it_clean_img'+nm_ind+'\*.sav',count=count)
;  data_sav_LR='.\sav\'+data_namet+'\it_clean_img'+nm_ind+'_sav_LR\'
;endif
;if sys eq 1 then begin
;  data_namet=data_name[7]
;  filename_savR=file_search('./sav/'+data_namet+'/it_clean_img'+nm_ind+'/*.sav',count=count)
;  data_namet=data_name[6]
;  filename_savL=file_search('./sav/'+data_namet+'/it_clean_img'+nm_ind+'/*.sav',count=count)
;  data_sav_LR='./sav/'+data_namet+'/it_clean_img'+nm_ind+'_sav_LR/'
;endif
;
;file_mkdir,data_sav_LR
;
;;;;============
;nfile=count
;
;for cc=0,nfile-1 do begin;count-1
;  restore,filename_savL[cc];,/ver
;  dimgl=dimg
;  imgl=img
;  img_cll=img_cl
; 
;  tmp0=rstrpos(filename_savL[cc],'clean_'+nm_ind+'_')
;  file=strmid(filename_savL[cc],tmp0,60)
;  print,file
;  print,cc,minmax(img_cl)
;  
;  restore,filename_savR[cc];,/ver
;  dimgr=dimg
;  imgr=img
;  img_clr=img_cl
;  print,cc,minmax(img_cl)
;  
;  dimg=dimgl+dimgr
;  img=imgl+imgr
;  img_cl=img_cll+img_clr
;  print,cc,minmax(img_cl)
;  
;  save,dimg,img,img_cl,gps_time_cl,ch,filename=data_sav_LR+file
;  
;endfor
;
;print,'.sav has saved to file :',data_sav_LR


;;-------------------------step 1: generate fits file-------------------------
data_namet=data_name[6]
if sys eq 0 then begin
  filename_sav=file_search('.\sav\'+data_namet+'\it_clean_img'+nm_ind+'_sav_LR\*.sav',count=count)
  data_fits_LR='.\sav\'+data_namet+'\it_clean_img'+nm_ind+'_fits_LR\'
  data_raw_fits_LR='.\sav\'+data_namet+'\it_clean_img'+nm_ind+'_fits_raw_LR\'
endif
if sys eq 1 then begin
  filename_sav=file_search('./sav/'+data_namet+'/it_clean_img'+nm_ind+'_sav_LR/*.sav',count=count)
  data_fits_LR='./sav/'+data_namet+'/it_clean_img'+nm_ind+'_fits_LR/'
  data_raw_fits_LR='./sav/'+data_namet+'/it_clean_img'+nm_ind+'_fits_raw_LR/'
endif

;file_mkdir,data_fits_LR
file_mkdir,data_raw_fits_LR

;;;============
nfile=count

for cc=0,nfile-1 do begin;count-1
  ;muser2fits, filename_sav[cc], data_fits_LR,1024,/fovaia
  muser2fits, filename_sav[cc], data_raw_fits_LR,512, /lr, /band4;,/fovaia
endfor

print,'.sav has saved to file :',data_fits_LR

;-------------------------file name for linux-------------------------
endfor

print,'it is ok.............'
end