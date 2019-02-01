;.r rm_int_cor_data17.pro
;.compile package_muser.pro
; 
; Name :
; rm_int_cor_data
; 
; Purpose :
; Do the rm data and do the integration for 1 min and for one file.
; 
; Explanation :
; Preparing rm data and it data.
;
; History :
; Xingyao Chen, 31 Aug 2017
; 
;-------------------------begins-------------------------
Device, RETAIN=2
!P.FONT = 0

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

for ooo=7,7 do begin
data_namet=data_name[ooo]

;sys=0 ;for windows system
sys=1 ;for linux   system

nm_ind='05'
nm_ind_sav='17'
;;ch=0
for ch=4,4 do begin
nm_ind_sf='sf'+strcompress(ch,/remove_all)

num=2400l

;;;-------------------------file name for windows-------------------------
if sys eq 0 then filename=file_search('.\'+data_namet+'\prep\*.sav',count=count)
if sys eq 0 then filename_sf=file_search('.\'+data_namet+'\prepsf'+strcompress(ch,/remove_all)+'\*.sav',count=count)
if sys eq 0 then filename_rm=file_search('.\sav\'+data_namet+'\rm_time'+nm_ind+'\*.sav',count=count)
if sys eq 0 then filename_rm_sf=file_search('.\sav\'+data_namet+'\rm_time'+nm_ind_sf+'\*.sav',count=count)

;;-------------------------file name for linux-------------------------
if sys eq 0 then filename=file_search('./'+data_namet+'/prep/*.sav',count=count)
if sys eq 1 then filename_sf=file_search('./'+data_namet+'/prepsf'+strcompress(ch,/remove_all)+'/*.sav',count=count)
if sys eq 1 then filename_rm=file_search('./sav/'+data_namet+'/rm_time'+nm_ind+'/rm_time*.sav',count=count)
if sys eq 1 then filename_rm_sf=file_search('./sav/'+data_namet+'/rm_time'+nm_ind_sf+'/rm_time*.sav',count=count)

;-------------------------basic para-------------------------
if sys eq 0 then data_sav='.\sav\'+data_namet+'\rm_cor_data17\'
if sys eq 1 then data_sav='./sav/'+data_namet+'/rm_cor_data17/'
file_mkdir,data_sav

print,count
nfile=count
;nfile=1
for ii=0,nfile-1 do begin
  restore,filename_rm[ii];,/ver
  tmp0=rstrpos(filename_rm[ii],'rm_time_')
  file=strmid(filename_rm[ii],tmp0,80)
  print,file
  help,rm_time
  array1=rm_time
  restore,filename_rm_sf[ii];,/ver
  help,rm_time
  array2=rm_time
  temp = [array1, array2]
  rm_time = temp(uniq(temp, sort(temp)))
  help,rm_time
  save,rm_time,filename='./sav/'+data_namet+'/rm_time17/'+file 
endfor

;-------------------------ending-------------------------
endfor
endfor
print,'OK'
end