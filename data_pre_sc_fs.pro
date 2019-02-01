;.r data_pre_sc_fs.pro
;.compile package_muser.pro
;
; History :
; Xingyao Chen, 16 Dec 2017
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

for ooo=6,6 do begin
data_namet=data_name[ooo]

sys=0 ;for windows system
sys=1 ;for linux   system

nm_ind='17'

num=2400l
;;;-------------------------file name for windows-------------------------
;if sys eq 0 then filename=file_search('.\'+data_namet+'\fs\*.sav',count=count)
;if sys eq 0 then filename_rm=file_search('.\sav\'+data_namet+'\rm_time'+nm_ind+'\*.sav',count=count)
;;-------------------------file name for linux-------------------------
if sys eq 1 then filename=file_search('./'+data_namet+'/fs/*.sav',count=count)
if sys eq 1 then filename_rm=file_search('./sav/'+data_namet+'/rm_time17/*.sav',count=count);'+nm_ind+'

;-------------------------basic para-------------------------
if sys eq 0 then data_sav='.\'+data_namet+'\fs'+nm_ind+'\'
if sys eq 1 then data_sav='./'+data_namet+'/fs'+nm_ind+'/'
file_mkdir,data_sav

print,count
nfile=count
nfile=1

for ii=15,15+20-1 do begin ;nfile-1 ; for which file

  tmp0=rstrpos(filename[ii],'cor_data_')
  file=strmid(filename[ii],tmp0,80)
  print,file
  restore,filename[ii]
  restore,filename_rm[ii];,/ver
  
  ;;============for case of rm_time
  sz=num-(size(rm_time))[1]
  cor_data_rm=dblarr(sz,40,40,2)
  gps_time_rm=intarr(sz,9)
  ind_temp=indgen(num)
  remove,rm_time,ind_temp
  for rr=0,sz-1 do begin
    cor_data_rm[rr,*,*,*]=cor_data_fs[ind_temp[rr],*,*,*]
    gps_time_rm[rr,*]=gps_time_fs[ind_temp[rr],*]
  endfor
  help,cor_data_rm
  
  ;;============integration
  int=6
  cor_data_it=cor_int(cor_data_rm,int)
  gps_time_it=gps_time_int(gps_time_rm,int)
  help,cor_data_it
  
  ;;============
  cor_data_fs=dblarr(1,40,40,2)
  gps_time_fs=intarr(1,9)
  for mm=0,int-1 do begin
    cor_data_fs=cor_data_it[mm,*,*,*]
    gps_time_fs=gps_time_it[mm,*]
    time_obs=time_gps_to_string(gps_time_fs)
    save,cor_data_fs,gps_time_fs,filename=data_sav+'cor_data_fs_it'+nm_ind+'_'+time_obs+'.sav'
  endfor
 
endfor

;-------------------------ending-------------------------
endfor
print,'OK'
end