;.r rm_int_cor_data.pro
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

for ooo=6,6 do begin
data_namet=data_name[ooo]

sys=0 ;for windows system
sys=1 ;for linux   system

nm_ind='17'
nm_ind='17_test'

;;ch=0
;for ch=4,4 do begin
;nm_ind='17sf'+strcompress(ch,/remove_all)

num=2400l
;
;;;;-------------------------file name for windows-------------------------
;;if sys eq 0 then filename=file_search('.\'+data_namet+'\prep\*.sav',count=count)
;;if sys eq 0 then filename=file_search('.\'+data_namet+'\prepsf'+strcompress(ch,/remove_all)+'\*.sav',count=count)
;;if sys eq 0 then filename_rm=file_search('.\sav\'+data_namet+'\rm_time'+nm_ind+'\*.sav',count=count)
;
;;;-------------------------file name for linux-------------------------
;if sys eq 1 then filename=file_search('./'+data_namet+'/prep/*.sav',count=count)
;;if sys eq 1 then filename=file_search('./'+data_namet+'/prepsf'+strcompress(ch,/remove_all)+'/*.sav',count=count)
;if sys eq 1 then filename_rm=file_search('./sav/'+data_namet+'/rm_time17/rm_time*.sav',count=count);'+nm_ind+'
;
;;-------------------------basic para-------------------------
;if sys eq 0 then data_sav='.\sav\'+data_namet+'\rm_cor_data'+nm_ind+'\'
;if sys eq 1 then data_sav='./sav/'+data_namet+'/rm_cor_data'+nm_ind+'/'
;file_mkdir,data_sav
;
;print,count
;nfile=count
;;nfile=1
;
;for ii=0,nfile-1 do begin ;nfile-1 ; for which file
;
;;-------------------------for praparing rm data-------------------------
;  tmp0=rstrpos(filename[ii],'cor_data_')
;  file=strmid(filename[ii],tmp0,80)
;  print,file
;  restore,filename[ii]
;  restore,filename_rm[ii];,/ver
;  
;  ;;============for case of rm_time
;  sz=num-(size(rm_time))[1]
;  cor_data_rm=dblarr(sz,40,40,2)
;  gps_time_rm=intarr(sz,9)
;  ind_temp=indgen(num)
;  remove,rm_time,ind_temp
;  for rr=0,sz-1 do begin
;    cor_data_rm[rr,*,*,*]=cor_data_prep[ind_temp[rr],*,*,*]
;    gps_time_rm[rr,*]=gps_time_prep[ind_temp[rr],*]
;  endfor
;  help,cor_data_rm
;  
;;  ;;============for case of rm_time1
;;  sz1=sz-(size(rm_time1))[1]
;;  cor_data_rm1=dblarr(sz1,40,40,2)
;;  gps_time_rm1=intarr(sz1,9)
;;  ind_temp=indgen(sz)
;;  remove,rm_time1,ind_temp
;;  for rr=0,sz1-1 do begin
;;    cor_data_rm1[rr,*,*,*]=cor_data_rm[ind_temp[rr],*,*,*]
;;    gps_time_rm1[rr,*]=gps_time_rm[ind_temp[rr],*]
;;  endfor
;;  help,cor_data_rm1
;;  cor_data_rm=cor_data_rm1
;;  gps_time_rm=gps_time_rm1
;  
;  ;;============
;  save,cor_data_rm,gps_time_rm,filename=data_sav+'cor_data_rm_'+strmid(file,14,50)
;    
;endfor
;print,'rm data has prepared...',data_sav


;;;-------------------------integration for 1 min-------------------------
;;;for rm data
;if sys eq 0 then filename=file_search('.\sav\'+data_namet+'\rm_cor_data'+nm_ind+'\*.sav',count=count)
;if sys eq 1 then filename=file_search('./sav/'+data_namet+'/rm_cor_data'+nm_ind+'/*.sav',count=count)
;
;if sys eq 0 then data_sav='.\sav\'+data_namet+'\it_cor_data'+nm_ind+'\'
;if sys eq 1 then data_sav='./sav/'+data_namet+'/it_cor_data'+nm_ind+'/'
;file_mkdir,data_sav
;
;print,count
;nfile=count
;int=60;;file number of 1 min, int=6 means 10s/file
;
;for kk=0,nfile-1 do begin ;count/kk_int-1
;  restore,filename[kk];,/ver
;  tmp0=rstrpos(filename[kk],'cor_data_')
;  file=strmid(filename[kk],tmp0,80)
;  print,file
;  cor_data_it=cor_int(cor_data_rm,int)
;  gps_time_it=gps_time_int(gps_time_rm,int)
;  
;  save,cor_data_it,gps_time_it,filename=data_sav+'cor_data_it_'+file
;
;endfor
;print,'it data has prepared...',data_sav


;;-------------------------deal with the keep_table-------------------------
;;for keep_table
if sys eq 0 then filename_kp=file_search('.\sav\'+data_namet+'\keep_table'+nm_ind+'\*.sav',count=count) & $
  keep_sav='.\sav\'+data_namet+'\'
if sys eq 1 then filename_kp=file_search('./sav/'+data_namet+'/keep_table'+nm_ind+'/*.sav',count=count) & $
  keep_sav='./sav/'+data_namet+'/'

keep_table=crt_keep_table(filename_kp)
keep_table_sub=crt_keep_table(filename_kp,/sub)

save,keep_table,keep_table_sub,filename=keep_sav+'keep_table'+nm_ind+'.sav'
help,keep_table

;restore,'./sav/1.6-2.0_Lnum2400/keep_table17.sav',/ver
;help,keep_table

;;;-------------------------integration for one file-------------------------
;;;for it data
;if sys eq 0 then filename=file_search('.\sav\'+data_namet+'\it_cor_data'+nm_ind+'\*.sav',count=count)
;if sys eq 1 then filename=file_search('./sav/'+data_namet+'/it_cor_data'+nm_ind+'/*.sav',count=count)
;
;if sys eq 0 then data_sav='.\sav\'+data_namet+'\it_cor_data00\'
;if sys eq 1 then data_sav='./sav/'+data_namet+'/it_cor_data00/'
;file_mkdir,data_sav
;
;tmp0=rstrpos(filename[0],'cor_data_')
;file=strmid(filename[0],tmp0,80)
;
;restore,filename[0],/ver
;
;nfile=15
;int=(size(cor_data_it))[1]
;
;cor_data=make_array(int*nfile,40,40,2)
;gps_time=intarr(int*nfile,9)
;for ff=0,nfile-1 do begin
;  restore,filename[ff];,/ver
;  print,filename[ff]
;  cor_data[(ff*int):((ff+1)*int-1),*,*,*]=cor_data_it[*,*,*,*]
;  gps_time[(ff*int):((ff+1)*int-1),*]=gps_time_it[*,*]
;  
;  if ff eq 0 then st_time=time_gps_to_string(gps_time_it[0,*])
;  if ff eq nfile-1 then ed_time=time_gps_to_string(gps_time_it[0,*])
;endfor
;cor_data_it=cor_data
;gps_time_it=gps_time
;
;save,cor_data_it,gps_time_it,filename=data_sav+'cor_data_it'+nm_ind+'_'+data_namet+ $
;  '_'+strmid(st_time,0,15)+'-'+strmid(ed_time,9,6)+'.sav'
;print,'it data has prepared...',data_sav
;
;if sys eq 0 then data_savqq='.\sav\'+data_namet+'\quiet\'
;if sys eq 1 then data_savqq='./sav/'+data_namet+'/quiet/'
;file_mkdir,data_savqq
;
;cor_data_it11=cor_int(cor_data_it[0:(int*nfile-1),*,*,*],1)
;gps_time_it11=gps_time_int(gps_time_it[0:(int*nfile-1),*],1)
;cor_data_it=cor_data_it11
;gps_time_it=gps_time_it11
;save,cor_data_it,gps_time_it,filename='./sav/'+data_namet+'/quiet/cor_data_it'+nm_ind+'_quiet.sav'


;-------------------------ending-------------------------
;endfor
endfor
print,'OK'
end