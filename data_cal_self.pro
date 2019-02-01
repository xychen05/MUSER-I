;.r data_cal_self.pro
;.compile package_muser.pro
;
; Purpose :
; Prepare for the self-correlation calibration.
; 
; History :
; Wei Wang, DATA_PRE_PROCESS.PRO
; Xingyao Chen, 24 Nov 2017
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

sys=0 ;for windows system
;sys=1 ;for linux   system

if sys eq 0 then begin  
  data_namet=data_name[2]
  file_mkdir,'.\'+data_namet+'\calsf\'
  file_pwr=file_search('I:\muserdata_20141217\muser_spec\sav\spec_data_L_122855_123-123855_100.sav')
  restore,file_pwr,/ver
  for ch=0,15 do begin
    pwrt=pwr_specl[*,16+ch,0]
    indp=7610;where(pwrt eq max(pwrt))
    print,gps_time_specl[indp,*]
    file_call=file_search('.\'+data_namet+'\fssf'+strcompress(ch,/remove_all)+'\cor_data_fs_sf'+strcompress(ch,/remove_all)+'_'+strmid(data_namet,0,8)+'L_int01_num2400_1217-1230.sav')
    restore,file_call,/ver
    tmp0=rstrpos(file_call,'cor_data_fs_sf')
    filet=strmid(file_call,tmp0+12,100)
    cor_data_cal=cor_data_fs[(indp mod 2400),*,*,*]
    gps_time_cal=gps_time_fs[(indp mod 2400),*]
    print,gps_time_fs[(indp mod 2400),*]
    save,cor_data_cal,gps_time_cal,filename='.\'+data_namet+'\calsf\cor_data_cal_'+filet
  endfor
  
  data_namet=data_name[3]
  file_mkdir,'.\'+data_namet+'\calsf\'
  file_pwr=file_search('I:\muserdata_20141217\muser_spec\sav\spec_data_R_122855_127-123855_103.sav')
  restore,file_pwr,/ver
  for ch=0,15 do begin
    pwrt=pwr_specr[*,16+ch,0]
    indp=7928;where(pwrt eq max(pwrt))
    print,gps_time_specr[indp,*]
    file_calr=file_search('.\'+data_namet+'\fssf'+strcompress(ch,/remove_all)+'\cor_data_fs_sf'+strcompress(ch,/remove_all)+'_'+strmid(data_namet,0,8)+'R_int01_num2400_1217-1230.sav')
    restore,file_calr,/ver
    tmp0=rstrpos(file_calr,'cor_data_fs_sf')
    filet=strmid(file_calr,tmp0+12,100)
    cor_data_cal=cor_data_fs[(indp mod 2400),*,*,*]
    gps_time_cal=gps_time_fs[(indp mod 2400),*]
    print,gps_time_fs[(indp mod 2400),*]
    save,cor_data_cal,gps_time_cal,filename='.\'+data_namet+'\calsf\cor_data_cal_'+filet
  endfor 
endif

;if sys eq 1 then begin
;  data_namet=data_name[4]
;  file_mkdir,'./'+data_namet+'/calsf/'
;  file_pwr=file_search('/Volumes/Seagate_cxy/muserdata_20141217/muser_spec/sav/spec_data_L_122855_123-123855_100.sav')
;  restore,file_pwr,/ver
;  for ch=12,15 do begin
;    pwrt=pwr_specl[*,32+ch,0]
;    indp=where(pwrt eq max(pwrt));;7610  2014-12-17T12:32:05.374
;    print,gps_time_specl[indp,*]
;    file_call=file_search('./'+data_namet+'/fssf'+strcompress(ch,/remove_all)+'/cor_data_fs_sf'+strcompress(ch,/remove_all)+'_'+strmid(data_namet,0,8)+'L_int01_num2400_1217-1230.sav')
;    restore,file_call,/ver
;    tmp0=rstrpos(file_call,'cor_data_fs_sf')
;    filet=strmid(file_call,tmp0+12,100)
;    cor_data_cal=cor_data_fs[(indp mod 2400),*,*,*]
;    gps_time_cal=gps_time_fs[(indp mod 2400),*]
;    print,gps_time_fs[(indp mod 2400),*]
;    save,cor_data_cal,gps_time_cal,filename='./'+data_namet+'/calsf/cor_data_cal_'+filet
;  endfor
;
;  data_namet=data_name[5]
;  file_mkdir,'./'+data_namet+'/calsf/'
;  file_pwr=file_search('/Volumes/Seagate_cxy/muserdata_20141217/muser_spec/sav/spec_data_R_122855_127-123855_103.sav')
;  restore,file_pwr,/ver
;  for ch=14,15 do begin
;    pwrt=pwr_specr[*,32+ch,0]
;    indp=where(pwrt eq max(pwrt));;8520  2014-12-17T12:32:28.127
;    
;    file_calr=file_search('./'+data_namet+'/fssf'+strcompress(ch,/remove_all)+'/cor_data_fs_sf'+strcompress(ch,/remove_all)+'_'+strmid(data_namet,0,8)+'R_int01_num2400_1217-1230.sav')
;    restore,file_calr,/ver
;    tmp0=rstrpos(file_calr,'cor_data_fs_sf')
;    filet=strmid(file_calr,tmp0+12,100)
;    cor_data_cal=cor_data_fs[(indp mod 2400),*,*,*]
;    gps_time_cal=gps_time_fs[(indp mod 2400),*]
;    print,gps_time_specr[indp,*]
;    print,gps_time_fs[(indp mod 2400),*]
;    save,cor_data_cal,gps_time_cal,filename='./'+data_namet+'/calsf/cor_data_cal_'+filet
;  endfor
;endif


print,'hahaha.........'
END
