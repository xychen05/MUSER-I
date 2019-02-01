;.r data_pre_process.pro
;.compile package_muser.pro
;
; Name :
; data_pre_process
;
; Purpose :
; Do the phase calibration before image processing of dirty map.
; 
; Inputs :
; Original data both the data for solar and satellite after fringe stoping.
;
; Outputs :
; prep data.
;
; Notes :
; Select good frames for the satellite data.
;
; History :
; Wei Wang, DATA_PRE_PROCESS.PRO
; Xingyao Chen, June 2017
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

data_namet=data_name[7]

;sys=0 ;for windows system
sys=1 ;for linux   system

;;;-------------------------file name for windows-------------------------
if sys eq 0 then begin
  file=file_search('.\'+data_namet+'\fs\*.sav',count=count)
  data_nametc=strmid(data_namet,0,9)+'_int01_'+strmid(data_namet,9,9)
  file_cal=file_search('.\cal\cor_data_cal_'+data_nametc+'*.sav')
  restore,file_cal,/ver
  
  data_sav='.\'+data_namet+'\prep\'
  file_mkdir,data_sav
  file_mkdir,'.\sav\'+data_namet
  file_mkdir,'.\imagetest\'+data_namet
endif
;;;-------------------------file name for linux-------------------------
if sys eq 1 then begin
  file=file_search('./'+data_namet+'/fs/*.sav',count=count)
  data_nametc=strmid(data_namet,0,9)+'_int01_'+strmid(data_namet,9,9)
  file_cal=file_search('./cal/cor_data_cal_'+data_nametc+'*.sav')
  
  restore,file_cal,/ver
  data_sav='./'+data_namet+'/prep/'
  file_mkdir,data_sav
  file_mkdir,'./sav/'+data_namet
  file_mkdir,'./imagetest/'+data_namet
endif

szcal=(size(cor_data_cal))[1]/6
;-------------------------begins-------------------------
for cc=0,count-1 do begin;count-1
  
  RESTORE,file[cc],/ver
  tmp0=rstrpos(file[cc],'cor_data_')
  file_fs=strmid(file[cc],tmp0,60)
  print,file_fs
  
  nsize=size(gps_time_fs)
  cadence=1
  ;cadence=10
  n=nsize[1]/cadence
  CH=4;for frequency 1.725GHz in band 1.6-2.0

  CALI_DATA=DBLARR(40,40,2)
  for kk=0,39 do begin
    for jj=0,39 do begin
      CALI_DATA[kk,jj,0]=mean(cor_data_cal[szcal:szcal*5,kk,jj,0]) ;select good frames.
      CALI_DATA[kk,jj,1]=mean(cor_data_cal[szcal:szcal*5,kk,jj,1]) ;select good frames.
    endfor
  endfor

  cor_data_prep=DBLARR(N,40,40,2)
  AMP=DBLARR(N,40,40)
  PHASE=DBLARR(N,40,40)
  PHASE_CALI=DBLARR(40,40)
  gps_time_prep=INTARR(N,9)

  FOR I=0,N-1 DO BEGIN
    cor_data_prep[I,*,*,*]=cor_data_fs[I*CADENCE,*,*,*]
    gps_time_prep[I,*]=gps_time_fs[I*CADENCE,*]
  ENDFOR

  AMP[*,*,*]=SQRT(cor_data_prep[*,*,*,0]^(2D)+cor_data_prep[*,*,*,1]^(2D))
  PHASE[*,*,*]=ATAN(cor_data_prep[*,*,*,1],cor_data_prep[*,*,*,0])

  PHASE_CALI[*,*]=ATAN(CALI_DATA[*,*,1],CALI_DATA[*,*,0])

  FOR I=0,N-1 DO BEGIN
    cor_data_prep[I,*,*,0]=AMP[I,*,*]*COS(PHASE[I,*,*]-PHASE_CALI[*,*])
    cor_data_prep[I,*,*,1]=AMP[I,*,*]*SIN(PHASE[I,*,*]-PHASE_CALI[*,*]);PHASE_CALI[*,*] from 'CALIBRATOR.SAV'
  ENDFOR

  save,cor_data_prep,gps_time_prep,filename=data_sav+'cor_data_prep_'+strmid(file_fs,12,60)
  
  print,'data has saved .....'
 
endfor

;restore,'./prep/'+'cor_data_prep_'+strmid(file_fs,12,30),/ver
print,'hahaha.........'
END
