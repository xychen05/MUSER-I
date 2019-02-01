;.r fringe_stoping.pro
; 
; Name :
; fringe_stoping
;
; Purpose :
; Do the fringe stoping of the original data, COR_DATA
; 
; History :
; Wei Wang, FRINGE_STOPING.PRO
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

temp=replicate('num0060',8)
data_name=data_name+temp

data_namet=data_name[6]

;sys=0 ;for windows system
sys=1 ;for linux   system

;;;;============move cal data to the cal file
;if sys eq 0 then file_cal=file_search('.\'+data_namet+'\*_1217-1730.sav',count=count)
;if sys eq 1 then file_cal=file_search('./'+data_namet+'/*_1217-1730.sav',count=count)
;if sys eq 0 then data_sav='.\cal\'
;if sys eq 1 then data_sav='./cal/'
;file_mkdir,data_sav
;if file_cal ne '' then file_move,[file_cal],data_sav
;;;;============for fs data
;if sys eq 0 then file=file_search('.\'+data_namet+'\*.sav',count=count)
;if sys eq 1 then file=file_search('./'+data_namet+'/*.sav',count=count)
;
;if sys eq 0 then data_sav='.\'+data_namet+'\fs\'
;if sys eq 1 then data_sav='./'+data_namet+'/fs/'
;file_mkdir,data_sav


;;;============for cal data
if sys eq 0 then file=file_search('.\cal\*.sav',count=count)
if sys eq 1 then file=file_search('./cal/*.sav',count=count)
if sys eq 0 then data_sav='.\cal\'
if sys eq 1 then data_sav='./cal/'

nfile=count

for ff=0,nfile-1 do begin;count-1
  restore,file[ff],/ver
;-------------------------hi-------------------------
  ;for int file
  ;COR_DATA=COR_DATA_INT
  ;GPS_TIME=GPS_TIME_INT
  
  DELAY_NS=[-56,  0, 48,921, 13, 59,  -3, 460, 49, 69,-675,-157,  363, -65,$
    30, 42, 51, 121, -2,  73,  35, 26,  74, 35, -3,   47, -71,$
    75,343, 56, 32,313,  678,  12,-30, 48,-18, 20,   10, -1666,$
    0,  0,  0,  0]
    
  DELAY_NS_ADD=[0,0,-30,-20, 10, 9,  -15, -1, -48, -33, 0, 186, 0, 0,$  ;20140120
    0,  -26, 0, 0, 0, 0, -24, -1, -57, 0, 0, 3, 66,$
    -12, 0, -20, -1,  0, -22, 0, 7, 10, 2, 1,  0, 0,$
    0,  0,  0,  0]
    
  DELAY_NS_ADD1=[0,0,0,0, 56,56,56,56,0,0,0,0,56,56,$  ;20140120
    56,56, 0, 0, 0, 0,56,56,56,56, 0, 0, 0,$
    0,56,56,56,56, 0,0,0,0,56,56,56,56,$
    0,  0,  0,  0]
  DELAY_NS_ADD2=[0,0,-102,-50,-51,-47,-48,-50,-62,-64,0,0,0,0,$  ;20141015
    -57,-57,-60,-60,-60,-60,-60,-69,-70,-54,-37,-47,0,$
    -60,-62,-52,-61,-59,-59,-49,-63,-65,-59,-60,-67,0,$
    0,0,0,0]
  DELAY_NS=DELAY_NS+DELAY_NS_ADD+DELAY_NS_ADD1+DELAY_NS_ADD2
  
  nn=size(gps_time)
  n_file=nn[1]
  gps_time_fs=make_array(nn[1],nn[2])
  gps_time_fs=gps_time
  
  delay=dblarr(n_file,44)
  
  for i=2,n_file-3 do begin
    delay[i,*]=par_delay[i,*]-delay_ns[*]
  endfor
  
  delay1=fix(delay*10000d,type=3)/10000d
  
  cor_data_fs=dblarr(n_file,40,40,2)
  
  pi=3.1415926535897932d
  ch=4;just for 4th frequency in 1.6-2.0_R (16 in total)
  for i=0,39 do begin
    for j=i+1,39 do begin
      for tt=0,n_file-1 do begin
        Frf=(1600+ch*25+12.5d)/1000.  ;12.5改为2.5   ;frf=1.7125  ;(400+ch*25+12.5d)/1000. ;for 0.4-0.8 channel
        Fif=(ch*25+12.5+50)/1000d     ;fif=0.1625
        tg=delay1[tt,j]-delay1[tt,i]
        tg0=fix(delay1[tt,j])-fix(delay1[tt,i])
        phai=2*pi*(Frf*tg-Fif*tg0)
        ;print,i,j,ch,frf,fif,tg,tg0,phai,format='(f20.6)'
        cor_data_fs[tt,i,j,0]=cor_data[tt,ch,i,j,0]*cos(phai)+cor_data[tt,ch,i,j,1]*sin(phai)
        cor_data_fs[tt,i,j,1]=cor_data[tt,ch,i,j,1]*cos(phai)-cor_data[tt,ch,i,j,0]*sin(phai)
        
      endfor
    endfor
  endfor
  print,'god',string(ff)
  
  
  tmp0=rstrpos(file[ff],'COR_data_')
;  ;;-------------------------for fs data-------------------------
;  ;;for fs file 
;  save,cor_data_fs,gps_time_fs,delay,filename=data_sav+'cor_data_fs_'+strmid(file[ff],tmp0+9,33)+'.sav'
;  print,'data has saved .....'

  ;;-------------------------for cal data-------------------------
  ;;for cal file
  cor_data_cal=cor_data_fs
  gps_time_cal=gps_time_fs
  save,cor_data_cal,gps_time_cal,filename=data_sav+'cor_data_cal_'+strmid(file[ff],tmp0+9,33)+'.sav'
  print,'data has saved .....'

endfor

print,'ready.......................'
end