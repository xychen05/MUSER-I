;.r plt_data_automode_new.pro
;.compile package_muser.pro
;
; Name :
; plt_data_automode_new_1min
; 
; Purpose :
; Transform the MUSER raw data into .sav format.
; 
; Inputs :
; MUSER raw data.
;
; Outputs :
; .sav file, which contains COR_DATA, PWR, GPS_TIME, PAR_DELAY.
; 
; Notes :
; This procedure is for the new vision of MUSER raw data ....
; 
; History :
; Wei Wang, 
; Xingyao Chen, June 2017
;
;-------------------------begins-------------------------
systim = SYSTIME(1)

;sys=0 ;for windows system
sys=1 ;for linux   system

if sys eq 0 then FILEt=FILE_SEARCH('I:\muserdata_20170906\MUSERraw\','*',count=count)
if sys eq 1 then FILEt=FILE_SEARCH('/Volumes/Seagate_cxy/muserdata_20170906/MUSERraw/','*',count=count)

nfile=count
nfile=1

FOR ii=0,nfile-1 DO BEGIN ;count-1

  tmp0=rstrpos(FILEt[ii],'CSRH_')
  file_temp=strmid(FILEt[ii],tmp0,40);'CSRH_20141217-123418_13352702'
  print,file_temp
  
  ;-------------------------adding-------------------------

if sys eq 0 then FILE=FILE_SEARCH('I:\muserdata_20170906\MUSERraw\',string(file_temp))
if sys eq 1 then FILE=FILE_SEARCH('/Volumes/Seagate_cxy/muserdata_20170906/MUSERraw/',string(file_temp))

FILE=FILE(SORT(FILE))
N_FILE=N_ELEMENTS(FILE)

N_FILE=1

NUM=2400L
gpst=bytarr(8)
data_2ch=bytarr(12)
pwr_4ant=bytarr(32)
bk0=bytarr(40)
bk1=bytarr(32)
gtime=INTarr(9)
par=bytarr(4)
par_bk=bytarr(128)
FREQ_CODE=BYTARR(1)

cor_data=lonarr(N_FILE*NUM,44,44,16,2)
pwr=lonarr(N_FILE*NUM,44,16)

INT_TIME=1

cor_data_INT=lonarr(44,44,16,2,INT_TIME)
pwr_INT=lonarr(44,16,INT_TIME)

GPS_TIME=INTARR(N_FILE*NUM,9)
PAR_DELAY=dblARR(N_FILE*NUM,44)

OFFSET=100000LL
OFFSET0=0ll

data_name=['0.4-0.8_L',$
           '0.4-0.8_R',$ 
           '0.8-1.2_L',$
           '0.8-1.2_R',$
           '1.2-1.6_L',$
           '1.2-1.6_R',$
           '1.6-2.0_L',$
           '1.6-2.0_R']
;-------------------------adding------------------------- 
data_name_sav=data_name

if num lt 100 then data_name_sav=data_name+replicate('num00',8)+replicate(strmid(NUM,10,2),8)
if num ge 100 and num lt 1000 then data_name_sav=data_name+replicate('num0',8)+replicate(strmid(NUM,10,2),8)
if num ge 1000 then data_name_sav=data_name+replicate('num0',8)+replicate(strmid(NUM,10,2),8)

if num lt 100 then data_name_add='_int01_num00'+replicate(strmid(NUM,10,2),8)+'_'+strmid(file_temp,9,9)
if num ge 100 and num lt 1000 then data_name_add='_int01_num0'+replicate(strmid(NUM,9,3),8)+'_'+strmid(file_temp,9,9)
if num ge 1000 then data_name_add='_int01_num'+replicate(strmid(NUM,8,4),8)+'_'+strmid(file_temp,9,9)

data_name=data_name+data_name_add

;-------------------------adding------------------------- 

for aa=0,7 do begin
file_mkdir,data_name_sav[aa]


OPENw,LUN2,'log_'+data_name[aa]+'.txt',/get_lun
AFF=0
FOR N_FF=0,N_FILE-1 DO BEGIN
OPENR, 1, FILE[N_FF+AFF]
   N_F=N_FF
   FOR KK=0,NUM-1 DO BEGIN                                                      ;READ 3MS EVERY 6 SECOND
 
 ;read title
    
;   IF N_F EQ 1 AND KK EQ 0 THEN BEGIN
;
;   readu,1,par                                                                ;32bits are for switch of center frequency
;    if par[0] eq 51 then print, 'Open'
;    if par[0] eq 204 then print,'Close'
;
;   readu,1,par_bk                                                             ;64*16bits are for backup
;
;   readu,1,par                                                                ;32bits are for bandwidth of signal
;    if par[0] eq 51 then print,'25MHz'
;    if par[0] eq 119 then print,'12.5MHz'
;    if par[0] eq 136 then print,'6.25MHz'
;    if par[0] eq 187 then print,'3.125MHz'
;    if par[0] eq 204 then print,'1.5625MHz'
;
;   readu,1,par                                                                ;32 bits are for quantanzition
;    if par[0] eq 34 then print,'2'
;    if par[0] eq 51 then print,'3'
;    if par[0] eq 68 then print,'4'
;    if par[2] eq 51 then print,'P_V'
;    if par[2] eq 204 then print,'M_V'
;
;   readu,1,par                                                                ;32 bits are for switch of delay
;    if par[0] eq 51 then print,'Delay On'
;    if par[0] eq 204 then print,'Delay Off'
;
;   readu,1,par                                                                ;32 bits are for switch of fringe
;    if par[0] eq 51 then print,'Fringe On'
;    if par[0] eq 204 then print,'Fringe Off'
;
;   readu,1,par                                                                ;32 bits are to display polarization and frequency band
;    if par[0] eq 51 then print,'Left Pol'
;    if par[0] eq 204 then print,'Right Pol'
;    if par[2] eq 51 then print,'0.4-0.8GHz'
;    if par[2] eq 119 then print,'0.8-1.2GHz'
;    if par[2] eq 187 then print,'1.2-1.6GHz'
;    if par[2] eq 204 then print,'1.6-2.0GHz'
;
;    ENDIF

bb=0
freq_code=100
while (freq_code ne AA ) do begin
  point_lun,1,208LL+KK*320LL*100000LL*60/NUM+OFFSET0+OFFSET*bb
  READU,1,FREQ_CODE
   bb=bb+1
endwhile

   point_lun,1,32+KK*320LL*100000LL*60/NUM+OFFSET0+OFFSET*(bb-1)                       ;read GPS time,100000L:DATA LENGTH OF EACH 3MS, 320LL:there are 320 sector in 1 second data.
   readu,1, gpst
    readtime,(gpst),gpstime
    GPS_TIME[N_F*NUM+KK,*]=GPSTIME
    printf,lun2,'Time',gpstime
    printf,lun2,data_name[aa],freq_code,bb-1

   FOR INT=0,INT_TIME-1 DO BEGIN

   point_lun,1,2944lL+KK*100000LL*320LL*60/NUM+OFFSET0+OFFSET*INT*8+OFFSET*(bb-1)

 for ch=0,7 do begin                                                         ;there are 16chs in each antenna,i,j is the number of antenna
   FOR I=0,43 DO BEGIN
     FOR J=I+1,43 DO BEGIN
       READU,1,DATA_2CH
       A=DATA_2CH[11]*256L*256L+DATA_2CH[7]*256L+DATA_2CH[3]
       B=DATA_2CH[10]*256L*256L+DATA_2CH[6]*256L+DATA_2CH[2]
       C=DATA_2CH[9]*256L*256L+DATA_2CH[5]*256L+DATA_2CH[1]
       D=DATA_2CH[8]*256L*256L+DATA_2CH[4]*256L+DATA_2CH[0]
       SIGNED,A
       SIGNED,B
       SIGNED,C
       SIGNED,D
       COR_DATA_INT[I,J,0+CH*2,0,INT]=A
       COR_DATA_INT[I,J,0+CH*2,1,INT]=B
       COR_DATA_INT[I,J,1+CH*2,0,INT]=C
       COR_DATA_INT[I,J,1+CH*2,1,INT]=D
     ENDFOR
   ENDFOR

   readu,1,bk0
   FOR I=0,10 DO BEGIN
     READU,1,PWR_4ANT
     PWR_INT[0+I*4,0+CH*2,INT]=PWR_4ANT[3]*256L*256L*256L+PWR_4ANT[2]*256L*256L+PWR_4ANT[1]*256L+PWR_4ANT[0]
     PWR_INT[0+I*4,1+CH*2,INT]=PWR_4ANT[7]*256L*256L*256L+PWR_4ANT[6]*256L*256L+PWR_4ANT[5]*256L+PWR_4ANT[4]

     PWR_INT[1+I*4,0+CH*2,INT]=PWR_4ANT[11]*256L*256L*256L+PWR_4ANT[10]*256L*256L+PWR_4ANT[9]*256L+PWR_4ANT[8]
     PWR_INT[1+I*4,1+CH*2,INT]=PWR_4ANT[15]*256L*256L*256L+PWR_4ANT[14]*256L*256L+PWR_4ANT[13]*256L+PWR_4ANT[12]

     PWR_INT[2+I*4,0+CH*2,INT]=PWR_4ANT[19]*256L*256L*256L+PWR_4ANT[18]*256L*256L+PWR_4ANT[17]*256L+PWR_4ANT[16]
     PWR_INT[2+I*4,1+CH*2,INT]=PWR_4ANT[23]*256L*256L*256L+PWR_4ANT[22]*256L*256L+PWR_4ANT[21]*256L+PWR_4ANT[20]

     PWR_INT[3+I*4,0+CH*2,INT]=PWR_4ANT[27]*256L*256L*256L+PWR_4ANT[26]*256L*256L+PWR_4ANT[25]*256L+PWR_4ANT[24]
     PWR_INT[3+I*4,1+CH*2,INT]=PWR_4ANT[31]*256L*256L*256L+PWR_4ANT[30]*256L*256L+PWR_4ANT[29]*256L+PWR_4ANT[28]

     READU,1,BK1
     IF CH EQ 7 THEN BEGIN
      PAR_DELAY[N_F*NUM+KK,0+I*4]=BK1[17]*256+BK1[16]+(BK1[9]*256+BK1[8])/10000.
      PAR_DELAY[N_F*NUM+KK,1+I*4]=BK1[19]*256+BK1[18]+(BK1[11]*256+BK1[10])/10000.
      PAR_DELAY[N_F*NUM+KK,2+I*4]=BK1[21]*256+BK1[20]+(BK1[13]*256+BK1[12])/10000.
      PAR_DELAY[N_F*NUM+KK,3+I*4]=BK1[23]*256+BK1[22]+(BK1[15]*256+BK1[14])/10000.
     endif

   ENDFOR
    readu,1,bk1
 endfor  ;END OF CH

ENDFOR  ;END OF INT

;===========================
FOR I=0,43 DO BEGIN
 FOR K=0,15 DO BEGIN
  FOR J=0,43 DO BEGIN
   FOR L=0,1 DO BEGIN
    COR_DATA[N_F*NUM+KK,I,J,K,L]=MEAN(COR_DATA_INT[I,J,K,L,*])
    ENDFOR
  ENDFOR
  PWR[N_F*NUM+KK,I,k]=MEAN(PWR_INT[I,k,*])
 ENDFOR
ENDFOR
;===============================

ENDFOR ;END OF KK,LOOP
FREE_LUN,1
ENDFOR   ; end of nf

free_lun,lun2

SAVE,COR_DATA,PWR,GPS_TIME,PAR_DELAY,FILENAME='COR_data_'+data_name[aa]+'.sav'

file_move,['log_'+data_name[aa]+'.txt','COR_data_'+data_name[aa]+'.sav'],data_name_sav[aa]
;file_move,'COR_data_'+data_name[aa]+'.sav'],data_name_sav[aa]

endfor

;-------------------------adding-------------------------
ENDFOR
print,'fine'
print, 'File written in '+strtrim(SYSTIME(1) - systim,1)+' seconds'
END

;=========================================================================
