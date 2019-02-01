;.r clean_int_data_prep_LR.pro
;.compile package_muser.pro
;
; Name :
; clean_int_data_prep_LR
;
; Purpose :
; Prepare the data.
;
; History :
; Xingyao Chen, 10 July 2017
; Xingyao Chen, 15 Aug 2017
;
;-------------------------begins-------------------------

;;-------------------------file name for windows-------------------------

;filename=file_search('I:\muserdata_20141217\muser-n\sa2\1.6-2.0_Rnum600\mv_img_beam_08\*.sav',count=count)
;file=strmid(filename[0],65,70)

;filename=file_search('I:\muserdata_20141217\muser-n\sa2\1.6-2.0_Rnum600\mv_int_img_beam_07\*.sav',count=count)
;file=strmid(filename[0],69,70)

;filename=file_search('I:\muserdata_20141217\muser-n\sa2\cleanLR_mv_1h_05.sav',count=count)
;file=strmid(filename[0],34,15)
;RESTORE,filename[0],/ver

;-------------------------file name for linux-------------------------
;
;filename=file_search('/Volumes/Seagate_cxy/muserdata_20141217/muser-n/sa2/1.6-2.0_Lnum600/rm_img_beam_02/*.sav',count=count)
;file=strmid(filename[0],65+18,53)
;print,count

;filename=file_search('/Volumes/Seagate_cxy/muserdata_20141217/muser-n/sa2/1.6-2.0_Lnum600/mv_int_img_beam_02/*.sav',count=count)
;file=strmid(filename[0],65+18,53)

;filename=file_search('/Volumes/Seagate_cxy/muserdata_20141217/muser-n/sa2/cleanR_mv_1h_02.sav',count=count)
;file=strmid(filename[0],34+18,15)
;RESTORE,filename[0],/ver

;;;-------------------------begins-------------------------
;PIX=512
;nfile=count
;;
;;;;;for int mv data
;img_1h=make_array(nfile,PIX,PIX)
;beam_1h=make_array(nfile,PIX*2,PIX*2)
;gps_time_1h=make_array(nfile,9)
;;
;for ii=0,nfile-1 do begin
; 
;;;  ;;;;-------------------------for non-int mv data-------------------------
;;;  restore,filename[ii],/ver
;;;  img_mv=abs(img_mv)
;;;  beam_mv=abs(beam_mv)
;;;  sz=(size(img_mv))[1]
;;;  
;;;  img_int=make_array(PIX,PIX)
;;;  for cc=0,PIX-1 do begin
;;;    for ee=0,PIX-1 do begin
;;;      img_int[cc,ee]=mean(img_mv[0:(sz-1),cc,ee])
;;;    endfor
;;;  endfor
;;;  print,ii,minmax(img_int)
;;;  print,ii,minmax(img_mv)
;;;
;;;  beam_int=make_array(PIX*2,PIX*2)
;;;  for cc=0,2*PIX-1 do begin
;;;    for ee=0,2*PIX-1 do begin
;;;      beam_int[cc,ee]=mean(beam_mv[0:(sz-1),cc,ee])
;;;    endfor
;;;  endfor
;;;  print,ii,minmax(beam_mv)  
;;;  print,ii,minmax(beam_int)  
;;;  
;;;  gps_time_int=gps_time_mv[sz/2,*]
;;;
;;;  save,img_int,beam_int,gps_time_int,filename='I:\muserdata_20141217\muser-n\sa2\1.6-2.0_Lnum600\img_beam_int_'+strmid(file,9,70)
;;
;;
; 
;;;;;-------------------------for int mv data-------------------------
;  restore,filename[ii],/ver
;  img_int=abs(img_int)
;  beam_int=abs(beam_int)
;  img_1h[ii,*,*]=img_int[*,*]
;  beam_1h[ii,*,*]=beam_int[*,*]  
;  gps_time_1h[ii,*]=gps_time_int[0,*]
;  print,'ii'
;  
;endfor

;;;-------------------------for int mv data-------------------------
;img_1h_mean=make_array(PIX,PIX)
;beam_1h_mean=make_array(PIX*2,PIX*2)
;gps_time_1h_mean=make_array(9)
;
;for xx=0,pix-1 do begin
;  for yy=0,pix-1 do begin
;    img_1h_mean[xx,yy]=mean(img_1h[*,xx,yy])
;  endfor
;endfor
;for xx=0,2*pix-1 do begin
;  for yy=0,2*pix-1 do begin
;    beam_1h_mean[xx,yy]=mean(beam_1h[*,xx,yy])
;  endfor
;endfor
;
;save,img_1h,beam_1h,img_1h_mean,beam_1h_mean,gps_time_1h,filename='I:\muserdata_20141217\muser-n\sa2\cleanR_mv_1h_07.sav'

filename=file_search('I:\muserdata_20141217\muser-n\sa2\cleanR_mv_1h_07.sav',count=count)
file=strmid(filename[0],34,15)
RESTORE,filename[0],/ver
a1=img_1h_mean & b1=beam_1h_mean
filename=file_search('I:\muserdata_20141217\muser-n\sa2\cleanL_mv_1h_07.sav',count=count)
file=strmid(filename[0],34,15)
RESTORE,filename[0],/ver
a2=img_1h_mean & b2=beam_1h_mean
img_1h_mean=(a1+a2)/2.
beam_1h_mean=(b1+b2)/2.
save,img_1h_mean,beam_1h_mean,gps_time_1h,filename='I:\muserdata_20141217\muser-n\sa2\cleanLR_mv_1h_07.sav'


PRINT,'It is ok............'
;DEVICE,/CLOSE
;SET_PLOT,'X'
END
