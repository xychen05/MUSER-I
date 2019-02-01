;.r main.pro
;.compile package_muser.pro
; 
; Name :
; main
; 
; Purpose :
; Display dirty images.
; 
; Explanation :
; This procedure is for plotting dirty images from original data without any processing.
; For the prep data.
; You should give a flag table.
; 
; History :
; Wei Wang, MAIN.PRO
; Xingyao Chen, June 2017
; 
;-------------------------begins-------------------------
Device, RETAIN=2
!P.FONT = 0
!P.Color=0
!P.Background=255

WINDOW,XS=512,YS=512

pos1 = [0.1,0.1,0.9,0.9]
posc1 = [pos1[2]-0.03,pos1[3]-0.3,pos1[2],pos1[3]-0.02]

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

;;;-------------------------file name for windows-------------------------
if sys eq 0 then filename=file_search('.\'+data_namet+'\prep\*.sav',count=count) & $
  image_sav='.\imagetest\'+data_namet+'\prep_dirty_img\' & data_sav='.\sav\'+data_namet+'\prep_img_beam\'

;;;-------------------------file name for linux-------------------------
if sys eq 1 then filename=file_search('./'+data_namet+'/prep/*.sav',count=count) & $
  image_sav='./imagetest/'+data_namet+'/prep_dirty_img/' & data_sav='./sav/'+data_namet+'/prep_img_beam/'

file_mkdir,image_sav
file_mkdir,data_sav

print,count
nfile=count
nfile=1
num=2400

;-------------------------loop 1-------------------------
for ii=140,140 do begin ; for which file ;;nfile-1
  
  tmp0=rstrpos(filename[ii],'cor_data_')
  file=strmid(filename[ii],tmp0,60)
  print,file
  restore,filename[ii],/ver
  
  ;FLAG_TABLE=[10,11,12,13,14,15,16,17,18,19,25,26,39]
  ;FLAG_TABLE=[10,11,12,13,15,16,17,18,19,24,25,26,27,38,39]   ;flag_table 1
  FLAG_TABLE=[10,11,12,13,14,15,16,17,18,19,24,25,26,27,38,39];flag_table 2
  ;FLAG_TABLE=[10,11,12,13,14,15,16,17,18,19,24,25,26,38,39]  ;flag_table 3
  ;FLAG_TABLE=[10,11,12,13,14,15,16,17,18,19,25,26,38,39]  ;flag_table 4

  SUN_ANGLE=22.47/180*!PI
  nsize=size(gps_time_prep)
  cadence=1
  n=nsize[1]/cadence
  
  cor_data1=dblarr(40,40,2)
  LOADCT,39;13
;-------------------------loop 2-------------------------  
  FOR tt=1,5 DO BEGIN ;;for which frame
    
    I=tt*10
    DAY=gps_time_prep[0,2]
    BJT=gps_time_prep[I*CADENCE,3]+gps_time_prep[I*CADENCE,4]/60.+gps_time_prep[I*CADENCE,5]/3600.
    CAL_UVW,DAY,BJT,U,V,W;print,bjt,h,ra,dec
    
    cor_data1[*,*,*]=cor_data_prep[I,*,*,*]
    
    FLAG_DATA,FLAG_TABLE,cor_data1,VIS,U,V
    ;mov_pos,vis,u,v,-(319-512/2)/512.,(357-512/2)/512. ;for Left polarization
    ;mov_pos,vis,u,v,-(303-512/2)/512.,(362-512/2)/512. ;for Right polarization 
    
    ;SUN_MODEL_DFT,U,V,VIS_MODEL
    ; PLOT,U,V,PSYM=4,XR=[-10000,10000],YR=[-10000,10000]
    ; ;OPLOT,-U,-V,PSYM=4,COLOR=255
    ; WRITE_JPEG,'UV'+STRCOMPRESS(I,/REMOVE_ALL)+'.JPG',TVRD(TRUE=3),TRUE=3
    ; PLOT,SQRT(U^2.+V^2.),ABS(VIS),PSYM=4
    ; WRITE_JPEG,'AMP_UV'+STRCOMPRESS(I,/REMOVE_ALL)+'.JPG',TVRD(TRUE=3),TRUE=3
    FFT_IMAGE,U,V,VIS,BEAM,IMG
    ;IMG=ROT(IMG,SUN_ANGLE)
;-------------------------basic para-------------------------    
    plot_image,abs(img),XTICKLEN=-0.01,YTICKLEN=-0.01,title='MUSER @1.725GHz '+ $
      strmid(gps_time_prep[I*CADENCE,3],6,6)+':'+strmid(gps_time_prep[I*CADENCE,4],6,6)+':'+ $
      strmid(gps_time_prep[I*CADENCE,5],6,6)+'.'+strmid(gps_time_prep[I*CADENCE,6],5,6),charsize=chars1,position=pos1, $
      drange=[min(abs(img)),max(abs(img))-3]
    cgcolorbar,range=[min(abs(img)),max(abs(img))],yticks = 2.,/ver,charsize=chars1,position=posc1,color=255
    
    write_jpeg,image_sav+'dirtyimg'+strmid(file,22,24)+'_'+strcompress(i,/remove_all)+'.jpg',tvrd(true=3),true=3
    
  endfor
  
  save,img,beam,filename=data_sav+'img_beam_'+file
  
endfor
;-------------------------ending-------------------------
print,'dirty image has prepared...',image_sav
print,'img data has prepared...',data_sav
print,'OK'
end