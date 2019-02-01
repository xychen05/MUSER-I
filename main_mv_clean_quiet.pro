;.r main_mv_clean_quiet.pro
;.compile package_muser.pro
;
; Name :
; main_mv_clean_quiet
; 
; Purpose :
; Display dirty images after moving to the image center and clean map.
; For quiet sun.
; 
; Explanation :
; This procedure is for plotting dirty images after removing bad frames and moving to the image center.
; The rm data has prepared 
; The left dirty image is the single frame data while the right dirty is for the intergration of the data you select.
; Both of the dirty has overlapped with the fitted solar disk, which is used the ideal solar disk model
;
; History :
; Wei Wang, MAIN.PRO & CLEAN.PRO
; Xingyao Chen, June 2017
; Xingyao Chen, 24 Nov 2017
; 
;-------------------------begins-------------------------
;Device, RETAIN=2
LOADCT,39
;WINDOW,2,XS=700*3/2.*1.0,YS=700*1.0,title='main_mv_clean'
WINDOW,2,XS=700*3/2.*1.3,YS=700*1.3,title='main_mv_clean'
; x-0.3, y-0.45
pos1 = [0.03,0.52,0.33,0.97]
pos2 = [0.36,0.52,0.66,0.97]
pos3 = [0.69,0.52,0.99,0.97]

pos4 = [0.03,0.02,0.33,0.47]
pos5 = [0.36,0.02,0.66,0.47]
pos6 = [0.69,0.02,0.99,0.47]

posc1 = [pos1[2]-0.01,pos1[3]-0.15,pos1[2],pos1[3]-0.01]
posc2 = [pos2[2]-0.01,pos2[3]-0.15,pos2[2],pos2[3]-0.01]
posc3 = [pos3[2]-0.01,pos3[3]-0.15,pos3[2],pos3[3]-0.01]

posc4 = [pos4[2]-0.01,pos4[3]-0.15,pos4[2],pos4[3]-0.01]
posc5 = [pos5[2]-0.01,pos5[3]-0.15,pos5[2],pos5[3]-0.01]
posc6 = [pos6[2]-0.01,pos6[3]-0.15,pos6[2],pos6[3]-0.01]

chars1=1.2/1.

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

for ooo=4,4 do begin
data_namet=data_name[ooo]

sys=0 ;for windows system
;sys=1 ;for linux   system

;nm_ind='05'

;ch=0
for ch=0,15 do begin
nm_ind='sf'+strcompress(ch,/remove_all)

;;;-------------------------file name-------------------------
;;;for prep data, original data
;if sys eq 0 then filename=file_search('.\'+data_namet+'\prep\*.sav',count=count) & $
;  image_sav='.\imagetest\'+data_namet+'\prep_clean_img\' & data_sav='.\sav\'+data_namet+'\prep_img_beam\'
;if sys eq 1 then filename=file_search('./'+data_namet+'/prep/*.sav',count=count) & $
;  image_sav='./imagetest/'+data_namet+'/prep_clean_img/' & data_sav='./sav/'+data_namet+'/prep_img_beam/'
;
;;;for rm data, after removing bad frames
;if sys eq 0 then filename=file_search('.\sav\'+data_namet+'\rm_cor_data'+nm_ind+'\*.sav',count=count) & $
;  image_sav='.\imagetest\'+data_namet+'\rm_clean_img'+nm_ind+'\' & data_sav='.\sav\'+data_namet+'\rm_img_beam'+nm_ind+'\'
;if sys eq 1 then filename=file_search('./sav/'+data_namet+'/rm_cor_data'+nm_ind+'/*.sav',count=count) & $
;  image_sav='./imagetest/'+data_namet+'/rm_clean_img'+nm_ind+'/' & data_sav='./sav/'+data_namet+'/rm_img_beam'+nm_ind+'/'
;
;;;for mv data after removing bad frames and moving to the image center
;if sys eq 0 then filename=file_search('.\sav\'+data_namet+'\mv_cor_data'+nm_ind+'\*.sav',count=count) & $
;  image_sav='.\imagetest\'+data_namet+'\mv_clean_img'+nm_ind+'\' & data_sav='.\sav\'+data_namet+'\mv_img_beam'+nm_ind+'\'
;if sys eq 1 then filename=file_search('./sav/'+data_namet+'/mv_cor_data'+nm_ind+'/*.sav',count=count) & $
;  image_sav='./imagetest/'+data_namet+'/mv_clean_img'+nm_ind+'/' & data_sav='./sav/'+data_namet+'/mv_img_beam'+nm_ind+'/'
;
;;;for pm data, after doing the phase modulation
;if sys eq 0 then filename=file_search('.\sav\'+data_namet+'\pm_cor_data'+nm_ind+'\*.sav',count=count) & $
;  image_sav='.\imagetest\'+data_namet+'\pm_clean_img'+nm_ind+'\' & data_sav='.\sav\'+data_namet+'\pm_img_beam'+nm_ind+'\'
;if sys eq 1 then filename=file_search('./sav/'+data_namet+'/pm_cor_data'+nm_ind+'/*.sav',count=count) & $
;  image_sav='./imagetest/'+data_namet+'/pm_clean_img'+nm_ind+'/' & data_sav='./sav/'+data_namet+'/pm_img_beam'+nm_ind+'/'

;;;for it data, after integration
;if sys eq 0 then filename=file_search('.\sav\'+data_namet+'\it_cor_data00\cor_data_it_'+'*.sav',count=count) & $
;  image_sav='.\imagetest\'+data_namet+'\it_clean_img00\' & data_sav='.\sav\'+data_namet+'\it_img_beam00\'
;if sys eq 1 then filename=file_search('./sav/'+data_namet+'/it_cor_data00/cor_data_it_'+'*.sav',count=count) & $
;  image_sav='./imagetest/'+data_namet+'/it_clean_img00/' & data_sav='./sav/'+data_namet+'/it_img_beam00/'

;;for it data, after integration
if sys eq 0 then filename=file_search('.\sav\'+data_namet+'\quiet\cor_data_it'+nm_ind+'_quiet.sav',count=count) & $
  image_sav='.\imagetest\'+data_namet+'\it_clean_img_quiet'+nm_ind+'\' & data_sav='.\sav\'+data_namet+'\quiet_clean_img\'
if sys eq 1 then filename=file_search('./sav/'+data_namet+'/quiet/cor_data_it'+nm_ind+'_quiet.sav',count=count) & $
  image_sav='./imagetest/'+data_namet+'/it_clean_img_quiet'+nm_ind+'/' & data_sav='./sav/'+data_namet+'/quiet_clean_img/'

;-------------------------basic para-------------------------
print,count

;if sys eq 0 then restore,'.\sav\keep_table.sav',/ver
;if sys eq 1 then restore,'./sav/keep_table.sav',/ver

if sys eq 0 then filename_kp=file_search('.\sav\'+data_namet+'\keep_table'+nm_ind+'.sav')
if sys eq 1 then filename_kp=file_search('./sav/'+data_namet+'/keep_table'+nm_ind+'.sav')

restore,filename_kp[0],/ver
help,keep_table

file_mkdir,image_sav
file_mkdir,data_sav

PIX=512
num=1

for ii=0,0 do begin ;count-1 ; for which file
;  ;;for prep data
;  tmp0=rstrpos(filename[ii],'cor_data_')
;  file=strmid(filename[ii],tmp0,80)
;  print,file
;  restore,filename[ii],/ver
;  cor_data_rm=cor_data_prep
;  gps_time_rm=gps_time_prep
;  cor_data=cor_data_prep

;  ;;for rm data
;  tmp0=rstrpos(filename[ii],'cor_data_')
;  file=strmid(filename[ii],tmp0,80)
;  print,file
;  restore,filename[ii],/ver
;  ;cor_data_rm, gps_time_rm, cor_data ;;cor_data_rm=cor_data
; 
;  ;;for mv data
;  tmp0=rstrpos(filename[ii],'cor_data_')
;  file=strmid(filename[ii],tmp0,80)
;  print,file
;  restore,filename[ii],/ver
;  cor_data_rm=cor_data_mv
;  gps_time_rm=gps_time_mv
;  
;  ;;;for pm data
;  tmp0=rstrpos(filename[ii],'cor_data_')
;  file=strmid(filename[ii],tmp0,80)
;  cor_data_rm=cor_data_pm[ii*500:((ii+1)*500-1),*,*,*]
;  gps_time_rm=gps_time_pm[ii*500:((ii+1)*500-1),*]
  
  ;;;for it data
  restore,filename[0],/ver
  tmp0=rstrpos(filename[0],'cor_data_')
  file=strmid(filename[0],tmp0,80)
  cor_data_rm=cor_data_it
  gps_time_rm=gps_time_it
  
  ;;;begins
  cor_data=cor_data_rm
  gps_time_mv=gps_time_rm
  gps_time_obs=gps_time_rm
  gps_time_obs[*,3]=gps_time_obs[*,3];-8
 
  ;FLAG_TABLE=[10,11,12,13,14,15,16,17,18,19,25,26,39]
  ;FLAG_TABLE=[10,11,12,13,15,16,17,18,19,24,25,26,27,38,39]   ;flag_table 1
  ;FLAG_TABLE=[10,11,12,13,14,15,16,17,18,19,24,25,26,27,38,39];flag_table 2
  ;FLAG_TABLE=[10,11,12,13,14,15,16,17,18,19,24,25,26,38,39]  ;flag_table 3
  ;FLAG_TABLE=[10,11,12,13,14,15,16,17,18,19,25,26,38,39]  ;flag_table 4
  ;FLAG_TABLE=[10,11,12,13,14,16,17,18,19,24,25,26,27,38,39]   ;flag_table 5

  SUN_ANGLE=22.47/180*!PI
  nsize=size(cor_data_rm)
  cadence=1
  n=nsize[1]/cadence
  print,'size(cor_data_rm)',n
  ;n=1
  
  cor_data1=dblarr(40,40,2)

  i_st=0 ;;which frame to begin in file
  i_ed=0 ;;n-1 ;which frame to end in file
  
  img_mv=complexarr((i_ed-i_st+1),PIX,PIX)
  beam_mv=complexarr((i_ed-i_st+1),PIX*2,PIX*2)

  LOADCT,39
;-------------------------basic para-------------------------
  FOR I=i_st,i_ed DO BEGIN

    DAY=gps_time_rm[0,2]
    BJT=gps_time_rm[I*CADENCE,3]+gps_time_rm[I*CADENCE,4]/60.+gps_time_rm[I*CADENCE,5]/3600.
    CAL_UVW,DAY,BJT,U,V,W;print,bjt,h,ra,dec

    cor_data1[*,*,*]=cor_data[I*CADENCE,*,*,*]

    ;FLAG_DATA,FLAG_TABLE,cor_data1,VIS,U,V
    flag_baseline,keep_table,cor_data1,vis,u,v
    
    ;mov_pos,vis,u,v,-(319-512/2)/512.,(357-512/2)/512. ;for Left polarization
    ;mov_pos,vis,u,v,-(303-512/2)/512.,(362-512/2)/512. ;for Right polarization 
    ;mov_pos,vis,u,v,-(319+25-512/2)/512.,(357+25-512/2)/512. ;for Left polarization
    ;mov_pos,vis,u,v,-(319-512/2)/512.,(357-512/2)/512. ;for Left polarization
    
    ;FFT_IMAGE,U,V,VIS,BEAM,IMG
    FFT_IMAGE,U,V,VIS,IMG
    FFT_BEAM,U,V,BEAM
    
    ;-------------------------clean beginning-------------------------
    beam=real_part(beam)
    img=real_part(img)
    BM=beam

    IMG=IMG
    BM=BM/MAX(BM)
    ;BM=REVERSE(BM,2)
    ;IMG=REVERSE(IMG,2)

    DIMG=IMG

    ;;-------------------------for cleaning dirty beam-------------------------
    bm_abs=abs(bm[(PIX-50):(PIX+50),(PIX-50):(PIX+50)])
    ind_bm=where(bm_abs ge 0.5*max(bm_abs))
    szbm=size(bm_abs);szbm[1],szbm[2]
    bm_sl=make_array(szbm[1],szbm[2])
    bm_sl[ind_bm]=bm_abs[ind_bm]
    bm_fit=gauss2dfit(bm_sl,A)
    print,a[2]*3600/512,' arcsec',a[3]*3600/512,' arcsec',a[0];width in the X and Y direction
    X=FINDGEN(PIX)-PIX/2
    Y=FINDGEN(PIX)-PIX/2
    UU=FLTARR(PIX,PIX)
    XX=X*COS(A[6])-Y*SIN(A[6])
    YY=X*SIN(A[6])-Y*COS(A[6])
    FOR II=0,PIX-1 DO BEGIN
      FOR JJ=0,PIX-1 DO BEGIN
        UU[II,JJ]=(XX[II]/A[2])^2+(YY[JJ]/A[3])^2
      ENDFOR
    ENDFOR

    C_BM=A[1]*EXP(-UU/2);+A[0]
    MODEL=FLTARR(PIX,PIX)
    
    gps_time_cl=gps_time_mv[i,*]
    time_obs=time_gps_to_string(gps_time_cl)
    
    ;-------------------------clean to get model, residual, clean map-------------------------
    TIMES=200;200;set the times for cleaning
    GAIN=0.1 ;always keep 0.1

    ratio_sl=1.5 ;select ratio_sl*r_sun region
    center_sl=[222,224] ;need to set a center first
    ratio_rms=3.0
    ;ratio_rms=2.7+ff*0.01 ;max(residual) gt ratio_rms*stddev(residual)

    r_sun=pix*979.684d/60d/60d ;in pixels
    region_sl=make_array(pix,pix)
;    ;for a circle
;    for xx=0,pix-1 do begin
;      for yy=0, pix-1 do begin
;        if sqrt((xx-center_sl[0])^2+(yy-center_sl[1])^2) le r_sun*ratio_sl then region_sl[xx,yy]=1
;      endfor
;    endfor
    ;for a square
    region_sl[(center_sl[0]-155):(center_sl[0]+155),(center_sl[1]-155):(center_sl[1]+155)]=1

    region_sl_index=where(region_sl eq 1)
    img_sl=make_array(pix,pix)
    img_sl[region_sl_index]=img[region_sl_index]

    iteration=0
    ;while (max(img_sl) gt ratio_rms*stddev(img_sl)) do begin
    while iteration lt 200 do begin
      img_sl=make_array(pix,pix)
      img_sl[region_sl_index]=img[region_sl_index]

      loc=where(img_sl eq max(img_sl))
      IF (N_ELEMENTS(LOC) GT 1) THEN BEGIN
        LOCXY=LOC[0]
      ENDIF ELSE BEGIN
        LOCXY=LOC
      ENDELSE

      LOCY=LOCXY/512
      LOCX=LOCXY MOD 512
      GAIN=MAX(IMG)*0.1

      MODEL[LOCX,LOCY]=MODEL[LOCX,LOCY]+GAIN
      IMG=IMG-BM[255-LOCX+256:255+511-LOCX+256,255-LOCY+256:255+511-LOCY+256]*GAIN
      C_IMG=CONVOL(MODEL,C_BM[256-30:256+30,256-30:256+30]);+abs(IMG)
      
;      ;;============
;      plot_image,beam,title='d_beam '+time_obs+' times:'+strcompress(iteration,/remove_all) $
;        ,charsize=chars1,position=pos1
;      ;;============
;      plot_image,c_bm,title='c_beam'+nm_ind,charsize=chars1,position=pos2,/noerase
;      ;;============
;      inputimg=dimg
;      plot_image,dimg,title='dirty img ',charsize=chars1,position=pos3,/noerase
;      cgcolorbar,range=[min(dimg),max(dimg)],yticks = 2.,/ver,charsize=chars1,position=posc3,color=255
;      ;;============
;      plot_image,model,title='model',charsize=chars1,position=pos4,/noerase
;      ;;============
;      inputimg=img
;      center_fit,inputimg,centx,centy
;      inputimg=img
;      add_disk,inputimg,0.01,[centx,centy],solar_disk
;      plot_image,img,title='residual img '+'['+strcompress(centx,/remove_all)+','+ $
;        strcompress(centy,/remove_all)+']',charsize=chars1,position=pos5,/noerase
;      contour,solar_disk,levels=[0.5]*max(solar_disk),/overplot,color=255
;      plots,[centx,centx],[centy,centy],psym=1,symsize=3,color=255
;      cgcolorbar,range=[min(img),max(img)],yticks = 2.,/ver,charsize=chars1,position=posc5,color=255
;      ;;============
;      inputimg=c_img
;      add_disk,inputimg,0.01,[centx,centy],solar_disk
;      plot_image,c_img,title='clean img '+'['+strcompress(centx,/remove_all)+','+ $
;        strcompress(centy,/remove_all)+']',charsize=chars1,position=pos6,/noerase
;      contour,solar_disk,levels=[0.5]*max(solar_disk),/overplot,color=255
;      plots,[centx,centx],[centy,centy],psym=1,symsize=3,color=255
;      cgcolorbar,range=[min(c_img),max(c_img)],yticks = 2.,/ver,charsize=chars1,position=posc6,color=255
;      
;      ;;============
      iteration=iteration+1
      print,iteration
      if iteration gt times then begin
        goto,nafik
      endif
 
      ;WRITE_JPEG,image_sav+'clean_img_mv_'+file+'times'+strcompress(iteration,/remove_all)+'_01.jpg',tvrd(true=3),true=3

    endwhile
    nafik:
      
      ;;============
      plot_image,beam,title='d_beam '+time_obs+' times:'+strcompress(iteration,/remove_all) $
        ,charsize=chars1,position=pos1
      ;;============
      plot_image,c_bm,title='c_beam'+nm_ind,charsize=chars1,position=pos2,/noerase
      ;;============
      inputimg=dimg
      plot_image,dimg,title='dirty img ',charsize=chars1,position=pos3,/noerase
      cgcolorbar,range=[min(dimg),max(dimg)],yticks = 2.,/ver,charsize=chars1,position=posc3,color=255
      ;;============
      plot_image,model,title='model',charsize=chars1,position=pos4,/noerase
      ;;============
      inputimg=img
      center_fit,inputimg,centx,centy
      inputimg=img
      add_disk,inputimg,0.01,[centx,centy],solar_disk
      plot_image,img,title='residual img '+'['+strcompress(centx,/remove_all)+','+ $
        strcompress(centy,/remove_all)+']',charsize=chars1,position=pos5,/noerase
      contour,solar_disk,levels=[0.5]*max(solar_disk),/overplot,color=255
      plots,[centx,centx],[centy,centy],psym=1,symsize=3,color=255
      cgcolorbar,range=[min(img),max(img)],yticks = 2.,/ver,charsize=chars1,position=posc5,color=255
      ;;============
      inputimg=c_img
      add_disk,inputimg,0.01,[centx,centy],solar_disk
      plot_image,c_img,title='clean img '+'['+strcompress(centx,/remove_all)+','+ $
        strcompress(centy,/remove_all)+']',charsize=chars1,position=pos6,/noerase
      contour,solar_disk,levels=[0.5]*max(solar_disk),/overplot,color=255
      plots,[centx,centx],[centy,centy],psym=1,symsize=3,color=255
      cgcolorbar,range=[min(c_img),max(c_img)],yticks = 2.,/ver,charsize=chars1,position=posc6,color=255
      
      ;;============    
    WRITE_JPEG,image_sav+'clean_img_mv_'+file+'times'+strcompress(iteration,/remove_all)+'_01.jpg',tvrd(true=3),true=3
    img_cl=c_img
    save,dimg,img,img_cl,gps_time_cl,ch,beam,beam_cl,filename=data_sav+'clean_'+nm_ind+'_'+file+'_times'+strcompress(iteration,/remove_all)+'.sav'

  endfor
endfor

;-------------------------ending-------------------------
endfor
endfor
print,'OK'
end