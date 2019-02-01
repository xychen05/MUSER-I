;.r main_mv_clean.pro
;.compile package_muser.pro
;
; Name :
; main_mv_clean
; 
; Purpose :
; Display dirty images after moving to the image center and clean map.
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
; 
;-------------------------begins-------------------------
Device, RETAIN=2
!P.FONT = 0
LOADCT,39

WINDOW,1,XS=512,YS=512
WINDOW,5,XS=512,YS=512

WINDOW,2,XS=700*3/2.*1.0,YS=700*1.0,title='main_mv_clean'
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

chars1=1.2

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

nm_ind='17'

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

;;;for it data, after integration
;if sys eq 0 then filename=file_search('.\sav\'+data_namet+'\cor_data_it_quiet.sav',count=count) & $
;  image_sav='.\imagetest\'+data_namet+'\it_clean_img00\' & data_sav='.\sav\'+data_namet+'\it_img_beam00\'
;if sys eq 1 then filename=file_search('./sav/'+data_namet+'/cor_data_it_quiet.sav',count=count) & $
;  image_sav='./imagetest/'+data_namet+'/it_clean_img00/' & data_sav='./sav/'+data_namet+'/it_img_beam00/'

;;for it data, after integration
if sys eq 0 then filename=file_search('.\sav\'+data_namet+'\it_cor_data'+nm_ind+'\cor_data_it_'+'*.sav',count=count) & $
  image_sav='.\imagetest\'+data_namet+'\it_clean_img'+nm_ind+'\' & data_sav='.\sav\'+data_namet+'\it_clean_img'+nm_ind+'\'
if sys eq 1 then filename=file_search('./sav/'+data_namet+'/it_cor_data'+nm_ind+'/cor_data_it_'+'*.sav',count=count) & $
  image_sav='./imagetest/'+data_namet+'/it_clean_img'+nm_ind+'/' & data_sav='./sav/'+data_namet+'/it_clean_img'+nm_ind+'/'

;-------------------------basic para-------------------------
print,count

if sys eq 0 then filename_kp=file_search('.\sav\'+data_namet+'\keep_table'+nm_ind+'.sav')
if sys eq 1 then filename_kp=file_search('./sav/'+data_namet+'/keep_table17.sav');;'+nm_ind+'
help,keep_table

restore,filename_kp[0]

file_mkdir,image_sav
file_mkdir,data_sav

PIX=512
num=60

for ii=15,15 do begin ;count-1 ; for which file
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
  restore,filename[ii],/ver
  tmp0=rstrpos(filename[ii],'cor_data_')
  file=strmid(filename[ii],tmp0,80)
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
  
  img_mv=make_array((i_ed-i_st+1),PIX,PIX)
  beam_mv=make_array((i_ed-i_st+1),PIX*2,PIX*2)
  LOADCT,39
;-------------------------basic para-------------------------
  FOR I=i_st,i_ed DO BEGIN

    DAY=gps_time_rm[0,2]
    BJT=gps_time_rm[I*CADENCE,3]+gps_time_rm[I*CADENCE,4]/60.+gps_time_rm[I*CADENCE,5]/3600.
    CAL_UVW,DAY,BJT,U,V,W;print,bjt,h,ra,dec

    cor_data1[*,*,*]=cor_data[I,*,*,*]

    ;FLAG_DATA,FLAG_TABLE,cor_data1,VIS,U,V
    flag_baseline,keep_table,cor_data1,vis,u,v
    
    wset,1
    plot,[0,0],xrange=[-17000.,17000],yrange=[-17000.,17000],xtitle='U',ytitle='V',position=[0.15,0.15,0.95,0.95],xstyle=1,ystyle=1
    for xy1=0,39 do begin
      for xy2=0,39 do begin
        plots,u[xy1,xy2],v[xy1,xy2],thick=1.,psym=2,color=255,symsize=1
      endfor
    endfor
    
    ;mov_pos,vis,u,v,-(319-512/2)/512.,(357-512/2)/512. ;for Left polarization
    ;mov_pos,vis,u,v,-(303-512/2)/512.,(362-512/2)/512. ;for Right polarization
    ;mov_pos,vis,u,v,-(319+25-512/2)/512.,(357+25-512/2)/512. ;for Left polarization

    FFT_IMAGE,U,V,VIS,BEAM,IMG
    
    if i eq i_st then time_string=time_gps_to_string(gps_time_obs[I,*])
    time_obs=time_gps_to_obs(gps_time_obs[I,*])
    inputimg=abs(img)
    ;;center_fit,inputimg,centx,centy
    ;;add_disk,inputimg,0.1,[centx,centy],solar_disk
    ;;print,'center:',centx,centy
    
    wset,5
    plot_image,abs(img),XTICKLEN=-0.01,YTICKLEN=-0.01,min=min(abs(img)),max=max(abs(img)), $
      charsize=chars1,position=[0.15,0.15,0.95,0.95], $
      title='MUSER dirty map '+time_obs;+' ['+strcompress(centx,/remove_all)+','+strcompress(centy,/remove_all)+']'
    ;contour,solar_disk,levels=[0.5]*max(solar_disk),/overplot,color=255
    ;plots,[centx,centx],[centy,centy],psym=1,symsize=3,color=255
    ;cgcolorbar,range=[min(abs(img)),max(abs(img))],yticks = 2.,/ver,charsize=chars1,position=posc1,color=255

    ;mov_pos,vis,u,v,-(280-512/2)/512.,(380-512/2)/512. ;for 1.7GHz
    
    FFT_IMAGE,U,V,VIS,BEAM,IMG
    
    ;-------------------------clean beginning-------------------------
    
    img_org=img
    ;;;for abs method
    BM=abs(beam)
    IMG=abs(IMG)
    
;    ;;for zero method
;    bm=beam
;    ind_temp=where(bm lt 0)
;    bm[ind_temp]=0
;    bm=abs(bm)
;    img=img
;    ind_temp=where(img lt 0)
;    img[ind_temp]=0
;    img=abs(img)
;    beam=abs(beam)
;    img=abs(img)
;    BM=beam
;    
    ;;for -min(img) method
    
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
    FOR II1=0,PIX-1 DO BEGIN
      FOR JJ1=0,PIX-1 DO BEGIN
        UU[II1,JJ1]=(XX[II1]/A[2])^2+(YY[JJ1]/A[3])^2
      ENDFOR
    ENDFOR

    C_BM=A[1]*EXP(-UU/2);+A[0]
    MODEL=FLTARR(PIX,PIX)
    
    gps_time_cl=gps_time_mv[i,*]
    time_obs=time_gps_to_string(gps_time_cl)
    
    ;-------------------------clean to get model, residual, clean map-------------------------
    TIMES=100;200;set the times for cleaning
    GAIN=0.1 ;always keep 0.1

    ratio_sl=1.3 ;select ratio_sl*r_sun region
    center_sl=[319+25-56,357];[256,256] ;need to set a center first
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
    region_sl[(center_sl[0]-155):(center_sl[0]+155),(center_sl[1]-155):(center_sl[1]+154)]=1
    
;    ;for radio source
;    img_rgsl=make_array(pix,pix)
;    img_rgsl[(center_sl[0]-155):(center_sl[0]+155),(center_sl[1]-155):(center_sl[1]+155)]=img[(center_sl[0]-155):(center_sl[0]+155),(center_sl[1]-155):(center_sl[1]+155)]
;    ind_rgsl=where(img_rgsl gt 0.5*max(img_rgsl))
;    help,ind_rgsl
;    ind_rgsl_x=ind_rgsl / pix
;    ind_rgsl_y=ind_rgsl mod pix
;    for vv=0,(size(ind_rgsl_x))[1]-1 do begin
;      region_sl[ind_rgsl_x[vv],ind_rgsl_y[vv]]=1
;    endfor
;    
    region_sl_index=where(region_sl eq 1)
    img_sl=make_array(pix,pix)
    img_sl[region_sl_index]=img[region_sl_index]

    iteration=0
    while (max(img_sl) gt ratio_rms*stddev(img_sl)) do begin
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
      
      ;;============
      iteration=iteration+1
      print,ii,iteration
      if iteration gt times then begin
        goto,nafik
      endif

    endwhile
    nafik:
    ;-------------------------ploting-------------------------
    zz=i
    gps_time_cl=gps_time_mv[zz,*]
    gps_timett=gps_time_cl
    gps_timett[*,3]=gps_time_cl[*,3];-8
    time_obs=time_gps_to_string(gps_timett)
    centx=pix/2 & centy=pix/2
    
    wset,2
    ;;============
    plot_image,abs(beam),title='d_beam '+time_obs+' times:'+strcompress(iteration,/remove_all) $
      ,charsize=chars1,position=pos1
    ;;============
    plot_image,abs(c_bm),title='c_beam',charsize=chars1,position=pos2,/noerase
    ;;============
    inputimg=abs(dimg)
    ;center_fit,inputimg,centx,centy
    inputimg=abs(img)
    add_disk,inputimg,0.01,[centx,centy],solar_disk
    plot_image,abs(dimg),title='dirty img '+'['+strcompress(centx,/remove_all)+','+ $
      strcompress(centy,/remove_all)+']',charsize=chars1,position=pos3,/noerase
    contour,solar_disk,levels=[0.5]*max(solar_disk),/overplot,color=255
    plots,[centx,centx],[centy,centy],psym=1,symsize=3,color=255
    cgcolorbar,range=[min(abs(dimg)),max(abs(dimg))],yticks = 2.,/ver,charsize=chars1,position=posc3,color=255
    ;;============
    plot_image,abs(model),title='model',charsize=chars1,position=pos4,/noerase
    ;;============
    inputimg=abs(img)
    center_fit,inputimg,centx,centy
    inputimg=abs(img)
    add_disk,inputimg,0.01,[centx,centy],solar_disk
    plot_image,abs(img),title='residual img '+'['+strcompress(centx,/remove_all)+','+ $
      strcompress(centy,/remove_all)+']',charsize=chars1,position=pos5,/noerase
    contour,solar_disk,levels=[0.5]*max(solar_disk),/overplot,color=255
    plots,[centx,centx],[centy,centy],psym=1,symsize=3,color=255
    cgcolorbar,range=[min(abs(img)),max(abs(img))],yticks = 2.,/ver,charsize=chars1,position=posc5,color=255
    ;;============
    inputimg=abs(c_img)
    add_disk,inputimg,0.01,[centx,centy],solar_disk
    plot_image,abs(c_img),title='clean img '+'['+strcompress(centx,/remove_all)+','+ $
      strcompress(centy,/remove_all)+']',charsize=chars1,position=pos6,/noerase
    contour,solar_disk,levels=[0.5]*max(solar_disk),/overplot,color=255
    plots,[centx,centx],[centy,centy],psym=1,symsize=3,color=255
    cgcolorbar,range=[min(abs(c_img)),max(abs(c_img))],yticks = 2.,/ver,charsize=chars1,position=posc6,color=255
    
    write_jpeg,image_sav+'clean_'+nm_ind+'_'+time_obs+'_times'+strcompress(iteration,/remove_all)+'.jpg',tvrd(true=3),true=3
    
    print,'dimg',minmax(abs(dimg))
    print,'img',minmax(abs(img))
    print,'c_img',minmax(abs(c_img))
    
    ;write_jpeg,image_sav+'clean_'+file+'_times'+strcompress(iteration,/remove_all)+'.jpg',tvrd(true=3),true=3
    img_cl=c_img
    beam_cl=c_bm
    ch=4
    ;save,dimg,img,img_cl,gps_time_cl,ch,filename=data_sav+'clean_'+time_obs+'_times'+strcompress(iteration,/remove_all)+'.sav'
    save,dimg,img,img_cl,gps_time_cl,ch,beam,beam_cl,filename=data_sav+'clean_'+nm_ind+'_'+time_obs+'_times'+strcompress(iteration,/remove_all)+'.sav'

  endfor
endfor

;-------------------------ending-------------------------
print,'OK'
end