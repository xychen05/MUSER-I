;.r clean.pro
;.compile package_muser.pro
;PRO CLEAN
;
; Name :
; clean
;
; Purpose :
; Clean the dirty image of the img_beam data for single frame.
; Also fit the solar center using the ideal solar disk model.
; Display dirty beam, clean beam, dirty image, model, residual image and clean image.
; Select the maximum piont from the whole image.

; Inputs :
; dirty image and dirty beam.
;
; Outputs :
; cleaned image.
;
; Notes :
; Different iteration times give different cleaned images.
; Please avoid the over-cleaning.
;
; History :
; Wei Wang, CLEAN.PRO
; Xingyao Chen, June 2017
;
;-------------------------begins-------------------------
Device, RETAIN=2
LOADCT,39
!P.FONT = 0
!P.Color=0
!P.Background=255

WINDOW,2,XS=700*3/2.,YS=700,title='clean'
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

chars1=0.7

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

nm_ind='05'

;;;-------------------------file name-------------------------
;;;for prep data
;if sys eq 0 then filename=file_search('.\sav\'+data_namet+'\img_beam'+nm_ind+'\*.sav',count=count) & $
;  image_sav='.\imagetest\'+data_namet+'\prep_clean_img'+nm_ind+'\' & data_sav='.\sav\'+data_namet+'\prep_clean_img_beam'+nm_ind+'\'
;if sys eq 1 then filename=file_search('./sav/'+data_namet+'/img_beam'+nm_ind+'/*.sav',count=count) & $
;  image_sav='./imagetest/'+data_namet+'/prep_clean_img'+nm_ind+'/' & data_sav='./sav/'+data_namet+'/prep_clean_img_beam'+nm_ind+'/'
;
;;;for rm data
;if sys eq 0 then filename=file_search('.\sav\'+data_namet+'\rm_img_beam'+nm_ind+'\*.sav',count=count) & $
;  image_sav='.\imagetest\'+data_namet+'\rm_clean_img'+nm_ind+'\' & data_sav='.\sav\'+data_namet+'\rm_clean_img_beam'+nm_ind+'\'
;if sys eq 1 then filename=file_search('./sav/'+data_namet+'/rm_img_beam'+nm_ind+'/*.sav',count=count) & $
;  image_sav='./imagetest/'+data_namet+'/rm_clean_img'+nm_ind+'/' & data_sav='./sav/'+data_namet+'/rm_clean_img_beam'+nm_ind+'/'
;
;;;for mv data
;if sys eq 0 then filename=file_search('.\sav\'+data_namet+'\mv_img_beam'+nm_ind+'\*.sav',count=count) & $
;  image_sav='.\imagetest\'+data_namet+'\mv_clean_img'+nm_ind+'\' & data_sav='.\sav\'+data_namet+'\mv_clean_img_beam'+nm_ind+'\'
;if sys eq 1 then filename=file_search('./sav/'+data_namet+'/mv_img_beam'+nm_ind+'/*.sav',count=count) & $
;  image_sav='./imagetest/'+data_namet+'/mv_clean_img'+nm_ind+'/' & data_sav='./sav/'+data_namet+'/mv_clean_img_beam'+nm_ind+'/'
;
;;;for pm data
;if sys eq 0 then filename=file_search('.\sav\'+data_namet+'\pm_img_beam'+nm_ind+'\*.sav',count=count) & $
;  image_sav='.\imagetest\'+data_namet+'\pm_clean_img'+nm_ind+'\' & data_sav='.\sav\'+data_namet+'\pm_clean_img_beam'+nm_ind+'\'
;if sys eq 1 then filename=file_search('./sav/'+data_namet+'/pm_img_beam'+nm_ind+'/*.sav',count=count) & $
;  image_sav='./imagetest/'+data_namet+'/pm_clean_img'+nm_ind+'/' & data_sav='./sav/'+data_namet+'/pm_clean_img_beam'+nm_ind+'/'

;;for it data
if sys eq 0 then filename=file_search('.\sav\'+data_namet+'\it_img_beam00\*.sav',count=count) & $
  image_sav='.\imagetest\'+data_namet+'\it_clean_img00\' & data_sav='.\sav\'+data_namet+'\it_clean_img_beam00\'
if sys eq 1 then filename=file_search('./sav/'+data_namet+'/it_img_beam00/*.sav',count=count) & $
  image_sav='./imagetest/'+data_namet+'/it_clean_img00/' & data_sav='./sav/'+data_namet+'/it_clean_img_beam00/'

;-------------------------file ending-------------------------
print,count

file_mkdir,image_sav
file_mkdir,data_sav

pix=512

for ff=0,0 do begin ;count-1 ;for which file
;-------------------------for single frame of img_beam data beginning-------------------------  
  
  tmp0=rstrpos(filename[ff],'img_beam_')
  file=strmid(filename[ff],tmp0,80)
  print,file
  restore,filename[ff],/ver
  
  for jj=0,0 do begin ;;(size(img_mv))[1]-1
    
    zz=round(jj*(size(img_mv))[1]/60.)
    img=make_array(pix,pix)
    beam=make_array(pix*2,pix*2)
    img[*,*]=img_mv[zz,*,*]
    beam[*,*]=beam_mv[zz,*,*]
    
    ;;-------------------------clean begins-------------------------
    
    BM=beam
    PIX=512
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
    ;plot_image,bm_abs
    ;contour,bm_abs,levels = [0.3,0.5,0.7]*max(bm_abs),/overplot
    ;;;original method to clean beam
    ;CEN_PIX=10
    ;CEN_BM=BM[512-CEN_PIX/2:512+CEN_PIX/2,512-CEN_PIX/2:512+CEN_PIX/2]
    ;CLEAN_BM=GAUSS2DFIT(CEN_BM,A)
    ;PRINT,A[2]*3600/512,' Arcsec',A[3]*3600/512,' Arcsec',A[0]
    X=FINDGEN(PIX)-PIX/2
    Y=FINDGEN(PIX)-PIX/2
    UU=FLTARR(PIX,PIX)
    
    XX=X*COS(A[6])-Y*SIN(A[6])
    YY=X*SIN(A[6])-Y*COS(A[6])
    
    FOR I=0,PIX-1 DO BEGIN
      FOR J=0,PIX-1 DO BEGIN
        UU[I,J]=(XX[I]/A[2])^2+(YY[J]/A[3])^2
      ENDFOR
    ENDFOR
    C_BM=A[1]*EXP(-UU/2);+A[0]  
    MODEL=FLTARR(PIX,PIX)
    
    gps_time_cl=gps_time_mv[zz,*]
    time_obs=time_gps_to_string(gps_time_cl)
    
    ;-------------------------clean to get model, residual, clean map-------------------------
    TIMES=2 ;200 ;set the times for cleaning
    GAIN=0.1 ;always keep 0.1
    
    FOR I=0,TIMES DO BEGIN
    
      LOC=WHERE(IMG EQ MAX(IMG))
      IF (N_ELEMENTS(LOC) GT 1) THEN BEGIN
        LOCXY=LOC[0]
      ENDIF ELSE BEGIN
        LOCXY=LOC
      ENDELSE
      
      LOCY=LOCXY/512
      LOCX=LOCXY MOD 512
      GAIN=MAX(IMG)*0.1
      ;PRINT,[LOCX,LOCY],GAIN,MAX(IMG),MIN(IMG)
      
      MODEL[LOCX,LOCY]=MODEL[LOCX,LOCY]+GAIN
      IMG=IMG-BM[255-LOCX+256:255+511-LOCX+256,255-LOCY+256:255+511-LOCY+256]*GAIN
      C_IMG=CONVOL(MODEL,C_BM[256-30:256+30,256-30:256+30])+abs(IMG)
      
      ;-------------------------ploting-------------------------
      ;;============
      plot_image,abs(beam),title='d_beam '+time_obs+' times:'+strcompress(times,/remove_all) $
        ,charsize=chars1,position=pos1
      ;;============
      plot_image,abs(c_bm),title='c_beam',charsize=chars1,position=pos2,/noerase
      ;;============
      inputimg=abs(dimg)
      center_fit,inputimg,centx,centy
      inputimg=abs(dimg)
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
      center_fit,inputimg,centx,centy
      inputimg=abs(c_img)
      add_disk,inputimg,0.01,[centx,centy],solar_disk
      plot_image,abs(c_img),title='clean img '+'['+strcompress(centx,/remove_all)+','+ $
        strcompress(centy,/remove_all)+']',charsize=chars1,position=pos6,/noerase
      contour,solar_disk,levels=[0.5]*max(solar_disk),/overplot,color=255
      plots,[centx,centx],[centy,centy],psym=1,symsize=3,color=255
      cgcolorbar,range=[min(abs(c_img)),max(abs(c_img))],yticks = 2.,/ver,charsize=chars1,position=posc6,color=255
      ;;============

      print,'dimg',minmax(abs(dimg))
      print,'img',minmax(abs(img))
      print,'c_img',minmax(abs(c_img))
    ENDFOR
  endfor
  
  img_cl=c_img
  write_jpeg,image_sav+'clean_'+file+'_times'+strcompress(times,/remove_all)+'.jpg',tvrd(true=3),true=3
  save,dimg,img,img_cl,gps_time_cl,filename=data_sav+'clean_'+file+'_times'+strcompress(times,/remove_all)+'.sav'
  
endfor

  
PRINT,'It is ok............'
;DEVICE,/CLOSE
;SET_PLOT,'X'
END