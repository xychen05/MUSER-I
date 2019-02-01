;.r main_mv_keep_table.pro
;.compile package_muser.pro
;
; Name :
; main_mv
; 
; Purpose :
; Display dirty images after moving to the image center.
; 
; Explanation :
; This procedure is for plotting dirty images after removing bad frames and moving to the image center.
; The rm data has prepared 
; The left dirty image is the single frame data while the right dirty is for the intergration of the data you select.
; Both of the dirty has overlapped with the fitted solar disk, which is used the ideal solar disk model
;
; History :
; Wei Wang, MAIN.PRO
; Xingyao Chen, June 2017
; 
;-------------------------begins-------------------------
Device, RETAIN=2
!P.FONT = 0

;WINDOW,1,XS=512,YS=512
WINDOW,5,XS=512,YS=512

pos1 = [0.07,0.1,0.47,0.9]
pos2 = [0.57,0.1,0.97,0.9]
posc1 = [pos1[2]-0.015,pos1[3]-0.25,pos1[2],pos1[3]-0.02]
posc2 = [pos2[2]-0.015,pos2[3]-0.25,pos2[2],pos2[3]-0.02]

chars1=0.9

pos1=[0.1,0.1,0.9,0.9]
posc1 = [pos1[2]-0.015,pos1[3]-0.25,pos1[2],pos1[3]-0.02]

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

ktkt=['_test_20_100_1500','_test_20_100_2000', $
  '_test_20_180_1500','_test_20_180_2000', $
  '_test_30_100_1500','_test_30_100_2000', $
  '_test_30_180_1500','_test_30_180_2000', $
  '_test_30','_test_40']

nm_ind='17'

for kk=0,0 do begin;;;for different keep_table

;ch=0
for ch=4,4 do begin
;nm_ind='17sf'+strcompress(ch,/remove_all)
;nm_ind='sf'+strcompress(ch,/remove_all)
;;;-------------------------file name-------------------------
;;;for prep data, original data
;if sys eq 0 then filename=file_search('.\'+data_namet+'\prep\*.sav',count=count) & $
;  image_sav='.\imagetest\'+data_namet+'\prep_dirty_img\' & data_sav='.\sav\'+data_namet+'\prep_img_beam\'
;if sys eq 1 then filename=file_search('./'+data_namet+'/prep/*.sav',count=count) & $
;  image_sav='./imagetest/'+data_namet+'/prep_dirty_img/' & data_sav='./sav/'+data_namet+'/prep_img_beam/'
;
;;;for prepsf data, original data
;if sys eq 0 then filename=file_search('.\'+data_namet+'\prepsf'+strcompress(ch,/remove_all)+'\*.sav',count=count) & $
;  image_sav='.\imagetest\'+data_namet+'\prep_dirty_img'+nm_ind+'\' & data_sav='.\sav\'+data_namet+'\prep_img_beam'+nm_ind+'\'
;if sys eq 1 then filename=file_search('./'+data_namet+'/prepsf'+strcompress(ch,/remove_all)+'/*.sav',count=count) & $
;  image_sav='./imagetest/'+data_namet+'/prep_dirty_img'+nm_ind+'/' & data_sav='./sav/'+data_namet+'/prep_img_beam'+nm_ind+'/'

;;;for rm data, after removing bad frames
;if sys eq 0 then filename=file_search('.\sav\'+data_namet+'\rm_cor_data'+nm_ind+'\*.sav',count=count) & $
;  image_sav='.\imagetest\'+data_namet+'\rm_dirty_img'+nm_ind+'\' & data_sav='.\sav\'+data_namet+'\rm_img_beam'+nm_ind+'\'
;if sys eq 1 then filename=file_search('./sav/'+data_namet+'/rm_cor_data'+nm_ind+'/*.sav',count=count) & $
;  image_sav='./imagetest/'+data_namet+'/rm_dirty_img'+nm_ind+'/' & data_sav='./sav/'+data_namet+'/rm_img_beam'+nm_ind+'/'
;;
;;;for mv data after removing bad frames and moving to the image center
;if sys eq 0 then filename=file_search('.\sav\'+data_namet+'\mv_cor_data'+nm_ind+'\*.sav',count=count) & $
;  image_sav='.\imagetest\'+data_namet+'\mv_dirty_img'+nm_ind+'\' & data_sav='.\sav\'+data_namet+'\mv_img_beam'+nm_ind+'\'
;if sys eq 1 then filename=file_search('./sav/'+data_namet+'/mv_cor_data'+nm_ind+'/*.sav',count=count) & $
;  image_sav='./imagetest/'+data_namet+'/mv_dirty_img'+nm_ind+'/' & data_sav='./sav/'+data_namet+'/mv_img_beam'+nm_ind+'/'
;
;;;for pm data, after doing the phase modulation
;if sys eq 0 then filename=file_search('.\sav\'+data_namet+'\pm_cor_data'+nm_ind+'\*.sav',count=count) & $
;  image_sav='.\imagetest\'+data_namet+'\pm_dirty_img'+nm_ind+'\' & data_sav='.\sav\'+data_namet+'\pm_img_beam'+nm_ind+'\'
;if sys eq 1 then filename=file_search('./sav/'+data_namet+'/pm_cor_data'+nm_ind+'/*.sav',count=count) & $
;  image_sav='./imagetest/'+data_namet+'/pm_dirty_img'+nm_ind+'/' & data_sav='./sav/'+data_namet+'/pm_img_beam'+nm_ind+'/'

;;for it data, after integration
if sys eq 0 then filename=file_search('.\sav\'+data_namet+'\it_cor_data'+nm_ind+'\cor_data_it_'+'*.sav',count=count) & $
  image_sav='.\imagetest\'+data_namet+'\it_dirty_img'+nm_ind+ktkt[kk]+'\' & data_sav='.\sav\'+data_namet+'\it_img_beam'+nm_ind+ktkt[kk]+'\'
if sys eq 1 then filename=file_search('./sav/'+data_namet+'/it_cor_data'+nm_ind+'/cor_data_it_'+'*.sav',count=count) & $
  image_sav='./imagetest/'+data_namet+'/it_dirty_img'+nm_ind+ktkt[kk]+'/' & data_sav='./sav/'+data_namet+'/it_img_beam'+nm_ind+ktkt[kk]+'/'



;-------------------------basic para-------------------------
print,count

PIX=512
nfile=count
;nfile=1

;num=600
;;;;just for pm data
;file=strmid(filename_pm[0],65,70)
;restore,filename_pm[0],/ver
;num=(size(cor_data_pm))[1]
;count=num/500

;-------------------------keep table-------------------------


;if sys eq 0 then filename_kp=file_search('.\sav\'+data_namet+'\keep_table'+nm_ind+'\*.sav')
;if sys eq 1 then filename_kp=file_search('./sav/'+data_namet+'/keep_table'+nm_ind+'/*.sav')

if sys eq 0 then filename_kp=file_search('.\sav\'+data_namet+'\keep_table'+nm_ind+ktkt[kk]+'.sav')
if sys eq 1 then filename_kp=file_search('./sav/'+data_namet+'/keep_table17.sav');;'+nm_ind+'
if sys eq 1 then filename_kp=file_search('./sav/'+data_namet+'/keep_table'+nm_ind+ktkt[kk]+'.sav')


;if sys eq 0 then filename_kp=file_search('.\sav\keep_table'+nm_ind+'.sav')
;if sys eq 1 then filename_kp=file_search('./sav/keep_table'+nm_ind+'.sav')

restore,filename_kp[0]
help,keep_table

file_mkdir,image_sav
file_mkdir,data_sav

;-------------------------begins-------------------------
for ii=15,15 do begin ;count-1 ; for which file;+20-1
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
;  cor_data_rm=cor_data_pm[ii*500:((ii+1)*500-1),*,*,*]
;  gps_time_rm=gps_time_pm[ii*500:((ii+1)*500-1),*]
  
  ;;;for it data
  restore,filename[ii],/ver
  tmp0=rstrpos(filename[ii],'cor_data_')
  file=strmid(filename[ii],tmp0,80)
  cor_data_rm=cor_data_it
  gps_time_rm=gps_time_it
  
  ;;;============begins
  cor_data=cor_data_rm
  
  ;gps_time_mv=gps_time_rm
  gps_time_obs=gps_time_rm
  gps_time_obs[*,3]=gps_time_obs[*,3]
 
  ;FLAG_TABLE=[10,11,12,13,14,15,16,17,18,19,25,26,39]
  ;FLAG_TABLE=[10,11,12,13,15,16,17,18,19,24,25,26,27,38,39]   ;flag_table 1
  ;FLAG_TABLE=[10,11,12,13,14,15,16,17,18,19,24,25,26,27,38,39];flag_table 2
  ;FLAG_TABLE=[10,11,12,13,14,15,16,17,18,19,24,25,26,38,39]  ;flag_table 3
  ;FLAG_TABLE=[10,11,12,13,14,15,16,17,18,19,25,26,38,39]  ;flag_table 4
  ;FLAG_TABLE=[10,11,12,13,14,16,17,18,19,24,25,26,27,38,39]   ;flag_table 5
   
  SUN_ANGLE=22.47/180*!PI
  nsize=size(cor_data_rm)
  cadence=12
  n=nsize[1]/cadence
  print,'size(cor_data_rm)',n
  ;n=1
  
  cor_data1=dblarr(40,40,2)

  i_st=0 ;;which frame to begin in file
  i_ed=n-1 ;;n-1 ;which frame to end in file
  
  img_mv=complexarr((i_ed-i_st+1),PIX,PIX)
  beam_mv=complexarr((i_ed-i_st+1),PIX*2,PIX*2)
  gps_time_mv=intarr((i_ed-i_st+1),9)
  LOADCT,39
;-------------------------basic para-------------------------
  FOR I=i_st,i_ed DO BEGIN
      
    DAY=gps_time_rm[0,2]
    BJT=gps_time_rm[I*CADENCE,3]+gps_time_rm[I*CADENCE,4]/60.+gps_time_rm[I*CADENCE,5]/3600.
    CAL_UVW,DAY,BJT,U,V,W;print,bjt,h,ra,dec

    cor_data1[*,*,*]=cor_data[I*CADENCE,*,*,*]

    ;FLAG_DATA,FLAG_TABLE,cor_data1,VIS,U,V      ;;;flag antenna
    flag_baseline,keep_table,cor_data1,vis,u,v   ;;;flag baseline
    
    ;;WEIGHTING,U,V,VIS
    
    ;;mov_pos,vis,u,v,-(319+25-512/2)/512.,(357+25-512/2)/512. 
    ;mov_pos,vis,u,v,-(319-512/2)/512.,(357-512/2)/512. ;for Left polarization
    ;mov_pos,vis,u,v,-(319+25-56-512/2)/512.,(357-512/2)/512. ;for 1.7GHz
    ;mov_pos,vis,u,v,-(319+25-56-512/2)/512.,(357+25-512/2)/512. ;for Left polarization
    mov_pos,vis,u,v,-(319-46-512/2)/512.,(357+25-512/2)/512. ;for Left polarization

;    wset,1
;    plot,[0,0],xrange=[-17000.,17000],yrange=[-17000.,17000],xtitle='U',ytitle='V',position=[0.15,0.15,0.95,0.95],xstyle=1,ystyle=1
;    for xy1=0,39 do begin
;      for xy2=0,39 do begin
;        plots,u[xy1,xy2],v[xy1,xy2],thick=1.,psym=2,color=255,symsize=1
;      endfor
;    endfor
    
    ;FFT_IMAGE,U,V,VIS,BEAM,IMG
    FFT_IMAGE,U,V,VIS,IMG
    FFT_BEAM,U,V,BEAM
    
;;-------------------------for one frame dirty image-------------------------    
    if i eq i_st then time_string=time_gps_to_string(gps_time_obs[I*CADENCE,*])
    time_obs=time_gps_to_obs(gps_time_obs[I*CADENCE,*])
    inputimg=abs(img)
    ;;center_fit,inputimg,centx,centy
    centx=256 & centy=256
    add_disk,inputimg,0.1,[centx,centy],solar_disk
    ;;print,'center:',centx,centy
    
    wset,5
    plot_image,real_part(img),XTICKLEN=-0.01,YTICKLEN=-0.01,min=min(real_part(img)),max=max(real_part(img)), $
      charsize=chars1,position=pos1, $
      title='MUSER dirty map(self '+strcompress(ch,/remove_all)+') -real Part- '+time_obs
    contour,solar_disk,levels=[0.5]*max(solar_disk),/overplot,color=255
    plots,[centx,centx],[centy,centy],psym=1,symsize=3,color=255  
    cgcolorbar,range=[min(real_part(img)),max(real_part(img))],yticks = 2.,/ver,charsize=chars1,position=posc1,color=255  

;    plot_image,abs(img),XTICKLEN=-0.01,YTICKLEN=-0.01,min=min(abs(img)),max=max(abs(img)), $
;      charsize=chars1,position=pos1, $
;      title='MUSER dirty map '+time_obs;+' ['+strcompress(centx,/remove_all)+','+strcompress(centy,/remove_all)+']'
;    ;contour,solar_disk,levels=[0.5]*max(solar_disk),/overplot,color=255
;    ;plots,[centx,centx],[centy,centy],psym=1,symsize=3,color=255  
;    cgcolorbar,range=[min(abs(img)),max(abs(img))],yticks = 2.,/ver,charsize=chars1,position=posc1,color=255
;      
    img_mv[i-i_st,*,*]=img[*,*];(abs(img))[*,*]
    beam_mv[i-i_st,*,*]=beam[*,*];(abs(beam))[*,*]
    gps_time_mv[i-i_st,*]=gps_time_rm[I*CADENCE,*]
;;;-------------------------intergration and fitting center of dirty image-------------------------

;    if i eq i_ed then begin
;      img_int=make_array(PIX,PIX)
;      for cc=0,PIX-1 do begin
;        for ee=0,PIX-1 do begin
;          img_int[cc,ee]=mean(img_mv[0:(i-i_st),cc,ee])
;        endfor
;      endfor
;      print,i,minmax(abs(img_mv))
;      print,i,minmax(img_int)
;      beam_int=make_array(PIX*2,PIX*2)
;      for cc=0,2*PIX-1 do begin
;        for ee=0,2*PIX-1 do begin
;          beam_int[cc,ee]=mean(beam_mv[0:(i-i_st),cc,ee])
;        endfor
;      endfor
;      print,i,minmax(abs(beam_mv))
;      print,i,minmax(beam_int) 
;    endif
;
;;    inputimg=img_int
;;    center_fit,inputimg,centx,centy
;;    add_disk,inputimg,0.1,[centx,centy],solar_disk
;;    print,'center:',centx,centy
; 
;    plot_image,img_int,XTICKLEN=-0.01,YTICKLEN=-0.01,min=min(img_int),max=max(img_int), $
;      title='Int: MUSER dirty map @1.725GHz '+time_obs+ $
;      ' ['+strcompress(centx,/remove_all)+','+strcompress(centy,/remove_all)+']',charsize=chars1,position=pos2,/noerase
;    ;contour,solar_disk,levels=[0.5]*max(solar_disk),/overplot,color=255
;    ;plots,[centx,centx],[centy,centy],psym=1,symsize=3,color=255
;    cgcolorbar,range=[min(img_int),max(img_int)],yticks = 2.,/ver,charsize=chars1,position=posc2,color=255

;;;-------------------------save images-------------------------       
    
    if ((i-i_st) mod (100)) eq 0 then $
      write_jpeg,image_sav+'dirtyimg_mv_'+file+'_'+strcompress(i,/remove_all)+'.jpg',tvrd(true=3),true=3
      
    if i eq i_ed then $
      write_jpeg,image_sav+'dirtyimg_mv_'+file+'_'+strcompress(i,/remove_all)+'.jpg',tvrd(true=3),true=3

  endfor
  ;;ch=4
  save,img_mv,beam_mv,gps_time_mv,ch,filename=data_sav+'img_beam_'+file
  ;gps_time_int=gps_time_mv[(i_ed-i_st)/2,*]
  ;save,img_int,beam_int,gps_time_int,filename=data_sav+'img_beam_int'+file
  ;;+time_string+
  

endfor

;-------------------------ending-------------------------
endfor
endfor
print,'OK'
end