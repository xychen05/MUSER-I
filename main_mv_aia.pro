;.r main_mv_aia.pro
;.compile package_muser.pro
;
; Name :
; main_mv_aia
; 
; Purpose :
; Display dirty images.
; 
; Explanation :
; This procedure is for the center fitting of 'AIA or NORH dirty map' and MUSER dirty map.
; All the dirty images have removed bad frames and moved to the image center.
; AIA image has been transformed into dirty map through the convolution with the dirty beam of MUSER.
; The left dirty image is 'AIA or NORH dirty map' while the right dirty MUSER dirty map.
; You can also use the ideal solar disk model to fit the solar center.
;
; Notes :
; Rotate after moving the solar center to the image center.
; 
; History :
; Wei Wang, MAIN.PRO
; Xingyao Chen, June 2017
; 
;-------------------------begins-------------------------
Device, RETAIN=2
!P.FONT = 0
;!P.multi=[0,2,1]

WINDOW,1,XS=512*2.,YS=512

pos1 = [0.07,0.1,0.47,0.9]
pos2 = [0.57,0.1,0.97,0.9]

posc1 = [pos1[2]-0.015,pos1[3]-0.25,pos1[2],pos1[3]-0.02]
posc2 = [pos2[2]-0.015,pos2[3]-0.25,pos2[2],pos2[3]-0.02]

chars1=0.9
pix=512

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

num=60
pix=512

;;;-------------------------file name-------------------------
;;;;for prep data
;if sys eq 0 then filename=file_search('.\'+data_namet+'\prep\*.sav',count=count) & $
;  image_sav='.\imagetest\'+data_namet+'\prep_dirty_img\' & data_sav='.\sav\'+data_namet+'\prep_img_beam\'
;if sys eq 1 then filename=file_search('./'+data_namet+'/prep/*.sav',count=count) & $
;  image_sav='./imagetest/'+data_namet+'/prep_dirty_img/' & data_sav='./sav/'+data_namet+'/prep_img_beam/'

;;;for rm data
;if sys eq 0 then filename=file_search('.\sav\'+data_namet+'\rm_cor_data'+nm_ind+'\*.sav',count=count) & $
;  image_sav='.\imagetest\'+data_namet+'\rm_dirty_img'+nm_ind+'\' & data_sav='.\sav\'+data_namet+'\rm_img_beam'+nm_ind+'\'
;if sys eq 1 then filename=file_search('./sav/'+data_namet+'/rm_cor_data'+nm_ind+'/*.sav',count=count) & $
;  image_sav='./imagetest/'+data_namet+'/rm_dirty_img'+nm_ind+'/' & data_sav='./sav/'+data_namet+'/rm_img_beam'+nm_ind+'/'

;;for it data, after integration
if sys eq 0 then filename=file_search('.\sav\'+data_namet+'\it_cor_data00\cor_data_it_'+'*.sav',count=count) & $
  image_sav='.\imagetest\'+data_namet+'\it_dirty_img00\' & data_sav='.\sav\'+data_namet+'\it_img_beam00\'
if sys eq 1 then filename=file_search('./sav/'+data_namet+'/it_cor_data00/cor_data_it_'+'*.sav',count=count) & $
  image_sav='./imagetest/'+data_namet+'/it_dirty_img00/' & data_sav='./sav/'+data_namet+'/it_img_beam00/'

;;
if sys eq 0 then filenamea='I:\muserdata_20141217\AIA\AIA20141217_030011_0171.fits'
if sys eq 0 then filenamen='I:\muserdata_20141217\AIA\ifa141217_024433'
if sys eq 1 then filenamea='/Volumes/Seagate_cxy/muserdata_20141217/AIA/AIA20141217_030011_0171.fits'
if sys eq 1 then filenamen='/Volumes/Seagate_cxy/muserdata_20141217/AIA/ifa141217_024433'

;aia_data_prep,filenamea,pix,1.5,imgaia
norh_data_prep,filenamen,pix,1.5,imgnorh

print,count

file_mkdir,image_sav
file_mkdir,data_sav

if sys eq 0 then restore,'.\sav\keep_table.sav',/ver
if sys eq 1 then restore,'./sav/keep_table.sav',/ver
help,keep_table

;-------------------------basic para-------------------------

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
  
  ;;;for it data
  restore,filename[0],/ver
  tmp0=rstrpos(filename[0],'cor_data_')
  file=strmid(filename[0],tmp0,80)
  cor_data_rm=cor_data_it
  gps_time_rm=gps_time_it
  
  ;;============  
  cor_data=cor_data_rm
  
  gps_time_mv=gps_time_rm
  gps_time_obs=gps_time_rm
  gps_time_obs[*,3]=gps_time_obs[*,3];-8
  
  ;% Flag_Ant =[3 5 10 11 12 13 14 16 17 18 19 20 25 26 27 28 39 40];  % 20141217Solar Observation Absence
  ;% Flag_Ant =[11 12 13 14 16 17 18 19 20 25 26 27 39 40];  % 20141217Solar Observation Absence modified
  ;% Flag_Ant =[10 11 12 13 16 17 18 19 25 26 39 40];  % 20141217Solar Observation Absence
  
  ;FLAG_TABLE=[10,11,12,13,15,16,17,18,19,24,25,26,27,38,39]   ;flag_table 1
  FLAG_TABLE=[10,11,12,13,14,15,16,17,18,19,24,25,26,27,38,39];flag_table 2
  ;FLAG_TABLE=[10,11,12,13,14,15,16,17,18,19,24,25,26,38,39]  ;flag_table 3
  ;FLAG_TABLE=[10,11,12,13,14,15,16,17,18,19,25,26,38,39]  ;flag_table 4
  
  SUN_ANGLE=22.47/180*!PI ;0.392175
  nsize=size(cor_data_rm)
  cadence=1
  n=nsize[1]/cadence
  print,'size(cor_data_rm)',n
  
  cor_data1=dblarr(40,40,2) ;set for one time point

  i_st=0 ;;which frame to begin in file
  i_ed=0 ;;n-1 ;which frame to end in file
  
  img_mv=make_array((i_ed-i_st+1),PIX,PIX) ;make for integration
  beam_mv=make_array((i_ed-i_st+1),PIX*2,PIX*2) ;make for integration
  LOADCT,39
  ;-------------------------basic para-------------------------
  FOR I=i_st,i_ed DO BEGIN

    DAY=gps_time_rm[0,2]
    BJT=gps_time_rm[I*CADENCE,3]+gps_time_rm[I*CADENCE,4]/60.+gps_time_rm[I*CADENCE,5]/3600.
    CAL_UVW,DAY,BJT,U,V,W ;print,bjt,h,ra,dec

    cor_data1[*,*,*]=cor_data[I,*,*,*]

    FLAG_DATA,FLAG_TABLE,cor_data1,VIS,U,V
    ;;flag_baseline,keep_table,cor_data1,vis,u,v   ;;;flag baseline
 
    mov_pos,vis,u,v,-(319-512/2)/512.,(357-512/2)/512. ;for Left polarization
    ;mov_pos,vis,u,v,-(303-512/2)/512.,(362-512/2)/512. ;for Right polarization
    
    ;mov_pos,vis,u,v,-(238-512/2)/512.,(310-512/2)/512. ;for ideal disk
    ;mov_pos,vis,u,v,-(229-512/2)/512.,(316-512/2)/512. ;for aia
    ;mov_pos,vis,u,v,-(230-512/2)/512.,(308-512/2)/512. ;for norh
    
    FFT_IMAGE,U,V,VIS,BEAM,IMG
    ;;angle = 9.35 for day 1217 & angle = 8.89 for day 1218
    angle=9.35+(8.89-9.35)/24.*3+(8.89-9.35)/24./60.*(ii+2) ;1 min for rotating
    ;img=polax_rot(abs(img), Angle) ; rotate after moving the solar center to the image center
    
    if i eq i_st then time_string=time_gps_to_string(gps_time_obs[I,*])
    time_obs=time_gps_to_obs(gps_time_obs[I,*])
        
    img_mv[i-i_st,*,*]=(abs(img))[*,*]
    beam_mv[i-i_st,*,*]=(abs(beam))[*,*]
    
    img_int=make_array(pix,pix)
    for cc=0,PIX-1 do begin
      for ee=0,PIX-1 do begin
        img_int[cc,ee]=mean(img_mv[0:(i-i_st),cc,ee])
      endfor 
    endfor
    print,i,minmax(abs(img))
    print,i,minmax(img_int)
    
    ;;-------------------------for aia or norh-------------------------
    dbeam=abs(beam)
    ;img2dirty,imgaia,dbeam,imgconv
    img2dirty,imgnorh,dbeam,imgconv
    inputimg=img
    center_fit_d2d,inputimg,imgconv,centx,centy  
    ;centx=512/2 & centy=512/2
    print,'center',centx,centy  
    add_disk,inputimg,0.1,[centx,centy],solar_disk ;to design the contour of solar disk
     
    ;plot_image,imgconv,title='AIA 171 to dirty map',charsize=chars1,position=pos1
    plot_image,imgconv,title='NORH 17GHz to dirty map',charsize=chars1,position=pos1
    cgcolorbar,range=[min(imgconv),max(imgconv)],yticks = 2.,/ver,charsize=chars1,position=posc1,color=255
    
;    ;;;-------------------------for ideal solar disk-------------------------
;    inputimg=abs(img)
;    center_fit,inputimg,centx,centy
;    ;centx=512/2 & centy=512/2
;    add_disk,inputimg,0.1,[centx,centy],solar_disk
;    
;    plot_image,solar_disk,title='ideal solar disk',charsize=chars1,position=pos1
;    cgcolorbar,range=[min(solar_disk),max(solar_disk)],yticks = 2.,/ver,charsize=chars1,position=posc1,color=255

    ;;-------------------------ploting-------------------------
    plot_image,abs(img),XTICKLEN=-0.01,YTICKLEN=-0.01,min=min(abs(img)),max=max(abs(img)), $
      charsize=chars1,position=pos2, $
      title='MUSER dirty map @1.725GHz '+time_obs,/noerase
    contour,solar_disk,levels=[0.5]*max(solar_disk),/overplot,color=255
    plots,[centx,centx],[centy,centy],psym=1,symsize=3,color=255
    ;plots,[512/2,512/2],[512/2,512/2],psym=1,symsize=3 ,color=255
    cgcolorbar,range=[min(abs(img)),max(abs(img))],yticks = 2.,/ver,charsize=chars1,position=posc2,color=255

    ;;-------------------------image saving-------------------------    
 
    if ((i-i_st) mod (100)) eq 0 then $
      write_jpeg,image_sav+'aia_dirtyimg_'+file+'_'+strcompress(i,/remove_all)+'.jpg',tvrd(true=3),true=3
      
    if i eq i_ed then $
      write_jpeg,image_sav+'aia_dirtyimg_'+file+'_'+strcompress(i,/remove_all)+'.jpg',tvrd(true=3),true=3

  endfor
  ;save,img_mv,beam_mv,gps_time_mv,filename=data_sav+'aia_img_beam_'+file
  
endfor
;-------------------------ending-------------------------
print,'It is OK..................'
end