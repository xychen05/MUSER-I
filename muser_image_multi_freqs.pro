;.r muser_image_multi_freqs.pro
;.compile package_muser.pro 
;.compile package_xy.pro
;
; Purpose :
; plot image of MUSER clean image.
; 
; History :
; Xingyao Chen, 31 Oct 2017.
; Xingyao Chen, 28 Nov 2017.
;-------------------------begins-------------------------
;WINDOW,0,xsize=700, ysize=700
;Device,RETAIN=2
;!P.FONT = 0

ENTRY_DEVICE = !D.NAME
SET_PLOT, 'PS'
!P.FONT = 0
xsz=8.0 & ysz=6.0;region 5

charsz=0.4 & charth=1.2 & th_ind =1

x1 = -1000. & x2 = 200 & y1 = -1200. & y2 = 0. ;region 2
x1 = -420. & x2 = 220 & y1 = -620. & y2 = -140. ;region 5
;x1 = -370. & x2 = 130 & y1 = -575. & y2 = -200. ;region 5
;x1 = -750. & x2 = 250 & y1 = -700. & y2 = 50. ;region 5

x11 = x1 & x21 = x2 & y11 = y1 & y21 = y2 ;region 2

xrange = [x1,x2] & yrange = [y1,y2]

pos1 = [0.12,0.12,0.92,0.92]
posc1 = [pos1[2]+0.07,pos1[3]-0.2,pos1[2]+0.09,pos1[2]]

systim = SYSTIME(1)

;;-------------------------171-------------------------
;nm_ind='05'
;filename = file_search('./sav/1.6-2.0_Lnum2400/it_clean_img'+nm_ind+'_fits_raw_LR/*.fits',count=count)

;-------------------------aia image-------------------------
filenamea='/Volumes/Seagate_cxy/muserdata_20141217/AIA/AIA20141217_043207_0304.fits';;;nn=725
;filenamea='/Volumes/Seagate_cxy/muserdata_20141217/AIA/AIA20141217_042407_0304.fits';;;nn=255

filenamea1 = file_search('/Volumes/TOSHIBA_cxy/data_20141217/AIA/prep/171/AIA20141217*.fits')
filenamea2 = file_search('/Volumes/TOSHIBA_cxy/data_20141217/AIA/prep/131/AIA20141217*.fits')
filenamea3 = file_search('/Volumes/TOSHIBA_cxy/data_20141217/AIA/prep/94/AIA20141217*.fits')

wavelength = ['94','171','193','211','304','335','131','1600']
mag_scal = [1.5,1.5,3,8.,3.,2.,1.5,5.]
rs_ratio = 1.0d
xllp=0.
nxp=4096.
yllp=0.
nyp=4096.

;wl=4
;wavelnth = wavelength[wl]
;aia_lct, rr, gg, bb, wavelnth=wavelnth, /load
;read_sdo, filenamea, indexa, dataa, /UNCOMP_DELETE,/noshell,/use_shared_lib,/HIDE
;dataa=dataa/indexa.EXPTIME
;imga = sdo_aia_scale_hdr(dataa/mag_scal[wl],indexa,xllp,yllp,nxp,nyp,wavelnth=wavelnth,rs_ratio=rs_ratio,/no_imgrscale)
;index2map,indexa,imga,mapa

;-------------------------basic parameters-------------------------
pixel=512

nt_ind=5;
file_st=1030;+600
nfile=count
nfile=1

nt_indm=6
file_stm=0

;-------------------------Part I: select points and regions-------------------------
for nn=1,19 do begin;;99
  
  imga = intarr(4096,4096,3)
  
  tmp0=rstrpos(filenamea2[file_st+nt_ind*nn],'AIA20141217_')
  file=strmid(filenamea2[file_st+nt_ind*nn],tmp0+3,15)
  print,file

  filename_eps='/Users/xychen/Desktop/movie/vs03/image_left_multi_freq_dfro_rsdalign/'+'muser_rd_'+file+'_5e6.eps';;guass_;;'+file+'
  device,file=filename_eps,xsize=xsz,ysize=ysz,BITS_PER_PIXEL=8,$
    SET_FONT='Helvetica',/ENCAPSULATED,scale_factor=2,/color
    
  wl=1;171
  wavelnth = wavelength[wl]
  read_sdo, filenamea1[file_st+nt_ind*nn], indexa1, dataa1, /UNCOMP_DELETE,/noshell,/use_shared_lib,/HIDE
  dataa1=dataa1/indexa1.EXPTIME
  img1 = sdo_aia_scale_hdr(dataa1/mag_scal[wl],indexa1,xllp,yllp,nxp,nyp,wavelnth=wavelnth,rs_ratio=rs_ratio,/no_imgrscale)
  index2map,indexa1,img1,map
  
  wl=6;304
  wavelnth = wavelength[wl]
  read_sdo, filenamea2[file_st+nt_ind*nn-1], indexa2, dataa2, /UNCOMP_DELETE,/noshell,/use_shared_lib,/HIDE
  dataa2=dataa2/indexa2.EXPTIME
  img2 = sdo_aia_scale_hdr(dataa2/mag_scal[wl],indexa2,xllp,yllp,nxp,nyp,wavelnth=wavelnth,rs_ratio=rs_ratio,/no_imgrscale)
  
  wl=0;304
  wavelnth = wavelength[wl]
  read_sdo, filenamea3[file_st+nt_ind*nn], indexa3, dataa3, /UNCOMP_DELETE,/noshell,/use_shared_lib,/HIDE
  dataa2=dataa2/indexa2.EXPTIME
  img3 = sdo_aia_scale_hdr(dataa3/mag_scal[wl],indexa3,xllp,yllp,nxp,nyp,wavelnth=wavelnth,rs_ratio=rs_ratio,/no_imgrscale)
  
  imga[*,*,0]=img1
  imga[*,*,1]=img2
  imga[*,*,2]=img3
  indexlzero=where(imga lt 0)
  imga[indexlzero]=0
  loadct,0
  x1pix=round(indexa1.crpix1+x1/indexa1.cdelt1-1)
  x2pix=round(indexa1.crpix1+x2/indexa1.cdelt1-1)
  y1pix=round(indexa1.crpix2+y1/indexa1.cdelt2-1)
  y2pix=round(indexa1.crpix2+y2/indexa1.cdelt2-1)
  
  plot_image,sqrt(imga[x1pix:x2pix,y1pix:y2pix,*]),origin=[x11,y11], $
    ;-[(indexa1.crpix1-1)*indexa1.cdelt1, $ ;,min=sqrt(4), max=sqrt(1600),true=3
    ;(indexa1.crpix2-1)*indexa1.cdelt2], $
    scale=[indexa1.cdelt1,indexa1.cdelt2],title=' ', $;;Quiet- guassian rotate
    xrange=xrange,yrange=yrange, $
    xtitle='X-postion (arcseconds)', $
    ytitle='Y-position (arcseconds)', $
    position = pos1,xthick=charth*th_ind,ythick=charth*th_ind,charthick=charth,charsize=charsz,/NOADJUST;,/normal
  plot_map,map,grid=15,/no_data,/noerase,color=255,title=' ', $
    xrange = xrange,yrange = yrange,position = pos1,/NOAXES,/NOLABELS,/NOXTICKS,/NOYTICKS
  
  xyouts,pos1[0]+0.45,pos1[1]+0.02, textoidl('SDO/AIA-')+strmid(string(indexa1.WAVELNTH),9,4)+' '+ $
    strmid(string(indexa2.WAVELNTH),9,4)+' '+strmid(string(indexa3.WAVELNTH),9,4)+' '+ $
    strmid(indexa1.date_OBS,11,8)+'(background)',/normal,color=255,align=0,charsize=charsz/1.2,charthick=charth
    
  print,strmid(filenamea1[file_st+nt_ind*nn],48,23)
  print,strmid(filenamea2[file_st+nt_ind*nn],48,23)
  print,strmid(filenamea3[file_st+nt_ind*nn],48,23)

  
  ;-------------------------muser-------------------------
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
  
  loadct,39
  for ott=0,1 do begin
    ooo=1-ott
    data_namet=data_name[ooo*2+4]
    
    for chtt=0,15 do begin
      ch=15-chtt
      freqch=(1200+400*ooo+ch*25+25)
      freqchstr=strcompress(freqch,/remove_all)
      nm_ind='sf'+strcompress(ch,/remove_all)
      filenamem = file_search('./sav/'+data_namet+'/it_clean_img'+nm_ind+'_fits_raw/'+'*.fits',count=count)
      ;'/quiet_clean_img_fits/'+'*'+freqchstr+'*.fits')
      ;'/it_clean_img'+nm_ind+'_fits_raw/'+'*.fits',count=count)
      
      tmp0=rstrpos(filenamem[file_stm+nt_indm*nn],'MUSER_raw_')
      tmp1=rstrpos(filenamem[file_stm+nt_indm*nn],'.fits')
      file=strmid(filenamem[file_stm+nt_indm*nn],tmp0,tmp1-tmp0)
      print,file
      
      datam = readfits( filenamem[file_stm+nt_indm*nn], indexm, /NOSCALE, /SILENT)
      
;      ;;============for lift-right plot 
;      dataml = datam
;      data_namemr=data_name[ooo*2+5]
;      filenamemr = file_search('./sav/'+data_namemr+'/it_clean_img'+nm_ind+'_fits_raw/'+'*.fits',count=count)
;      filemr=strmid(filenamemr[file_stm+nt_indm*nn],tmp0,tmp1-tmp0)
;      print,filemr
;      datamr = readfits( filenamemr[file_stm+nt_indm*nn], indexmr, /NOSCALE, /SILENT)
;      datam=datam-datamr
;  
;      ;;============data_align
;      data_align=file_search('./sav/'+data_namet+'/quiet_clean_img_fits_align/'+'*'+freqchstr+'.sav')
;      print,data_align
;      restore,data_align;,/ver
;      sz=size(alignm)
;      diff3=alignm[*,*,*,3]
;      temp=min(diff3,loca3)
;      inf3=array_indices(diff3, loca3)
;      print,inf3
;      print,alignm[inf3[0],inf3[1],inf3[2],3]
;      xx=inf3[0] & yy=inf3[1] & zz=inf3[2];inf3[2]
;      xx0=sz[1] & yy0=sz[2] & zz0=sz[3]
;      print,xx,yy,zz
;      ;;xx=15 & yy=15 & zz=5
;      imgmt=mv_img(datam,-16+xx-xx0/2,-52+yy-yy0/2)
;      imgm=polax_rot(imgmt,19+zz-zz0/2);19+zz-5
;      
      ;;============data_align of residual image
      data_align=file_search('./sav/'+data_namet+'/quiet_clean_img_fits_residual_align/'+'*.sav')
      restore,data_align
      imgmt1=mv_img(datam,-11,-47)  ;;;for 1.7 GHz left
      imgm1=shift_img(imgmt1, disp[*,ch+1])
      
      imgm=polax_rot(imgm1,18) ;;for  left
      ;imgmt1=mv_img(datam,-11,-47)  ;;;for 1.7 GHz right
      ;imgm=polax_rot(imgm1,17)  ;;for  right
      
;      ;;============original
;      datamt=mv_img(datam,-16,-52)
;      angle=19 ;1 min for rotating
;      datam=datamt
;      imgm=polax_rot(datam, angle) ; rotate after moving the solar center to the image center
;           
;      ;dmin=5 & dmax=15 ;;for ch=0
      ;imgm=datam
      ;index2map,indexm,imgm,mapm;
      
      ;;============guassian fit==========method 1
      index1=where(imgm ge 0.75*max(imgm));10.0 for flare;;3 for quiet ;;1 for left-right
      print,'pixels of radio source:',(size(index1))[1]
      locx=index1 mod pixel
      locy=index1 / pixel
      ind_rg=[min(locx),max(locx),min(locy),max(locy)]
      ind1=min(ind_rg)-10 & ind2=max(ind_rg);+10
      ind1=170 & ind2=270
      if ind1 lt 0 then ind1=min(ind_rg)
      print,ind1,ind2
      imgm_rg=imgm[ind1:ind2,ind1:ind2]
      imgm_rg_fit=GAUSS2DFIT(imgm_rg, coeff, /TILT)
      gcenx=coeff[4]+ind1 & gceny=coeff[5]+ind1
;      imgm_fit=imgm
;      imgm_fit[ind1:ind2,ind1:ind2]=imgm_rg_fit
;
;      index2map,indexm,imgm_fit,mapm;

;      ;;============guassian fit==========method 2
;      index1=where(imgm ge 0.75*max(imgm))
;      print,'pixels of radio source:',(size(index1))[1]
;      locx=index1 mod pixel
;      locy=index1 / pixel
;      ind_rg=[min(locx),max(locx),min(locy),max(locy)]
;      ind1=min(ind_rg)-10 & ind2=max(ind_rg);+10
;      if ind1 lt 0 then ind1=min(ind_rg)
;      print,ind1,ind2
;      imgm_rg=imgm[ind1:ind2,ind1:ind2]
;      imgm_rg_fit=GAUSS2DFIT(imgm_rg, coeff, /TILT)
;      gcenx=coeff[4]+ind1 & gceny=coeff[5]+ind1
  
      ;;============maximum brightness
      indmm=where(imgm eq max(imgm))
      mcenx=indmm[0] mod pixel
      mceny=indmm[0] / pixel

      ;;============
      imgm=smooth(imgm,2)
      index2map,indexm,imgm,mapm;
      
;      ;;============polarization calibration
;      polar=(mean(dataml[190:270,170:250])-mean(datamr[190:270,170:250]))/(mean(dataml[190:270,170:250])+mean(datamr[190:270,170:250]))
;      print,freqchstr,' MHz ',polar
;      
      ;;============ploting
      plot_map,mapm,position = pos1,/notitle, $
        xrange = xrange,yrange = yrange,xstyle=1,ystyle=1, $ ;,/NOLABELS,/NOAXES,/NOXTICKS,/NOYTICKS,/NOTITLE
        thick=charth*th_ind,charthick=charth,charsize=charsz,LEVELS=[5e6],/cont,color=250-(8*ch+ooo*127),/over
      ;0.5*max(imgm)
      plots,(gcenx-256)*60.*60./512.,(gceny-256)*60.*60./512.,psym=1,thick=charth*th_ind,color=250-(8*ch+ooo*127),symsize=charsz*2
      
      xyouts,pos1[2]-0.1,pos1[1]+0.02*(ch+2.5)+0.32*ooo, freqchstr+' MHz' ,$
        /normal,color=250-(8*ch+ooo*127),align=0,charsize=charsz/1.3,charthick=charth        
      
;      xyouts,pos1[0]+0.01+0.08*ooo,pos1[1]+0.02*(ch+2), freqchstr+' MHz' ,$
;        /normal,color=8*ch+ooo*127,align=0,charsize=charsz,charthick=charth
        
    endfor 
  endfor
  xyouts,pos1[2]-0.1-0.02,0.82, '  MUSER-I-L' ,$;+strmid(data_namet,8,1)
    /normal,color=255,align=0,charsize=charsz/1.3,charthick=charth
  
endfor

;-------------------------plot-------------------------
device,/CLOSE
set_plot,'X'
print,'It is all ready..........................'
end