;.r muser_image_multi_freqs_movie_spec_flux_image.pro
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
xsz=12 & ysz=6.0;region 5

charsz=0.4 & charth=1.2 & th_ind =1

charsz=0.4 & charth=1.2 & th_ind =1

x1 = -1000. & x2 = 200 & y1 = -1200. & y2 = 0. ;region 2
x1 = -420. & x2 = 220 & y1 = -620. & y2 = -140. ;region 5
x1 = -390. & x2 = 250 & y1 = -620. & y2 = -140. ;region 5
;x1 = -370. & x2 = 130 & y1 = -575. & y2 = -200. ;region 5
;x1 = -750. & x2 = 250 & y1 = -700. & y2 = 50. ;region 5
x1 = -310. & x2 = 210 & y1 = -620. & y2 = -200. ;region 5


x11 = x1 & x21 = x2 & y11 = y1 & y21 = y2 ;region 2

xrange = [x1,x2] & yrange = [y1,y2]

pos1 = [0.12,0.12,0.92,0.92]
;pos1 = [0.,0.,1.,1.]
;posc1 = [pos1[2]+0.07,pos1[3]-0.2,pos1[2]+0.09,pos1[2]]

pos1 = [0.13,0.05,0.93,0.45]
pos2= [0.13,0.5,0.93,0.72]
pos3= [0.13,0.72,0.93,0.98]
framedims   = [800*2.5,1200*2.5]

pos1 = [0.13,0.95-6./9.*0.8,0.93,0.95]
pos3 = [0.13,0.06,0.93,0.35]
framedims   = [800*2.5,900*2.5]

pos1 = [0.13,0.98-6./10.*0.8,0.93,0.98]
pos3 = [0.13,0.1,0.93,0.43]
framedims   = [800*2.5,1000*2.5]

pos1 = [0.035,0.10,0.035+0.86/42*52./2.,0.96]
pos3 = [0.625,0.10,0.985,0.50]
pos2 = [0.625,0.56,0.985,0.96]

systim = SYSTIME(1)

;;-------------------------171-------------------------
;nm_ind='05'
;filename = file_search('./sav/1.6-2.0_Lnum2400/it_clean_img'+nm_ind+'_fits_raw_LR/*.fits',count=count)

;-------------------------aia image-------------------------
filenamea='/Volumes/Seagate_cxy/muserdata_20141217/AIA/AIA20141217_043207_0304.fits';;;nn=725
;filenamea='/Volumes/Seagate_cxy/muserdata_20141217/AIA/AIA20141217_042407_0304.fits';;;nn=255

;xy1='94'
;xy2='193'
;xy3='335'

xx=['94','171','193','211','304','335','131']

for xyi=0,0 do begin
  for xyj=0,0 do begin
    for xyk=0,0 do begin
      
      xy1=xx[xyi]
      xy2=xx[xyj]
      xy3=xx[xyk]
      
      xy1='94'
      xy2='193'
      xy3='171'
    
filenamea1 = file_search('/Volumes/TOSHIBA_cxy/data_20141217/AIA/prep/'+xy1+'/AIA20141217*.fits')
filenamea2 = file_search('/Volumes/TOSHIBA_cxy/data_20141217/AIA/prep/'+xy2+'/AIA20141217*.fits')
filenamea3 = file_search('/Volumes/TOSHIBA_cxy/data_20141217/AIA/prep/'+xy3+'/AIA20141217*.fits')

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

nt_ind=1;
file_st=1030;+600
nfile=count
nfile=1

nt_indm=1
file_stm=0

;-------------------------for spectrum plot and flux plot-------------------------
restore,'/Volumes/Seagate_cxy/muserdata_20141217/muser_spec/sav/spec_time_dataLR_041955-044555_A1.sav',/ver
timerange=[anytim('17-Dec-14 04:25:00.000'),anytim('17-Dec-14 04:46:00.000')]
;timerange=[time[0],max(time)]
;-------------------------Part I: select points and regions-------------------------
for nn=30,30 do begin;;99;49

  filename_eps='/Users/xychen/Desktop/'+'muser_l-r_rd'+string(nn)+'-'+xy1+'-'+xy2+'-'+xy3+'.eps';;guass_;;'
  device,file=filename_eps,xsize=xsz,ysize=ysz,BITS_PER_PIXEL=8,$
    SET_FONT='Helvetica',/ENCAPSULATED,scale_factor=2,/color
  
  imga = intarr(4096,4096,3)
  
  tmp0=rstrpos(filenamea2[file_st+nt_ind*nn],'AIA20141217_')
  file=strmid(filenamea2[file_st+nt_ind*nn],tmp0+3,15)
  print,file
; ;;image_l-r_multi_freq_dfro_rsdalign/
;  filename_eps='/Users/xychen/Desktop/movie/'+'muser_l_rd_'+file+'.eps';;guass_;;'+file+'
;  device,file=filename_eps,xsize=xsz,ysize=ysz,BITS_PER_PIXEL=8,$
;    SET_FONT='Helvetica',/ENCAPSULATED,scale_factor=2,/color
    
  wl=where(wavelength eq xy1);171
  wl=wl[0]
  wavelnth = wavelength[wl]
  read_sdo, filenamea1[file_st+nt_ind*nn], indexa1, dataa1, /UNCOMP_DELETE,/noshell,/use_shared_lib,/HIDE
  dataa1=dataa1/indexa1.EXPTIME
  img1 = sdo_aia_scale_hdr(dataa1/mag_scal[wl],indexa1,xllp,yllp,nxp,nyp,wavelnth=wavelnth,rs_ratio=rs_ratio,/no_imgrscale)
  index2map,indexa1,img1,map
  
  wl=where(wavelength eq xy2);304
  wl=wl[0]
  wavelnth = wavelength[wl]
  read_sdo, filenamea2[file_st+nt_ind*nn-1], indexa2, dataa2, /UNCOMP_DELETE,/noshell,/use_shared_lib,/HIDE
  dataa2=dataa2/indexa2.EXPTIME
  img2 = sdo_aia_scale_hdr(dataa2/mag_scal[wl],indexa2,xllp,yllp,nxp,nyp,wavelnth=wavelnth,rs_ratio=rs_ratio,/no_imgrscale)
  
  wl=where(wavelength eq xy3);304
  wl=wl[0]
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
    position = pos1,xthick=charth*th_ind,ythick=charth*th_ind,charthick=charth,charsize=charsz,/NOADJUST,color=-1
  plot_map,map,grid=15,/no_data,/noerase,color=255,title=' ', $
    xrange = xrange,yrange = yrange,position = pos1,/NOAXES,/NOLABELS,/NOXTICKS,/NOYTICKS
  
  xyouts,pos1[0]+0.27,pos1[1]+0.01, textoidl('SDO/AIA-')+strmid(string(indexa1.WAVELNTH),9,4)+' '+ $
    strmid(string(indexa2.WAVELNTH),9,4)+' '+strmid(string(indexa3.WAVELNTH),9,4)+' '+ $
    strmid(indexa1.date_OBS,11,8)+'(background)',/normal,color=255,align=0,charsize=charsz/1.1,charthick=charth
    
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
      
      ;;============guassian fit
      index1=where(imgm ge 1.4);10.0 for flare;;3 for quiet ;;1 for left-right
      print,'pixels of radio source:',(size(index1))[1]
      locx=index1/pixel
      locy=index1 mod pixel
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
      
      ;;============
      imgm=smooth(imgm,2)
      index2map,indexm,imgm,mapm;
      
;      ;;============polarization calibration
;      polar=(mean(dataml[190:270,170:250])-mean(datamr[190:270,170:250]))/(mean(dataml[190:270,170:250])+mean(datamr[190:270,170:250]))
;      print,freqchstr,' MHz ',polar
;      
      
      if (nn eq 37) and (freqchstr eq '1550') then goto,nafik
      if (nn eq 35) and (freqchstr eq '1550') then goto,nafik
      ;;============ploting
      plot_map,mapm,position = pos1,/notitle, $
        xrange = xrange,yrange = yrange,xstyle=1,ystyle=1, $ ;,/NOLABELS,/NOAXES,/NOXTICKS,/NOYTICKS,/NOTITLE
        thick=charth*th_ind*1.2,charthick=charth,charsize=charsz,LEVELS=[12],/cont,color=250-(8*ch+ooo*127),/over
      ;0.6*max(imgm);;0.7*max(imgm)
      plots,(gcenx-256)*60.*60./512.,(gceny-256)*60.*60./512.,psym=1,thick=charth*th_ind,color=250-(8*ch+ooo*127),symsize=charsz*2
      nafik: 
      xyouts,pos1[2]-0.06,pos1[1]+0.022*(ch+2.+ooo*16), freqchstr+' MHz' ,$
        /normal,color=250-(8*ch+ooo*127),align=0,charsize=charsz/1.3,charthick=charth        
      
;      xyouts,pos1[0]+0.01+0.08*ooo,pos1[1]+0.02*(ch+2), freqchstr+' MHz' ,$
;        /normal,color=8*ch+ooo*127,align=0,charsize=charsz,charthick=charth
       
    endfor 
  endfor
  xyouts,pos1[2]-0.075,pos1[3]-0.1, 'MUSER-I-'+strmid(data_namet,8,1) ,$;+strmid(data_namet,8,1);(L-R)
    /normal,color=255,align=0,charsize=charsz/1.1,charthick=charth
  
  ;-------------------------for goes flux plot-------------------------
  restore,'/Users/xychen/Desktop/mpro/goes/goes201412170406.sav',/VERBOSE
  ;17-Dec-14 04:25:00.000 17-Dec-14 04:51:00.000 17-Dec-14 05:20:00.000
  timegoes=TARRAY+UTBASE
  timerange_goes=[anytim('2014-12-17T04:00:00'),anytim('2014-12-17T05:35:00')]
  st_time_goes=anytim('17-Dec-14 04:25:00.000')
  st_index_goes=where(abs(timegoes-st_time_goes) lt 1)
  max_time_goes=anytim('17-Dec-14 04:51:00.000')
  max_index_goes=where(abs(timegoes-max_time_goes) lt 1)
  end_time_goes=anytim('17-Dec-14 05:20:00.000')
  end_index_goes=where(abs(timegoes-end_time_goes) lt 1)


  utplot,timegoes,YCLEAN[*,0],timerange=timerange_goes,XSTYLE=1,YSTYLE=1,thick=charth*th_ind,charthick=charth,charsize=charsz $
    ,position=pos2,/nodata $
    ,/noerase,yrange=[1e-8,1e-4],/ylog,xtitle=' ',xtickname = replicate(' ',20),ytitle=' ',ytickname = replicate(' ',20)
  
  coord = CONVERT_COORD( [anytim('17-Dec-14 04:25:00.000'),anytim('17-Dec-14 04:46:00.000')],[1e-4,1e-8],/data,/to_normal)
  X_polyf = [coord[0,0],pos2[2],pos2[2],coord[0,1],coord[0,0]]
  Y_polyf = [coord[1,0],pos2[0],pos2[1],coord[1,1],coord[1,0]]
  
  X_polyf = [coord[0,0],coord[0,1],coord[0,1],coord[0,0],coord[0,0]]
  Y_polyf = [coord[1,1],coord[1,1],coord[1,0],coord[1,0],coord[1,1]]
  LOADCT,0
  POLYFILL, X_polyf, Y_polyf, COLOR=230,/NORMAL
  POLYFILL, [X_polyf[0],X_polyf[1],pos3[2],pos3[0],X_polyf[0]], [Y_polyf[0],Y_polyf[1],pos3[3],pos3[3],Y_polyf[0]], COLOR=230,/NORMAL

  utplot,timegoes,YCLEAN[*,0],timerange=timerange_goes,XSTYLE=1,YSTYLE=1,thick=charth*th_ind,charthick=charth,charsize=charsz $
    ,position=pos2,ytitle=textoidl('watts m^{-2}'),/nodata $
    ,/noerase,yrange=[1e-8,1e-4],/ylog,xtitle=' '

  oplot,[min(timegoes),max(timegoes)],[1e-5,1e-5],linestyle=1
  oplot,[min(timegoes),max(timegoes)],[1e-6,1e-6],linestyle=1
  oplot,[min(timegoes),max(timegoes)],[1e-7,1e-7],linestyle=1
  outplot,timegoes,YCLEAN[*,0],thick=charth*th_ind,color=cgColor('red')
  outplot,timegoes,YCLEAN[*,1],thick=charth*th_ind,color=cgColor('blue')
  oplot,[timegoes[st_index_goes],timegoes[st_index_goes]],[1e-8,1e-4]
  oplot,[timegoes[max_index_goes],timegoes[max_index_goes]],[1e-8,1e-4]
  oplot,[timegoes[end_index_goes],timegoes[end_index_goes]],[1e-8,1e-4]
  ;axis,yaxis=1,charthick=charth,charsize=charsz,ytickname=['A','B','C','M','X']
  xyouts,pos2[0]+0.007,pos2[3]-0.03,'GOES15 1.0-8.0 A' $;;pos2[0]+0.007
    ,/normal,color=cgColor('red'),align=0,charsize=charsz/1.2,charthick=charth/1.2
  xyouts,pos2[0]+0.007,pos2[3]-0.06,'GOES15 0.5-4.0 A' $
    ,/normal,color=cgColor('blue'),align=0,charsize=charsz/1.2,charthick=charth/1.2
;  xyouts,1.1*pos2[0],(1-0.012*4)*pos2[3],'Start: 04:25:00' $
;    ,/normal,color=cgColor('black'),align=0,charsize=charsz,charthick=charth
;  xyouts,1.1*pos2[0],(1-0.012*5)*pos2[3],'Peak: 04:51:00' $
;    ,/normal,color=cgColor('black'),align=0,charsize=charsz,charthick=charth
;  xyouts,1.1*pos2[0],(1-0.012*6)*pos2[3],'End  : 05:20:00' $
;    ,/normal,color=cgColor('black'),align=0,charsize=charsz,charthick=charth
  print,'goes..........'
 
  loadct,3
  ;-------------------------for spectrum plot and flux plot-------------------------
  spectro_plot,sp,time,freqs,position=pos3,xstyle=1+4,ystyle=1+4 $
    ,yrange=[400,2000],drange=[1.0,2.6],xrange=timerange $; $; $;;;
    ,ytitle='Frequency [MHz]',title=' ',xtitle=' ' $
    ,thick=charth*th_ind,charthick=charth,charsize=charsz,/noerase,color=-1
  axis,yaxis=0,yrange=[0.4,2],ystyle=1,charthick=charth,charsize=charsz,ytitle=textoidl('Frequency [GHz]'),color=0,ythick=charth*th_ind
  axis,yaxis=1,yrange=[0.4,2],charthick=charth,charsize=charsz,color=0,ythick=charth*th_ind,ytickname=replicate(' ',8)
  ;axis,xaxis=1,charthick=charth,charsize=charsz,xtitle=' ',color=255,xthick=charth*th_ind*1.2,xtickname=replicate(' ',8)
  ;cgcolorbar,range=[1.0,2.6],yticks = 2.,/ver,position=posc4,color=255,charsize=charsize,charthick=charthick
  ;outplot,[anytim(indexa1.date_OBS),anytim(indexa1.date_OBS)],[0,10],thick=charth*th_ind*1.2,color=255;,linestyle=1

  utplot,time,sp[*,48+4],position=pos3,xstyle=1,ystyle=1+4, $
    timerange=timerange,title=' ',ytitle='Flux [Arbitrary]',/noerase,color=255, $
    thick=charth*th_ind,charthick=charth,charsize=charsz,yrange=[0.8,2.2]
  ;axis,yaxis=0,yrange=[0.8,2.2],ystyle=1,charthick=charth,charsize=charsz,ytitle='Flux [Arbitrary]',color=255,ythick=charth*th_ind*1.2
  ;axis,yaxis=1,yrange=[0.8,2.2],charthick=charth,charsize=charsz,color=255,ythick=charth*th_ind*1.2,ytickname=replicate(' ',8)
  ;axis,xaxis=1,charthick=charth,charsize=charsz,xtitle=' ',color=255,xthick=charth*th_ind*1.2,xtickname=replicate(' ',8)
  outplot,[anytim(indexa1.date_OBS),anytim(indexa1.date_OBS)],[0.8,2.2],thick=charth*th_ind,color=255;,linestyle=1
  
  xyouts,pos3[2]-0.1,pos3[3]-0.06, 'Flux- 1.7 GHz' ,$;+strmid(data_namet,8,1);(L-R)
    /normal,color=255,align=0,charsize=charsz/1.2,charthick=charth/1.2
  xyouts,pos3[2]-0.1,pos3[3]-0.03, 'Spectrum-MUSER',$;+strmid(data_namet,8,1);(L-R)
    /normal,color=255,align=0,charsize=charsz/1.2,charthick=charth/1.2
  
  xyouts,pos1[0]+0.007,pos1[3]-0.04, '(a)',/normal,color=255,align=0,charsize=charsz/0.8,charthick=charth/1.
  xyouts,pos2[2]-0.025,pos2[3]-0.04, '(b)',/normal,color=0,align=0,charsize=charsz/0.8,charthick=charth/1.
  xyouts,pos3[0]+0.007,pos3[3]-0.04, '(c)',/normal,color=255,align=0,charsize=charsz/0.8,charthick=charth/1.

endfor

endfor
endfor
endfor

;-------------------------plot-------------------------
device,/CLOSE
set_plot,'X'
print,'It is all ready..........................'
end