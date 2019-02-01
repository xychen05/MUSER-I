;.r muser_movie_image_flux_spec_aia_wl3_norh.pro
;.compile package_muser.pro 
;.compile package_xy.pro
;.compile packageaia.pro
;
; Purpose :
; make movie of MUSER clean image changing with the flux and the spectrum.
;
; History :
; Xingyao Chen, 26 July 2017.
;-------------------------begins-------------------------

nm_ind='05'
;filename = file_search('./sav/1.6-2.0_Lnum2400/it_clean_img'+nm_ind+'_fits_raw/*.fits',count=count)
;filename = file_search('./sav/1.6-2.0_Lnum2400/it_clean_img'+nm_ind+'_fits_raw/*.fits',count=count)
filename = file_search('./sav/1.6-2.0_Lnum2400/it_clean_img'+nm_ind+'_fits_raw_LR/*.fits',count=count)

rrr='_diffrg01_raw_';;'_'
rrr='_diffrg01_'
rrr='01_';;'_'
rrr='01_raw_'

ttt='_thred2'

filepoints = file_search('./sav/1.6-2.0_Lnum2400/flux_image_dataprep'+nm_ind+'/points'+rrr+'*_thred2.sav')
fileflux = file_search('./sav/1.6-2.0_Lnum2400/flux_image_dataprep'+nm_ind+'/fluximg'+rrr+'*_thred2.sav')
restore,filepoints,/ver
restore,fileflux,/ver

restore,'/Volumes/Seagate_cxy/muserdata_20141217/muser_spec/sav/spec_time_dataLR_041955-044555_A1.sav',/ver

outfile = '/Users/xychen/Desktop/muser_movie_image_flux_spec_aia_wl3_norh.mp4'

pixel=512

;-------------------------basic parameters-------------------------
filenamea1 = file_search('/Volumes/TOSHIBA_cxy/data_20141217/AIA/prep/171/AIA20141217*.fits')
filenamea2 = file_search('/Volumes/TOSHIBA_cxy/data_20141217/AIA/prep/131/AIA20141217*.fits')
filenamea3 = file_search('/Volumes/TOSHIBA_cxy/data_20141217/AIA/prep/94/AIA20141217*.fits')
filename17 = file_search('/Users/xychen/Desktop/mpro/norh/data/ipa/*')
filename34 = file_search('/Users/xychen/Desktop/mpro/norh/data/ipz/*')

;-------------------------basic parameters-------------------------
x1 = -300. & x2 = 250 & y1 = -600. & y2 = -200. ;region 5
x11 = -300. & x21 = 250 & y11 = -600. & y21 = -200. ;region 5

xrange = [x11,x21] & yrange = [y11,y21]

pos1 = [0.10,0.1/2.,0.90,0.9/2.]
pos3= [0.10,0.73,0.90,0.98]
pos4= [0.10,0.48,0.90,0.73]

;xysz=2./3.
;pos1 = [0.12,0.15*xysz,0.92,0.55*xysz]
;posc1 = [pos1[2]+0.045,pos1[3]-0.2,pos1[2]+0.06,pos1[3]]
;pos3 = [0.12,pos1[3]+0.06,0.92,pos1[3]+0.06+0.26]
;pos4 = [0.12,pos3[3],0.92,pos3[3]+0.26]
;posc4 = [pos4[2]+0.045,pos4[3]-0.2,pos4[2]+0.06,pos4[3]]
;;framedims   = [1000,1000/xysz]

charsizep=0.9*2.5
charsizex=0.75*2.5
CHARTHICKp=1.5*2.5
CHARTHICKx=1.5*2.5
thick=1.5*2.5

systim = SYSTIME(1)
wavelength = ['94','171','193','211','304','335','131']
mag_scal = [2.5,1.5,5.,8.,3.,2.,3.0]
rs_ratio = 1.0d
xllp=0.
nxp=4096.
yllp=0.
nyp=4096.
;-------------------------movie prepared-------------------------
video_file = outfile
video = idlffvideowrite(video_file)
framerate = 12

framedims   = [550*2.8,800*2.8]

stream = video.addvideostream(framedims[0], framedims[1], framerate,BIT_RATE=5.0e6)
ENTRY_DEVICE = !D.NAME
set_plot,'z',/copy
device, set_resolution=framedims, set_pixel_depth=24, decomposed=0


;-------------------------plot-------------------------
;for norp
restore,'/Users/xychen/Desktop/mpro/norp/norp20141217_0432.xdr',/VERBOSE
;colorxy=['Red','green','Yellow','black','Violet','Pink','Brown'];,'Olive','Sky Blue'
freqs_np=['1GHz','2GHz','4GHz','9GHz','17GHz','35GHz','80GHz']
fluxnp=FI
sznp=size(fluxnp)
;17-Dec-14 04:41:28.000 17-Dec-14 04:42:58.000 17-Dec-14 05:33:36.000
st_timenp=anytim('17-Dec-14 04:03:11.000')
timesnp=st_timenp+findgen(sznp[2])*0.1

;-------------------------plot-------------------------
nt_ind=1;5*12s
file_st=900+26*5;+600
nfile=100
;nfile=2

timerange=[times[0],max(times)]

nt_ind_n=1;nt_ind*10s
file_st_n=137;6*23-1

nt_ind_m=1;nt_ind*10s
file_st_m=365

;-------------------------for aia-------------------------
nn=0
wl=1;171
wavelnth = wavelength[wl]
read_sdo, filenamea1[file_st+nt_ind*nn-3], indexa1, dataa1, /UNCOMP_DELETE,/noshell,/use_shared_lib,/HIDE
dataa1=dataa1/indexa1.EXPTIME
img1 = sdo_aia_scale_hdr(dataa1/mag_scal[wl],indexa1,xllp,yllp,nxp,nyp,wavelnth=wavelnth,rs_ratio=rs_ratio,/no_imgrscale)
index2map,indexa1,img1,map
x1pix=round(indexa1.crpix1+x11/indexa1.cdelt1-1)
x2pix=round(indexa1.crpix1+x21/indexa1.cdelt1-1)
y1pix=round(indexa1.crpix2+y11/indexa1.cdelt2-1)
y2pix=round(indexa1.crpix2+y21/indexa1.cdelt2-1)

colorxy=['yellow','red','green','black','Violet','Pink','Brown','Olive','Sky Blue']
timerange=[anytim('2014-12-17T04:26:00'),anytim('2014-12-17T04:45:00')]

;-------------------------plot-------------------------
for mm=0,nfile-1 do begin;nfile;count-1
  
  print,'For loop:',mm
  ;-------------------------nn image of AIA-------------------------
  nn=mm
  ;wavelength = ['94','171','193','211','304','335','131']
  imga = intarr(4096,4096,3)
  wl=1;171
  wavelnth = wavelength[wl]
  read_sdo, filenamea1[file_st+nt_ind*nn], indexa1, dataa1, /UNCOMP_DELETE,/noshell,/use_shared_lib,/HIDE
  dataa1=dataa1/indexa1.EXPTIME
  img1 = sdo_aia_scale_hdr(dataa1/mag_scal[wl],indexa1,xllp,yllp,nxp,nyp,wavelnth=wavelnth,rs_ratio=rs_ratio,/no_imgrscale)
  index2map,indexa1,img1,map
  
  wl=6;304
  wavelnth = wavelength[wl]
  read_sdo, filenamea2[file_st+nt_ind*nn], indexa2, dataa2, /UNCOMP_DELETE,/noshell,/use_shared_lib,/HIDE
  dataa2=dataa2/indexa2.EXPTIME
  img2 = sdo_aia_scale_hdr(dataa2/mag_scal[wl],indexa2,xllp,yllp,nxp,nyp,wavelnth=wavelnth,rs_ratio=rs_ratio,/no_imgrscale)
  
  wl=0;304
  wavelnth = wavelength[wl]
  read_sdo, filenamea3[file_st+nt_ind*nn+1], indexa3, dataa3, /UNCOMP_DELETE,/noshell,/use_shared_lib,/HIDE
  dataa3=dataa3/indexa3.EXPTIME
  img3 = sdo_aia_scale_hdr(dataa3/mag_scal[wl],indexa3,xllp,yllp,nxp,nyp,wavelnth=wavelnth,rs_ratio=rs_ratio,/no_imgrscale)
  
  ;for aia prep data
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
    scale=[indexa1.cdelt1,indexa1.cdelt2], $
    xrange=xrange,yrange=yrange, $
    xtitle='X-postion (arcseconds)', $
    ytitle='Y-position (arcseconds)', $
    position = pos1,xthick=thick,ythick=thick,CHARTHICK=CHARTHICKp,charsize=charsizep,/NOADJUST,title=' ';,/normal
  plot_map,map,grid=15,/no_data,/noerase,color=255,title=' ',CHARTHICK=CHARTHICKp, $
    xrange = xrange,yrange = yrange,position = pos1,/NOAXES,/NOLABELS,/NOXTICKS,/NOYTICKS,thick=thick+1
  xyouts,1.1*pos1[0],1.15*pos1[1], textoidl('SDO/AIA-')+strmid(string(indexa1.WAVELNTH),9,4)+'/'+ $
    strmid(string(indexa2.WAVELNTH),9,4)+'/'+strmid(string(indexa3.WAVELNTH),9,4)+' '+ $
    strmid(string(indexa2.date_OBS),11,8)+'(background)',/normal,color=255,align=0,charsize=charsizex,CHARTHICK=CHARTHICKx
  
  print,mm
  print,strmid(filenamea1[file_st+nt_ind*nn],48,23)
  print,strmid(filenamea2[file_st+nt_ind*nn],48,23)
  print,strmid(filenamea3[file_st+nt_ind*nn],48,23)
  
  ;-------------------------image of NORH-------------------------
  
  if (mm ge 0) and (mm lt 304/2.) then begin
    ii=round(mm*12./10)
    read_sdo, filename17[file_st_n+nt_ind_n*ii], index17, data17, /UNCOMP_DELETE,/noshell,/use_shared_lib,/HIDE
    index2map,index17,data17,map17
    plot_map,map17,/cont,/noerase,xstyle=1+4,ystyle=1+4, $
      xrange = xrange,yrange = yrange,title=' ', $;,/NOAXES,/NOLABELS,/NOXTICKS,/NOYTICKS
      position = pos1,LEVELS = [0.5,0.7,0.95]*max(data17),c_color =cgColor('red'),c_thick=5,thick=5.
    xyouts,1.1*pos1[0],1.75*pos1[1], 'NORH- 17GHz- R+L '+strmid(filename17[file_st_n+nt_ind_n*ii],51,2) $
      +':'+strmid(filename17[file_st_n+nt_ind_n*ii],53,2)+':'+strmid(filename17[file_st_n+nt_ind_n*ii],55,2), $
      /normal,color=cgColor('red'),align=0,charsize=charsizex,CHARTHICK=CHARTHICKx
    print,strmid(filename17[file_st_n+nt_ind_n*ii],41,20)
  endif
  if (mm ge 10.) and (mm lt 304/2.) then begin
    ii=round(mm*12./10.)
    read_sdo, filename34[file_st_n+nt_ind_n*ii], index34, data34, /UNCOMP_DELETE,/noshell,/use_shared_lib,/HIDE
    index2map,index34,data34,map34
    
    plot_map,map34,/cont,/noerase,xstyle=1+4,ystyle=1+4, $
      xrange = xrange,yrange = yrange,title=' ', $;,/NOAXES,/NOLABELS,/NOXTICKS,/NOYTICKS
      position = pos1,LEVELS = [0.5,0.7,0.95]*max(data34),c_color =cgColor('blue'),c_thick=5,thick=5.
    xyouts,1.1*pos1[0],2.05*pos1[1], 'NORH- 34GHz- R+L '+strmid(filename34[file_st_n+nt_ind_n*ii],51,2) $
      +':'+strmid(filename34[file_st_n+nt_ind_n*ii],53,2)+':'+strmid(filename34[file_st_n+nt_ind_n*ii],55,2)+' (blue contour)', $
      /normal,color=cgColor('white'),align=0,charsize=charsizex,CHARTHICK=CHARTHICKx
      
    print,strmid(filename34[file_st_n+nt_ind_n*ii],41,20)
    
  endif
  
  ;-------------------------image of MUSER-------------------------
  hh=12*mm
  read_sdo, filename[file_st_m+nt_ind_m*hh], indexm, datam, /UNCOMP_DELETE,/noshell,/use_shared_lib,/HIDE
  
  angle=9.35+(8.89-9.35)/24.*4+(8.89-9.35)/24./60.*25+(8.89-9.35)/24./60./60*(nn) ;1 min for rotating
  imgm=polax_rot(datam, angle) ; rotate after moving the solar center to the image center
  
  datamt=mv_img(datam,17,11)
  datam=datamt
  angle=19 ;1 min for rotating
  imgm=polax_rot(datam, angle) ; rotate after moving the solar center to the image center
  
  ;  thred=0.3
  ;  index_temp=where(imgm lt thred*max(imgm))
  ;  imgm[index_temp]=thred*max(imgm)
  
  index2map,indexm,imgm,mapm;
  loadct,39
  
  plot_map,mapm,/cont,/noerase,xstyle=1+4,ystyle=1+4, $
    xrange = xrange,yrange = yrange,title=' ', $;,/NOAXES,/NOLABELS,/NOXTICKS,/NOYTICKS
    position = pos1,LEVELS = [18,28,38],c_color =cgColor('green'),c_thick=5,thick=5.;LEVELS = [0.5,0.75,0.9]*max(datam)
  xyouts,1.1*pos1[0],1.45*pos1[1], 'MUSER-1.7GHz- R+L '+strmid((indexm.T_OBS),11,8) $
    ,/normal,color=cgColor('green'),align=0,charsize=charsizex,CHARTHICK=CHARTHICKx
  
 
  ;-------------------------for lines-------------------------
  ;timerange=[anytim('2014-12-17T04:00:00'),anytim('2014-12-17T05:32:00')]
  timeline=anytim(indexa2.date_OBS);timerange[0]+mm*12.*60
  ;-------------------------for goes-------------------------
  ;for goes
  restore,'/Users/xychen/Desktop/mpro/goes/goes201412170453.sav',/VERBOSE
  ;17-Dec-14 04:25:00.000 17-Dec-14 04:51:00.000 17-Dec-14 05:20:00.000
  times=TARRAY+UTBASE
  st_time=anytim('17-Dec-14 04:25:00.000')
  st_index=where(abs(times-st_time) lt 1)
  max_time=anytim('17-Dec-14 04:51:00.000')
  max_index=where(abs(times-max_time) lt 1)
  end_time=anytim('17-Dec-14 05:20:00.000')
  end_index=where(abs(times-end_time) lt 1)
  utplot,times,YCLEAN[*,0],timerange=timerange,XSTYLE=1,YSTYLE=1+8,thick=thick,charthick=CHARTHICKp,charsize=charsizep $
    ,position=pos3,ytitle=textoidl('watts m^{-2}'),/nodata,xthick=thick,ythick=thick $
    ,/noerase,yrange=[1e-8,1e-3],/ylog,xtitle=' ',xtickname = replicate(' ',20);
  oplot,[min(times),max(times)],[1e-4,1e-4],linestyle=1
  oplot,[min(times),max(times)],[1e-5,1e-5],linestyle=1
  oplot,[min(times),max(times)],[1e-6,1e-6],linestyle=1
  oplot,[min(times),max(times)],[1e-7,1e-7],linestyle=1
  outplot,times,YCLEAN[*,0],thick=thick,color=cgColor('white')
  outplot,times,YCLEAN[*,1],thick=thick,color=cgColor('sky blue')
  oplot,[times[st_index],times[st_index]],[1e-8,1e-3],linestyle=1,thick=thick*1.5
  oplot,[times[max_index],times[max_index]],[1e-8,1e-3],linestyle=1,thick=thick*1.5
  oplot,[times[end_index],times[end_index]],[1e-8,1e-3],linestyle=1,thick=thick*1.5
  ;axis,yaxis=1,charthick=1.,charsize=0.9,ytickname=['A','B','C','M','X',' ']
  xyouts,1.1*pos3[0],(1-0.013*1)*pos3[3],'GOES15 1.0-8.0 A' $
    ,/normal,color=cgColor('white'),align=0,charsize=charsizex,charthick=CHARTHICKx
  xyouts,1.1*pos3[0],(1-0.013*2)*pos3[3],'GOES15 0.5-4.0 A' $
    ,/normal,color=cgColor('sky blue'),align=0,charsize=charsizex,charthick=CHARTHICKx

  print,'goes..........'
  
  oplot,[timeline,timeline],[1e-8,1e-3],thick=thick
  
  
  ;-------------------------for norh------------------------
  ;for norp
  restore,'/Users/xychen/Desktop/mpro/norp/norp20141217_0432.xdr',/VERBOSE
  ;colorxy=['Red','green','Yellow','black','Violet','Pink','Brown'];,'Olive','Sky Blue'
  freqs=['1GHz','2GHz','4GHz','9GHz','17GHz','35GHz','80GHz']
  fluxp=FI
  sz=size(fluxp)
  ;17-Dec-14 04:41:28.000 17-Dec-14 04:42:58.000 17-Dec-14 05:33:36.000
  st_time=anytim('17-Dec-14 04:03:11.000')
  times=st_time+findgen(sz[2])*0.1
  st_index=where(abs(times-st_time) lt 1)
  max_time=anytim('17-Dec-14 04:32:08.000')
  max_index=where(abs(times-max_time) lt 1)
  end_time=anytim('17-Dec-14 05:30:02.000')
  end_index=where(abs(times-end_time) lt 1)
  
  ii=1
  utplot,times,fluxp[0,*]-fluxp[0,0],XSTYLE=1+4,YSTYLE=4,thick=thick,charthick=CHARTHICKp,charsize=charsizep $
    ,position=pos3,yrange=[-50.,450.] $
    ,timerange=timerange,/noerase,/nodata,ytitle=' ',xtitle=' ',/noaxis,color=cgColor('blue')
  axis,yaxis=1,charthick=charthickp,charsize=charsizep,ytitle=textoidl('Flux [Arbitrary]'),color=255,ythick=thick
  outplot,times,fluxp[2,*]-fluxp[2,0],thick=thick,color=cgColor(colorxy[ii])
  outplot,times,fluxp[1,*]-fluxp[1,0],thick=thick*0.7,color=cgColor('blue')
  
  oplot,[times[max_index],times[max_index]],[-100.,500.],linestyle=1,thick=thick*1.5
  xyouts,1.1*pos3[0],(1-0.013*3)*pos3[3],'NORP  - 4 GHz' $
    ,/normal,color=cgColor(colorxy[ii]),align=0,charsize=charsizex,charthick=CHARTHICKx
  xyouts,1.1*pos3[0],(1-0.013*4)*pos3[3],'NORP  - 2 GHz' $
    ,/normal,color=cgColor('blue'),align=0,charsize=charsizex,charthick=CHARTHICKx

  
  ;-------------------------for muser------------------------
  restore,'/Volumes/Seagate_cxy/muserdata_20141217/muser_spec/sav/spec_time_dataLR_041955-044555_A1.sav',/ver

  utplot,time,sp[*,48+4],position=pos3,xstyle=1+4,ystyle=1+4, $
    timerange=timerange,title=' ',ytitle=' ',/noerase,color=cgcolor('green'), $
    thick=thick/2.,xthick=charthickx,ythick=charthickx,charsize=charsizex,charthick=CHARTHICKx
  xyouts,1.1*pos3[0],(1-0.013*5)*pos3[3],'MUSER - 1.7 GHz' $
    ,/normal,color=cgcolor('green'),align=0,charsize=charsizex,charthick=CHARTHICKx
  outplot,[times[file_st+nt_ind*nn],times[file_st+nt_ind*nn]],[0,10],thick=charthickx,color=255;,linestyle=1

  ;-------------------------for spectrum------------------------
  loadct,39
  spectro_plot,sp,time,freqs,position=pos4,ytitle='Frequency [MHz]',yrange=[400,2000] $
    ,xstyle=1,ystyle=1,charsize=charsizep,thick=thick,drange=[1.0,2.6] $;,xthick=thick,ythick=thick
    ,/noerase,xrange=timerange,xtitle=' ',CHARTHICK=CHARTHICKp
  oplot,[timeline,timeline],[400,2000],thick=thick,color=255


  
  timestamp = video.put(stream, tvrd(true=1))
  ;print,filename[nn]
  print,minmax(datam)

endfor

;-------------------------ending-------------------------
close,/all
device, /close
video.cleanup
set_plot,'X'

print,'It is all ready..........................'

end