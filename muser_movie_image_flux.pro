;.r muser_movie_image_flux.pro
;.compile package_muser.pro 
;.compile package_xy.pro
;
; Purpose :
; make movie of MUSER clean image changing with the flux.
;
; History :
; Xingyao Chen, 26 July 2017.
;-------------------------begins-------------------------

nm_ind='05'
filename = file_search('./sav/1.6-2.0_Lnum2400/it_clean_img'+nm_ind+'_fits_raw/*.fits',count=count)
;filename = file_search('./sav/1.6-2.0_Lnum2400/it_clean_img'+nm_ind+'_fits_raw/*.fits',count=count)
;filename = file_search('./sav/1.6-2.0_Lnum2400/it_clean_img'+nm_ind+'_fits_raw_LR/*.fits',count=count)

rrr='_diffrg01_raw_';;'_'
rrr='_diffrg01_'
rrr='01_';;'_'
rrr='01_raw_'

filepoints = file_search('./sav/1.6-2.0_Lnum2400/flux_image_dataprep'+nm_ind+'_L/points'+rrr+'*.sav')
fileflux = file_search('./sav/1.6-2.0_Lnum2400/flux_image_dataprep'+nm_ind+'_L/fluximg'+rrr+'*.sav')

restore,filepoints,/ver
restore,fileflux,/ver  

outfile = '/Users/xychen/Desktop/muser1217_image_flux_L+R_test01.mp4'

pixel=512

;-------------------------basic parameters-------------------------
x1 = -1200. & x2 = 0 & y1 = -1200. & y2 = 0. ;region 2

xrange = [x1,x2] & yrange = [y1,y2]

xysz=3./4.

pos1 = [0.12,0.12*xysz,0.90,0.90*xysz]
posc1 = [pos1[2]+0.07,pos1[3]-0.2,pos1[2]+0.09,pos1[3]]

pos2 = [0.12,pos1[3]+0.06,0.90,0.95]

charsize=1.7
charthick=2.8

systim = SYSTIME(1)

;-------------------------movie prepared-------------------------
video_file = outfile
video = idlffvideowrite(video_file)
framerate = 12

framedims   = [1000,1000/xysz]

stream = video.addvideostream(framedims[0], framedims[1], framerate,BIT_RATE=5.0e6)
ENTRY_DEVICE = !D.NAME
set_plot,'z',/copy
device, set_resolution=framedims, set_pixel_depth=24, decomposed=0

;-------------------------plot-------------------------
nt_ind=1;5*12s
file_st=0;+600

timerange=[times[0],max(times)]

loadct,39
for nn=0,count-2 do begin;nfile;count-1
  
  ;-------------------------flux-------------------------
  utplot,times,flux_total,timerange=timerange,xstyle=1,ystyle=1,thick=charthick,charthick=charthick,charsize=charsize $
    ,position=pos2,ytitle=textoidl('Flux [Arbitrary]'),yrange=[min(flux_rg1)-5,max(flux_rg1)+10]
  outplot,times,flux_rg1,thick=charthick,color=cgColor('green')
  outplot,times,flux_rg2,thick=charthick,color=cgColor('red')
  
  outplot,[times[nn],times[nn]],[min(flux_rg1)-5,max(flux_rg1)+10],thick=charthick,color=255;,linestyle=1
  
  ;-------------------------image-------------------------
  read_sdo, filename[file_st+nt_ind*nn], indexm, datam, /UNCOMP_DELETE,/noshell,/use_shared_lib,/HIDE
  
  angle=9.35+(8.89-9.35)/24.*4+(8.89-9.35)/24./60.*25+(8.89-9.35)/24./60./60*(nn) ;1 min for rotating
  imgm=polax_rot(datam, angle) ; rotate after moving the solar center to the image center
  
  thred=0.3
  index_temp=where(imgm lt thred*max(imgm))
  imgm[index_temp]=thred*max(imgm)
  
  index2map,indexm,imgm,map;
  
  plot_map,map,position = pos1,xrange = xrange,yrange = yrange,xstyle=1,ystyle=1,$
    thick=charthick,charthick=charthick,charsize=charsize,/limb,title=' ',/noerase;,dmin =2,dmax =25.;,/LOG_SCALE;
  xyouts,pos1[0]+0.01,pos1[1]+0.02, textoidl('MUSER-L+R-1.725GHz ')+(indexm.T_OBS) ,$
    /normal,color=255,align=0,charsize=charsize+0.2,charthick=charthick-0.4
  cgColorbar,range=[2.0,max(imgm)],yticks=3,/ver,color=255,charsize=charsize,position=posc1,charthick=charthick
  
  plot_map,map,/cont,position = pos1,xrange = xrange,yrange = yrange, $
    xstyle=1+4,ystyle=1+4,title=' ',/noerase, $;,/NOAXES,/NOLABELS,/NOXTICKS,/NOYTICKS
    LEVELS = [0.5,0.7,0.9]*max(imgm),c_color =0,c_thick=charthick,thick=charthick;

;  ;-------------------------gaussian fit of selected region-------------------------
;  index1=where(imgm ne min(imgm))
;  locx=index1/1024
;  locy=index1 mod 1024
;  ind_rg=[min(locx),max(locx),min(locy),max(locy)]
;  ind1=min(ind_rg) & ind2=max(ind_rg)
;  imgm_rg=imgm[ind1:ind2,ind1:ind2]
;  imgm_rg_fit=GAUSS2DFIT(imgm_rg, coeff, /TILT)
;  imgm_fit=imgm
;  imgm_fit[*,*]=min(imgm_rg_fit)
;  imgm_fit[ind1:ind2,ind1:ind2]=imgm_rg_fit
;  
;  ;-------------------------definition of 6 points-------------------------
;  width=2.0
;  gauss_points, coeff,width,/xdir,lx1,lx2,ly1,ly2,lx3,lx4,ly3,ly4,lx5,lx6,ly5,ly6
;  
;  lx1=lx1+ind1 & ly1=ly1+ind1 & lx2=lx2+ind1 & ly2=ly2+ind1
;  lx3=lx3+ind1 & ly3=ly3+ind1 & lx4=lx4+ind1 & ly4=ly4+ind1
;  lx5=lx5+ind1 & ly5=ly5+ind1 & lx6=lx6+ind1 & ly6=ly6+ind1
;  
  
  plots,[lx1,lx2]*(map.dx)-pixel*(map.dx)/2.,[ly1,ly2]*(map.dy)-pixel*(map.dy)/2.,symsize=3,color=255,thick=3
  plots,[lx3,lx4]*(map.dx)-pixel*(map.dx)/2.,[ly3,ly4]*(map.dy)-pixel*(map.dy)/2.,symsize=3,color=255,thick=3
  plots,[lx5,lx6]*(map.dx)-pixel*(map.dx)/2.,[ly5,ly6]*(map.dy)-pixel*(map.dy)/2.,symsize=3,color=255,thick=3
  
  plots,[lx3,lx5]*(map.dx)-pixel*(map.dx)/2.,[ly3,ly5]*(map.dy)-pixel*(map.dy)/2.,symsize=3,color=255,thick=3
  plots,[lx4,lx6]*(map.dx)-pixel*(map.dx)/2.,[ly4,ly6]*(map.dy)-pixel*(map.dy)/2.,symsize=3,color=255,thick=3

  
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