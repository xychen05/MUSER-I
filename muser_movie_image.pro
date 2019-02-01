;.r muser_movie_image.pro
;.compile package_muser.pro 
;.compile package_xy.pro
;
; Purpose :
; make movie of MUSER clean image.
;
; History :
; Xingyao Chen, July 2017.
;-------------------------begins-------------------------

nm_ind='17'
ch=4
for ch=4,4 do begin
nm_ind='17sf'+strcompress(ch,/remove_all)

;;;============
filename = file_search('./sav/1.6-2.0_Rnum2400/it_clean_img'+nm_ind+'_fits_raw/*.fits',count=count)
;filename = file_search('./sav/1.6-2.0_Rnum2400/it_clean_img'+nm_ind+'_fits_raw/*.fits',count=count)
;filename = file_search('./sav/1.6-2.0_Lnum2400/it_clean_img'+nm_ind+'_fits_raw_LR/*.fits',count=count)
;;;============

freqch=(1600+ch*25+25)
freqchstr=strcompress(freqch,/remove_all)

;filename = file_search('./sav/1.6-2.0_Lnum2400/quiet_clean_img_fits/'+'*'+'_'+freqchstr+'.fits',count=count)
;print,freqchstr


rrr='_diffrg01_raw_';;'_'
rrr='_diffrg01_'
rrr='01_';;'_'
rrr='01_raw_'

ttt='_thred2'

;filepoints = file_search('./sav/1.6-2.0_Lnum2400/flux_image_dataprep'+nm_ind+'/points'+rrr+'*_thred2.sav')
;fileflux = file_search('./sav/1.6-2.0_Lnum2400/flux_image_dataprep'+nm_ind+'/fluximg'+rrr+'*_thred2.sav')
;restore,filepoints,/ver
;restore,fileflux,/ver  

outfile = '/Users/xychen/Desktop/muser_image_R_'+freqchstr+'_fix_selfcal.mp4';movie/;;selfcal;satellite

pixel=512

;-------------------------aia image-------------------------
filenamea='/Volumes/Seagate_cxy/muserdata_20141217/AIA/AIA20141217_041955_0304.fits'
wavelength = ['94','171','193','211','304','335','131','1600']
mag_scal = [1.5,1.5,3,8.,3.,2.,1.5,5.]
rs_ratio = 1.0d
xllp=0.
nxp=4096.
yllp=0.
nyp=4096.
wl=4
wavelnth = wavelength[wl]

;read_sdo, filenamea, indexa, dataa, /UNCOMP_DELETE,/noshell,/use_shared_lib,/HIDE
;dataa=dataa/indexa.EXPTIME
;imga = sdo_aia_scale_hdr(dataa/mag_scal[wl],indexa,xllp,yllp,nxp,nyp,wavelnth=wavelnth,rs_ratio=rs_ratio,/no_imgrscale)
;index2map,indexa,imga,mapa

;-------------------------basic parameters-------------------------
x1 = -1000. & x2 = 200 & y1 = -1200. & y2 = 0. ;region 2
x1 = -1200. & x2 = 1200 & y1 = -1200. & y2 = 1200. ;region 2

xrange = [x1,x2] & yrange = [y1,y2]

pos1 = [0.15,0.15,0.90,0.90]
posc1 = [pos1[2]+0.07,pos1[3]-0.2,pos1[2]+0.09,pos1[2]]

charsize=1.8
charthick=1.5

systim = SYSTIME(1)

;-------------------------movie prepared-------------------------
video_file = outfile
video = idlffvideowrite(video_file)
framerate = 12

framedims   = [1000,1000]

stream = video.addvideostream(framedims[0], framedims[1], framerate,BIT_RATE=5.0e6)
ENTRY_DEVICE = !D.NAME
set_plot,'z',/copy
device, set_resolution=framedims, set_pixel_depth=24, decomposed=0

;-------------------------plot-------------------------
nt_ind=1;5*12s
file_st=0;+600
nfile=count
;nfile=1
;
;-------------------------plot-------------------------
loadct,39
for nn=0,nfile-1 do begin;nfile;count-1

  ;read_sdo, filename[file_st+nt_ind*nn], indexm, datam, /UNCOMP_DELETE,/noshell,/use_shared_lib,/HIDE
  datam = readfits( filename[file_st+nt_ind*nn], indexm, /NOSCALE, /SILENT)
  print,filename[file_st+nt_ind*nn]
  
  angle=9.35+(8.89-9.35)/24.*4+(8.89-9.35)/24./60.*25+(8.89-9.35)/24./60./60*(nn) ;1 min for rotating
  imgm=polax_rot(datam, angle) ; rotate after moving the solar center to the image center

;  datamt=mv_img(datam,26,16)  ;;for 1.7 GHz nm_ind=05
;;  datamt=mv_img(datam,17,11)
;  angle=19 ;1 min for rotating

;  datamt=mv_img(datam,15,-30)  ;;for 1.7 GHz nm_ind=17
;  angle=19 ;1 min for rotating
; 
  datamt=mv_img(datam,-15,-42)  
  angle=9.3 ;1 min for rotating
  datamt=mv_img(datam,-16,-52)  ;;for self calibration
  angle=19 ;1 min for rotating
  
  dmin=4 & dmax=15 ;;for ch=0
  datam=datamt
  imgm=polax_rot(datam, angle) ; rotate after moving the solar center to the image center

;  thred=0.3
;  index_temp=where(imgm lt thred*max(imgm))
;  imgm[index_temp]=thred*max(imgm)
  
;  index_temp=where(imgm lt 5)
;  imgm[index_temp]=5
  
  index2map,indexm,imgm,map;bytscl((imgm),min=min(imgm))
  
  plot_map,map,position = pos1,xrange = xrange,yrange = yrange,$
    thick=2.4,charthick=charthick,charsize=charsize,/limb,dmin =dmin,dmax=dmax;,dmin =5,dmax =65.;,/LOG_SCALE;
  
;  plot_map,mapa,position = pos1,/notitle, $
;    xrange = xrange,yrange = yrange,xstyle=1,ystyle=1, $ ;,/NOLABELS,/NOAXES,/NOXTICKS,/NOYTICKS,/NOTITLE
;    thick=charthick,charthick=charthick,charsize=charsize,LEVELS=[200],/cont,color=255,/over 
  xyouts,pos1[0]+0.01,pos1[1]+0.02, 'MUSER-'+strcompress(indexm[13],/remove_all)+freqchstr+' MHz', $+textoidl('AIA-304-white contour') ,$
    /normal,color=255,align=0,charsize=charsize,charthick=charthick
  cgColorbar,range=[dmin,dmax],yticks=3,/ver,color=255,charsize=charsize,position=posc1,charthick=charsize
   
;  plots,[lx1,lx2]*(map.dx)-pixel*(map.dx)/2.,[ly1,ly2]*(map.dy)-pixel*(map.dy)/2.,symsize=3,color=255,thick=3
;  plots,[lx3,lx4]*(map.dx)-pixel*(map.dx)/2.,[ly3,ly4]*(map.dy)-pixel*(map.dy)/2.,symsize=3,color=255,thick=3
;  plots,[lx5,lx6]*(map.dx)-pixel*(map.dx)/2.,[ly5,ly6]*(map.dy)-pixel*(map.dy)/2.,symsize=3,color=255,thick=3
;  
;  plots,[lx3,lx5]*(map.dx)-pixel*(map.dx)/2.,[ly3,ly5]*(map.dy)-pixel*(map.dy)/2.,symsize=3,color=255,thick=3
;  plots,[lx4,lx6]*(map.dx)-pixel*(map.dx)/2.,[ly4,ly6]*(map.dy)-pixel*(map.dy)/2.,symsize=3,color=255,thick=3

  timestamp = video.put(stream, tvrd(true=1))
    
  ;print,filename[nn]
  print,minmax(datam)

endfor

;-------------------------ending-------------------------
close,/all
device, /close
video.cleanup
set_plot,'X'

endfor

print,'It is all ready..........................'

end