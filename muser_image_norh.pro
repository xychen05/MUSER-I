;.r muser_image_norh.pro
;.compile package_muser.pro 
;.compile package_xy.pro
;
; Purpose :
; plot image of MUSER clean image.
; 
; History :
; Xingyao Chen, 31 Oct 2017.
;-------------------------begins-------------------------
WINDOW,0,xsize=700, ysize=700
Device,RETAIN=2
;!P.FONT = 0

x1 = -800. & x2 = 200 & y1 = -1000. & y2 = 0. ;region 2

xrange = [x1,x2] & yrange = [y1,y2]

pos1 = [0.12,0.12,0.90,0.90]
posc1 = [pos1[2]+0.07,pos1[3]-0.2,pos1[2]+0.09,pos1[2]]

charsize=1.5
charthick=1.3

systim = SYSTIME(1)

;-------------------------171-------------------------
nm_ind='05'
filename = file_search('./sav/1.6-2.0_Lnum2400/it_clean_img'+nm_ind+'_fits_raw_LR/*.fits',count=count)

filename17 = file_search('/Users/xychen/Desktop/mpro/norh/data/ipa/*')
filename34 = file_search('/Users/xychen/Desktop/mpro/norh/data/ipz/*')

;-------------------------aia image-------------------------
filenamea='/Volumes/Seagate_cxy/muserdata_20141217/AIA/AIA20141217_043207_0304.fits';;;nn=725
;filenamea='/Volumes/Seagate_cxy/muserdata_20141217/AIA/AIA20141217_042407_0304.fits';;;nn=255
wavelength = ['94','171','193','211','304','335','131','1600']
mag_scal = [1.5,1.5,3,8.,3.,2.,1.5,5.]
rs_ratio = 1.0d
xllp=0.
nxp=4096.
yllp=0.
nyp=4096.
wl=4
wavelnth = wavelength[wl]

aia_lct, rr, gg, bb, wavelnth=wavelnth, /load
read_sdo, filenamea, indexa, dataa, /UNCOMP_DELETE,/noshell,/use_shared_lib,/HIDE
dataa=dataa/indexa.EXPTIME
imga = sdo_aia_scale_hdr(dataa/mag_scal[wl],indexa,xllp,yllp,nxp,nyp,wavelnth=wavelnth,rs_ratio=rs_ratio,/no_imgrscale)
index2map,indexa,imga,mapa

;-------------------------basic parameters-------------------------
pixel=512

nt_ind=1;
file_st=0;+600
nfile=count

;-------------------------Part I: select points and regions-------------------------
for nn=725,725 do begin
  
  loadct,3
  ;read_sdo, filename[file_st+nt_ind*nn], indexm, datam, /UNCOMP_DELETE,/noshell,/use_shared_lib,/HIDE
  datam = readfits( filename[file_st+nt_ind*nn], indexm, /NOSCALE, /SILENT)
  print,minmax(datam)
  
  angle=9.35+(8.89-9.35)/24.*4+(8.89-9.35)/24./60.*20+(8.89-9.35)/24./60./60*(nn) ;1 min for rotating
  imgm=polax_rot(datam, angle) ; rotate after moving the solar center to the image center

;  ;datamt=mv_img(datam,17,11);;;nn=255
;  datamt=mv_img(datam,22,12);;;nn=725
;  datam=datamt
;  angle=19 ;1 min for rotating
;  imgm=polax_rot(datam, angle) ; rotate after moving the solar center to the image center

  thred=0.15
  index_temp=where(imgm lt thred*max(imgm))
  imgm[index_temp]=thred*max(imgm)
  
  index2map,indexm,imgm,map;
  
  plot_map,map,position = pos1,xrange = xrange,yrange = yrange,$
    thick=charthick,charthick=charthick,charsize=charsize,/limb;,title=string(pixel)+'*'+string(pixel);,dmin =2,dmax =25.;,/LOG_SCALE;
  
  plot_map,mapa,position = pos1,/notitle, $
    xrange = xrange,yrange = yrange,xstyle=1,ystyle=1, $ ;,/NOLABELS,/NOAXES,/NOXTICKS,/NOYTICKS,/NOTITLE
    thick=charthick,charthick=charthick,charsize=charsize,LEVELS=[200],/cont,color=255,/over
  
  ;;-------------------------image of NORH-------------------------
  nt_ind_n=1;nt_ind*10s
  file_st_n=101
  ii=round((file_st+nt_ind*nn)/10.)
  read_sdo, filename17[file_st_n+nt_ind_n*ii], index17, data17, /UNCOMP_DELETE,/noshell,/use_shared_lib,/HIDE
  index2map,index17,data17,map17
  plot_map,map17,position = pos1,xrange = xrange,yrange = yrange,xstyle=1+4,ystyle=1+4, $
    /noerase,title=' ',/cont,LEVELS = [0.5,0.7,0.95]*max(data17),c_color =cgColor('green'), $
    thick=charthick,charthick=charthick,charsize=charsize  
  
  read_sdo, filename34[file_st_n+nt_ind_n*ii], index34, data34, /UNCOMP_DELETE,/noshell,/use_shared_lib,/HIDE
  index2map,index34,data34,map34  
  plot_map,map34,position = pos1,xrange = xrange,yrange = yrange,xstyle=1+4,ystyle=1+4, $
    /noerase,title=' ',/cont,LEVELS = [0.5,0.7,0.95]*max(data34),c_color =cgColor('blue'), $
    thick=charthick,charthick=charthick,charsize=charsize
 
  print,strmid(filename17[file_st_n+nt_ind_n*ii],41,20)
  print,strmid(filename34[file_st_n+nt_ind_n*ii],41,20)

  xyouts,pos1[0]+0.01,pos1[1]+0.06, 'NORH- 17GHz',/normal,color=cgColor('green'),align=0,charsize=charsize,charthick=charthick
  xyouts,pos1[0]+0.01,pos1[1]+0.04, 'NORH- 34GHz',/normal,color=cgColor('blue'),align=0,charsize=charsize,charthick=charthick
  xyouts,pos1[0]+0.01,pos1[1]+0.02, textoidl('AIA-304-white contour  MUSER ')+strcompress(indexm[13],/remove_all) ,$
    /normal,color=255,align=0,charsize=charsize,charthick=charthick
  loadct,3
;  xyouts,pos1[0]+0.01,pos1[1]+0.02, textoidl('MUSER-1.725GHz-')+(indexm.POLARIZA)+' '+(indexm.T_OBS) ,$
;    /normal,color=255,align=0,charsize=charsize,charthick=charthick
  cgColorbar,range=[0.5*max(imgm),max(imgm)],yticks=3,/ver,color=255,charsize=charsize-0.4,position=posc1,charthick=charthick
  
  
endfor

;-------------------------plot-------------------------
print,'It is all ready..........................'

end