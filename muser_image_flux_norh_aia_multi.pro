;.r muser_image_flux_norh_aia_multi.pro
;.compile packageaia.pro
;.compile package_xy.pro
;.compile package_muser.pro
;set_plot,'X'
;window,0,xsize=500,ysize=500
ENTRY_DEVICE = !D.NAME
set_plot,'PS'
!P.FONT = 0
xsz=(3.0+3.0+3.0)*1.5 & ysz=(3.0+5.0)*1.5
filename_eps='/Users/xychen/Desktop/muser_image_flux_norh_aia_multi-01'+'.eps'
device,file=filename_eps,xsize=xsz,ysize=ysz,BITS_PER_PIXEL=8,$
  SET_FONT='Helvetica',/ENCAPSULATED,scale_factor=2,/color
;-------------------------basic parameters-------------------------
x11 = -300. & x21 = 0 & y11 = -490. & y21 = -190.
xrange0 = [x11,x21] & yrange0 = [y11,y21]
;-------------------------3row 5 image xsz=3+3+3 ysz3+4+4-------------------------

pos1 = [0.08, 0.98-0.27*9./8.,0.35,0.98]
pos2 = [0.40, 0.98-0.27*9./8.,0.67,0.98]
pos3 = [0.72, 0.98-0.27*9./8.,0.99,0.98]
pos4 = [0.08,0.08,0.99,0.60]

charsz=0.5 & charth=0.8
th_ind =2.

colorxy=['Black','Red','blue','Violet','purple','Yellow','green','Sky Blue','Pink','Olive','orange','white']
timerange=['17-Dec-14 04:00:00.000','17-Dec-14 04:46:00.000']


;;;-------------------------for AIA------------------------
filenamea='/Volumes/Seagate_cxy/muserdata_20141217/AIA/AIA20141217_043207_0304.fits'
filenamea_fx='/Users/xychen/Desktop/mpro/work08/sav1/aia_flux_region_multi7_norotate_304_rg01.sav'

wl=4
wavelength = ['94','171','193','211','304','335','131','1600']
mag_scal = [2.5,1.5,5.,8.,3.,2.,3.0,0.5]
wavelnth = wavelength[wl]
rs_ratio = 1.0d
xllp=0.
nxp=4096.
yllp=0.
nyp=4096.

;;;;for ribbons in detail
x1_1 = -231. & x2_1 = -220. & y1_1 = -380. & y2_1 = -365.
x1_2 = -223. & x2_2 = -210. & y1_2 = -400. & y2_2 = -380.  ;;_rg02
x1_3 = -214. & x2_3 = -203. & y1_3 = -415. & y2_3 = -400.
x1_10 = -185. & x2_10 = -172. & y1_10 = -424. & y2_10 = -410.
x1_12 = -172. & x2_12 = -154. & y1_12 = -427. & y2_12 = -414.
x1_13 = -110. & x2_13 = -95. & y1_13 = -411. & y2_13 = -390.
x1_22 = -95. & x2_22 = -85. & y1_22 = -405. & y2_22 = -385.

rg1 = [x1_1, x2_1, y1_1, y2_1] & rg1c=colorxy[1] ;sl region1
rg2 = [x1_2, x2_2, y1_2, y2_2] & rg2c=colorxy[2] ;sl region1
rg3 = [x1_3, x2_3, y1_3, y2_3] & rg3c=colorxy[3] ;sl region2
rg4 = [x1_10, x2_10, y1_10, y2_10] &   rg4c=colorxy[4] ;sl region4
rg5 = [x1_12, x2_12, y1_12, y2_12] & rg5c=colorxy[5] ;sl region5
rg6 = [x1_13, x2_13, y1_13, y2_13] & rg6c=colorxy[6] ;sl region1
rg7 = [x1_22, x2_22, y1_22, y2_22] & rg7c=colorxy[7] ;sl region1

aia_lct, rr, gg, bb, wavelnth=wavelnth, /load
read_sdo, filenamea, indexa, dataa, /UNCOMP_DELETE,/noshell,/use_shared_lib,/HIDE
dataa=dataa/indexa.EXPTIME
imga = sdo_aia_scale_hdr(dataa/mag_scal[wl],indexa,xllp,yllp,nxp,nyp,wavelnth=wavelnth,rs_ratio=rs_ratio,/no_imgrscale)
index2map,indexa,imga,mapa

plot_map,mapa,position = pos1, $
  xrange = xrange0,yrange = yrange0,$;,/NOAXES,/NOLABELS,/NOXTICKS,/NOYTICKS,xstyle=1+4,ystyle=1+4,$
  thick=charth*th_ind,charsize=charsz,charthick=charth*th_ind,/notitle,/noerase
plots,[rg1[0],rg1[1],rg1[1],rg1[0],rg1[0]],[rg1[2],rg1[2],rg1[3],rg1[3],rg1[2]],color=cgColor(rg1c),thick=charth*th_ind
;plots,[rg2[0],rg2[1],rg2[1],rg2[0],rg2[0]],[rg2[2],rg2[2],rg2[3],rg2[3],rg2[2]],color=cgColor(rg2c),thick=charth*th_ind
;plots,[rg3[0],rg3[1],rg3[1],rg3[0],rg3[0]],[rg3[2],rg3[2],rg3[3],rg3[3],rg3[2]],color=cgColor(rg3c),thick=charth*th_ind
;plots,[rg4[0],rg4[1],rg4[1],rg4[0],rg4[0]],[rg4[2],rg4[2],rg4[3],rg4[3],rg4[2]],color=cgColor(rg4c),thick=charth*th_ind
;plots,[rg5[0],rg5[1],rg5[1],rg5[0],rg5[0]],[rg5[2],rg5[2],rg5[3],rg5[3],rg5[2]],color=cgColor(rg5c),thick=charth*th_ind
plots,[rg6[0],rg6[1],rg6[1],rg6[0],rg6[0]],[rg6[2],rg6[2],rg6[3],rg6[3],rg6[2]],color=cgColor(rg6c),thick=charth*th_ind
;plots,[rg7[0],rg7[1],rg7[1],rg7[0],rg7[0]],[rg7[2],rg7[2],rg7[3],rg7[3],rg7[2]],color=cgColor(rg7c),thick=charth*th_ind

restore,filenamea_fx,/ver
flux_plot=flux_rg1nr
yrange=[min(flux_plot)-4,max(flux_plot)+4]
utplot,times,flux_plot,/nodata,yrange=yrange,ystyle=1,xstyle=1,position=pos4,timerange=timerange, $
  ytitle='Brightness (arbitrary)',thick=charth*th_ind,charthick=charth*th_ind,charsize=charsz,/noerase
  
utplot,times,flux_plot,yrange=yrange,ystyle=1+4,xstyle=1+4,position=pos4,timerange=timerange, $
  ytitle=' ',thick=charth*th_ind,charthick=charth*th_ind,charsize=charsz,/noerase,xtitle=' ',xtickname=replicate(' ',10),color=cgColor(rg1c)
outplot,[anytim('17-Dec-14 04:25:00.000'),anytim('17-Dec-14 04:25:00.000')],[yrange[0],yrange[1]],thick=charth,color=0
outplot,[anytim('17-Dec-14 04:51:00.000'),anytim('17-Dec-14 04:51:00.000')],[yrange[0],yrange[1]],thick=charth,color=0

;flux_plot=flux_rg2nr
;yrange=[min(flux_plot)-4,max(flux_plot)+4]
;utplot,times,flux_plot,yrange=yrange,ystyle=1+4,xstyle=1+4,position=pos4,timerange=timerange, $
;  ytitle=' ',thick=charth*th_ind,charthick=charth*th_ind,charsize=charsz,/noerase,xtitle=' ',xtickname=replicate(' ',10),color=cgColor(rg2c)

;flux_plot=flux_rg3nr
;yrange=[min(flux_plot)-4,max(flux_plot)+4]
;utplot,times,flux_plot,yrange=yrange,ystyle=1+4,xstyle=1+4,position=pos4,timerange=timerange, $
;  ytitle=' ',thick=charth*th_ind,charthick=charth*th_ind,charsize=charsz,/noerase,xtitle=' ',xtickname=replicate(' ',10),color=cgColor(rg3c)

;flux_plot=flux_rg4nr
;yrange=[min(flux_plot)-4,max(flux_plot)+4]
;utplot,times,flux_plot,yrange=yrange,ystyle=1+4,xstyle=1+4,position=pos4,timerange=timerange, $
;  ytitle=' ',thick=charth*th_ind,charthick=charth*th_ind,charsize=charsz,/noerase,xtitle=' ',xtickname=replicate(' ',10),color=cgColor(rg4c)

flux_plot=flux_rg6nr
yrange=[min(flux_plot)-4,max(flux_plot)+4]
utplot,times,flux_plot,yrange=yrange,ystyle=1+4,xstyle=1+4,position=pos4,timerange=timerange, $
  ytitle=' ',thick=charth*th_ind,charthick=charth*th_ind,charsize=charsz,/noerase,color=cgColor(rg6c)

;flux_plot=flux_rg7nr
;yrange=[min(flux_plot)-4,max(flux_plot)+4]
;utplot,times,flux_plot,yrange=yrange,ystyle=1+4,xstyle=1+4,position=pos4,timerange=timerange, $
;  ytitle=' ',thick=charth*th_ind,charthick=charth*th_ind,charsize=charsz,/noerase,color=cgColor(rg7c)

;;;-------------------------for NORH------------------------
filename17=file_search('/Users/xychen/Desktop/mpro/norh/data/ipa/*',count=count)
pixel=64

filename17_fx=file_search('./sav/flux_image/norh_flux01_17G_'+'*.sav')
filename34_fx=file_search('./sav/flux_image/norh_flux01_34G_'+'*.sav')

read_sdo, filename17[176], index17, data17, /UNCOMP_DELETE,/noshell,/use_shared_lib,/HIDE
img17=data17

thred=0.3
index_temp=where(img17 lt thred*max(img17))
img17[index_temp]=thred*max(img17)

index2map,index17,img17,map17
loadct,39
plot_map,map17,position = pos2,xrange = xrange0,yrange = yrange0,$
  thick=charth*th_ind,charsize=charsz,charthick=charth*th_ind,/limb,title=' ',ytitle=' ',/noerase

img17_fit=GAUSS2DFIT(img17, coeff, /TILT)
width=3.0
gauss_points, coeff,width,/ydir,lx1,lx2,ly1,ly2,lx3,lx4,ly3,ly4,lx5,lx6,ly5,ly6

plots,[lx1,lx2]*(map17.dx)+(index17.CRVAL1)-pixel*(map17.dx)/2.,[ly1,ly2]*(map17.dy)+(index17.CRVAL2)-pixel*(map17.dy)/2.,symsize=3,color=255
plots,[lx3,lx4]*(map17.dx)+(index17.CRVAL1)-pixel*(map17.dx)/2.,[ly3,ly4]*(map17.dy)+(index17.CRVAL2)-pixel*(map17.dy)/2.,symsize=3,color=255
plots,[lx5,lx6]*(map17.dx)+(index17.CRVAL1)-pixel*(map17.dx)/2.,[ly5,ly6]*(map17.dy)+(index17.CRVAL2)-pixel*(map17.dy)/2.,symsize=3,color=255

plots,[lx3,lx5]*(map17.dx)+(index17.CRVAL1)-pixel*(map17.dx)/2.,[ly3,ly5]*(map17.dy)+(index17.CRVAL2)-pixel*(map17.dy)/2.,symsize=3,color=255
plots,[lx4,lx6]*(map17.dx)+(index17.CRVAL1)-pixel*(map17.dx)/2.,[ly4,ly6]*(map17.dy)+(index17.CRVAL2)-pixel*(map17.dy)/2.,symsize=3,color=255

restore,filename17_fx[0],/ver
yrange=[min(flux_total_data)-1e5,max(flux_total_data)+4e5]
utplot,times,flux_total_data,yrange=yrange,ystyle=1+4,xstyle=1+4,position=pos4,timerange=timerange, $
  ytitle=' ',thick=charth*th_ind,charthick=charth*th_ind,charsize=charsz,/noerase,color=cgColor('orange')


;;;-------------------------for MUSER------------------------
filenamem = file_search('./sav/1.6-2.0_Lnum2400/it_clean_img05_fits_raw/*.fits',count=count)
pixel=512

filenamem_fx=file_search('./sav/1.6-2.0_Lnum2400/flux_image_dataprep05_LR/fluximg01.sav')

loadct,39
datam = readfits( filenamem[725], indexm, /NOSCALE, /SILENT)
print,minmax(datam)

;datamt=mv_img(datam,17,11);;;nn=255
datamt=mv_img(datam,22,12);;;nn=725
datam=datamt
angle=19 ;1 min for rotating
imgm=polax_rot(datam, angle) ; rotate after moving the solar center to the image center

thred=0.15
index_temp=where(imgm lt thred*max(imgm))
imgm[index_temp]=thred*max(imgm)

index2map,indexm,imgm,map;

plot_map,map,position = pos3,xrange = xrange0,yrange = yrange0,$
  thick=charth*th_ind,charsize=charsz,charthick=charth*th_ind,/limb,/noerase,title=' ',ytitle=' '


index1=where(imgm ne min(imgm))
print,'pixels of radio source:',(size(index1))[1]
locx=index1/pixel
locy=index1 mod pixel
ind_rg=[min(locx),max(locx),min(locy),max(locy)]
ind1=min(ind_rg) & ind2=max(ind_rg)
imgm_rg=imgm[ind1:ind2,ind1:ind2]
imgm_rg_fit=GAUSS2DFIT(imgm_rg, coeff, /TILT)
imgm_fit=imgm
imgm_fit[*,*]=min(imgm_rg_fit)
imgm_fit[ind1:ind2,ind1:ind2]=imgm_rg_fit

width=2.0
gauss_points, coeff,width,/xdir,lx1,lx2,ly1,ly2,lx3,lx4,ly3,ly4,lx5,lx6,ly5,ly6

lx1=lx1+ind1 & ly1=ly1+ind1 & lx2=lx2+ind1 & ly2=ly2+ind1
lx3=lx3+ind1 & ly3=ly3+ind1 & lx4=lx4+ind1 & ly4=ly4+ind1
lx5=lx5+ind1 & ly5=ly5+ind1 & lx6=lx6+ind1 & ly6=ly6+ind1

plots,[lx1,lx2]*(map.dx)-pixel*(map.dx)/2.,[ly1,ly2]*(map.dy)-pixel*(map.dy)/2.,symsize=3,color=255
plots,[lx3,lx4]*(map.dx)-pixel*(map.dx)/2.,[ly3,ly4]*(map.dy)-pixel*(map.dy)/2.,symsize=3,color=255
plots,[lx5,lx6]*(map.dx)-pixel*(map.dx)/2.,[ly5,ly6]*(map.dy)-pixel*(map.dy)/2.,symsize=3,color=255

plots,[lx3,lx5]*(map.dx)-pixel*(map.dx)/2.,[ly3,ly5]*(map.dy)-pixel*(map.dy)/2.,symsize=3,color=255
plots,[lx4,lx6]*(map.dx)-pixel*(map.dx)/2.,[ly4,ly6]*(map.dy)-pixel*(map.dy)/2.,symsize=3,color=255

restore,filenamem_fx,/ver
yrange=[min(flux_total_data)-0.5,max(flux_total_data)+2.2]
utplot,times,flux_total_data,yrange=yrange,ystyle=1+4,xstyle=1+4,position=pos4,timerange=timerange, $
  ytitle=' ',thick=charth*th_ind,charthick=charth*th_ind,charsize=charsz,/noerase,color=cgColor('blue')


xyouts,pos1[0]+0.008, pos1[1]+0.01, textoidl('SDO/AIA-304 04:32'),$
  /normal,color=255,align=0,charsize=charsz,charthick=charth*th_ind
xyouts,pos2[0]+0.008, pos2[1]+0.01, textoidl('NORH-17GHz 04:32'),$
  /normal,color=0,align=0,charsize=charsz,charthick=charth*th_ind
xyouts,pos3[0]+0.008, pos3[1]+0.01, textoidl('MUSER-1.7GHz 04:32'),$
  /normal,color=255,align=0,charsize=charsz,charthick=charth*th_ind

xyouts,pos4[0]+0.008, pos4[3]-0.025, textoidl('AIA- region 01'),$
  /normal,align=0,charsize=charsz*1.3,charthick=charth*th_ind,color=cgColor(rg1c)
xyouts,pos4[0]+0.008, pos4[3]-0.05, textoidl('AIA- region 02'),$
  /normal,align=0,charsize=charsz*1.3,charthick=charth*th_ind,color=cgColor(rg6c)
xyouts,pos4[0]+0.008, pos4[3]-0.075, textoidl('NORH-  17GHz'),$
  /normal,align=0,charsize=charsz*1.3,charthick=charth*th_ind,color=cgColor('brown')
xyouts,pos4[0]+0.008, pos4[3]-0.1, textoidl('MUSER-1.7GHz'),$
  /normal,align=0,charsize=charsz*1.3,charthick=charth*th_ind,color=cgColor('blue')
   
;-------------------------ending flux for 34GHz-------------------------
device,/CLOSE
set_plot,'X'
print,'It is all ready..........................'
end