;.r muser_movie_image_flux_spec.pro
;.compile package_muser.pro 
;.compile package_xy.pro
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

outfile = '/Users/xychen/Desktop/muser_image_flux_spec_L+R_test01.mp4'

pixel=512

;-------------------------basic parameters-------------------------

x1 = -800. & x2 = 400 & y1 = -650. & y2 = -50. ;region 2

xrange = [x1,x2] & yrange = [y1,y2]

;xysz=2./3.
;pos1 = [0.08,0.1*xysz,0.48,0.5*xysz]
;;posc1 = [pos1[2]+0.02,pos1[3]-0.2,pos1[2]+0.03,pos1[3]]
;pos2 = [0.56,pos1[1],0.96,pos1[3]]
;posc2 = [pos2[2]+0.02,pos2[3]-0.2,pos2[2]+0.03,pos2[3]]
;pos3 = [0.08,pos1[3]+0.06,0.96,pos1[3]+0.06+0.26]
;pos4 = [0.08,pos3[3],0.96,pos3[3]+0.26]

xysz=2./3.
pos1 = [0.12,0.15*xysz,0.92,0.55*xysz]
posc1 = [pos1[2]+0.045,pos1[3]-0.2,pos1[2]+0.06,pos1[3]]
pos3 = [0.12,pos1[3]+0.06,0.92,pos1[3]+0.06+0.26]
pos4 = [0.12,pos3[3],0.92,pos3[3]+0.26]
posc4 = [pos4[2]+0.045,pos4[3]-0.2,pos4[2]+0.06,pos4[3]]

charsize=1.8
charthick=2.6

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
nt_ind=2;5*12s
file_st=365;+600
nfile=580
;nfile=2
timerange=[times[0],max(times)]

;-------------------------plot-------------------------
for nn=0,nfile-1 do begin;nfile;count-1
  loadct,39
  ;-------------------------spectrum-------------------------
  spectro_plot,sp,time,freqs,position=pos4,xstyle=1+4,ystyle=1+4 $
    ,yrange=[400,2000],drange=[1.0,2.6],xrange=timerange $; $; $;;;
    ,ytitle='Frequency [MHz]',title=' ',xtitle=' ' $
    ,charsize=charthick,thick=charthick,xthick=charthick,ythick=charthick
  axis,yaxis=0,yrange=[0.4,2],ystyle=1,charthick=charthick,charsize=charsize,ytitle=textoidl('Frequency [GHz]'),color=255,ythick=charthick
  axis,yaxis=1,yrange=[0.4,2],charthick=charthick,charsize=charsize,color=255,ythick=charthick,ytickname=replicate(' ',8)
  axis,xaxis=1,charthick=charthick,charsize=charsize,xtitle=' ',color=255,xthick=charthick,xtickname=replicate(' ',8)
  cgcolorbar,range=[1.0,2.6],yticks = 2.,/ver,position=posc4,color=255,charsize=charsize,charthick=charthick
  
  utplot,time,sp[*,48+4],position=pos3,xstyle=1,ystyle=1, $
    timerange=timerange,title=' ',ytitle='Flux [Arbitrary]',/noerase,color=255, $
    charsize=charthick/1.3,thick=charthick,xthick=charthick,ythick=charthick
  xyouts,1.1*pos3[0],(1-0.01*3)*pos3[3],'MUSER - 1.7 GHz' $
    ,/normal,color=255,align=0,charsize=charsize,charthick=charthick  
  outplot,[times[file_st+nt_ind*nn],times[file_st+nt_ind*nn]],[0,10],thick=charthick,color=255;,linestyle=1
  
;  ;;;-------------------------flux------------------------- 
;  utplot,times,flux_total,timerange=timerange,xstyle=1+4,ystyle=1+4,thick=charthick,charthick=charthick,charsize=charsize $
;    ,position=pos3,ytitle=' ',xtitle=' ',timerange=timerange,yrange=[min(flux_rg1)-5,max(flux_rg1)+10],/noerase,color=cgColor('blue')
;;  outplot,times,flux_rg1,thick=charthick,color=cgColor('green')
;;  outplot,times,flux_rg2,thick=charthick,color=cgColor('red')
;;  
;;  outplot,[times[nn],times[nn]],[min(flux_rg1)-5,max(flux_rg1)+10],thick=charthick,color=255;,linestyle=1
;;  
  ;-------------------------flux of NORP-------------------------
  utplot,timesnp,fluxnp[0,*]-fluxnp[0,0],xstyle=1+4,ystyle=1+4,thick=charthick,charthick=charthick,charsize=charsize $
    ,position=pos3,yrange=[-50.,450.],timerange=timerange,/nodata,ytitle=' ',xtitle=' ',/noaxis,color=cgColor('blue'),/noerase
  outplot,timesnp,fluxnp[2,*]-fluxnp[2,0],thick=charthick,color=cgColor('pink')
  outplot,timesnp,fluxnp[1,*]-fluxnp[1,0],thick=charthick,color=cgColor('yellow')
  outplot,[times[file_st+nt_ind*nn],times[file_st+nt_ind*nn]],[-100,1e3],thick=charthick,color=cgColor('green');,linestyle=1
  
  xyouts,1.1*pos3[0],(1-0.01*7)*pos3[3],'NORP  - 4 GHz' $
    ,/normal,color=cgColor('pink'),align=0,charsize=charsize,charthick=charthick
  xyouts,1.1*pos3[0],(1-0.01*5)*pos3[3],'NORP  - 2 GHz' $
    ,/normal,color=cgColor('yellow'),align=0,charsize=charsize,charthick=charthick

  ;-------------------------image-------------------------

  read_sdo, filename[file_st+nt_ind*nn], indexm, datam, /UNCOMP_DELETE,/noshell,/use_shared_lib,/HIDE
  
  angle=9.35+(8.89-9.35)/24.*4+(8.89-9.35)/24./60.*25+(8.89-9.35)/24./60./60*(nn) ;1 min for rotating
  imgm=polax_rot(datam, angle) ; rotate after moving the solar center to the image center
  
  datamt=mv_img(datam,17,11)
  datam=datamt
  angle=19 ;1 min for rotating
  imgm=polax_rot(datam, angle) ; rotate after moving the solar center to the image center
  
;  thred=0.3
;  index_temp=where(imgm lt thred*max(imgm))
;  imgm[index_temp]=thred*max(imgm)
  
  index2map,indexm,imgm,map;
  loadct,39
  plot_map,map,position = pos1,xrange = xrange,yrange = yrange,xstyle=1,ystyle=1,$
    thick=charthick,charthick=charthick,charsize=charsize,/limb,title=' ',/noerase,dmin =10,dmax=50
  xyouts,pos1[0]+0.01,pos1[1]+0.02, textoidl('MUSER-L+R-1.725GHz ')+(indexm.T_OBS) ,$
    /normal,color=255,align=0,charsize=charsize+0.2,charthick=charthick-0.4
  cgColorbar,range=[10,50],yticks=3,/ver,color=255,charsize=charsize,position=posc1,charthick=charthick;;max(imgm)
  
;  plot_map,map,/cont,position = pos1,xrange = xrange,yrange = yrange, $
;    xstyle=1+4,ystyle=1+4,title=' ',/noerase, $;,/NOAXES,/NOLABELS,/NOXTICKS,/NOYTICKS
;    LEVELS = [0.5,0.7,0.9]*max(imgm),c_color =0,c_thick=charthick,thick=charthick;

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

print,'It is all ready..........................'

end