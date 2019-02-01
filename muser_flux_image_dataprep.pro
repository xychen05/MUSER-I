;.r muser_flux_image_dataprep.pro
;.compile package_muser.pro 
;.compile package_xy.pro
;
; Purpose :
; make image of MUSER clean image.
; preparation for the flux of selected region.
;
; History :
; Xingyao Chen, 25 July 2017.
;-------------------------begins-------------------------
ENTRY_DEVICE = !D.NAME
set_plot,'X'
Device,RETAIN=2

WINDOW,0,xsize=700, ysize=700
WINDOW,2,xsize=1100, ysize=700

x1 = -1000. & x2 = 200 & y1 = -1200. & y2 = 0. ;region 2

xrange = [x1,x2] & yrange = [y1,y2]

pos1 = [0.15,0.15,0.90,0.90]
posc1 = [pos1[2]+0.07,pos1[3]-0.2,pos1[2]+0.09,pos1[2]]

charsize=1.8
charthick=1.5

systim = SYSTIME(1)

;-------------------------171-------------------------
nm_ind='17'
filename = file_search('./sav/1.6-2.0_Lnum2400/it_clean_img'+nm_ind+'_fits_raw/*.fits',count=count)
data_sav='./sav/1.6-2.0_Lnum2400/flux_image_dataprep'+nm_ind+'/'
image_sav='./imagetest/1.6-2.0_Lnum2400/flux_image_dataprep'+nm_ind+'/'

;filename = file_search('./sav/1.6-2.0_Rnum2400/it_clean_img'+nm_ind+'_fits_raw/*.fits',count=count)
;data_sav='./sav/1.6-2.0_Rnum2400/flux_image_dataprep'+nm_ind+'/'
;image_sav='./imagetest/1.6-2.0_Rnum2400/flux_image_dataprep'+nm_ind+'/'

;filename = file_search('./sav/1.6-2.0_Lnum2400/it_clean_img'+nm_ind+'_fits_raw_LR/*.fits',count=count)
;data_sav='./sav/1.6-2.0_Lnum2400/flux_image_dataprep'+nm_ind+'_LR/'
;image_sav='./imagetest/1.6-2.0_Lnum2400/flux_image_dataprep'+nm_ind+'_LR/'


file_mkdir,data_sav
file_mkdir,image_sav

rrr='01_';;'_'
rrr='01_raw_'

ttt='_thred2'

pixel=512

;-------------------------basic parameters-------------------------
nt_ind=1;
file_st=0;+600
nfile=count

flux_rg1=make_array(nfile)
flux_rg2=make_array(nfile)
flux_total=make_array(nfile)

flux_rg1_data=make_array(nfile)
flux_rg2_data=make_array(nfile)
flux_total_data=make_array(nfile)

wset,0
wshow,0
loadct,39

;-------------------------Part I: select points and regions-------------------------
for nn=1,1 do begin
  
  ;read_sdo, filename[file_st+nt_ind*nn], indexm, datam, /UNCOMP_DELETE,/noshell,/use_shared_lib,/HIDE
  
  datam = readfits( filename[file_st+nt_ind*nn], indexm, /NOSCALE, /SILENT)
  ;tmp0=rstrpos(filename[file_st+nt_ind*nn],'MUSER_')
  ;file=strmid(filename[file_st+nt_ind*nn],tmp0+12,19)
  
  print,minmax(datam)
  
;  angle=9.35+(8.89-9.35)/24.*4+(8.89-9.35)/24./60.*25+(8.89-9.35)/24./60./60*(nn) ;1 min for rotating
;  imgm=polax_rot(datam, angle) ; rotate after moving the solar center to the image center
;  
;  ;datamt=mv_img(datam,26,16)
;  datamt=mv_img(datam,17,11)
;  datam=datamt
;  angle=19 ;1 min for rotating
;  imgm=polax_rot(datam, angle) ; rotate after moving the solar center to the image center

  angle=15. ;1 min for rotating
  imgm=polax_rot(datam, angle) ; rotate after moving the solar center to the image center


  thred=0.0
  index_temp=where(imgm lt thred*max(imgm))
  imgm[index_temp]=thred*max(imgm)
  
  index2map,indexm,imgm,map;
  
  plot_map,map,position = pos1,xrange = xrange,yrange = yrange,$
    thick=charthick,charthick=charthick,charsize=charsize,/limb,title=' ';,dmin =2,dmax =25.;,/LOG_SCALE;
  xyouts,pos1[0]+0.01,pos1[1]+0.02, textoidl('MUSER-1.725GHz ')+(map.TIME)+' '+strcompress(indexm[13],/remove_all) ,$
    /normal,color=255,align=0,charsize=charsize/1.25,charthick=charthick
  cgColorbar,range=[0.5*max(imgm),max(imgm)],yticks=3,/ver,color=255,charsize=charsize-0.4,position=posc1,charthick=charthick
  
  ;-------------------------gaussian fit of selected region-------------------------
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
  
  ;-------------------------definition of 6 points-------------------------
  width=2.0
  gauss_points, coeff,width,/xdir,lx1,lx2,ly1,ly2,lx3,lx4,ly3,ly4,lx5,lx6,ly5,ly6

;  theta=coeff[6]
;  width=2. ;;you can change here.
;  widthx=width*coeff[3]*sin(theta)
;  widthy=width*coeff[3]*cos(theta)
;  lx1=coeff[4]-widthx & ly1=coeff[5]-widthy
;  lx2=coeff[4]+widthx & ly2=coeff[5]+widthy
;  lx3=coeff[4]-width*coeff[2]*cos(-theta)+widthx & ly3=coeff[5]-width*coeff[2]*sin(-theta)+widthy
;  lx4=coeff[4]+width*coeff[2]*cos(-theta)+widthx & ly4=coeff[5]+width*coeff[2]*sin(-theta)+widthy
;  lx5=coeff[4]-width*coeff[2]*cos(-theta)-widthx & ly5=coeff[5]-width*coeff[2]*sin(-theta)-widthy
;  lx6=coeff[4]+width*coeff[2]*cos(-theta)-widthx & ly6=coeff[5]+width*coeff[2]*sin(-theta)-widthy

  lx1=lx1+ind1 & ly1=ly1+ind1 & lx2=lx2+ind1 & ly2=ly2+ind1
  lx3=lx3+ind1 & ly3=ly3+ind1 & lx4=lx4+ind1 & ly4=ly4+ind1
  lx5=lx5+ind1 & ly5=ly5+ind1 & lx6=lx6+ind1 & ly6=ly6+ind1
  
  plots,[lx1,lx2]*(map.dx)-pixel*(map.dx)/2.,[ly1,ly2]*(map.dy)-pixel*(map.dy)/2.,symsize=3,color=255
  plots,[lx3,lx4]*(map.dx)-pixel*(map.dx)/2.,[ly3,ly4]*(map.dy)-pixel*(map.dy)/2.,symsize=3,color=255
  plots,[lx5,lx6]*(map.dx)-pixel*(map.dx)/2.,[ly5,ly6]*(map.dy)-pixel*(map.dy)/2.,symsize=3,color=255
  
  plots,[lx3,lx5]*(map.dx)-pixel*(map.dx)/2.,[ly3,ly5]*(map.dy)-pixel*(map.dy)/2.,symsize=3,color=255
  plots,[lx4,lx6]*(map.dx)-pixel*(map.dx)/2.,[ly4,ly6]*(map.dy)-pixel*(map.dy)/2.,symsize=3,color=255

;  plot_image,imgm_fit;,xrange=[200,500],yrange=[200,500]
;  plots,[lx1,lx2],[ly1,ly2],symsize=3,color=255
;  plots,[lx3,lx4],[ly3,ly4],symsize=3,color=255
;  plots,[lx5,lx6],[ly5,ly6],symsize=3,color=255
  
  ;-------------------------2 regions on the left and right-------------------------
  region1=sl_region_index([lx1,ly1],[lx2,ly2],[lx3,ly3],[lx5,ly5])
  region2=sl_region_index([lx1,ly1],[lx2,ly2],[lx4,ly4],[lx6,ly6])
  regiont=sl_region_index([lx3,ly3],[lx5,ly5],[lx4,ly4],[lx6,ly6])
;  for ii=0,(size(region1))[2]-1 do begin
;    plots,[region1[0,ii]],[region1[1,ii]],psym=1,symsize=2,color=255
;  endfor
  
  ;time_string=indexm.T_STR
  time_string=time_obs_to_string(map.TIME)
  write_jpeg,image_sav+'img'+rrr+time_string+ttt+'.jpg',tvrd(true=3),true=3
  
endfor
save,lx1,lx2,ly1,ly2,lx3,lx4,ly3,ly4,lx5,lx6,ly5,ly6,filename=data_sav+'points'+rrr+time_string+ttt+'.sav'
save,region1,region2,regiont,filename=data_sav+'region'+rrr+time_string+ttt+'.sav'


;-------------------------Part II: flux data from the image in the two regions-------------------------
;restore,'./sav/'+'region_raw_20141217_043300_527.sav'
;nfile=10
times=dblarr(nfile)
print,nfile
loadct,39
for nn=0,nfile-1 do begin;nfile;count-1
  ;read_sdo, filename[file_st+nt_ind*nn], indexm, datam, /UNCOMP_DELETE,/noshell,/use_shared_lib,/HIDE
  datam = readfits( filename[file_st+nt_ind*nn], indexm, /NOSCALE, /SILENT)
  
;  angle=9.35+(8.89-9.35)/24.*4+(8.89-9.35)/24./60.*25+(8.89-9.35)/24./60./60*(nn) ;1 min for rotating
;  imgm=polax_rot(datam, angle) ; rotate after moving the solar center to the image center
;  
;  ;datamt=mv_img(datam,26,16)
;  datamt=mv_img(datam,17,11)
;  datam=datamt
;  angle=19 ;1 min for rotating
;  imgm=polax_rot(datam, angle) ; rotate after moving the solar center to the image center
;  datam=imgm
  
  angle=15. ;1 min for rotating
  imgm=polax_rot(datam, angle) ; rotate after moving the solar center to the image center  
  
  thred=0.0
  index_temp=where(imgm lt thred*max(imgm))
  imgm[index_temp]=thred*max(imgm)

  index2map,indexm,imgm,map;
  
;  plot_map,map,position = pos1,xrange = xrange,yrange = yrange,$
;    thick=charthick,charthick=charthick,charsize=charsize,/limb,title=' ';,dmin =2,dmax =25.;,/LOG_SCALE;
;  xyouts,pos1[0]+0.01,pos1[1]+0.02, textoidl('MUSER-L+R-1.725GHz ')+(indexm.T_OBS) ,$
;    /normal,color=255,align=0,charsize=chsize,charthick=charthick
;  cgColorbar,range=[0.5*max(imgm),max(imgm)],yticks=3,/ver,color=255,charsize=1.4,position=posc1,charthick=charthick
;  
  
  ;-------------------------flux data from the image-------------------------
  
  flux_rg1_total=make_array((size(region1))[2])
  flux_rg1_total_data=make_array((size(region1))[2])
  for ii=0,(size(region1))[2]-1 do begin
    flux_rg1_total[ii]=imgm[region1[0,ii],region1[1,ii]]
    flux_rg1_total_data[ii]=datam[region1[0,ii],region1[1,ii]]
    ;plots,[region1[0,ii]],[region1[1,ii]],psym=1,symsize=2,color=255
  endfor
  flux_rg2_total=make_array((size(region2))[2])
  flux_rg2_total_data=make_array((size(region2))[2])
  for ii=0,(size(region2))[2]-1 do begin
    flux_rg2_total[ii]=imgm[region2[0,ii],region2[1,ii]]
    flux_rg2_total_data[ii]=datam[region2[0,ii],region2[1,ii]]
    ;plots,[region2[0,ii]],[region2[1,ii]],psym=1,symsize=2,color=255
  endfor
  flux_tl=make_array((size(regiont))[2])
  flux_tl_data=make_array((size(regiont))[2])
  for ii=0,(size(regiont))[2]-1 do begin
    flux_tl[ii]=imgm[regiont[0,ii],regiont[1,ii]]
    flux_tl_data[ii]=datam[regiont[0,ii],regiont[1,ii]]
    ;plots,[regiont[0,ii]],[regiont[1,ii]],psym=1,symsize=2,color=255
  endfor
  
  flux_rg1[nn]=mean(flux_rg1_total)
  flux_rg2[nn]=mean(flux_rg2_total)
  
  temp=where(flux_tl ne 0)
  flux_total[nn]=mean(flux_tl[temp])
  
  flux_rg1_data[nn]=mean(flux_rg1_total_data)
  flux_rg2_data[nn]=mean(flux_rg2_total_data)
  flux_total_data[nn]=mean(flux_tl_data)
  
  times[nn]=anytim(map.TIME)
  
  if nn eq 0 then time_st=time_obs_to_string(map.TIME)
  if nn eq nfile-1 then time_ed=time_obs_to_string(map.TIME)
  
endfor

print,time_st
print,time_ed
save,flux_rg1,flux_rg2,flux_total,flux_rg1_data,flux_rg2_data,flux_total_data,times $
  ,filename=data_sav+'fluximg'+rrr+time_string+'-'+time_st+'-'+time_ed+ttt+'.sav'

;-------------------------plot flux curve-------------------------
timerange=[times[0],times[nfile-1]]

flux_rg1_data=flux_rg1
flux_rg2_data=flux_rg2
flux_total_data=flux_total

wset,2
wshow,2
rgmax=max([max(flux_rg1_data),max(flux_rg2_data),max(flux_total_data)])
rgmin=min([min(flux_rg1_data),min(flux_rg2_data),min(flux_total_data)])
;;restore,'./'+'sav/'+'flux_'+rrr+time_st+'-'+time_ed+'.sav'
utplot,times,flux_total_data,timerange=timerange,xstyle=1,ystyle=1,thick=chthick,charthick=charthick,charsize=charsize $
  ,position=pos1,ytitle=textoidl('Average Flux [Arbitrary]'),yrange=[rgmin,rgmax],title=' '
outplot,times,flux_rg1_data,thick=charthick,color=cgColor('green')
outplot,times,flux_rg2_data,thick=charthick,color=cgColor('red')
;
;outplot,times,flux_total_data,thick=charthick,color=cgColor('yellow')
;outplot,times,flux_rg1_data,thick=charthick,color=cgColor('blue')
;outplot,times,flux_rg2_data,thick=charthick,color=cgColor('pink')

write_jpeg,image_sav+'fluximg'+rrr+time_st+'-'+time_ed+ttt+'.jpg',tvrd(true=3),true=3


xyouts,pos1[0]+0.01,pos1[3]-0.03,'Flux_total',/normal,align=0,charsize=charsize*1.3,color=255,charthick=charthick
xyouts,pos1[0]+0.01,pos1[3]-0.06,'Flux_left',/normal,align=0,charsize=charsize*1.3,color=cgColor('green'),charthick=charthick
xyouts,pos1[0]+0.01,pos1[3]-0.09,'Flux_right',/normal,align=0,charsize=charsize*1.3,color=cgColor('red'),charthick=charthick

;-------------------------plot-------------------------

print,'It is all ready..........................'
;DEVICE,/CLOSE
SET_PLOT,'X'
end