;.r norh_flux_image_dataprep.pro
;.compile package_muser.pro 
; Name :
; norh_flux_image_dataprep
;
; Purpose :
; make image of NORH 17, 34 GHz image.
; preparation for the flux of selected region.
; 
; History :
; Xingyao Chen, 30 July 2017.
;-------------------------begins-------------------------
WINDOW,0,xsize=700, ysize=700
WINDOW,2,xsize=1100, ysize=700

Device,RETAIN=2

x1 = -300. & x2 = 0 & y1 = -490. & y2 = -190. ;region 2

xrange = [x1,x2] & yrange = [y1,y2]

pos1 = [0.15,0.15,0.90,0.90]
posc1 = [pos1[2]+0.07,pos1[3]-0.2,pos1[2]+0.09,pos1[2]]

charsize=1.8
charthick=1.8

;;;;-------------------------for 17 GHz-------------------------
;filename17=file_search('/Users/xychen/Desktop/mpro/norh/data/ipa/*',count=count)
;
;data_sav='./sav/flux_image/'
;file_mkdir,data_sav
;
;rrr='01_17G_';;'_'
;pixel=64

;;;-------------------------for 34 GHz-------------------------
filename17=file_search('/Users/xychen/Desktop/mpro/norh/data/ipz/*',count=count)
data_sav='./sav/flux_image/'
file_mkdir,data_sav

rrr='01_34G_';;'_'
pixel=64*2

;-------------------------start-------------------------
nt_ind=1;5*12s
file_st=101;+600
nfile=count-file_st
nfile=150

flux_rg1=make_array(nfile)
flux_rg2=make_array(nfile)
flux_total=make_array(nfile)

flux_rg1_data=make_array(nfile)
flux_rg2_data=make_array(nfile)
flux_total_data=make_array(nfile)

wset,0
wshow,0
loadct,39
;-------------------------plot-------------------------
for nn=75,75 do begin
  
  read_sdo, filename17[file_st+nt_ind*nn], index17, data17, /UNCOMP_DELETE,/noshell,/use_shared_lib,/HIDE
  img17=data17
  
  thred=0.3
  index_temp=where(img17 lt thred*max(img17))
  img17[index_temp]=thred*max(img17)
  
  index2map,index17,img17,map17
  
  plot_map,map17,position = pos1,xrange = xrange,yrange = yrange,$
    thick=charthick,charthick=charthick,charsize=charsize,/limb,title='NORH-'+string(strmid(filename17[file_st+nt_ind*nn],41,20))
  cgColorbar,range=[thred*max(img17),max(img17)],yticks=3,/ver,color=255,charsize=charsize-0.4,position=posc1,charthick=charthick
  
  ;-------------------------gaussian fit of selected region-------------------------
  img17_fit=GAUSS2DFIT(img17, coeff, /TILT)
  width=3.0
  gauss_points, coeff,width,/ydir,lx1,lx2,ly1,ly2,lx3,lx4,ly3,ly4,lx5,lx6,ly5,ly6
  
  ;-------------------------definition of 6 points-------------------------
  plots,[lx1,lx2]*(map17.dx)+(index17.CRVAL1)-pixel*(map17.dx)/2.,[ly1,ly2]*(map17.dy)+(index17.CRVAL2)-pixel*(map17.dy)/2.,symsize=3,color=255
  plots,[lx3,lx4]*(map17.dx)+(index17.CRVAL1)-pixel*(map17.dx)/2.,[ly3,ly4]*(map17.dy)+(index17.CRVAL2)-pixel*(map17.dy)/2.,symsize=3,color=255
  plots,[lx5,lx6]*(map17.dx)+(index17.CRVAL1)-pixel*(map17.dx)/2.,[ly5,ly6]*(map17.dy)+(index17.CRVAL2)-pixel*(map17.dy)/2.,symsize=3,color=255
  
  plots,[lx3,lx5]*(map17.dx)+(index17.CRVAL1)-pixel*(map17.dx)/2.,[ly3,ly5]*(map17.dy)+(index17.CRVAL2)-pixel*(map17.dy)/2.,symsize=3,color=255
  plots,[lx4,lx6]*(map17.dx)+(index17.CRVAL1)-pixel*(map17.dx)/2.,[ly4,ly6]*(map17.dy)+(index17.CRVAL2)-pixel*(map17.dy)/2.,symsize=3,color=255

  ;-------------------------definition of 6 points-------------------------
  region1=sl_region_index([lx1,ly1],[lx2,ly2],[lx3,ly3],[lx5,ly5])
  region2=sl_region_index([lx1,ly1],[lx2,ly2],[lx4,ly4],[lx6,ly6])
  regiont=sl_region_index([lx3,ly3],[lx4,ly4],[lx5,ly5],[lx6,ly6])
  
  time_string=strmid(filename17[file_st+nt_ind*nn],41,20)
  
endfor
save,lx1,lx2,ly1,ly2,lx3,lx4,ly3,ly4,lx5,lx6,ly5,ly6,filename=data_sav+'norh_points'+rrr+time_string+'.sav'
save,region1,region2,regiont,filename=data_sav+'norh_region'+rrr+time_string+'.sav'

;-------------------------for 17 GHz-------------------------
times=dblarr(nfile)
print,nfile
loadct,39
for nn=0,nfile-1 do begin;nfile;count-1
  
  read_sdo, filename17[file_st+nt_ind*nn], index17, data17, /UNCOMP_DELETE,/noshell,/use_shared_lib,/HIDE
  
  img17=data17
  
  thred=0.3
  index_temp=where(img17 lt thred*max(img17))
  img17[index_temp]=thred*max(img17)
  
  ;;-------------------------flux data from the image-------------------------
  
  flux_rg1_total=make_array((size(region1))[2])
  flux_rg1_total_data=make_array((size(region1))[2])
  for ii=0,(size(region1))[2]-1 do begin
    flux_rg1_total[ii]=img17[region1[0,ii],region1[1,ii]]
    flux_rg1_total_data[ii]=data17[region1[0,ii],region1[1,ii]]
    ;plots,[region1[0,ii]],[region1[1,ii]],psym=1,symsize=2,color=255
  endfor
  flux_rg2_total=make_array((size(region2))[2])
  flux_rg2_total_data=make_array((size(region2))[2])
  for ii=0,(size(region2))[2]-1 do begin
    flux_rg2_total[ii]=img17[region2[0,ii],region2[1,ii]]
    flux_rg2_total_data[ii]=data17[region2[0,ii],region2[1,ii]]
    ;plots,[region2[0,ii]],[region2[1,ii]],psym=1,symsize=2,color=255
  endfor
  flux_tl=make_array((size(regiont))[2])
  flux_tl_data=make_array((size(regiont))[2])
  for ii=0,(size(regiont))[2]-1 do begin
    flux_tl[ii]=img17[regiont[0,ii],regiont[1,ii]]
    flux_tl_data[ii]=data17[regiont[0,ii],regiont[1,ii]]
    ;plots,[regiont[0,ii]],[regiont[1,ii]],psym=1,symsize=2,color=255
  endfor
  
  flux_rg1[nn]=mean(flux_rg1_total)
  flux_rg2[nn]=mean(flux_rg2_total)
  flux_total[nn]=mean(flux_tl)
  
  flux_rg1_data[nn]=mean(flux_rg1_total_data)
  flux_rg2_data[nn]=mean(flux_rg2_total_data)
  flux_total_data[nn]=mean(flux_tl_data)
  
  times[nn]=anytim(index17.DATE_OBS)
  if nn eq 0 then time_st=time_obs_to_string(index17.DATE_OBS)
  if nn eq nfile-1 then time_ed=time_obs_to_string(index17.DATE_OBS)
  
endfor

save,flux_rg1,flux_rg2,flux_total,flux_rg1_data,flux_rg2_data,flux_total_data,times $
  ,filename=data_sav+'norh_flux'+rrr+time_string+'_'+time_st+'-'+time_ed+'.sav'
  
timerange=[times[0],times[nfile-1]]
;-------------------------plot-------------------------
wset,2
wshow,2
rgmax=max([max(flux_rg1_data),max(flux_rg2_data),max(flux_total_data)])
rgmin=min([min(flux_rg1_data),min(flux_rg2_data),min(flux_total_data)])

utplot,times,flux_total_data,timerange=timerange,xstyle=1,ystyle=1,thick=chthick,charthick=charthick,charsize=charsize $
  ,position=pos1,ytitle=textoidl('Average Flux [Arbitrary]'),title='NORH'+rrr+'GHz',yrange=[rgmin,rgmax]
  ;,yrange=[min(flux_rg1)-5,max(flux_rg1)+10],title=string(pixel)+'*'+string(pixel)
outplot,times,flux_rg1_data,thick=charthick,color=cgColor('green')
outplot,times,flux_rg2_data,thick=charthick,color=cgColor('red')
;
;outplot,times,flux_total_data,thick=charthick,color=cgColor('yellow')
;outplot,times,flux_rg1_data,thick=charthick,color=cgColor('blue')
;outplot,times,flux_rg2_data,thick=charthick,color=cgColor('pink')
;

;xyouts,pos1[0]+0.01,pos1[3]-0.03,'flux_total',/normal,align=0,charsize=charsize*1.3,color=255,charthick=charthick
;xyouts,pos1[0]+0.01,pos1[3]-0.06,'flux_left',/normal,align=0,charsize=charsize*1.3,color=cgColor('green'),charthick=charthick
;xyouts,pos1[0]+0.01,pos1[3]-0.09,'flux_right',/normal,align=0,charsize=charsize*1.3,color=cgColor('red'),charthick=charthick

;-------------------------plot-------------------------
print,'It is all ready..........................'
end