;.r muser_movie_rhessi_norh.pro
;.compile package_muser.pro 
; 
; Name :
; muser_movie_rhessi_norh
;
; Purpose :
; make movie of MUSER clean image.
; also the movie of rhessi, and norh
;
; History :
; Xingyao Chen, July 2017.
;-------------------------begins-------------------------

nm_ind='05'
filename = file_search('./sav/1.6-2.0_Lnum2400/it_clean_img'+nm_ind+'_fits_raw/*.fits',count=count)
;filename = file_search('./sav/1.6-2.0_Lnum2400/it_clean_img'+nm_ind+'_fits_raw/*.fits',count=count)
;filename = file_search('./sav/it_clean_img'+nm_ind+'_fits_raw_LR/*.fits',count=count)

filename17 = file_search('/Users/xychen/Desktop/mpro/norh/data/ipa/*')
filename34 = file_search('/Users/xychen/Desktop/mpro/norh/data/ipz/*')

filename_r=file_search('/Users/xychen/Desktop/mpro/rhessi/data/hsi_imagecube_60tx9e_20141217_044130.fits',count=nfile_total)

filenamea='/Users/xychen/Desktop/mpro/hmi/AIA20141217_045207_0304.fits'
wavelength = ['94','171','193','211','304','335','131','1600']
mag_scal = [1.5,1.5,3,8.,3.,2.,1.5,5.]
rs_ratio = 1.0d
xllp=0.
nxp=4096.
yllp=0.
nyp=4096.
wl=4
wavelnth = wavelength[wl]

outfile = '/Users/xychen/Desktop/muser_rhessi_norh_L+R_01.mp4'

;-------------------------basic parameters-------------------------
x1 = -1200. & x2 = 0 & y1 = -1200. & y2 = 0. ;region 2

xrange = [x1,x2] & yrange = [y1,y2]

pos1 = [0.12,0.12,0.90,0.90]
posc1 = [pos1[2]+0.07,pos1[3]-0.2,pos1[2]+0.09,pos1[2]]

systim = SYSTIME(1)

;-------------------------171-------------------------
video_file = outfile
video = idlffvideowrite(video_file)
framerate = 12

framedims   = [1000,1000]

stream = video.addvideostream(framedims[0], framedims[1], framerate,BIT_RATE=5.0e6)
ENTRY_DEVICE = !D.NAME
set_plot,'z',/copy
device, set_resolution=framedims, set_pixel_depth=24, decomposed=0

;-------------------------aia data-------------------------
aia_lct, rr, gg, bb, wavelnth=wavelnth, /load
read_sdo, filenamea, indexa, dataa, /UNCOMP_DELETE,/noshell,/use_shared_lib,/HIDE
dataa=dataa/indexa.EXPTIME
imga = sdo_aia_scale_hdr(dataa/mag_scal[wl],indexa,xllp,yllp,nxp,nyp,wavelnth=wavelnth,rs_ratio=rs_ratio,/no_imgrscale)
index2map,indexa,imga,mapa

;-------------------------plot-------------------------
nt_ind=1;5*12s
file_st=0;+600

;-------------------------plot-------------------------
for nn=0,count-1 do begin;nfile
  ;;-------------------------image of MUSER-------------------------
  read_sdo, filename[file_st+nt_ind*nn], indexm, datam, /UNCOMP_DELETE,/noshell,/use_shared_lib,/HIDE
  
  ;  index_temp=where(datam gt 0.5)
  ;  img_temp=make_array((size(datam))[1],(size(datam))[2])
  ;  img_temp[index_temp]=datam[index_temp]
  ;  imgm=GAUSS2DFIT(img_temp, coeff, /TILT)
  
  angle=9.35+(8.89-9.35)/24.*4+(8.89-9.35)/24./60.*25+(8.89-9.35)/24./60./60*(nn) ;1 min for rotating
  imgm=polax_rot(datam, angle) ; rotate after moving the solar center to the image center
  
  index2map,indexm,bytscl((imgm),min=2.),map
  
  loadct,39
  plot_map,map,position = pos1,xrange = xrange,yrange = yrange,$
    thick=3.,charthick=2.7,charsize=1.4,/limb,title=' ',/noerase
  xyouts,1.1*pos1[0],1.2*pos1[1], textoidl('MUSER-1.725GHz -L ')+(indexm.T_OBS)+' AIA/304' ,$
    /normal,color=255,align=0,charsize=1.5,charthick=2.
  ;cgColorbar,range=[2.0,max(imgm)],yticks=3,/ver,color=255,charsize=1.4,position=posc1,charthick=2.4,/noerase

  ;;-------------------------image of AIA-------------------------
  plot_map,mapa,position = pos1,/notitle, $
    xrange = xrange,yrange = yrange,xstyle=1,ystyle=1, $ ;,/NOLABELS,/NOAXES,/NOXTICKS,/NOYTICKS,/NOTITLE
    thick=1.4,charthick=1.2,charsize=0.3,LEVELS=[200],/cont,color=255,/over

  ;;-------------------------image of NORH-------------------------
  
  nt_ind_n=1;nt_ind*10s
  file_st_n=130
  ii=round(nn/10.)
  ;loadct,3
  read_sdo, filename17[file_st_n+nt_ind_n*ii], index17, data17, /UNCOMP_DELETE,/noshell,/use_shared_lib,/HIDE
;  ;aldata17=make_array(1024,1024)
  ;index_temp=where(data17 lt 0.5*max(data17))
  ;data17[index_temp]=0
  index2map,index17,data17,map17
  plot_map,map17,/cont,/noerase,xstyle=1+4,ystyle=1+4, $
    xrange = xrange,yrange = yrange,title=' ', $;,/NOAXES,/NOLABELS,/NOXTICKS,/NOYTICKS
    position = pos1,LEVELS = [0.5,0.7,0.95]*max(data17),c_color =cgColor('red'),c_thick=3.2,thick=2.
  xyouts,1.1*pos1[0],1.4*pos1[1], 'NORH- 17GHz- R+L 2014-12-17T'+strmid(filename17[file_st_n+nt_ind_n*ii],51,2) $
    +':'+strmid(filename17[file_st_n+nt_ind_n*ii],53,2)+':'+strmid(filename17[file_st_n+nt_ind_n*ii],55,2), $
    /normal,color=cgColor('red'),align=0,charsize=1.5,charthick=2.0
  
  ;loadct,8
  read_sdo, filename34[file_st_n+nt_ind_n*ii], index34, data34, /UNCOMP_DELETE,/noshell,/use_shared_lib,/HIDE  
;  data34=congrid(data34,1024,1024)
  ;index_temp=where(data34 lt 0.5*max(data34))
  ;data34[index_temp]=0
  index2map,index34,data34,map34
  plot_map,map34,/cont,/noerase,xstyle=1+4,ystyle=1+4, $
    xrange = xrange,yrange = yrange,title=' ', $;,/NOAXES,/NOLABELS,/NOXTICKS,/NOYTICKS
    position = pos1,LEVELS = [0.5,0.7,0.95]*max(data34),c_color =cgColor('green'),c_thick=3.2,thick=2.
  xyouts,1.1*pos1[0],1.6*pos1[1], 'NORH- 34GHz- R+L 2014-12-17T'+strmid(filename34[file_st_n+nt_ind_n*ii],51,2) $
    +':'+strmid(filename34[file_st_n+nt_ind_n*ii],53,2)+':'+strmid(filename34[file_st_n+nt_ind_n*ii],55,2), $
    /normal,color=cgColor('green'),align=0,charsize=1.5,charthick=2.0
    
  print,strmid(filename17[file_st_n+nt_ind_n*ii],41,20)
  print,strmid(filename34[file_st_n+nt_ind_n*ii],41,20)
  
;  ;;-------------------------image of RHESSI-------------------------
;  ;filename_r=file_search('/Users/xychen/Downloads/hsi_image_20141217_044258.fits',count=nfile_total)
;  ;04:41:30-04:46:30 ;42,43,44,45,46min
;  hsi_fits2map, filename_r, hsi_map,/sep
;  nt_ind_r=1;nt_ind*10s
;  file_st_r=6
;  loadct,39
;  if jj ge 16*60 then begin
;    for rr=2,2 do begin;for energy band
;      hh=(jj-16*60)/5
;      plot_map,hsi_map[rr,file_st_r+hh*nt_ind_r],/cont,/noerase,xstyle=1+4,ystyle=1+4, $
;        xrange = xrange,yrange = yrange,title=' ', $;,/NOAXES,/NOLABELS,/NOXTICKS,/NOYTICKS
;        LEVELS = [0.5,0.7,0.95]*max((hsi_map[rr,file_st_r+hh*nt_ind_r]).DATA), $
;        position = position,c_color =cgColor('black'),c_thick=2.5
;      xyouts,1.1*position[0],1.70*position[1], string((hsi_map[rr,file_st_r+hh*nt_ind_r]).ID)+' 2014-12-17T' $
;        +strmid((hsi_map[rr,file_st_r+hh*nt_ind_r]).TIME,12,8)  $
;        ,/normal,color=cgColor('black'),align=0,charsize=1.2,charthick=2.0
;        
;      print,(hsi_map[rr,file_st_r+hh*nt_ind_r]).ID,(hsi_map[rr,file_st_r+hh*nt_ind_r]).TIME
;    endfor
;  endif
   
   
   
  ;print,filename[nn]
  print,minmax(datam)
  timestamp = video.put(stream, tvrd(true=1))

endfor

;-------------------------ending-------------------------
close,/all
device, /close
video.cleanup
set_plot,'X'

print,'It is all ready..........................'

end