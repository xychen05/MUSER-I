;.r muser_image_align.pro
;.compile package_muser.pro 
;.compile package_xy.pro
;
; Purpose :
; Alignment between MUSER clean image and AIA image.
; Return the offset.
; 
; History :
; Xingyao Chen, 30 Nov 2017.
;-------------------------begins-------------------------
;WINDOW,0,xsize=700, ysize=700
;Device,RETAIN=2
;!P.FONT = 0

ENTRY_DEVICE = !D.NAME
SET_PLOT, 'PS'
!P.FONT = 0
xsz=9.0 & ysz=6.0;region 5

charsz=0.3*1.8 & charth=1.2*1.8 & th_ind =1

x1 = -1000. & x2 = 200 & y1 = -1200. & y2 = 0. ;region 2
x1 = -350. & x2 = 50 & y1 = -500. & y2 = -200. ;region 5
x1 = -900. & x2 = 300 & y1 = -650. & y2 = 150. ;region 2
x11 = x1 & x21 = x2 & y11 = y1 & y21 = y2 ;region 2

xrange = [x1,x2] & yrange = [y1,y2]

pos1 = [0.12,0.12,0.90,0.90]
posc1 = [pos1[2]+0.07,pos1[3]-0.2,pos1[2]+0.09,pos1[2]]

systim = SYSTIME(1)

sys=0 ;for windows system
sys=1 ;for linux   system

;-------------------------AIA image-------------------------
if sys eq 0 then filenamea='I:\muserdata_20141217\AIA\AIA20141217_041955_0304.fits';;;nn=725
if sys eq 1 then filenamea='/Volumes/Seagate_cxy/muserdata_20141217/AIA/AIA20141217_041955_0304.fits';;;nn=725

;filenamea='/Volumes/Seagate_cxy/muserdata_20141217/AIA/AIA20141217_042407_0304.fits';;;nn=255
;filenamea1 = file_search('/Volumes/TOSHIBA_cxy/data_20141217/AIA/prep/171/AIA20141217*.fits')
;filenamea2 = file_search('/Volumes/TOSHIBA_cxy/data_20141217/AIA/prep/131/AIA20141217*.fits')
;filenamea3 = file_search('/Volumes/TOSHIBA_cxy/data_20141217/AIA/prep/94/AIA20141217*.fits')

wavelength = ['94','171','193','211','304','335','131','1600']
mag_scal = [1.5,1.5,3,8.,3.,2.,1.5,5.]
rs_ratio = 1.0d
xllp=0.
nxp=4096.
yllp=0.
nyp=4096.

wl=4
wavelnth = wavelength[wl]

read_sdo, filenamea, indexa, dataa, /UNCOMP_DELETE,/noshell,/use_shared_lib,/HIDE
dataa=dataa/indexa.EXPTIME
imga = sdo_aia_scale_hdr(dataa/mag_scal[wl],indexa,xllp,yllp,nxp,nyp,wavelnth=wavelnth,rs_ratio=rs_ratio,/no_imgrscale)
index2map,indexa,imga,mapa

;;;;;============aia fitting
;ind1=950 & ind2=1450 & ind3=1500 & ind4=2000
;imga_rg=imga[ind1:ind2,ind3:ind4]
;indt=where(imga_rg lt 140)
;imga_rg[indt]=0
;plot_image,imga_rg
;imga_rg_fit=GAUSS2DFIT(imga_rg, coeffa, /TILT)
;thetaa=coeffa[6]
;cenax=4096/2.-(coeffa[4]+ind1) & cenay=4096/2.-(coeffa[5]+ind3)
;cenax=cenax*(indexa.CDELT1) & cenay=cenay*(indexa.CDELT2)
;print,'distance to the center[arcsec]:'
;print,cenax,cenay
;aligna=[thetaa,cenax,cenay,sqrt(cenax^2+cenay^2)]


;;-------------------------MUSER image-------------------------
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

data_namet=data_name[5]

;;nm_ind='05'
for ch=0,15 do begin
nm_ind='sf'+strcompress(ch,/remove_all)

;;;============
if sys eq 0 then begin
  filenamem = file_search('.\sav\'+data_namet+'\quiet_clean_img_fits\'+'*.fits',count=count)
  data_sav='.\sav\'+data_namet+'\quiet_clean_img_fits_align\'
endif
if sys eq 1 then begin
  filenamem = file_search('./sav/'+data_namet+'/quiet_clean_img_fits/'+'*.fits',count=count)
  data_sav='./sav/'+data_namet+'/quiet_clean_img_fits_align/'
endif

file_mkdir,data_sav

freqch=(1200+ch*25+25)
freqchstr=strcompress(freqch,/remove_all)
pixel=512

;-------------------------guassian fit-------------------------
nt_ind=1;5*12s
file_st=0;+600
nfile=1

loadct,39

;nn=ch

for nn=ch,ch do begin
  
  tmp0=rstrpos(filenamem[file_st+nt_ind*nn],'MUSER_raw_')
  tmp1=rstrpos(filenamem[file_st+nt_ind*nn],'.fits')
  file=strmid(filenamem[file_st+nt_ind*nn],tmp0,tmp1-tmp0)
  print,file
  
  datam = readfits( filenamem[file_st+nt_ind*nn], indexm, /NOSCALE, /SILENT)
  print,filenamem[file_st+nt_ind*nn]  
  
;  ;;;============guassian fit
;  datamt=mv_img(datam,-16,-52)
;  angle=19 ;1 min for rotating
;  ;;dmin=0 & dmax=4. ;;for ch=0
;  ;;datam=datamt
;  imgmorg=polax_rot(datamt, 9.3) ; rotate after moving the solar center to the image center
;  
  xx0=10 & yy0=10 & zz0=10
  alignm=make_array(xx0+1,yy0+1,zz0+1,4)
  for xx=0,xx0 do begin
    for yy=0,yy0 do begin
      for zz=0,zz0 do begin
        imgm1=mv_img(datam,-16+xx-xx0/2,-52+yy-yy0/2)
        imgm=polax_rot(imgm1,19+zz-zz0/2);;19+zz-5
        ;;============subtract the radio source in flare region
        index1=where(imgm ge 3.0)
        ;print,'pixels of radio source:',(size(index1))[1]
        locx=index1/pixel
        locy=index1 mod pixel
        ind_rg=[min(locx),max(locx),min(locy),max(locy)]
        ind1=min(ind_rg) & ind2=max(ind_rg);+10
        ;print,ind1,ind2
        imgmt=imgm
        imgmt[(ind1+5):ind2,(ind1-8):ind2]=0
        ;plot_image,imgmt
        ;;============guassion fit
        index1=where(imgm eq max(imgm))
        ;print,'pixels of radio source:',(size(index1))[1]
        locx=index1[0]/pixel
        locy=index1[0] mod pixel
        ind1=locx-65 & ind2=locx+15
        ind3=locy-40 & ind4=locy+40
        ;print,'X',ind1,ind2
        ;print,'Y',ind3,ind4  
        imgm_rg=imgmt[ind1:ind2,ind3:ind4]
        imgm_rg_fit=GAUSS2DFIT(imgm_rg, coeffm, /TILT)
        thetam=coeffm[6]
        cenmx=512/2.-(coeffm[4]+ind1) & cenmy=512/2.-(coeffm[5]+ind3)
        cenmx=cenmx*60.*60./512. & cenmy=cenmy*60.*60./512.
        ;print,'distance to the center[arcsec]:'
        ;print,cenmx,cenmy  
        alignm[xx,yy,zz,0]=thetam
        alignm[xx,yy,zz,1]=cenmx
        alignm[xx,yy,zz,2]=cenmy
        alignm[xx,yy,zz,3]=sqrt((cenmx-aligna[1])^2+(cenmy-aligna[2])^2)
      endfor
    endfor
  endfor
  save,aligna,alignm,filename=data_sav+'muser_image_align_guass_'+file+'.sav'
 
  ;;;============plotting
  filename_eps='/Users/xychen/Desktop/movie/image_align/'+'muser_image_align_guass_'+file+'.eps';;/movie/image
  device,file=filename_eps,xsize=xsz,ysize=ysz,BITS_PER_PIXEL=8,$
    SET_FONT='Helvetica',/ENCAPSULATED,scale_factor=2,/color
  
  restore,data_sav+'muser_image_align_guass_'+file+'.sav',/ver
  
;  diff0=abs(alignm[*,*,*,0]-aligna[0])
;  temp=min(diff0,loca0)
;  inf0=array_indices(diff0, loca0)
;  print,inf0
;  print,alignm[inf0[0],inf0[1],inf0[2],0],aligna[0];inf0[2]
;  
;  diff1=abs(alignm[*,*,*,1]-aligna[1])
;  temp=min(diff1,loca1)
;  inf1=array_indices(diff1, loca1)
;  print,inf1
;  print,alignm[inf1[0],inf1[1],inf1[2],1],aligna[1]
; 
;  diff2=abs(alignm[*,*,*,2]-aligna[2])
;  temp=min(diff2,loca2)
;  inf2=array_indices(diff2, loca2)
;  print,inf2
;  print,alignm[inf2[0],inf2[1],inf2[2],2],aligna[2]

  sz=size(alignm)
  diff3=alignm[*,*,*,3]
  temp=min(diff3,loca3)
  inf3=array_indices(diff3, loca3)
  print,inf3
  print,alignm[inf3[0],inf3[1],inf3[2],3]
  xx=inf3[0] & yy=inf3[1] & zz=inf3[2];inf3[2]
  xx0=sz[1] & yy0=sz[2] & zz0=sz[3]
  print,xx,yy,zz
  ;;xx=15 & yy=15 & zz=5
  imgmt=mv_img(datam,-16+xx-xx0/2,-52+yy-yy0/2)
  imgm=polax_rot(imgmt,19+zz-zz0/2);19+zz-5
  
  ;imgmt=mv_img(datam,-20,-35)
  ;imgm=polax_rot(imgmt,9.3);19+zz-5
  index2map,indexm,imgm,mapm;
  
  ;;============plot
  dmin=0 & dmax=4
  
  plot_map,mapm,position = pos1,xrange = xrange,yrange = yrange,$
    thick=charth*th_ind,charthick=charth,charsize=charsz,/limb,dmin =dmin;,dmin =5,dmax =65.;,/LOG_SCALE;
    
  plot_map,mapa,position = pos1,/notitle, $
    xrange = xrange,yrange = yrange,xstyle=1,ystyle=1, $ ;,/NOLABELS,/NOAXES,/NOXTICKS,/NOYTICKS,/NOTITLE
    thick=charth*th_ind,charthick=charth,charsize=charsz,LEVELS=[200],/cont,color=255,/over
  xyouts,pos1[0]+0.01,pos1[1]+0.02, 'MUSER-'+strcompress(indexm[13],/remove_all)+freqchstr+' MHz', $+textoidl('AIA-304-white contour') ,$
    /normal,color=255,align=0,charthick=charth,charsize=charsz
  cgColorbar,range=[dmin,dmax],yticks=3,/ver,color=0,charsize=charsz,position=posc1,charthick=charsz

endfor

endfor
;-------------------------plot-------------------------
device,/CLOSE
set_plot,'X'
print,'It is all ready..........................'
end