;.r img2vis.pro
;.compile package_muser.pro 
; 
; Name :
; img2vis
;
; Purpose :
; transform the image to visibility.
;
; History :
; Xingyao Chen, 16 Dec 2017.
;-------------------------begins-------------------------
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

data_namet=data_name[6]

;sys=0 ;for windows system
sys=1 ;for linux   system

nm_ind='17010203'

nm_vis='04'

;;;============;;for clean_img data
if sys eq 0 then begin
  filename=file_search('.\sav\'+data_namet+'\it_clean_img'+nm_ind+'\*.sav',count=count)
  data_sav='.\sav\'+data_namet+'\it_clean_img'+nm_ind+'_vis'+nm_vis+'\'
endif
if sys eq 1 then begin
  filename=file_search('./sav/'+data_namet+'/it_clean_img'+nm_ind+'/*.sav',count=count)
  data_sav='./sav/'+data_namet+'/it_clean_img'+nm_ind+'_vis'+nm_vis+'/'
endif
file_mkdir,data_sav

nfile=count
nfile=1
pixel=512
pix=pixel
for ff=0,0 do begin ;nfile-1
  tmp0=rstrpos(filename[ff],'clean_'+nm_ind+'_')
  file=strmid(filename[ff],tmp0,100)
  print,file
  restore,filename[ff],/ver
  
  img_soc=make_array(pixel,pixel)
  temp=where(img_cl gt 0.8*max(img_cl))
  help,temp
  img_soc[temp]=img_cl[temp]
  ;window,1,xs=512,ys=512
  ;plot_image,img_soc,position=[0.1,0.1,0.9,0.9],title='img_source 0.8*max'
 
  img=img_soc
  ;img=rotate(img,5)
  ;img=shift(img,pix/2,pix/2)
  uv_img=fft(img)
  
  day=gps_time_cl[0,2]
  bjt=gps_time_cl[0,3]+gps_time_cl[0,4]/60.+gps_time_cl[0,5]/3600.
  cal_uvw,day,bjt,u,v,w;print,bjt,h,ra,dec
  
  maxu=512.*180./!pi & maxv=512.*180./!pi 
  vis=complexarr(40,40)
  cor_data_cl=fltarr(40,40,2)
  for i=0,39 do begin
    for j=0,39 do begin
      vis[i,j]=uv_img[u[i,j]/maxu*pix+pix/2,v[i,j]/maxv*pix+pix/2]
      cor_data_cl[i,j,0]=real_part(vis[i,j])
      cor_data_cl[i,j,1]=imaginary(vis[i,j])
    endfor
  endfor
  
  save,cor_data_cl,gps_time_cl,filename=data_sav+'vis'+nm_vis+'_'+file
endfor
;
;cor_data1=cor_data_cl
;restore,'/Volumes/Seagate_cxy/muserdata_20141217/muser-e/sav/1.6-2.0_Lnum2400/keep_table17.sav',/ver
;flag_baseline,keep_table,cor_data1,vis,u,v   ;;;flag baseline
;FFT_IMAGE,U,V,VIS,BEAM,IMG

;-------------------------end-------------------------
print,'it is ok.............'
end