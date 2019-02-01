;.r residual_align.pro
;.compile package_muser.pro 
;.compile package_xy.pro
;
; Purpose :
; Alignment in residual images of 1.7G and other frequencies.
; Return the offset and data cube.
; 
; History :
; Xingyao Chen, 05 Dec 2017.
;-------------------------begins-------------------------
;WINDOW,0,xsize=700, ysize=700
;Device,RETAIN=2
;!P.FONT = 0

ENTRY_DEVICE = !D.NAME
SET_PLOT, 'PS'
!P.FONT = 0
xsz=9.0 & ysz=6.0;region 5

charsz=0.3*1.8 & charth=1.2*1.8 & th_ind =1

pos1 = [0.12,0.12,0.90,0.90]
posc1 = [pos1[2]+0.07,pos1[3]-0.2,pos1[2]+0.09,pos1[2]]

systim = SYSTIME(1)

sys=0 ;for windows system
sys=1 ;for linux   system

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
pixel=512
;data_namet=data_name[5]
;;nm_ind='05'

nch=16
for ooo=0,0 do begin
  
  filenamem17=file_search('./sav/1.6-2.0_Lnum2400/quiet_clean_img_fits_residual/MUSER_rsd_L_20141217_041055_396_1725.fits')
  ;filenamem17=file_search('./sav/1.6-2.0_Rnum2400/quiet_clean_img_fits_residual/MUSER_rsd_R_20141217_041055_400_1725.fits')
  datam17 = readfits( filenamem17, indexm17, /NOSCALE, /SILENT)

  data=make_array(pixel,pixel,nch+1)
  data[*,*,0]=datam17[*,*] 
  index=strarr((size(indexm17))[1],nch+1)
  index[*,0]=indexm17[*] 
  
  data_namet=data_name[ooo*2+4]
  
  for ch=0,nch-1 do begin
    nm_ind='sf'+strcompress(ch,/remove_all)
    
    ;;;============
    if sys eq 0 then begin
      filenamem = file_search('.\sav\'+data_namet+'\quiet_clean_img_fits_residual\'+'*.fits',count=count)
      data_sav='.\sav\'+data_namet+'\quiet_clean_img_fits_residual_align\'
    endif
    if sys eq 1 then begin
      filenamem = file_search('./sav/'+data_namet+'/quiet_clean_img_fits_residual/'+'*.fits',count=count)
      data_sav='./sav/'+data_namet+'/quiet_clean_img_fits_residual_align/'
    endif
    
    file_mkdir,data_sav
    freqch=(1200+400*ooo+ch*25+25)
    freqchstr=strcompress(freqch,/remove_all)
    
    
    tmp0=rstrpos(filenamem[ch],'MUSER_raw_')
    tmp1=rstrpos(filenamem[ch],'.fits')
    file=strmid(filenamem[ch],tmp0,tmp1-tmp0)
    print,file
    
    datam = readfits( filenamem[ch], indexm, /NOSCALE, /SILENT)
    
    data[*,*,ch+1]=datam[*,*]
    index[*,ch+1]=indexm[*]
 
  endfor 
  
  
  datatt=data
  disp = tr_get_disp(datatt ,/debug,/nowindow)
  datas = shift_img(datatt, disp)
  print,disp
  save,data,datas,disp,filename=data_sav+'residual_align.sav'
  
;  ;;fg_rigidalign, index, data, index_out, data_out, nt=1, dx=400, dy=400, x0=70, y0=70
;  
;  temp=data[*,*,4]
;  temp1=max(temp,loc1)
;  ind11=array_indices(temp, loc1)
;  print,ind11
;  
;  temp=data[*,*,5]
;  temp1=max(temp,loc1)
;  ind12=array_indices(temp, loc1)
;  print,ind12
;  
;  temp=datas[*,*,4]
;  temp1=max(temp,loc1)
;  ind1=array_indices(temp, loc1)
;  print,ind1
;  temp=datas[*,*,5]
;  temp1=max(temp,loc1)
;  ind1=array_indices(temp, loc1)
;  print,ind1
;  
endfor




;-------------------------plot-------------------------
device,/CLOSE
set_plot,'X'
print,'It is all ready..........................'
end