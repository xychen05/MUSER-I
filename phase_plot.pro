;.r phase_plot.pro
;.compile package_muser.pro
; 
; Name :
; phase_plot
; 
; Purpose :
; Display the amplitude and phase plot of cross-correlation between two antennas.
; 
; Inputs :
; A file with cor_data and gps_time.
;
; Outputs :
; Amplitude and phase plot of one antenna cross-correlating with the other antennas.
; 
; Notes :
; You can choose plots for original data, fs data, calibration data and prep data separately.
; 
; History :
; Xingyao Chen, June 2017
;-------------------------begins-------------------------
DEVICE, DECOMPOSED = 0
Device, RETAIN=2
!P.FONT = 0
!P.Color=0
!P.Background=255

data_name=['0.4-0.8_L',$
  '0.4-0.8_R',$
  '0.8-1.2_L',$
  '0.8-1.2_R',$
  '1.2-1.6_L',$
  '1.2-1.6_R',$
  '1.6-2.0_L',$
  '1.6-2.0_R']
temp=replicate('num0060',8)
data_name=data_name+temp

num=60
data_namet=data_name[6]

;sys=0 ;for windows system
sys=1 ;for linux   system

nm_ind='05'

;;-------------------------file name for windows-------------------------
;;;for original data
;if sys eq 0 then filename=file_search('.\'+data_namet+'\*.sav',count=count)
;if sys eq 1 then filename=file_search('./'+data_namet+'/*.sav',count=count)

;;;for fs data
;if sys eq 0 then filename=file_search('.\'+data_namet+'\fs\*.sav',count=count)
;if sys eq 1 then filename=file_search('./'+data_namet+'/fs/*.sav',count=count)
;
;;;for calibration data
;if sys eq 0 then filename=file_search('.\cal\fs\*.sav',count=count)
;if sys eq 1 then filename=file_search('./cal/fs/*.sav',count=count)
;
;;;for prep data
if sys eq 0 then filename=file_search('.\'+data_namet+'\prep\*.sav',count=count)
if sys eq 1 then filename=file_search('./'+data_namet+'/prep/*.sav',count=count)

;;;for rm data
;if sys eq 0 then filename=file_search('.\sav\'+data_namet+'\rm_cor_data'+nm_ind+'\*.sav',count=count)
;if sys eq 1 then filename=file_search('./sav/'+data_namet+'/rm_cor_data'+nm_ind+'/*.sav',count=count)
;
;;;for pm data
;if sys eq 0 then filename=file_search('.\sav\'+data_namet+'\pm_cor_data'+nm_ind+'\*.sav',count=count)
;if sys eq 1 then filename=file_search('./sav/'+data_namet+'/pm_cor_data'+nm_ind+'/*.sav',count=count)

;;;for it data
;if sys eq 0 then filename=file_search('.\sav\'+data_namet+'\it_cor_data'+nm_ind+'\' + $
;  ;'cor_data_it_L20141217_020155_363-20141217_040955_371.sav',count=count)
;  ;'cor_data_it_L20141217_030155_365-20141217_040955_371.sav',count=count)
;  ;'cor_data_it_L20141217_013655_365-20141217_040955_371.sav',count=count)
;  ;'cor_data_it_L20141217_013755_364-20141217_040955_371.sav',count=count)
;  'cor_data_it_R20141217_013755_368-20141217_040955_374.sav',count=count)  
;if sys eq 1 then filename=file_search('./sav/'+data_namet+'/it_cor_data'+nm_ind+'/' + $
;  ;'cor_data_it_L20141217_020155_363-20141217_040955_371.sav',count=count)
;  ;'cor_data_it_L20141217_030155_365-20141217_040955_371.sav',count=count)
;  'cor_data_it_L20141217_013655_365-20141217_040955_371.sav',count=count)
;


if sys eq 0 then image_sav='.\imagetest\'+data_namet+'\phase_prep\'
if sys eq 1 then image_sav='./imagetest/'+data_namet+'/phase_prep/'
file_mkdir,image_sav

;;-------------------------begins-------------------------

kk_int=count ;for 'kk_int' min per phase image in total 
;num=fix(num/kk_int/1.)
num=long(num)

for kk_st=0,count/kk_int-1 do begin ;less than count/kk_int-1 ;for file
  
  ;cor_data1=make_array(num*kk_int,40,40,2) ;just for original data
  
  cor_data=make_array(num*kk_int,40,40,2) ;just for rm mv pm data
  gps_time=intarr(num*kk_int,9)
  
  for kk=0,kk_int-1 do begin ;for max(kk) duration
    
;   ;;for original data
;    ff=kk_st*kk_int+kk ;kk_st*kk_int,(kk_st+1)*kk_int-1
;    restore,filename[ff] ;,/ver
;    tmp0=rstrpos(filename[ff],'COR_data_')
;    file=strmid(filename[ff],tmp0,50)
;    print,file
;    sz=(size(cor_data))[1]
;    cor_ind=findgen(num)*kk_int
;    for z=0,39 do begin
;      for zz=0,39 do begin
;        for zzz=0,1 do begin
;          cor_data1[kk*num:((kk+1)*num-1),z,zz,zzz]=cor_data[cor_ind,4,z,zz,zzz]
;        endfor
;      endfor
;    endfor  
   
;    ;;for fs data
;    ff=kk_st*kk_int+kk ;kk_st*kk_int,(kk_st+1)*kk_int-1
;    restore,filename[ff] ;,/ver
;    tmp0=rstrpos(filename[ff],'cor_data_')
;    file=strmid(filename[ff],tmp0,50)
;    print,file
;    sz=(size(cor_data_fs))[1]
;    cor_ind=findgen(num);*kk_int   
;    cor_data[kk*num:((kk+1)*num-1),*,*,*]=cor_data_fs[cor_ind,*,*,*]
;    gps_time[kk*num:((kk+1)*num-1),*] = gps_time_fs[cor_ind,*]
    
;    ;;for calibration data
;    ff=kk_st*kk_int+kk ;kk_st*kk_int,(kk_st+1)*kk_int-1
;    restore,filename[ff] ;,/ver
;    tmp0=rstrpos(filename[ff],'cor_data_')
;    file=strmid(filename[ff],tmp0,50)
;    print,file
;    sz=(size(cor_data_cal))[1]
;    cor_ind=findgen(num)*kk_int
;    cor_data[kk*num:((kk+1)*num-1),*,*,*]=cor_data_cal[cor_ind,*,*,*]
;    gps_time[kk*num:((kk+1)*num-1),*] = gps_time_cal[cor_ind,*]

    ;;for prep data
    ff=kk_st*kk_int+kk ;kk_st*kk_int,(kk_st+1)*kk_int-1
    restore,filename[ff] ;,/ver
    tmp0=rstrpos(filename[ff],'cor_data_')
    file=strmid(filename[ff],tmp0,50)
    print,file
    sz=(size(cor_data_prep))[1]
    cor_ind=findgen(num);*kk_int
    cor_data[kk*num:((kk+1)*num-1),*,*,*]=cor_data_prep[cor_ind,*,*,*]
    gps_time[kk*num:((kk+1)*num-1),*] = gps_time_prep[cor_ind,*]
    
;    ;;for rm data
;    ff=kk_st*kk_int+kk ;kk_st*kk_int,(kk_st+1)*kk_int-1
;    restore,filename[ff] ;,/ver
;    tmp0=rstrpos(filename[ff],'cor_data_')
;    file=strmid(filename[ff],tmp0,50)
;    print,file
;    sz=(size(cor_data_rm))[1]
;    cor_ind=round(findgen(num)*sz/1./num)
;    print,sz
;    cor_data[kk*num:((kk+1)*num-1),*,*,*]=cor_data_rm[cor_ind,*,*,*]
;    gps_time[kk*num:((kk+1)*num-1),*] = gps_time_rm[cor_ind,*]
      
;    ;;for pm data
;    ff=kk_st*kk_int+kk ;kk_st*kk_int,(kk_st+1)*kk_int-1
;    restore,filename[ff] ;,/ver
;    tmp0=rstrpos(filename[ff],'cor_data_')
;    file=strmid(filename[ff],tmp0,50)
;    print,file
;    sz=(size(cor_data_pm))[1]
;    print,sz
;    cor_data= cor_data_pm
;    gps_time = gps_time_pm
    
;    ;;for it data
;    ff=kk_st*kk_int+kk ;kk_st*kk_int,(kk_st+1)*kk_int-1
;    restore,filename[ff] ;,/ver
;    tmp0=rstrpos(filename[ff],'cor_data_')
;    file=strmid(filename[ff],tmp0,80)
;    print,file
;    sz=(size(cor_data_it))[1]
;    print,sz
;    cor_data= cor_data_it
;    gps_time = gps_time_it
    
    ;;;============
    if kk eq 0 then start_time=gps_time[0,*] ;time_gps_to_spec(gps_time[0,*])
    if kk eq kk_int-1 then end_time=gps_time[sz-1,*]  ;time_gps_to_spec(gps_time[num*kk_int-1,*])

    if kk eq 0 then file0=file
    
  endfor
  ;cor_data=cor_data1 ;just for original data
  amp=sqrt((cor_data[*, *, *, 1])^2+(cor_data[*, *, *, 0])^2)
  phase=atan(cor_data[*, *, *, 1],cor_data[*, *, *, 0])*180./!pi

  timerg=strcompress(start_time[3],/remove_all)+strcompress(start_time[4],/remove_all) $
    +'_to_'+strcompress(end_time[3],/remove_all)+strcompress(end_time[4],/remove_all)
  
;;;;-------------------------remove bad points again-------------------------   
;  rm_time_ag=rm_time_index(phase,1,300,50,3,index=index)
;  help,rm_time_ag
;  
;  sz=(size(phase))[1]-(size(rm_time_ag))[1]
;  amp_ag=dblarr(sz,40,40)
;  phase_ag=dblarr(sz,40,40)
;  ind_temp=indgen((size(phase))[1])
;  remove,rm_time_ag,ind_temp
;  for rr=0,sz-1 do begin
;    amp_ag[rr,*,*]=amp[ind_temp[rr],*,*]
;    phase_ag[rr,*,*]=phase[ind_temp[rr],*,*]
;  endfor
;
;  amp=amp_ag
;  phase=phase_ag
;  
;;;;-------------------------phase plot after flag antennas-------------------------  
;  flag_table=[10,11,12,13,14,16,17,18,19,24,25,26,27,38,39]   ;flag_table 5
;
;  amp_flag=amp
;  phase_flag=phase
;  
;  for ff=0,(size(flag_table))[1]-1 do begin
;    for xx=0,39 do begin
;      if xx eq  flag_table[ff] then begin
;        amp_flag[*,xx,*]=0
;        phase_flag[*,xx,*]=0
;      endif
;      for yy=0,39 do begin
;        if yy eq  flag_table[ff] then begin
;          amp_flag[*,*,yy]=0
;          phase_flag[*,*,yy]=0
;        endif
;      endfor
;    endfor
;  endfor
;  
;  amp=amp_flag
;  phase=phase_flag

;;;;-------------------------phase plot after flag baselines-------------------------
;  keep_table=keep_baseline(phase, 20, phase_rms=phase_rms)  ;;phase_rms [40,40]
;  ;;save,keep_table,filename='.\sav\keep_table.sav'
;  save,keep_table,filename='./sav/keep_table.sav'
;
;  amp_flag=make_array((size(amp))[1],(size(amp))[2],(size(amp))[3])
;  phase_flag=make_array((size(phase))[1],(size(phase))[2],(size(phase))[3])
;  
;  for ff=0,(size(keep_table))[2]-1 do begin
;    for xx=0,39 do begin
;      for yy=xx+1,39 do begin
;        if (xx eq keep_table[0,ff]) and (yy eq keep_table[1,ff]) then begin 
;          amp_flag[*,xx,yy]=amp[*,xx,yy]
;          phase_flag[*,xx,yy]=phase[*,xx,yy]
;        endif      
;      endfor
;    endfor
;  endfor
;
;  amp=amp_flag
;  phase=phase_flag
;  
;;;-------------------------plot begin-------------------------
  loadct,0
  WINDOW,1,XS=512*3.2/1.5,YS=512*1.8/1.5
  ;!P.multi=[0,1,1];[0,8,5]
  !P.multi=[0,8,5];[0,8,5]
  
  ;mode=0 ;;plot amptitute and phase
  ;mode=1 ;;plot amptitute
  ;mode=2 ;;plot phase
  ;mode=3 ;;plot Real Part of COR_DATA
  ;mode=4 ;;plot Imaginary Part of COR_DATA
  mode=2
  
  pp=['amp_phase','amp','phase','cor_Re','cor_Im']
  pp=pp[mode]
  
  sz=(size(phase))[1]
  
  for jj=0,0 do begin
    for ii=0,39 do begin
      
      phase_plot=phase[*,jj,ii]
      amp_plot=amp[*,jj,ii]     
      
      amp_plot=alog(amp_plot)
      
      if mode eq 0 then begin       
        ;;;for both amptitute and phase plot
        plot,amp_plot $
          ,thick=1.5,ystyle=1+8,xstyle=1,ytitle='  ',yrange=[5,13],/nodata,color=cgcolor('red')
        plot,amp_plot $
          ,thick=2,ystyle=4,xstyle=4+1,ytitle='  ',yrange=[5,13],/noerase,color=cgcolor('red');,yrange=[0,max(amp[*,0,ii])];
        plot,phase_plot,title='A'+strcompress(jj,/remove_all)+' and A'+strcompress(ii,/remove_all)+' / '+timerg $
          ,thick=1.5,/noerase,ystyle=4,xstyle=1+4,color=cgcolor('black'),yr=[-200,200] ;;[-4,4]
        axis,yaxis=1,ytitle='A '+strcompress(ii,/remove_all),color=cgcolor('black');,charthick=1.5,charsize=0.35  
      endif
      
      if mode eq 1 then begin
        ;;just for amptitute plot
        plot,amp_plot,thick=1.5,ystyle=1,xstyle=1,ytitle='Amp',yrange=[5,13] $;,yrange=[1e3,3e4],/ylog
          ,title='A'+strcompress(jj,/remove_all)+' and A'+strcompress(ii,/remove_all)+' / '+timerg
      endif
      
      if mode eq 2 then begin
        ;;just for phase plot
        plot,phase_plot,thick=1.5,ystyle=1,xstyle=1,ytitle='Phase',yrange=[-200,200],psym=6 $;;[-4,4]
          ,title='A'+strcompress(jj,/remove_all)+' and A'+strcompress(ii,/remove_all)+' / '+timerg;+' / rms'+string(phase_rms[jj,ii])
        ;oplot,phase_fit,color=cgcolor('red')
        oplot,[0,sz-1],[-180,-180],color=cgcolor('red'),linestyle=1
        oplot,[0,sz-1],[180,180],color=cgcolor('red'),linestyle=1
      endif
      
      if mode eq 3 then begin
        ;;just for Real Part of COR_DATA
        plot,cor_data[*,jj,ii,0],thick=1.5,ystyle=1,xstyle=1,ytitle='cor_Re' $;,yrange=[1e3,3e4],/ylog
          ,title='A'+strcompress(jj,/remove_all)+' and A'+strcompress(ii,/remove_all)+' / '+timerg
      endif
      
      if mode eq 4 then begin
        ;;just for Imaginary Part of COR_DATA
        plot,cor_data[*,jj,ii,1],thick=1.5,ystyle=1,xstyle=1,ytitle='cor_Im' $;,yrange=[1e3,3e4],/ylog
          ,title='A'+strcompress(jj,/remove_all)+' and A'+strcompress(ii,/remove_all)+' / '+timerg
      endif
      
;      write_jpeg,image_sav+pp+'_A'+strcompress(jj,/remove_all)+'_A' $
;        +strcompress(ii,/remove_all)+'_'+strmid(file,0,46)+'.jpg',tvrd(true=3),true=3
  
    endfor
    
    write_jpeg,image_sav+pp+'_A'+strcompress(jj,/remove_all)+'_A' $
      +strcompress(ii,/remove_all)+'_'+strmid(file,0,46)+'.jpg',tvrd(true=3),true=3
    
  endfor

endfor

;-------------------------ending-------------------------

print,'OK'
end