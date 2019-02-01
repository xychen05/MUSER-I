;.r phase_rm_time_old.pro
;.compile package_muser.pro
;
; Name :
; phase_rm_time
;
; Purpose :
; select bad points according to the phase plot of cross-correlation between two antennas.
; save the bad points to rm_time to file rm_time_.sav.
; 
; Inputs :
; A file with cor_data_prep and gps_time_prep.
;
; Outputs :
; Amplitude and phase plot of one antenna cross-correlating with the other antennas.
;
; History :
; Xingyao Chen, June 2017
; 
;-------------------------begins-------------------------
DEVICE, DECOMPOSED = 0
Device, RETAIN=2
!P.FONT = 0
!P.Color=0
!P.Background=255
WINDOW,1,XS=512*3.2/1.5,YS=512*1.8/1.5
;WINDOW,1,XS=512*2.2,YS=512*1.2

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

num=2400l
data_namet=data_name[6]

;sys=0 ;for windows system
sys=1 ;for linux   system

nm_ind='05'

;-------------------------file name-------------------------
;;for prep data
if sys eq 0 then filename=file_search('.\'+data_namet+'\prep\*.sav',count=count)
if sys eq 1 then filename=file_search('./'+data_namet+'/prep/*.sav',count=count)

print,count
;count=51

if sys eq 0 then image_sav='.\imagetest\'+data_namet+'\phase_rm'+nm_ind+'\'
if sys eq 1 then image_sav='./imagetest/'+data_namet+'/phase_rm'+nm_ind+'/'
file_mkdir,image_sav

if sys eq 0 then data_sav='.\sav\'+data_namet+'\rm_time'+nm_ind+'\'
if sys eq 1 then data_sav='./sav/'+data_namet+'/rm_time'+nm_ind+'/'
file_mkdir,data_sav

;-------------------------plot for 1 min-------------------------
;kk_st=0 ;for which file to start
;kk_num=1

kk_int=1 ;always set to be 1. ;for 'kk_int' min per phase image in total
num=num/kk_int

for kk_st=0,0 do begin ;count/kk_int-1

  cor_data=make_array(num*kk_int,40,40,2)
  
  for kk=0,kk_int-1 do begin
    
    ff=kk_st*kk_int+kk;kk_st*kk_int,(kk_st+1)*kk_int-1
    restore,filename[ff];,/ver
    ;;for prep data
    tmp0=rstrpos(filename[ff],'cor_data_')
    file=strmid(filename[ff],tmp0,60)
    print,file
    cor_ind=findgen(num)*kk_int
    cor_data[kk*num:((kk+1)*num-1),*,*,*]=cor_data_prep[cor_ind,*,*,*]
    gps_time = gps_time_prep
    
    if kk eq 0 then start_time=gps_time[0,*]  ;time_gps_to_spec(gps_time[0,*])
    if kk eq kk_int-1 then end_time=gps_time[num*kk_int-1,*]  ;time_gps_to_spec(gps_time[num*kk_int-1,*])

  endfor
  
  amp=sqrt((cor_data[*, *, *, 1])^2+(cor_data[*, *, *, 0])^2)
  phase=atan(cor_data[*, *, *, 1],cor_data[*, *, *, 0])*180./!pi
  
  timerg=strcompress(start_time[3],/remove_all)+strcompress(start_time[4],/remove_all) $
    +'_to_'+strcompress(end_time[3],/remove_all)+strcompress(end_time[4],/remove_all)

  ;-------------------------selecting begins-------------------------
  loadct,0
  !P.multi=[0,8,5];[0,8,5]
  
  ;mode=0 ;;plot amptitute and phase
  ;mode=1 ;;plot amptitute
  mode=2 ;;plot phase
  
  pp=['amp_phase','amp','phase','no plot']
  pp=pp[mode]
  
  ;;for num=2400
  rm_be=10 ;remove points at the beginning and in the end of the original phase data
  rm_max=100 ;permit select max number of every phase image
  ind_sm=100 ;'smooth index'
  ind_std=3 ;'stddev index'
  
  rm_be=10 & rm_max=100 & ind_sm=50 & ind_std=3 ;;if num eq 2400 then 
  ;rm_be=10 & rm_max=10 & ind_sm=50 & ind_std=3 ;; if num eq 600 then 
  ;rm_be=1 & rm_max=4 & ind_sm=20 & ind_std=2 ;; ;if num eq 60 then 
  
  index=intarr(rm_max,40,40)

  for jj=0,0 do begin
  
    for ii=0,39 do begin
    
      ph=phase[rm_be:(num-1-rm_be),jj,ii]
      ph_sm=ph-smooth(ph,ind_sm) ;you can set 'smooth index'
      ph_sm=abs(ph_sm)
      
      ind_temp=where(ph_sm gt ind_std*stddev(ph_sm)) ;you can set 'stddev index'
      if (min(ind_temp) gt 0) and (n_elements(ind_temp) lt rm_max) then begin
        index[0:((size(ind_temp))[1]-1),jj,ii]=ind_temp[*]
        ;print,'[',jj,ii,']',n_elements(ind_temp),n_elements(where(index[*,jj,ii] ne 0))
      endif
      
      amp_plot=amp[rm_be:(num-1-rm_be),jj,ii]
      amp_plot=alog(amp_plot)
      
      if jj eq 0 then begin
        
        if mode eq 0 then begin
          ;;;for both amptitute and phase plot
          if ii eq 0 then print,'plot amptitute and phase'
          plot,amp_plot $
            ,thick=1.5,ystyle=1+8,xstyle=1,ytitle='  ',yrange=[5,13],/nodata,color=cgcolor('red')
          plot,amp_plot $
            ,thick=2,ystyle=4,xstyle=4+1,ytitle='  ',yrange=[5,13],/noerase,color=cgcolor('red');,yrange=[0,max(amp[*,0,ii])];
          plot,ph,title='A'+strcompress(jj,/remove_all)+' and A'+strcompress(ii,/remove_all)+' / '+timerg $
            ,thick=1.5,/noerase,ystyle=4,xstyle=1+4,color=cgcolor('black'),yr=[-200,200] ;;[-4,4]
          axis,yaxis=1,ytitle='A '+strcompress(ii,/remove_all),color=cgcolor('black');,charthick=1.5,charsize=0.35
        endif
        if mode eq 1 then begin
          ;;just for amptitute plot
          if ii eq 0 then print,'plot amptitute'
          plot,amp_plot,thick=1.5,ystyle=1,xstyle=1,ytitle='Amp',yrange=[5,13] $;,yrange=[1e3,3e4],/ylog
            ,title='A'+strcompress(jj,/remove_all)+' and A'+strcompress(ii,/remove_all)+' / '+timerg
        endif
        if mode eq 2 then begin
          ;;just for phase plot
          if ii eq 0 then print,'plot phase'
          plot,ph,thick=1.5,ystyle=1,xstyle=1,ytitle='Phase',yrange=[-200,200],psym=6 $;
            ,title='A'+strcompress(jj,/remove_all)+' and A'+strcompress(ii,/remove_all)+' / '+timerg;+' / rms'+string(phase_rms[jj,ii])
          oplot,ph_sm,color=cgcolor('red'),psym=6
          ;oplot,[0,sz-1],[-180,-180],color=cgcolor('red'),linestyle=1
          ;oplot,[0,sz-1],[180,180],color=cgcolor('red'),linestyle=1
          if n_elements(where(index[*,jj,ii] gt 0)) gt 0 then begin
            for kk=0,rm_max-1 do begin;(size(ind_temp))[1]-1
              plots,[index[kk,jj,ii],index[kk,jj,ii]],[-4,ph[index[kk,jj,ii]]],color=cgcolor('green')
            endfor
          endif
        endif
        if mode eq 3 then print,'no phase ploting........'
            
      endif
    endfor
    
    write_jpeg,image_sav+'phase_A'+STRCOMPRESS(jj,/REMOVE_ALL)+'_A'+STRCOMPRESS(ii,/REMOVE_ALL)+'_'+strmid(file,0,46)+'_01.JPG',TVRD(TRUE=3),TRUE=3
    
  endfor
  
  temp=where(index ne 0)
  ;size(temp)
  value_temp=index[temp]
  rm_time_temp=value_temp[sort(value_temp)]
  help,rm_time_temp
  
  sz=(size(rm_time_temp))[1]
  
  rm_time1=intarr(sz)
  rm_time1[0]=rm_time_temp[0]
  for tt=1, sz-1 do begin
  
    if rm_time_temp[tt] ne rm_time_temp[tt-1] then rm_time1[tt]=rm_time_temp[tt]
    
  endfor
  temp=where(rm_time1 ne 0)
  rm_time2=rm_time1[temp]
  help,rm_time2
  
  sz=(size(rm_time2))[1]
  
  rm_time=intarr(sz+2*rm_be)
  rm_time[0:(rm_be-1)]=indgen(rm_be)
  rm_time[rm_be:(sz+rm_be-1)]=rm_be+rm_time2
  rm_time[(sz+rm_be):(sz+2*rm_be-1)]=reverse(num-indgen(rm_be)-1)
  
  help,rm_time
  
  save,rm_time,filename=data_sav+'rm_time_'+file
  
  print,'rm_time has saved..........'

endfor

;-------------------------ending-------------------------

print,'OK'
end