;.r phase_rm_time.pro
;.compile package_muser.pro
;
; Name :
; phase_rm_time
;
; Purpose :
; select bad points according to the phase plot of cross-correlation between two antennas.
; determine keep_table for baselines according to the RMS of phase plot.
; save the keep_table.sav
; select the bad points and the 'rm_time' again.
; save the bad points to rm_time to file rm_time_.sav.
; 
; Inputs :
; A file with cor_data_prep and gps_time_prep.
;
; Outputs :
; Amplitude and phase plot of one antenna cross-correlating with the other antennas.
; rm_time.sav
; 
; notes:
; one .sav file per minute.
; keep_table is different per minute.
; 
; History :
; Xingyao Chen, 12 Aug 2017
; 
;-------------------------begins-------------------------
DEVICE, DECOMPOSED = 0
Device, RETAIN=2
!P.FONT = 0
!P.Color=0
!P.Background=255
WINDOW,1,XS=512*3.2/1.5,YS=512*1.8/1.5

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

num=2400
for ooo=6,6 do begin
data_namet=data_name[ooo]

;sys=0 ;for windows system
sys=1 ;for linux   system

nm_ind='05'
nm_ind='17_test'

ch=4
;for ch=0,15 do begin
;nm_ind='sf'+strcompress(ch,/remove_all)

;-------------------------file name-------------------------
;;for prep data
if sys eq 0 then filename=file_search('.\'+data_namet+'\prep\*.sav',count=count)
if sys eq 1 then filename=file_search('./'+data_namet+'/prep/*.sav',count=count)

;;;for prepsf data
;if sys eq 0 then filename=file_search('.\'+data_namet+'\prepsf'+strcompress(ch,/remove_all)+'\*.sav',count=count)
;if sys eq 1 then filename=file_search('./'+data_namet+'/prepsf'+strcompress(ch,/remove_all)+'/*.sav',count=count)

print,count
;count=1
if sys eq 0 then image_sav='.\imagetest\'+data_namet+'\phase_rm'+nm_ind+'\'
if sys eq 1 then image_sav='./imagetest/'+data_namet+'/phase_rm'+nm_ind+'/'

if sys eq 0 then data_sav='.\sav\'+data_namet+'\rm_time'+nm_ind+'\'
if sys eq 1 then data_sav='./sav/'+data_namet+'/rm_time'+nm_ind+'/'

if sys eq 0 then keep_sav='.\sav\'+data_namet+'\keep_table'+nm_ind+'\'
if sys eq 1 then keep_sav='./sav/'+data_namet+'/keep_table'+nm_ind+'/'

file_mkdir,image_sav
file_mkdir,data_sav
file_mkdir,keep_sav

;-------------------------plot for 1 min-------------------------
;kk_st=0 ;for which file to start
;kk_num=1

kk_int=1 ;always set to be 1. ;for 'kk_int' min per phase image in total
num=num/kk_int

for kk_st=0,count/kk_int-1 do begin ;count/kk_int-1

  cor_data=make_array(num*kk_int,40,40,2)
  
  for kk=0,kk_int-1 do begin
    
    ff=kk_st*kk_int+kk;kk_st*kk_int,(kk_st+1)*kk_int-1
    restore,filename[ff];,/ver
    ;;for prep data
    tmp0=rstrpos(filename[ff],'cor_data_')
    file=strmid(filename[ff],tmp0,80)
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
    +'-'+strcompress(end_time[3],/remove_all)+strcompress(end_time[4],/remove_all)

  ;-------------------------selecting begins-------------------------
  loadct,0
  !P.multi=[0,8,5];[0,8,5]
  
  ;;for num=600 setting
  rm_be=10
  rm_max=10 ;permit select max number of every phase image
  ind_sm=50 ;'smooth index'
  ind_std=3 ;'stddev index'
  
  rm_be=10 & rm_max=60 & ind_sm=50 & ind_std=3 ;;if num eq 2400 then 
  ;rm_be=10 & rm_max=10 & ind_sm=50 & ind_std=3 ;; if num eq 600 then 
  ;rm_be=1 & rm_max=4 & ind_sm=20 & ind_std=2 ;; ;if num eq 60 then 
  
  ;;;-------------------------rm for the first time-------------------------
  rm_time=rm_time_index(phase, rm_be, rm_max, ind_sm, ind_std, index=rm_index_temp)
  help,rm_time
  
  ind_temp=indgen((size(phase))[1])
  remove,rm_time,ind_temp
  ph_rm_temp=phase[ind_temp,*,*]
  help,ph_rm_temp
  ph_rm=ph_rm_temp
  ;;;-------------------------calculate rms of the phase after rm-------------------------
  ;keep_table=keep_baseline(ph_rm_temp,20,phase_rms=phase_rms)  ;;phase_rms [40,40]
  ;;;keep_baseline, phase, rms_crit, phase_rms=phase_rms
  
  ;keep_table=keep_baseline(ph_rm_temp,40,phase_rms=phase_rms)  ;;phase_rms [40,40]
  keep_table=keep_baseline_mod(ph_rm_temp,30,100,2000.,1725,phase_rms=phase_rms)
  
  sz_tb=(size(keep_table))[2]
  help,keep_table
  
  ph_tb=make_array((size(ph_rm_temp))[1],(size(ph_rm_temp))[2],(size(ph_rm_temp))[3])
  for jj=0,39 do begin
    for ii=jj+1,39 do begin
      for tt=0,sz_tb-1 do begin
        if (jj eq keep_table[0,tt]) and (ii eq keep_table[1,tt]) then ph_tb[*,jj,ii]=ph_rm_temp[*,jj,ii]  
      endfor
    endfor
  endfor
  help,ph_tb
  
  ;;;-------------------------rm for the second time------------------------- 
  rm_time1=rm_time_index(ph_tb, 1, rm_max, ind_sm, ind_std, index=rm_index)
  help,rm_time1
  
  ind_temp=indgen((size(ph_tb))[1])
  remove,rm_time1,ind_temp
  ph_rm=ph_tb[ind_temp,*,*]
  help,ph_rm
  
  ;;;-------------------------plot------------------------- 
  for jj=0,0 do begin
    for ii=0,39 do begin
      plot,phase[*,jj,ii],thick=1.5,yr=[-200,200],YSTYLE=1,XSTYLE=1,color=0, $
        title='A'+STRCOMPRESS(jj,/REMOVE_ALL)+'-A'+STRCOMPRESS(ii,/REMOVE_ALL)+ $
        '/'+timerg+'/RMS-'+STRCOMPRESS(round(phase_rms[jj,ii]),/REMOVE_ALL)
      ;axis,yaxis=1,ytitle='A'+STRCOMPRESS(ii+1,/REMOVE_ALL),color=cgColor('red');,charthick=1.5,charsize=0.35
      ;oplot,ph_rm[*,jj,ii],color=cgcolor('green')
      oplot,ph_tb[*,jj,ii],color=cgcolor('red')
    endfor
    
    write_jpeg,image_sav+'phase_A'+STRCOMPRESS(jj,/REMOVE_ALL)+'_A'+STRCOMPRESS(ii,/REMOVE_ALL)+'_'+strmid(file,0,80)+'_01.JPG',TVRD(TRUE=3),TRUE=3
      
  endfor
  
  save,rm_time,rm_time1,filename=data_sav+'rm_time_'+file
  save,keep_table,filename=keep_sav+'keep_table_'+file
  print,'rm_time & keep_table has saved..........'
endfor

;endfor
;-------------------------ending-------------------------
endfor
print,data_sav
print,'OK'
end