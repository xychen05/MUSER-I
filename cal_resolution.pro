;.r cal_resolution.pro
;.compile package_muser.pro
;
; Purpose :
; Prepare for the self-correlation calibration.
; 
; History :
; Wei Wang, DATA_PRE_PROCESS.PRO
; Xingyao Chen, 24 Nov 2017
;-------------------------begins-------------------------
;;============antenna position
ant_pos=fltarr(44,3)
openr,lun1,'ANT_POS.TXT',/get_lun
for i=0,39 do begin
  readf,lun1,x,y,z
  ant_pos[i,*]=[x,y,z]
endfor
free_lun,lun1

;;============MUSER
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

sys=0 ;for windows system
sys=1 ;for linux   system

ktkt1=['_test_20_100_1500','_test_20_100_2000', $
  '_test_20_180_1500','_test_20_180_2000', $
  '_test_30_100_1500','_test_30_100_2000', $
  '_test_30_180_1500','_test_30_180_2000', $
  '_test_30','_test_40','_test_20']
  
ktkt=ktkt1[10]

print,'For Left Polarization:'
;print,'For Right Polarization:'
for ooo=1,1 do begin
  data_namet=data_name[ooo*2+4]
  
  for ch=4,4 do begin;0,15
    freqch=(1200+400*ooo+ch*25+25)
    freqchstr=strcompress(freqch,/remove_all)
    nm_ind='sf'+strcompress(ch,/remove_all)
    
    nm_ind='17'+ktkt
    
    if sys eq 0 then filename_kp=file_search('.\sav\'+data_namet+'\keep_table'+nm_ind+'.sav')
    if sys eq 1 then filename_kp=file_search('./sav/'+data_namet+'/keep_table'+nm_ind+'.sav')
    restore,filename_kp[0];,/ver
    print,filename_kp[0]
    help,keep_table
    sz=(size(keep_table))[2]
    dists=make_array(sz)
    for dd=0,sz-1 do begin
      a1=keep_table[0,dd] & a2=keep_table[1,dd]
      dists[dd]=sqrt((ant_pos[a1,0]-ant_pos[a2,0])^2+(ant_pos[a1,1]-ant_pos[a2,1])^2+(ant_pos[a1,2]-ant_pos[a2,2])^2)
      
    endfor
    freq=(1200+400*ooo+ch*25+25)/1000.
    resolution=0.3/freq/max(dists)*57.3*3600.
    print,freq,' GHz ',max(dists),' m ',resolution
  endfor
  
endfor


print,'hahaha.........'
END
