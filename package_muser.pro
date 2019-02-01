;.compile package_muser.pro
; Name:
; package_muser
; Note :
; Please compile it before MUSER image and data processing.
; History :
; Chen X.Y, June 2017
;=========================================================================
;
;
;PRO weighting,u,v,vis
;PRO grid,Lin,Lout
;FUNCTION crt_keep_table, filename_kp, sub=sub
;FUNCTION gps_time_int, gps_time, int
;FUNCTION cor_int, cor_data, int
;FUNCTION keep_baseline, phase, rms_crit, phase_rms=phase_rms
;FUNCTION rm_time_index, phase, rm_be, rm_max, ind_sm, ind_std, index=index
;PRO gauss_points, coeff,width,xdir=xdir,ydir=ydir,lx1,lx2,ly1,ly2,lx3,lx4,ly3,ly4,lx5,lx6,ly5,ly6
;FUNCTION slice_pixel_index_smp, xs1, xs2, ys1, ys2
;FUNCTION sl_region_index, p1,p2,p3,p4
;FUNCTION Polax_rot, Inputimg, Angle, Xcim, Ycim, XCEN=Xcen, YCEN=Ycen
;PRO img_rot,inputimg,timeh,timem,outputimg
;PRO img2dirty,inputimg,dbeam,imgconv
;PRO center_fit,inputimg,centx,centy
;PRO center_fit_d2d,inputimg1,inputimg2,centx,centy
;PRO norh_data_prep,filenamen,pix,rsun_radio,imgnorh
;PRO Mov_pos,vis,u,v,Dra,Ddec
;PRO add_disk,inputimage,factor,solar_center,solar_disk
;PRO contour_trim,Input_image,trim,output_image,region
;PRO my_convol,inputimage,funct,outputimage
;FUNCTION time_obs_to_string,obs_time
;FUNCTION time_gps_to_obs,gps_time
;FUNCTION time_gps_to_string,gps_time
;FUNCTION time_gps_to_spec,gps_time
;PRO muser2fits, inputfile, outputfile, pixel, fovaia=fovaia, left=left, right=right, lr=lr
;PRO D2B,X,B
;PRO READTIME,GPST,GPSTIME
;PRO signed,a
;FUNCTION interpolate_bessel,a,b
;PRO Cal_uvw,day,bjt,u,v,w
;PRO flag_baseline,keep_table,cor_data1,vis,u,v
;PRO FLAG_DATA,FLAG_TABLE,COR_DATA1,VIS,U,V
;PRO FFT_image,U,V,VIS,IMG
;PRO FFT_beam,U,V,BEAM
;PRO FFT_IMAGE,U,V,VIS,BEAM,IMG
;FUNCTION sdo_aia_scale_hdr,
;PRO aia_data_prep,filenamea,pix,rsun_radio,imgaia
;
;=========================================================================
; Purpose :
; This procedure is prepared for MUSER image and data processing.
; Please compile this package_muser.pro before MUSER image and data processing.
;
; History :
; Xingyao Chen, June 2017
;
;
;
;
;
;

;=========================================================================



;=========================================================================


;=========================================================================


;=========================================================================


;=========================================================================


;=========================================================================


;=========================================================================
FUNCTION keep_baseline_mod, phase, rms_crit, rms_crit_l, dist_crit, freqs, phase_rms=phase_rms
  ; Name:
  ; keep_baseline_mod
  ;
  ; Purpose :
  ; This function can select good baselines from the phase plot.
  ; For short baseline and long baseline, respectively
  ; 
  ; Inputs :
  ; phase
  ; rms_crit, you can set 'the value of rms'
  ; rms_crit_l, for long baseline
  ; 
  ; Outputs :
  ; keep_table: [2,*]
  ; phase_rms: [40,40]
  ;
  ; Examples :
  ; keep_table=keep_baseline_mod(phase,20,90,1500.,1725,phase_rms=phase_rms)
  ;
  ; History :
  ; Xingyao Chen, 31 Jul 2018
  ;-------------------------begins-------------------------
  sz=(size(phase))[1]
  jsz=(size(phase))[2]
  isz=(size(phase))[3]
  
  keep_table1=intarr(2,720)
  phase_rms=make_array(jsz,isz)
  
  dists1=make_array(720)
  
  ;;============antenna position
  ant_pos=fltarr(44,3)
  openr,lun1,'ANT_POS.TXT',/get_lun
  for i=0,39 do begin
    readf,lun1,x,y,z
    ant_pos[i,*]=[x,y,z]
  endfor
  free_lun,lun1
  
  ;;============cal
  tt=0
  for jj=0,jsz-1 do begin;jsz-1
    for ii=0,isz-1 do begin
    
      temp=linfit(findgen(sz),phase[*,jj,ii],yfit=phase_fit)
      phase_rms[jj,ii]=stddev(phase[*,jj,ii]-phase_fit)
      
      a1=jj & a2=ii
      distij=sqrt((ant_pos[a1,0]-ant_pos[a2,0])^2+(ant_pos[a1,1]-ant_pos[a2,1])^2+(ant_pos[a1,2]-ant_pos[a2,2])^2)
  
      if (distij le dist_crit) and (phase_rms[jj,ii] le rms_crit) and (phase_rms[jj,ii] ne 0) then begin
        keep_table1[*,tt]=[jj,ii]
        dists1[tt]=distij
        tt=tt+1
      endif
      
      if (distij gt dist_crit) and (phase_rms[jj,ii] le rms_crit_l) and (phase_rms[jj,ii] ne 0) then begin
        keep_table1[*,tt]=[jj,ii]
        dists1[tt]=distij
        tt=tt+1
      endif 
 
    endfor
  endfor
  
  keep_table=keep_table1[*,0:(tt-1)]
  dists=dists1[0:(tt-1)]
  
  resolution=0.3/freqs/max(dists)*57.3*3600.
  print,freqs,' GHz ',max(dists),' m ',resolution
  
  help,keep_table
  return,keep_table
end

;=========================================================================
PRO weighting,u,v,vis
; History :
; Wei Wang, 2018
;-------------------------begins-------------------------
  NN=size(u)
  N=NN(1)
  w=fltarr(n,n)
  fov=1./180*!pi
  r=2*1/fov
  for i=0,n-1 do begin
    for j=0,n-1 do begin
      count=1
      for ii=0,n-1 do begin
        for jj=0,n-1 do begin
          if (i ne j) and (i ne ii) and (j ne jj) and (ABS(vis[i,j]) ne 0) and (ABS(vis[ii,jj]) ne 0) then begin
            d=sqrt((U[i,j]-U[ii,jj])^2+(V[i,j]-V[ii,jj])^2)
            if d lt r then begin
              count=count+1
              ;if count gt 10 then print,i,j,ii,jj,d,count
            endif
          endif
        endfor
      endfor
      w[i,j]=1;/(count)
    endfor
  endfor
  vis=vis*w
  ;PRINT,'OK'
END

;=========================================================================
PRO grid,Lin,Lout
; Name:
; grid
;
; Purpose :
; Do gridding while uvw calibration.
;
; History :
; Wei Wang, 2018
;-------------------------begins-------------------------
  select_func=1
  case(select_func) of
    0: begin ;pillbox
      cfunc=0
    end
    1: begin ; exponential
      w=1.
      alfa=2.
      Lout=exp(-(Lin/w)^alfa)
    end
    2: begin ; sinc
      w=1.
      Lout=sin(Lin/w)/Lin/w
    end
    3: begin ; exponential times sinc
      w1=1.;2.52
      w2=1.55
      alfa=2.
      Lout=exp(-abs(Lin/w1)^alfa)*sin(Lin/w2)/Lin/w2
      if Lin eq 0 then Lout=1.
    end
    else: begin
    end
  endcase
  
END

;=========================================================================
FUNCTION crt_keep_table, filename_kp, sub=sub
  ; Name:
  ; makeup_keep_table
  ;
  ; Purpose :
  ; This function can adjust multi baselines.
  ; Do summation and subtraction. 
  ;
  ; Examples :
  ; keep_table_sum=makeup_keep_table(filename_kp)
  ; keep_table_sub=makeup_keep_table(filename_kp,/sub)
  ; 
  ; History :
  ; Xingyao Chen, 30 Oct 2017
  ;-------------------------begins-------------------------
  szfile=(size(filename_kp))[1]
  if szfile le 1 then print,'Please input more than 2 file!'
  if szfile le 1 then GOTO,nafik
  
  nfile=szfile
  
  keep_table_org=intarr(2,780)
  keep_table_temp=intarr(2,780)
  restore,filename_kp[0];,/ver
  help,keep_table
  szkp0=(size(keep_table))[2]
  keep_table_org[*,0:(szkp0-1)]=keep_table[*,*]
  
  ind=0
  for ii=1,nfile-1 do begin ;nfile-1 ; for which file
    restore,filename_kp[ii]
    szkp=(size(keep_table))[2]
    ind=ind
    szkp11=szkp0+ind
    for mm=0,szkp-1 do begin
      mark='save'
      for nn=0,szkp11-1 do begin
        if (keep_table[0,mm] eq keep_table_org[0,nn]) and (keep_table[1,mm] eq keep_table_org[1,nn]) then mark='no save'
      endfor
      if mark eq 'save' then ind=ind+1
      if mark eq 'save' then keep_table_org[*,(szkp0-1+ind)]=keep_table[*,mm]
    endfor
  endfor
  keep_table_sum=keep_table_org[*,0:(szkp0-1+ind)]
  help,keep_table_sum
  keep_table_ss=keep_table_sum
  
  if keyword_set( sub) then begin
    keep_table_total=intarr(nfile,2,(szkp0+ind))
    for ii=0,nfile-1 do begin
      restore,filename_kp[ii]
      szkp=(size(keep_table))[2]
      keep_table_total[ii,0,0:(szkp-1)]=keep_table[0,*]
      keep_table_total[ii,1,0:(szkp-1)]=keep_table[1,*]
    endfor
    keep_table_sub=intarr(2,(szkp0+ind))
    indss=0
    for xx=0,(szkp0+ind)-1 do begin
      indss=indss
      indtt=0
      for yy=0,nfile-1 do begin
        for zz=0,(szkp0+ind)-1 do begin
          if (keep_table_sum[0,xx] eq keep_table_total[yy,0,zz]) and (keep_table_sum[1,xx] eq keep_table_total[yy,1,zz]) then indtt=indtt+1
        endfor
      endfor
      if indtt eq nfile then keep_table_sub[*,indss]=keep_table_sum[*,xx]
      if indtt eq nfile then indss=indss+1
    endfor
    keep_table_ss=keep_table_sub[*,0:(indss-1)]
    keep_table_sub=keep_table_sub
    help,keep_table_sub
  endif
  
  return,keep_table_ss
  nafik:
  
END
;=========================================================================
FUNCTION gps_time_int, gps_time, int
  ; Name:
  ; gps_time_int
  ;
  ; Purpose :
  ; Same with cor_int.
  ;
  ; History :
  ; Xingyao Chen, 09 Oct 2017
  
  sz=(size(gps_time))[1]/int
  szi=(size(gps_time))[2]
  gps_time_it=intarr(int,szi)
  
  for l=0,int-1 do begin
    gps_time_it[l,*]=gps_time[l*sz,*]
  endfor
  
  return,gps_time_it
end

;=========================================================================
FUNCTION cor_int, cor_data, int
  ; Name:
  ; cor_int
  ;
  ; Purpose :
  ; Do the Re and Im part intergration of cor_data.
  ; 
  ; Inputs :
  ; int : if it is set to be 60, it will be split into 60 for one file.
  ; 
  ; Examples :
  ; cor_data_it=cor_int(cor_data,60)  ; cor_data_it[60,*,*,*]
  ; 
  ; History :
  ; Xingyao Chen, 09 Oct 2017
  
  sz=(size(cor_data))[1]/int
  szi=(size(cor_data))[2]
  szj=(size(cor_data))[3]
  szk=(size(cor_data))[4]
  
  cor_data_it=dblarr(int,szi,szj,szk)
  
  for i=0,szi-1 do begin
    for j=0,szj-1 do begin
      for k=0,szk-1 do begin
        for l=0,int-1 do begin
          cor_data_it[l,i,j,k]=mean(cor_data[(l*sz):((l+1)*sz-1),i,j,k])
        endfor
      endfor
    endfor
  endfor

  return,cor_data_it
end

;=========================================================================
FUNCTION keep_baseline, phase, rms_crit, phase_rms=phase_rms
  ; Name:
  ; keep_baseline
  ;
  ; Purpose :
  ; This function can select good baselines from the phase plot.
  ;
  ; Inputs :
  ; phase
  ; rms_crit, you can set 'the value of rms'
  ;
  ; Outputs :
  ; keep_table: [2,*]
  ; phase_rms: [40,40]
  ;
  ; Examples :
  ; keep_table=keep_baseline(phase,20,phase_rms=phase_rms)
  ;
  ; History :
  ; Xingyao Chen, 12 Aug 2017
  ;-------------------------begins------------------------- 
  sz=(size(phase))[1]
  jsz=(size(phase))[2]
  isz=(size(phase))[3]
  
  keep_table1=intarr(2,720)
  phase_rms=make_array(jsz,isz)
  
  tt=0
  for jj=0,jsz-1 do begin;jsz-1
    for ii=0,isz-1 do begin
      
      temp=linfit(findgen(sz),phase[*,jj,ii],yfit=phase_fit)
      phase_rms[jj,ii]=stddev(phase[*,jj,ii]-phase_fit)
      
      if (phase_rms[jj,ii] le rms_crit) and (phase_rms[jj,ii] ne 0) then begin
        keep_table1[*,tt]=[jj,ii]
        tt=tt+1
      endif
      
    endfor
  endfor
  
  keep_table=keep_table1[*,0:(tt-1)]
  
  help,keep_table
  return,keep_table
end

;=========================================================================
FUNCTION rm_time_index, phase, rm_be, rm_max, ind_sm, ind_std, index=index
  ; Name:
  ; rm_time_index
  ;
  ; Purpose :
  ; This function can select bad points from the phase plot and remove them.
  ;
  ; Inputs :
  ; phase
  ; rm_be, remove points at the beginning and in the end of the original phase data
  ; rm_max, permit select max number of every phase image
  ; ind_sm, you can set 'smooth index'
  ; ind_std, you can set 'stddev index'
  ;
  ; Outputs :
  ; rm_time
  ;
  ; Examples :
  ; rm_time=rm_time_index(phase,10,20,30,3)
  ; rm_time=rm_time_index(phase, rm_be, rm_max, ind_sm, ind_std)
  ;
  ; History :
  ; Xingyao Chen, 25 July 2017
  ;-------------------------begins-------------------------
  
  rm_be=rm_be ;remove points at the beginning and in the end of the original phase data
  rm_max=rm_max ;permit select max number of every phase image
  ind_sm=ind_sm ;you can set 'smooth index'
  ind_std=ind_std ;you can set 'stddev index'
  
  num=(size(phase))[1]
  jsz=(size(phase))[2]
  isz=(size(phase))[3]
  
  index=intarr(rm_max,jsz,isz)
  
  for jj=0,jsz-1 do begin;jsz-1
    for ii=0,isz-1 do begin
    
      ph=phase[rm_be:(num-1-rm_be),jj,ii]
      
      temp=linfit(findgen(num-2*rm_be),ph,yfit=ph_fit)
      ph_sm=ph-ph_fit
      ph_sm=ph_sm-smooth(ph_sm,ind_sm) ;you can set 'smooth index'
      ph_sm=abs(ph_sm)
      
;      ph_sm=ph-smooth(ph,ind_sm) ;you can set 'smooth index'
;      ph_sm=abs(ph_sm)
      
      ind_temp=where(ph_sm gt ind_std*stddev(ph_sm)) ;you can set 'stddev index'
      
      if (min(ind_temp) gt 0) and (n_elements(ind_temp) lt rm_max) then begin
        index[0:((size(ind_temp))[1]-1),jj,ii]=ind_temp[*]
      endif
      
    endfor
  endfor
  
  temp=where(index ne 0)
  ;size(temp)
  value_temp=index[temp]
  rm_time_temp=value_temp[sort(value_temp)]
  ;help,rm_time_temp
  sz=(size(rm_time_temp))[1]
  
  rm_time1=intarr(sz)
  rm_time1[0]=rm_time_temp[0]
  for tt=1, sz-1 do begin
    if rm_time_temp[tt] ne rm_time_temp[tt-1] then rm_time1[tt]=rm_time_temp[tt]
  endfor
  temp=where(rm_time1 ne 0)
  rm_time2=rm_time1[temp]
  ;help,rm_time2
  
  sz=(size(rm_time2))[1]
  
  rm_time=intarr(sz+2*rm_be)
  rm_time[0:(rm_be-1)]=indgen(rm_be)
  rm_time[rm_be:(sz+rm_be-1)]=rm_time2+rm_be
  rm_time[(sz+rm_be):(sz+2*rm_be-1)]=reverse(num-indgen(rm_be)-1)
  
  print,'rm_time has prepared..........'
  return,rm_time
  
end

;=========================================================================
PRO gauss_points, coeff,width,xdir=xdir,ydir=ydir,lx1,lx2,ly1,ly2,lx3,lx4,ly3,ly4,lx5,lx6,ly5,ly6
  ; Name:
  ; gauss_points
  ;
  ; Purpose :
  ; This function gives 6 points.
  ;
  ; Inputs :
  ; coeff after function gauss2dfit.
  ;
  ; Outputs :
  ; pixel index of 6 points.
  ;
  ; Examples :
  ; gauss_points, coeff,width=2.0,/ydir
  ;
  ; History :
  ; Xingyao Chen, 30 July 2017
  ;-------------------------begins-------------------------
  if keyword_set( xdir) then begin
    print,'X direction!'
    theta=coeff[6]
    width=width ;;you can change here.
    coeff2=coeff[2]
    coeff3=coeff[3]
    coeff4=coeff[4]
    coeff5=coeff[5]
    widthx=width*coeff3*sin(theta)
    widthy=width*coeff3*cos(theta)
    lx1=coeff4-widthx & ly1=coeff5-widthy
    lx2=coeff4+widthx & ly2=coeff5+widthy
    print,lx1,lx2,ly1,ly2
    ;plots,[lx1,lx2],[ly1,ly2],symsize=3,color=255
    lx3=coeff4-width*coeff2*cos(-theta)+widthx & ly3=coeff5-width*coeff2*sin(-theta)+widthy
    lx4=coeff4+width*coeff2*cos(-theta)+widthx & ly4=coeff5+width*coeff2*sin(-theta)+widthy
    print,lx3,lx4,ly3,ly4
    ;plots,[lx3,lx4],[ly3,ly4],symsize=3,color=255
    lx5=coeff4-width*coeff2*cos(-theta)-widthx & ly5=coeff5-width*coeff2*sin(-theta)-widthy
    lx6=coeff4+width*coeff2*cos(-theta)-widthx & ly6=coeff5+width*coeff2*sin(-theta)-widthy
    print,lx5,lx6,ly5,ly6
    ;plots,[lx5,lx6],[ly5,ly6],symsize=3,color=255
  endif
  if keyword_set( ydir) then begin
    print,'Y direction!'
    theta=abs(coeff[6])
    width=width ;;you can change here.
    coeff3=coeff[2];;x width
    coeff2=coeff[3];;y
    coeff4=coeff[4];;x center
    coeff5=coeff[5];;y
    widthy=width*coeff3*sin(theta)
    widthx=width*coeff3*cos(theta)
    lx1=coeff4-widthx & ly1=coeff5-widthy
    lx2=coeff4+widthx & ly2=coeff5+widthy
    print,lx1,lx2,ly1,ly2
    ;plots,[lx1,lx2],[ly1,ly2],symsize=3,color=255
    lx3=coeff4+width*coeff2*sin(theta)+widthx & ly3=coeff5-width*coeff2*cos(theta)+widthy
    lx4=coeff4+width*coeff2*sin(theta)-widthx & ly4=coeff5-width*coeff2*cos(theta)-widthy
    print,lx3,lx4,ly3,ly4
    ;plots,[lx3,lx4],[ly3,ly4],symsize=3,color=255
    lx5=coeff4-width*coeff2*sin(theta)+widthx & ly5=coeff5+width*coeff2*cos(theta)+widthy
    lx6=coeff4-width*coeff2*sin(theta)-widthx & ly6=coeff5+width*coeff2*cos(theta)-widthy
    print,lx5,lx6,ly5,ly6
  endif
end

;=========================================================================
FUNCTION slice_pixel_index_smp, xs1, xs2, ys1, ys2
  ;from axis to pixel for straight lines, from y1 to y2 and x1 to x2
  xs1=float(xs1) & xs2=float(xs2) & ys1=float(ys1) & ys1=float(ys1)
  ;case: parallel to y-axis
  if xs2 ne xs1 then begin
    ks = (ys2-ys1)/(xs2-xs1)
  endif else begin
    num=ceil(abs(ys2-ys1+1))
    point=make_array(2,num)
    ys=round(ys1)
    for ii=0,num-1 do begin
      point[0,ii]=xs1
      point[1,ii]=ys+ii
    endfor
    GOTO,nafik
  endelse
  ;case: parallel to x-axis
  if ks eq 0 then begin
    num=ceil(abs(xs2-xs1+1))
    point=make_array(2,num)
    xs=round(xs1)
    for ii=0,num-1 do begin
      point[0,ii]=xs+ii
      point[1,ii]=ys1
    endfor
  endif
  
  ;case: x plays a leading role
  if ks gt 0 then begin
    num=ceil(abs(xs2-xs1+1))
    point=make_array(2,num)
    xs=round(xs1)
    for ii=0,num-1 do begin
      ys=ys1+ii*ks
      point[0,ii]=xs+ii
      point[1,ii]=round(ys)
    endfor
  endif
  
  ;case: x plays a leading role
  if ks lt 0 then begin
    num=ceil(abs(xs2-xs1+1))
    point=make_array(2,num)
    xs=round(xs1)
    for ii=0,num-1 do begin
      ys=ys1+ii*ks
      point[0,ii]=xs+ii
      point[1,ii]=round(ys)
    endfor
  endif
  nafik:
  RETURN, point
  
END

;-------------------------basic parameters-------------------------
FUNCTION sl_region_index, p1,p2,p3,p4
  ; Name:
  ; sl_region_index
  ;
  ; Purpose :
  ; This function gives the indexes of selected region, which is indecated by four points.
  ;
  ; Inputs :
  ; points1,points2,points3,points4
  ;
  ; Outputs :
  ; index.
  ;
  ; Examples :
  ; pixel_index=sl_region_index([x1,y1],[x2,y2],[x3,y3],[x4,y4])
  ;
  ; History :
  ; Xingyao Chen, 25 July 2017
  ;-------------------------begins-------------------------
  p1=round(p1) & p2=round(p2) & p3=round(p3) & p4=round(p4)
  
  p_temp=[[p1],[p2],[p3],[p4]]
  
  x1=p1[0] & x2=p2[0] & x3=p3[0] & x4=p4[0]
  
  temp=[x1,x2,x3,x4]
  temp1=sort(temp)
  
  p1=p_temp[*,temp1[0]] & p2=p_temp[*,temp1[1]] & p3=p_temp[*,temp1[2]] & p4=p_temp[*,temp1[3]]
  
  x1=p1[0] & x2=p2[0] & x3=p3[0] & x4=p4[0]
  
  y1=p1[1] & y2=p2[1] & y3=p3[1] & y4=p4[1]
  
  tempy=[y1,y2,y3,y4]
  
  print,x1,x2,x3,x4
  print,y1,y2,y3,y4
  ;plots, [p1[0],p2[0],p3[0],p4[0],p1[0]], [p1[1],p2[1],p3[1],p4[1],p1[1]]
  
  point12=slice_pixel_index_smp(x1, x2, y1, y2)
  point13=slice_pixel_index_smp(x1, x3, y1, y3)
  
  point24=slice_pixel_index_smp(x2, x4, y2, y4)
  point34=slice_pixel_index_smp(x3, x4, y3, y4)
  
  point23=slice_pixel_index_smp(x2, x3, y2, y3)
  
  ;for line 1 & 2
  point14=slice_pixel_index_smp(x1, x4, y1, y4)
  index=make_array(ceil(x4-x1+1),ceil(max(tempy)-min(tempy)+1))
  
  ;-------------------------4 cases-------------------------
  if ((y2-point14[1,(x2-x1)]) eq 0) or ((y3-point14[1,(x3-x1)]) eq 0) then begin
    print,'error: three points on one line!'
    GOTO,nafik
  endif
  
  ;;for y2 > the line point14 & y3 > the line point14
  if ((y2-point14[1,(x2-x1)]) gt 0) and ((y3-point14[1,(x3-x1)]) gt 0) then begin
    for xx=x1,x2-1 do begin
      for yy=point14[1,(xx-x1)],point12[1,(xx-x1)] do begin
        index[(xx-x1),(yy-point14[1,(xx-x1)])]=yy
      endfor
    endfor
    for xx=x2,x3-1 do begin
      for yy=point14[1,(xx-x1)],point23[1,(xx-x2)] do begin
        index[(xx-x1),(yy-point14[1,(xx-x1)])]=yy
      endfor
    endfor
    for xx=x3,x4-1 do begin
      for yy=point14[1,(xx-x1)],point34[1,(xx-x3)] do begin
        index[(xx-x1),(yy-point14[1,(xx-x1)])]=yy
      endfor
    endfor
    print,'case 1...........'
  endif
  ;;for y2 > the line point14 & y3 < the line point14
  if ((y2-point14[1,(x2-x1)]) gt 0) and ((y3-point14[1,(x3-x1)]) lt 0) then begin
    for xx=x1,x2-1 do begin
      for yy=point13[1,(xx-x1)],point12[1,(xx-x1)] do begin
        index[(xx-x1),(yy-point13[1,(xx-x1)])]=yy
      endfor
    endfor
    for xx=x2,x3-1 do begin
      for yy=point13[1,(xx-x1)],point24[1,(xx-x2)] do begin
        index[(xx-x1),(yy-point13[1,(xx-x1)])]=yy
      endfor
    endfor
    for xx=x3,x4-1 do begin
      for yy=point34[1,(xx-x3)],point24[1,(xx-x2)] do begin
        index[(xx-x1),(yy-point34[1,(xx-x3)])]=yy
      endfor
    endfor
    print,'case 2...........'
  endif
  ;;for y2 < the line point14 & y3 > the line point14
  if ((y2-point14[1,(x2-x1)]) lt 0) and ((y3-point14[1,(x3-x1)]) gt 0) then begin
    for xx=x1,x2-1 do begin
      for yy=point12[1,(xx-x1)],point13[1,(xx-x1)] do begin
        index[(xx-x1),(yy-point12[1,(xx-x1)])]=yy
      endfor
    endfor
    for xx=x2,x3-1 do begin
      for yy=point24[1,(xx-x2)],point13[1,(xx-x1)] do begin
        index[(xx-x1),(yy-point24[1,(xx-x2)])]=yy
      endfor
    endfor
    for xx=x3,x4-1 do begin
      for yy=point24[1,(xx-x2)],point34[1,(xx-x3)] do begin
        index[(xx-x1),(yy-point24[1,(xx-x2)])]=yy
      endfor
    endfor
    print,'case 3...........'
  endif
  ;;for y2 < the line point14 & y3 < the line point14
  if ((y2-point14[1,(x2-x1)]) lt 0) and ((y3-point14[1,(x3-x1)]) lt 0) then begin
    for xx=x1,x2-1 do begin
      for yy=point12[1,(xx-x1)],point14[1,(xx-x1)] do begin
        index[(xx-x1),(yy-point12[1,(xx-x1)])]=yy
      endfor
    endfor
    for xx=x2,x3-1 do begin
      for yy=point23[1,(xx-x2)],point14[1,(xx-x1)] do begin
        index[(xx-x1),(yy-point23[1,(xx-x2)])]=yy
      endfor
    endfor
    for xx=x3,x4-1 do begin
      for yy=point34[1,(xx-x3)],point14[1,(xx-x1)] do begin
        index[(xx-x1),(yy-point34[1,(xx-x3)])]=yy
      endfor
    endfor
    print,'case 4...........'
  endif
  
  ;-------------------------basic parameters-------------------------
  sz=size(index)
  index_total=intarr(2,(sz[1])*(sz[2]))
  for ii=0,sz[1]-1 do begin
    index_total[0,(ii*sz[2]):((ii+1)*sz[2]-1)]=x1+ii
    index_total[1,(ii*sz[2]):((ii+1)*sz[2]-1)]=index[ii,*]
  endfor
  
  vert_x=index_total[0,*]
  vert_y=index_total[1,*]
  print,minmax(vert_x)
  print,minmax(vert_y)
  
  ind_temp=where(vert_y eq 0)
  
  remove,ind_temp,vert_x
  remove,ind_temp,vert_y
  
  index_region=intarr(2,(size(vert_x))[1])
  index_region[0,*]=vert_x
  index_region[1,*]=vert_y
  help,index_region
  ;  for ii=0,(size(vert_x))[1]-1 do begin
  ;    plots,[index_region[0,ii]],[index_region[1,ii]],psym=1,symsize=2,color=255
  ;  endfor
  
  nafik:
  RETURN, index_region
  
end

;=========================================================================
FUNCTION Polax_rot, inputimg, angle, xcim, ycim, xcen=xcen, ycen=ycen
; Name:
; Polax_rot
;
; Purpose :
; This function rotates an image.
;
; Inputs :
; inputimg = original image
; angle = angle of rotation in degrees.
; xcim, ycim = coordinates of the center of the image after rotation if xcen and ycen have been declared.
;
; Outputs :
; Image after rotation.
;
; History :
; Xingyao Chen, June 2017
;
  SZ= SIZE(Inputimg) & DIM = FIX(SQRT(SZ(1)^2 + SZ(2)^2))
  IMAGE= FLTARR (Dim, Dim)
  
  IF NOT KEYWORD_SET(XCEN) THEN Xcen=Sz(1)/2
  IF NOT KEYWORD_SET(YCEN) THEN Ycen=Sz(2)/2
  
  Xl = Dim/2 - Sz[1]/2
  Yl = Dim/2 - Sz[2]/2
  Image[Xl:Xl+Sz[1]-1,Yl:Yl+Sz[2]-1] = Inputimg
  
  IMAGE = ROT(Image, Angle)
  
  Outputimg=Image[Xl:Xl+Sz[1]-1,Yl:Yl+Sz[2]-1]
  
  XCIM=Xcen + Xl
  YCIM=Ycen + Yl
  
  RETURN, Outputimg
END
;=========================================================================
PRO img_rot,inputimg,timeh,timem,outputimg
  ;the image need to be move to the center after fitting the solor center.
  ;1217 degrees = 9.35
  ;1218 degrees = 8.89
  deg17=9.35
  deg18=8.89  
  deg1h=(deg18-deg17)/24.
  deg1m=deg1h/60.  
  ;print,deg1h,deg1m
  ;-0.0191667  -0.000319444
  deg=deg17+deg1h*timeh+deg1m*timem  
  outputimg=rot(inputimg,deg)
END

;=========================================================================
PRO img2dirty,inputimg,dbeam,imgconv
; Name:
; img2dirty
;
; Purpose :
; This procedure is for transforming AIA or NORH image into 'AIA or NORH dirty map'.
;
; Inputs :
; inputimg = inputimg
; dbeam = dirty beam.
;
; Outputs :
; imgconv = 'AIA or NORH dirty map'.
;
; History :
; Xingyao Chen, June 2017
;
pix=(size(inputimg))[1]

imgconv_temp=convol(dbeam,inputimg)
imgconv=make_array(pix,pix)
imgconv=imgconv_temp[pix/2:(pix/2+pix-1),pix/2:(pix/2+pix-1)]
imgconv1=reverse(imgconv,2)
imgconv=reverse(imgconv1,1)

END
;=========================================================================
PRO center_fit,inputimg,centx,centy
; Name :
; center_fit_test
; 
; Purpose :
; For testing.
; Display dirty images, residual image, clean image with fitted center.
; The convolution between residual image and ideal solar disk.
; 
; Inputs :
; inputimg = inputimg.
;
; Outputs :
; centx = solar center in X axis.
; centy = solar center in Y axis.
; 
; History :
; Xingyao Chen, June 2017
;

  pix=512
  img_temp=abs(inputimg)
  add_disk,img_temp,0.01,[PIX/2,PIX/2],solar_disk
  ;add_disk,inputimage,factor,solar_center,solar_disk
  fft_solar_disk=fft(solar_disk)
  fft_img=fft(inputimg)
  dft_img=fft(fft_solar_disk*fft_img,/inverse)
  dft_img=shift(dft_img,pix/2,pix/2)
  
  dft_ind=where(dft_img eq max(dft_img))
  
  centx=(dft_ind mod PIX);+PIX/2
  centy=dft_ind/PIX;+PIX/2
  ;print,'center x:',centx
  ;print,'center x:',centy

END

;=========================================================================
PRO center_fit_d2d,inputimg1,inputimg2,centx,centy
; Name :
; center_fit_d2d
;
; Purpose :
; To fit the solar center using the model of 'AIA or NORH dirty map'.
; Display img, solar disk, fft_solar_disk, dft_img.
;
; Inputs :
; inputimg1 = inputimg of MUSER dirty map.
; inputimg2 = inputimg of 'AIA or NORH dirty map'.
;
; Outputs :
; centx = solar center in X axis.
; centy = solar center in Y axis.
;
; Notes :
; give a dirty beam and dirty img firstly.
; aia/norh_data_prep,filenamea,pix,1.5,imgaia
; img2dirty,imgaia/norh,dbeam,imgconv
; center_fit_d2d,img,imgconv,centx,centy
;
; History :
; Xingyao Chen, June 2017
;
  inputimg1=abs(inputimg1)
  pix=(size(inputimg1))[1]
  solar_disk=inputimg2
  ;add_disk,inputimage,factor,solar_center,solar_disk
  fft_solar_disk=fft(solar_disk)
  fft_img=fft(inputimg1)
  dft_img=fft(fft_solar_disk*fft_img,/inverse);,/inverse
  dft_img=shift(dft_img,pix/2,pix/2)
  
  dft_img=abs(dft_img)
  
  dft_ind=where(dft_img eq max(dft_img))
  
  centx=(dft_ind mod PIX);+PIX/2
  centy=dft_ind/PIX;+PIX/2
  ;print,'center x:',centx
  ;print,'center y:',centy
  ;add_disk,img_temp,0.01,[centx,centy],solar_disk1
END
;=========================================================================
PRO norh_data_prep,filenamen,pix,rsun_radio,imgnorh
; Name :
; norh_data_prep
;
; Purpose :
; the data prep for NORH, transform from original size into the same size with dirty image of MUSER.
;
; Inputs :
; filenamea = the fits file path of the NORH data.
; pix = the pixel size of MUSER image, 512 for MUSER-I and 1024 for MUSER-II normally.
; rsun_radio = 1 or larger, the value set to be 0 beyond rsun_radio*R_SUN.
;
; Outputs :
; imgnorh = transform into the [pix,pix] image.
;
; Note :
; This procedure also recorded in package_muser.pro
; Please compile package_muser.pro firstly.
;
; History :
; Xingyao Chen, June 2017.

r_sun=1625.0607/4096*512

read_sdo, filenamen, indexn, datan, /UNCOMP_DELETE,/noshell,/use_shared_lib,/HIDE

xysize=indexn.NAXIS1
solar_center=[indexn.CRPIX1,indexn.CRPIX2]

datan_temp=make_array(xysize,xysize)
for i=0,xysize-1 do begin
  for j=0,xysize-1 do begin  
    r=sqrt((i-solar_center[0])^2+(j-solar_center[1])^2);pixel
    ;r=sqrt((i-xsize/2)^2+(j-ysize/2)^2);pixel
    if r le ((r_sun)*rsun_radio) then datan_temp[i,j]=datan[i,j]    
  endfor
endfor
dataa=datan_temp

pixnorh=round(indexn.NAXIS1*indexn.CDELT1/3600.*PIX)
imagenorh=congrid(datan,pixnorh,pixnorh,/interp)

edge=(pix-pixnorh)/2

imgnorh=make_array(pix,pix)
for xx=edge,edge+pixnorh-1 do begin
  for yy=edge,edge+pixnorh-1 do begin
    imgnorh[xx,yy]=imagenorh[xx-edge,yy-edge]
  endfor
endfor

imgnorh=imgnorh/max(imgnorh)

END

;=========================================================================

;=========================================================================
PRO Mov_pos,vis,u,v,Dra,Ddec
; Name :
; Mov_pos
;
; Purpose :
; Add the phase to the visibility in order to Move the dirty map to the image center.
; 
; Examples:
; mov_pos,vis,u,v,-(319-512/2)/512.,(357-512/2)/512. ;for Left polarization
; move the dirty map center [319,357] to image center [512/2,512/2] .
; 
; History :
; Wei Wang.
;
  dra=dra/180d*!pi
  ddec=ddec/180d*!pi
  vis=vis*exp(complex(0,-2*!pi*(u*dra+v*ddec)))

END
;=========================================================================
PRO add_disk,inputimage,factor,solar_center,solar_disk
; Name :
; add_disk
; 
; Purpose :
; Add a solar disk to the inputimage.
; 
; Explanation :
; This procedure is to design a ideal disk and plot the solar disk edge overlaid on a image.
;
; Inputs :
; inputimage = inputimage.
; factor = 0.1 or large, to enhance the brightness of solar disk.
; solar_center = [ , ], position of solar center in pixels.
; 
; Outputs :
; solar_disk = [ , ], image with solar disk, same size with inputimage.
; 
; History :
; Wei Wang
; Xingyao Chen, June 2017.
;

factor=factor/1.
solar_center=solar_center/1.

nsize=size(inputimage)
xsize=nsize[1]
ysize=nsize[2]

pi=3.1415926D
Fov=60.0d;41.9d ;arcemin
d_solar=979.684d*2/60d; arcmin;32.656
onepix=fov/xsize

solar_disk=dblarr(xsize,ysize)
limb_pix=1
limb=gaussian_function(limb_pix)
for i=0,xsize-1 do begin
  for j=0,ysize-1 do begin
    r=sqrt((i-solar_center[0])^2+(j-solar_center[1])^2);pixel
    
    if r le (d_solar/2/onepix) then solar_disk[i,j]=1d
    if (r gt (d_solar/2/onepix)) and (r-(d_solar/2/onepix) lt 3*limb_pix) then solar_disk[i,j]=limb(r-(d_solar/2/onepix)+3*limb_pix)
  endfor
endfor
inputimage=abs(inputimage)+factor*solar_disk

END
;=========================================================================
PRO contour_trim,Input_image,trim,output_image,region

  temp=input_image
  Max_image=max(temp)
  temp=temp/max_image
  xysize=size(temp)
  xsize=xysize[1]
  ysize=xysize[2]
  output_image=dblarr(xsize,ysize)
  
  for i=0,xsize-1 do begin
    for j=0,ysize-1 do begin
      if (temp[i,j] ge trim) and (region[i,j] eq 1) then output_image[i,j]=temp[i,j]
    endfor
  endfor

END
;=========================================================================
PRO my_convol,inputimage,funct,outputimage
; Name :
; my_convol
;
; Purpose :
; Do the convolution between inputimage and funct.
;
; History :
; Wei Wang.

  nsize=size(inputimage)
  pix=nsize[1]
  
  in_ft=fft(inputimage,/center)
  funct_ft=fft(funct,/center)
  
  out_ft=in_ft*funct_ft
  
  outputimage=real_part(fft(out_ft,/inverse,/center))
  outputimage=shift(outputimage,pix/2,pix/2)*pix*pix
  ;outputimage=rotate(outputimage,2)
END

;=========================================================================
FUNCTION time_obs_to_string,obs_time
  ; Name :
  ; time_obs_to_string
  ;
  ; Purpose :
  ; 17-Dec-2014 04:32:59.199 to 20141217_043259_199.
  ;
  ; History :
  ; Xingyao Chen, 02 Nov 2017.
  ;
  
  tmp0=rstrpos(obs_time,'.')
  
  if tmp0 eq 20 then begin
    year=strmid(obs_time,7,4)
    mon=strmid(obs_time,3,3)
    day=strmid(obs_time,0,2)
    hour=strmid(obs_time,12,2)
    min=strmid(obs_time,15,2)
    sec=strmid(obs_time,18,2)
    msec=strmid(obs_time,21,3)
    if mon eq 'Jan' then mon='01'
    if mon eq 'Feb' then mon='02'
    if mon eq 'Mar' then mon='03'
    if mon eq 'Apr' then mon='04'
    if mon eq 'May' then mon='05'
    if mon eq 'Jun' then mon='06'
    if mon eq 'Jul' then mon='07'
    if mon eq 'Aug' then mon='08'
    if mon eq 'Sep' then mon='09'
    if mon eq 'Oct' then mon='10'
    if mon eq 'Nov' then mon='11'
    if mon eq 'Dec' then mon='12'    
  endif
  if tmp0 eq 19 then begin
    year=strmid(obs_time,0,4)
    mon=strmid(obs_time,5,2)
    day=strmid(obs_time,8,2)
    hour=strmid(obs_time,11,2)
    min=strmid(obs_time,14,2)
    sec=strmid(obs_time,17,2)
    msec=strmid(obs_time,20,3)
  endif
  if tmp0 ne 19 and tmp0 ne 20 then begin
    print,'time transformation from obs to string is failed!'
    GOTO,nafik
  endif
  
  str_time=year+mon+day+'_'+hour+min+sec+'_'+msec
   
  return,str_time
  nafik:
END

;=========================================================================
FUNCTION time_gps_to_obs,gps_time1
  ; Name :
  ; time_gps_to_obs
  ;
  ; Purpose :
  ; change the gps_time format to specified format, eg. 2017-12-17T04:28:01.098, same as T_OBS of AIA
  ; complement 0.
  ;
  ; History :
  ; Xingyao Chen, July 2017.
  ;
  gps_time=gps_time1
  gps_time[*,3]=gps_time[*,3]-8
  s_gps=size(gps_time)
  gps_string=strarr(s_gps[1])
  
  for i=0,0 do begin;s_gps[1]-1
  
    if gps_time[i,1] lt 10 then begin
      mon='0'+strmid(gps_time[i,1],7,1)
    endif else begin
      mon=strmid(gps_time[i,1],6,2)
    endelse
    
    if gps_time[i,2] lt 10 then begin
      day='0'+strmid(gps_time[i,2],7,1)
    endif else begin
      day=strmid(gps_time[i,2],6,2)
    endelse
    
    if gps_time[i,3] lt 10 then begin
      hour='0'+strmid(gps_time[i,3],7,1)
    endif else begin
      hour=strmid(gps_time[i,3],6,2)
    endelse
    
    if gps_time[i,4] lt 10 then begin
      min='0'+strmid(gps_time[i,4],7,1)
    endif else begin
      min=strmid(gps_time[i,4],6,2)
    endelse
    
    if gps_time[i,5] lt 10 then begin
      sec='0'+strmid(gps_time[i,5],7,1)
    endif else begin
      sec=strmid(gps_time[i,5],6,2)
    endelse
    
    if (gps_time[i,6] ge 100) then msec=strmid(gps_time[i,6],5,3)
    if (gps_time[i,6] lt 100) and (gps_time[i,6] ge 10) then msec='0'+strmid(gps_time[i,6],6,2)
    if (gps_time[i,6] lt 10) then msec='00'+strmid(gps_time[i,6],7,1)
    
    gps_string[i]='20'+strmid(gps_time[i,0],6,2)+'-'+mon+'-'+day+'T'+hour+':'+min+':'+sec+'.'+msec
    
  endfor
  
  return,gps_string
END
;=========================================================================
FUNCTION time_gps_to_string,gps_time1
  ; Name :
  ; time_gps_to_string
  ;
  ; Purpose :
  ; change the gps_time format to specified format, eg. 20171217_042801_098.
  ; complement 0.
  ;
  ; History :
  ; Xingyao Chen, July 2017.
  ;
  gps_time=gps_time1
  gps_time[*,3]=gps_time[*,3]-8
  s_gps=size(gps_time)
  gps_string=strarr(s_gps[1])
  
  for i=0,s_gps[1]-1 do begin;
  
    if gps_time[i,1] lt 10 then begin
      mon='0'+strmid(gps_time[i,1],7,1)
    endif else begin
      mon=strmid(gps_time[i,1],6,2)
    endelse
    
    if gps_time[i,2] lt 10 then begin
      day='0'+strmid(gps_time[i,2],7,1)
    endif else begin
      day=strmid(gps_time[i,2],6,2)
    endelse
    
    if gps_time[i,3] lt 10 then begin
      hour='0'+strmid(gps_time[i,3],7,1)
    endif else begin
      hour=strmid(gps_time[i,3],6,2)
    endelse
    
    if gps_time[i,4] lt 10 then begin
      min='0'+strmid(gps_time[i,4],7,1)
    endif else begin
      min=strmid(gps_time[i,4],6,2)
    endelse
    
    if gps_time[i,5] lt 10 then begin
      sec='0'+strmid(gps_time[i,5],7,1)
    endif else begin
      sec=strmid(gps_time[i,5],6,2)
    endelse
    
    if (gps_time[i,6] ge 100) then msec=strmid(gps_time[i,6],5,3)
    if (gps_time[i,6] lt 100) and (gps_time[i,6] ge 10) then msec='0'+strmid(gps_time[i,6],6,2)
    if (gps_time[i,6] lt 10) then msec='00'+strmid(gps_time[i,6],7,1)
    
    gps_string[i]='20'+strmid(gps_time[i,0],6,2)+mon+day+'_'+hour+min+sec+'_'+msec
    
  endfor
  
  return,gps_string
END
;=========================================================================
FUNCTION time_gps_to_spec,gps_time1
; Name :
; time_gps_to_spec
;
; Purpose :
; change the gps_time format in order to use 'spectro_plot'
;
; History :
; Xingyao Chen, June 2017.
;
  gps_time=gps_time1
  gps_time[*,3]=gps_time[*,3]-8
  s_gps=size(gps_time)
  gps_string=strarr(s_gps[1])
  spec_time=dblarr(s_gps[1])
  month=strarr(1)
  for i=0,s_gps[1]-1 do begin
    mon=gps_time[i,1]
    case mon of
      '1':month='Jan'
      '2':month='Feb'
      '3':month='Mar'
      '4':month='Apr'
      '5':month='May'
      '6':month='Jun'
      '7':month='Jul'
      '8':month='Aug'
      '9':month='Sep'
      '10':month='Oct'
      '11':month='Nov'
      '12':month='Dec'
    endcase
    
    if gps_time[i,2] lt 10 then begin
      day='0'+strmid(gps_time[i,2],7,1)
    endif else begin
      day=strmid(gps_time[i,2],6,2)
    endelse
    
    if gps_time[i,3] lt 10 then begin
      hour='0'+strmid(gps_time[i,3],7,1)
    endif else begin
      hour=strmid(gps_time[i,3],6,2)
    endelse
    
    if gps_time[i,4] lt 10 then begin
      min='0'+strmid(gps_time[i,4],7,1)
    endif else begin
      min=strmid(gps_time[i,4],6,2)
    endelse
    
    if gps_time[i,5] lt 10 then begin
      sec='0'+strmid(gps_time[i,5],7,1)
    endif else begin
      sec=strmid(gps_time[i,5],6,2)
    endelse
    
    if (gps_time[i,6] ge 100) then msec=strmid(gps_time[i,6],5,3)
    if (gps_time[i,6] lt 100) and (gps_time[i,6] ge 10) then msec='0'+strmid(gps_time[i,6],6,2)
    if (gps_time[i,6] lt 10) then msec='00'+strmid(gps_time[i,6],7,1)
    
    gps_string[i]=day+'-'+month+'-'+strmid(gps_time[i,0],6,2)+' '$
      +hour+':'+min+':'+sec+'.'+msec 
    ;  return,month
;    gps_string[i]=strmid(gps_time[i,2],6,2)+'-'+month+'-'+strmid(gps_time[i,0],6,2)+' '$
;      +strmid(gps_time[i,3],6,2)+':'+strmid(gps_time[i,4],6,2)+':'+strmid(gps_time[i,5],6,2)$
;      +'.'+strmid(gps_time[i,6],5,3)
  endfor
  
  spec_time=anytim(gps_string)
  return,spec_time
END

;=========================================================================
PRO muser2fits, inputfile, outputfile, pixel, fovaia=fovaia, left=left, right=right, $
  lr=lr, band1=band1, band2=band2, band3=band3, band4=band4, residual=residual
  ; Name :
  ; muser2fits
  ;
  ; Purpose :
  ; convert clean map to fits file.
  ;
  ; Inputs :
  ; inputfile = inputfile & pixel
  ;
  ; Outputs :
  ; outputfile = outputfile
  ;
  ; Keywords :
  ; /fovaia        : the field of view same as AIA image
  ; /left & /right : for the left & right polarization
  ; /lr            : for the left + right polarization
  ; /band          : frequency band for 0.4-0.8 & 0.8-1.2 & 1.2-1.6 & 1.6-2.0 GHz
  ; /residual      : for residual image after clean
  ;
  ; Examples:
  ; muser2fits, inputfile, outputfile,512,/fovaia, /band4
  ;
  ; History :
  ; Chen X.Y, 04 July 2017.
  ;-------------------------begins-------------------------
  restore,inputfile
  if keyword_set( residual) then img_cl=img
  ;;img_cl = Array[512, 512]
  ;;beam_cl = Array[1024, 1024]
  ;;gps_time_cl
  ;;ch
  freqch=1725
  if keyword_set( band1) then freqch=(400+ch*25+25)
  if keyword_set( band2) then freqch=(800+ch*25+25)
  if keyword_set( band3) then freqch=(1200+ch*25+25)
  if keyword_set( band4) then freqch=(1600+ch*25+25)
  freqchstr=strcompress(freqch,/remove_all)
  
  time_str=time_gps_to_string(gps_time_cl)
  time_obs=time_gps_to_obs(gps_time_cl)
  
  if keyword_set( fovaia) then begin
    pixset=pixel
    img1=congrid(img_cl,pixset,pixset)
    centx=pixset/2
    centy=pixset/2
    fov_aia=0.60*4096/60. ;40.96 arcmin
    pixsz=round(pixset*fov_aia/60./2) ;half pixel size of the cutting image
    imgcut=img1[(centx-pixsz):(centx+pixsz),(centy-pixsz):(centy+pixsz)]
    
    imgcut=congrid(imgcut,pixset,pixset)
    fov=strcompress(fov_aia,/remove_all)+' arcmin'
    cdelt1=0.600000023842*4096/pixel
    cdelt2=0.600000023842*4096/pixel
    
    filenamefits=outputfile+'MUSER'+time_str[0]+'_'+freqchstr+'.fits'
    if keyword_set( left) then filenamefits=outputfile+'MUSER_L_'+time_str[0]+'_'+freqchstr+'.fits'
    if keyword_set( right) then filenamefits=outputfile+'MUSER_R_'+time_str[0]+'_'+freqchstr+'.fits'
    if keyword_set( lr) then filenamefits=outputfile+'MUSER_LR_'+time_str[0]+'_'+freqchstr+'.fits'
    
  endif else begin
    pixset=pixel
    imgcut=congrid(img_cl,pixset,pixset)
    fov='60 arcmin'
    cdelt1=60.*60./pixset
    cdelt2=60.*60./pixset
    
    type='raw'
    if keyword_set( residual) then type='rsd'
    filenamefits=outputfile+'MUSER_'+type+'_'+time_str[0]+'_'+freqchstr+'.fits'
    if keyword_set( left) then filenamefits=outputfile+'MUSER_'+type+'_L_'+time_str[0]+'_'+freqchstr+'.fits'
    if keyword_set( right) then filenamefits=outputfile+'MUSER_'+type+'_R_'+time_str[0]+'_'+freqchstr+'.fits'
    if keyword_set( lr) then filenamefits=outputfile+'MUSER_'+type+'_LR_'+time_str[0]+'_'+freqchstr+'.fits'
  endelse
  
  writefits,filenamefits,imgcut,header1
  header=header1
  
  sxaddpar, header, 'TELESCOPE','MUSER','Name of Telescope',BEFORE='SIMPLE'
  sxaddpar, header, 'INSTRUME','MUSER-I',BEFORE='SIMPLE'
  sxaddpar, header, 'FREQS',freqchstr+' MHz',BEFORE='SIMPLE'
  sxaddpar, header, 'INT_TIME','1 s',BEFORE='SIMPLE'
  sxaddpar, header, 'T_OBS',time_obs[0],'Start time of observation',BEFORE='SIMPLE'
  sxaddpar, header, 'T_STR',time_str[0],'Start time of observation',BEFORE='SIMPLE'
  sxaddpar, header, 'FOV',string(fov),BEFORE='SIMPLE'
  sxaddpar, header, 'CDELT1',string(cdelt1),BEFORE='SIMPLE'
  sxaddpar, header, 'CDELT2',string(cdelt2),BEFORE='SIMPLE'
  sxaddpar, header, 'ROTATE','no',BEFORE='SIMPLE'
  sxaddpar, header, 'CLEAN','yes',BEFORE='SIMPLE'
  sxaddpar, header, 'CLE_MTH','hogbom clean',BEFORE='SIMPLE'
  sxaddpar, header, 'RSUN_OBS','969.033367000',BEFORE='SIMPLE'
  sxaddpar, header, 'DATE',string(time_obs[0])
  if keyword_set( left) then sxaddpar, header, 'POLARIZATION','left',BEFORE='SIMPLE'
  if keyword_set( right) then sxaddpar, header, 'POLARIZATION','right',BEFORE='SIMPLE'
  if keyword_set( lr) then sxaddpar, header, 'POLARIZATION','left+right',BEFORE='SIMPLE'
  sxdelpar, header, 'SIMPLE'
  writefits,filenamefits,imgcut,header
  
  print,'fits file has written to:'+filenamefits
  
END

;=========================================================================

PRO D2B,X,B
  B=INTARR(8)
  B[7]=X MOD 2
  B[6]=FIX(X/2) MOD 2
  B[5]=FIX(X/4) MOD 2
  B[4]=FIX(X/8) MOD 2
  B[3]=FIX(X/16) MOD 2
  B[2]=FIX(X/32) MOD 2
  B[1]=FIX(X/64) MOD 2
  B[0]=FIX(X/128) MOD 2
  
END

;========================================================================
PRO READTIME,GPST,GPSTIME

  TIME=GPST
  D2B,TIME[0],T0
  D2B,TIME[1],T1
  D2B,TIME[2],T2
  D2B,TIME[3],T3
  D2B,TIME[4],T4
  D2B,TIME[5],T5
  D2B,TIME[6],T6
  D2B,TIME[7],T7
  
  T=[T7,T6,T5,T4,T3,T2,T1,T0]
  ;T=REVERSE(T)
  Y=(T[0]*2048+T[1]*1024+T[2]*512+T[3]*256+T[4]*128+T[5]*64+T[6]*32+T[7]*16+T[8]*8+T[9]*4+T[10]*2+T[11]) ;12bits
  MN=(T[12]*8+T[13]*4+T[14]*2+T[15])                                                                     ;4bits
  D=(T[16]*16+T[17]*8+T[18]*4+T[19]*2+T[20])                                                             ;5bits
  H=(T[21]*16+T[22]*8+T[23]*4+T[24]*2+T[25])                                                             ;5bits
  M=(T[26]*32+T[27]*16+T[28]*8+T[29]*4+T[30]*2+T[31])                                                    ;6bits
  S=(T[32]*32+T[33]*16+T[34]*8+T[35]*4+T[36]*2+T[37])                                                    ;6bits
  MS=(T[38]*512+T[39]*256+T[40]*128+T[41]*64+T[42]*32+T[43]*16+T[44]*8+T[45]*4+T[46]*2+T[47])            ;10bits
  uS=(T[48]*512+T[49]*256+T[50]*128+T[51]*64+T[52]*32+T[53]*16+T[54]*8+T[55]*4+T[56]*2+T[57])            ;10bits
  nS=(T[58]*32+T[59]*16+T[60]*8+T[61]*4+T[62]*2+T[63])*20                                                ;6bits
  
  ;H=H-8.0
  ;IF H LT 0 THEN H=H+24.
  
  PRINT,'TIME',FIX(Y),'Y',FIX(MN),'M',FIX(D),'D',FIX(H),'H',FIX(M),'M',FIX(S),'S',FIX(MS),'MS',FIX(US),'uS',FIX(NS),'nS'
  GPSTIME=[FIX(Y),FIX(MN),FIX(D),FIX(H),FIX(M),FIX(S),FIX(MS),FIX(uS),FIX(nS)]
  
END
;=========================================================================
PRO signed,a
  if a ge 128l*256l*256l then begin
    a=a-256l*256l*256L
  endif
END

;=========================================================================

FUNCTION interpolate_bessel,a,b

  n_a=n_elements(a)
  a_1=dblarr(n_a)
  a_2=dblarr(n_a)
  a_3=dblarr(n_a)
  for i=1,n_a-1 do begin
    a_1[i]=a[i]-a[i-1]
  endfor
  for j=1,n_a-2 do begin
    a_2[j]=a[j+1]-a[j]
  endfor
  for k=1,n_a-2 do begin
    a_2[k]=a_1[k+1]-a_1[k]
  endfor
  for l=2,n_a-2 do begin
    a_3[l]=a_2[l+1]-a_2[l]
  endfor
  
  b1=fix(b)
  bn=b-b1
  
  b2=bn*(bn-1)/(2.*2.)
  b3=bn*(bn-1)(bn-0.5)/(2.*3.)
  
  c=double(a[b1]+bn*a_1[b1+1]+b2*(a_2[b1]+a_2[b1+1])+b3*a_2[b1+1])
  
  return,c
  
END
;=========================================================================

PRO Cal_uvw,day,bjt,u,v,w

ANT_POS=FLTARR(44,3)
BSL=FLTARR(44,44,3)
DXYZ=FLTARR(44,44,3)
W=FLTARR(40,40)
U=FLTARR(40,40)
V=FLTARR(40,40)

HEIGHT=1356.
LAT=42.+12.710/60
LAT=!PI*LAT/180.

LONT=115.+15.030/60.
LOC_TIME=LONT/15.

RA_SUN_TB=[[16,27,54.87],$
  [16,32,13.80],$
  [16,36,33.34],$
  [16,40,53.48],$;4
  
  [16,45,14.20],$
  [16,49,35.47],$
  [16,53,57.27],$
  [16,58,19.59],$
  [17,02,42.40],$;9
  
  [17,07,05.67],$
  [17,11,29.37],$
  [17,15,53.48],$
  [17,20,17.97],$
  [17,24,42.82],$;14
  
  [17,29,07.98],$
  [17,33,33.43],$
  [17,37,59.14],$;2014.12.17         ;[15,28,44.69]
  [17,42,25.07],$
  [17,46,51.20],$
  [17,51,17.47],$;20
  
  [17,55,43.87],$
  [18,00,25.07],$
  [18,04,36.84],$
  [18,09,03.34],$
  [18,13,29.79],$;25
  
  [18,17,56.15],$
  [18,22,22.39],$
  [18,26,48.47],$
  [18,31,14.36],$
  [18,35,40.02],$;30
  
  [18,40,05.43]]
  
DEC_SUN_TB= [[21,44,56.3],$
  [21,54,11.2],$
  [22,03,00.8],$
  [22,11,24.9],$;4
  
  [22,19,23.3],$
  [22,26,55.7],$
  [22,34,01.9],$
  [22,40,41.6],$
  [22,46,54.7],$;9
  
  [22,52,40.9],$
  [22,58,00.1],$
  [23,02,52.0],$
  [23,07,16.6],$
  [23,11,13.6],$;14
  
  [23,14,42.9],$
  [23,17,44.4],$
  [23,20,17.9],$;2014.12.17      [18,54,11.8]
  [23,22,23.5],$
  [23,24,01.0],$
  [23,25,10.4],$;20
  
  [23,25,51.5],$
  [23,26,04.5],$
  [23,25,49.3],$
  [23,25,05.8],$
  [23,23,54.1],$;25
  
  [23,22,14.3],$
  [23,20,06.3],$
  [23,17,30.1],$
  [23,14,26.0],$
  [23,10,53.9],$;30
  
  [23,06,54.0]]

;RA_SUN_TB=[[14,24,25.94],$
;           [14,28,20.91],$
;           [14,32,16.95],$
;           [14,36,13.19],$;4
;
;           [14,40,10.54],$
;           [14,44,08.71],$
;           [14,48,07.70],$
;           [14,52,07.54],$
;           [14,56,14.03],$;9
;
;           [15,00,09.77],$
;           [15,04,12.17],$
;           [15,08,15.43],$
;           [15,12,19.55],$
;           [15,16,24.55],$;14
;
;           [15,20,30.40],$
;           [17,33,33.43],$
;           [17,37,59.14],$;2014.12.17         ;[15,28,44.69]
;           [17,42,25.07],$
;           [15,37,02.39],$
;           [15,41,12.50],$;20
;                      
;           [15,45,23.44],$
;           [15,49,35.20],$
;           [15,53,47.76],$
;           [15,58,01.11],$
;           [16,02,15.22],$;25
;
;           [16,06,30.08],$
;           [16,10,45.66],$
;           [16,15,01.96],$
;           [16,19,18.93],$
;           [16,23,36.58],$;30
;           
;           [16,27,54.87]]
;
;DEC_SUN_TB= [[14,19,45.7],$
;             [14,38,55.5],$
;             [14,57,51.1],$
;             [15,16,32.0],$;4
;
;             [15,34,57.8],$
;             [15,53,08.1],$
;             [16,11,02.7],$
;             [16,28,41.0],$
;             [16,46,02.8],$;9
;
;             [17,03,07.5],$
;             [17,19,54.9],$
;             [17,36,24.5],$
;             [17,52,35.9],$
;             [18,08,28.7],$;14
;
;             [18,24,02.6],$
;             [23,17,44.4],$
;             [23,20,17.9],$;2014.12.17      [18,54,11.8]
;             [23,22,23.5],$
;             [19,23,00.4],$
;             [19,36,53.6],$;20
;             
;             [19,50,25.5],$
;             [20,03,50.0],$
;             [20,16,24.0],$
;             [20,28,49.9],$
;             [20,40,53.1],$;25
;
;             [20,52,33.2],$
;             [21,03,50.0],$
;             [21,14,43.0],$
;             [21,25,11.9],$
;             [21,35,16.4],$;30
;             
;             [21,44,56.3]]
;


DEC_SUN_TB=-DEC_SUN_TB

RA_SUN=DBLARR(N_ELEMENTS(RA_SUN_TB[0,*]))
DEC_SUN=DBLARR(N_ELEMENTS(RA_SUN_TB[0,*]))

RA_SUN[*]=RA_SUN_TB[0,*]+RA_SUN_TB[1,*]/(60.0D)+RA_SUN_TB[2,*]/(3600.0D)
DEC_SUN[*]=DEC_SUN_TB[0,*]+DEC_SUN_TB[1,*]/(60.0D)+DEC_SUN_TB[2,*]/(3600.0D)

GMST1=DBLARR(31)

GMST_TB=[[4,39,05.9194],$
  [4,43,02.4708],$
  [4,46,59.0261],$
  [4,50,55.5815],$;4
  
  [4,54,52.1369],$                  ; 5th
  [4,58,48.6922],$
  [5,02,45.2476],$
  [5,06,41.8030],$
  [5,10,38.3583],$;9
  
  [5,14,34.9137],$                  ; 10th
  [5,18,31.4691],$
  [5,22,28.0244],$
  [5,26,24.5798],$
  [5,30,21.1352],$;14
  
  [5,34,17.6905],$                  ; 15th
  [5,38,14.2456],$
  [5,42,10.8013],$;2014.12.17       ;[3,43,54.1402]
  [5,46,07.3566],$
  [5,50,03.9120],$
  [5,54,00.4674],$;20
  
  [5,57,57.0228],$                  ;21th
  [6,01,53.5781],$
  [6,05,50.1335],$
  [6,19,46.6889],$
  [6,13,43.2442],$;25
  
  [6,17,39.7996],$                 ;26th
  [6,21,36.3550],$
  [6,25,32.9103],$
  [6,29,29.4657],$
  [6,33,26.0211],$;30
  
  [6,37,22.5464]]

;GMST_TB=[[2,40,49.2543],$
;         [2,44,45.8097],$
;         [2,48,42.3651],$
;         [2,52,38.9204],$;4
;
;         [2,56,35.4758],$                  ; 5th
;         [3,00,32.0312],$
;         [3,04,28.5866],$
;         [3,08,25.1419],$
;         [3,12,21.6973],$;9
;
;         [3,16,18.2527],$                  ; 10th
;         [3,20,14.8080],$
;         [3,24,11.3634],$
;         [3,28,07.9188],$
;         [3,32,04.4741],$;14
;
;         [3,36,01.0295],$                  ; 15th
;         [3,39,57.5849],$
;         [5,42,10.8013],$;2014.12.17       ;[3,43,54.1402]
;         [3,47,50.6956],$
;         [3,51,47.2510],$
;         [3,55,43.8063],$;20
;
;         [3,59,40.3617],$                  ;21th
;         [4,03,36.9171],$
;         [4,07,33.4724],$
;         [4,11,30.0278],$
;         [4,15,26.5832],$;25
;
;         [4,19,23.1385],$                 ;26th
;         [4,23,19.6939],$
;         [4,27,16.2493],$
;         [4,31,12.8047],$
;         [4,35,09.3600],$;30
;         
;         [4,39,05.9154]]


GMST1[*]=GMST_TB[0,*]+GMST_TB[1,*]/(60.0D)+GMST_TB[2,*]/(3600.0D)

OPENR,LUN1,'ANT_POS.TXT',/GET_LUN
FOR I=0,39 DO BEGIN
READF,LUN1,X,Y,Z
ANT_POS[I,*]=[Y,X,Z]
ENDFOR
FREE_LUN,LUN1

FOR I=0,39 DO BEGIN
 FOR J=0,39 DO BEGIN
  BSL[I,J,*]=ANT_POS[I,*]-ANT_POS[J,*]
  DXYZ[I,J,0]=-SIN(LAT)*BSL[I,J,0]+COS(LAT)*BSL[I,J,2]
  DXYZ[I,J,1]=BSL[I,J,1]
  DXYZ[I,J,2]=COS(LAT)*BSL[I,J,0]+SIN(LAT)*BSL[I,J,2]
 ENDFOR
ENDFOR

       ;DAY=11
       HOUR_REAL=bjt-8.
       GMST=INTERPOLATE_BESSEL(GMST1,(DAY+(HOUR_REAL)/24.)-1)
       RA=INTERPOLATE_BESSEL(RA_SUN,(DAY+(HOUR_REAL+DOUBLE(67./3600.))/24.)-1)
       DEC=INTERPOLATE_BESSEL(DEC_SUN,(DAY+(HOUR_REAL+DOUBLE(67./3600.))/24.)-1)

       LOC_MST=(GMST+HOUR_REAL+LOC_TIME) MOD 24
       H=(LOC_MST-RA)*15D
       ;print,time_s,h

;       H=0.994141048299d*15d
;       DEC= 18.1132053319d
       print,bjt,h,ra,dec
       DEC=DEC/180.*!PI
       H=H/180.*!PI
       
       FOR I=0,39 DO BEGIN
        FOR J=I,39 DO BEGIN
         U[I,J]=DXYZ[I,J,0]*SIN(H)+DXYZ[I,J,1]*COS(H)
         V[I,J]=DXYZ[I,J,0]*(-SIN(DEC))*Cos(H)+DXYZ[I,J,1]*SIN(DEC)*Sin(H)+COS(DEC)*DXYZ[I,J,2]
         W[I,J]=DXYZ[I,J,0]*Cos(DEC)*Cos(H)+DXYZ[I,J,1]*(-Cos(DEC)*Sin(H))+Sin(DEC)*DXYZ[I,J,2]
        ENDFOR
       ENDFOR
;save,u,v,w,filename='uvw_'+strcompress(day,remove_all)+strcompress(hour_real,/remove_all)+'.sav'

END
;=========================================================================

PRO flag_baseline,keep_table,cor_data1,vis,u,v
; Name :
; flag_baseline
;
; Purpose :
; flag baseline.
;
; Inputs :
; keep_table = which baseline should be kept. array of [2,*]
; cor_data1 = [40,40,2]
; 
; Outputs :
; vis, u, v
;
; History :
; Wei Wang, flag_data
; Xingyao Chen, 10 Aug 2017

  vis=complexarr(40,40)
  
  for i=0,39 do begin
    for j=I+1,39 do begin
      for tt=0, (size(keep_table))[2]-1 do begin
        if (i eq keep_table[0,tt]) and (j eq keep_table[1,tt]) then begin

          vis[i,j]=complex(cor_data1[i,j,0],cor_data1[i,j,1])
          u[j,i]=-u[i,j]
          v[j,i]=-v[i,j]
          vis[j,i]=conj(vis[i,j])
          
        endif
      endfor
    endfor
  endfor
  
  freq=1.7025  ;ghz
  lamda=0.3/freq
  
  u=u/lamda
  v=v/lamda
  
END
;=========================================================================

PRO FLAG_DATA,FLAG_TABLE,COR_DATA1,VIS,U,V

  VIS=COMPLEXARR(40,40)
  FOR I=0,39 DO BEGIN
    FOR J=I+1,39 DO BEGIN
      VIS[I,J]=COMPLEX(COR_DATA1[I,J,0],COR_DATA1[I,J,1])
      U[J,I]=-U[I,J]
      V[J,I]=-V[I,J]
      VIS[J,I]=CONJ(VIS[I,J])
    ENDFOR
  ENDFOR
  ;;
  FOR I=0,39 DO BEGIN
    FOR J=0,39 DO BEGIN
      IF (WHERE(I EQ FLAG_TABLE) NE -1) OR (WHERE(J EQ FLAG_TABLE) NE -1) THEN BEGIN
        VIS[I,J]=COMPLEX(0,0)
        U[I,J]=0
        V[I,J]=0
      ENDIF
    ENDFOR
  ENDFOR
  
  FREQ=1.7025  ;GHZ
  LAMDA=0.3/FREQ
  
  U=U/LAMDA
  V=V/LAMDA
  
END
;;=========================================================================
;;;==========After Gradding
PRO FFT_IMAGE,U,V,VIS,IMG

  PIX=512.
  nant=40
  Fov=1.0/180*!pi
  maxUV=Pix/Fov
  
  IMG=FLTARR(PIX,PIX)
  UV_IMG=COMPLEXARR(PIX,PIX)
  
  FOR I=0,Nant-1 DO BEGIN
    FOR J=0,Nant-1 DO BEGIN
      X=U[I,J]/MAXUV*PIX
      Y=V[I,J]/MAXUV*PIX
      IF Round(abs(X)+pix/2) GE PIX-3 THEN X=0
      IF Round(abs(Y)+pix/2) GE PIX-3 THEN Y=0
      
      option=3
      case(option) of
        0: begin
          UV_IMG[round(X+PIX/2),round(Y+PIX/2)]=VIS[I,J]+UV_IMG[round(X+PIX/2),round(Y+PIX/2)]
        end
        1: begin
          ;  ;==============================,case 1, Assign to adjacent 4 points.
          Sec0=(X+PIX/2) mod 1
          Sec1=(Y+PIX/2) mod 1
          Sec2=1-Sec0
          Sec3=1-Sec1
          L0=sqrt(Sec0^2+Sec1^2)
          L1=sqrt(Sec0^2+Sec3^2)
          L2=sqrt(Sec2^2+Sec3^2)
          L3=sqrt(Sec1^2+Sec2^2)
          TOTAL_L=TOTAL(L0+L1+L2+L3)
          L0=L0/TOTAL_L
          L1=L1/TOTAL_L
          L2=L2/TOTAL_L
          L3=L3/TOTAL_L
          UV_IMG[fix(X+PIX/2),fix(Y+PIX/2)]=VIS[I,J]*(1-L0)+UV_IMG[fix(X+PIX/2),fix(Y+PIX/2)]
          UV_IMG[fix(X+PIX/2),fix(Y+PIX/2)+1]=VIS[I,J]*(1-L1)+UV_IMG[fix(X+PIX/2),fix(Y+PIX/2)+1]
          UV_IMG[fix(X+PIX/2)+1,fix(Y+PIX/2)+1]=VIS[I,J]*(1-L2)+UV_IMG[fix(X+PIX/2)+1,fix(Y+PIX/2)+1]
          UV_IMG[fix(X+PIX/2)+1,fix(Y+PIX/2)]=VIS[I,J]*(1-L3)+UV_IMG[fix(X+PIX/2)+1,fix(Y+PIX/2)]
        end
        2: begin
          ;  ;==============================,case 2, Assign to adjacent 16 points.
          Sec0=(X+PIX/2) mod 1
          Sec1=(Y+PIX/2) mod 1
          Sec2=1-Sec0
          Sec3=1-Sec1
          L0=sqrt(Sec0^2+Sec1^2)
          L1=sqrt(Sec0^2+Sec3^2)
          L2=sqrt(Sec2^2+Sec3^2)
          L3=sqrt(Sec1^2+Sec2^2)
          L4=sqrt((1+Sec0)^2+(1+Sec1)^2)
          L5=sqrt((1+Sec0)^2+(Sec1)^2)
          L6=sqrt((1+Sec0)^2+(Sec3)^2)
          L7=sqrt((1+Sec0)^2+(1+Sec3)^2)
          L8=sqrt((Sec0)^2+(1+Sec3)^2)
          L9=sqrt((Sec2)^2+(1+Sec3)^2)
          L10=sqrt((1+Sec2)^2+(1+Sec3)^2)
          L11=sqrt((1+Sec2)^2+(Sec3)^2)
          L12=sqrt((1+Sec2)^2+(Sec1)^2)
          L13=sqrt((1+Sec2)^2+(1+Sec1)^2)
          L14=sqrt((Sec2)^2+(1+Sec1)^2)
          L15=sqrt((Sec0)^2+(1+Sec1)^2)
          
          Grid,L0,Lw0
          Grid,L1,Lw1
          Grid,L2,Lw2
          Grid,L3,Lw3
          Grid,L4,Lw4
          Grid,L5,Lw5
          Grid,L6,Lw6
          Grid,L7,Lw7
          Grid,L8,Lw8
          Grid,L9,Lw9
          Grid,L10,Lw10
          Grid,L11,Lw11
          Grid,L12,Lw12
          Grid,L13,Lw13
          Grid,L14,Lw14
          Grid,L15,Lw15
          UV_IMG[fix(X+PIX/2),fix(Y+PIX/2)]=VIS[I,J]*Lw0+UV_IMG[fix(X+PIX/2),fix(Y+PIX/2)]
          UV_IMG[fix(X+PIX/2),fix(Y+PIX/2)+1]=VIS[I,J]*Lw1+UV_IMG[fix(X+PIX/2),fix(Y+PIX/2)+1]
          UV_IMG[fix(X+PIX/2)+1,fix(Y+PIX/2)+1]=VIS[I,J]*Lw2+UV_IMG[fix(X+PIX/2)+1,fix(Y+PIX/2)+1]
          UV_IMG[fix(X+PIX/2)+1,fix(Y+PIX/2)]=VIS[I,J]*Lw3+UV_IMG[fix(X+PIX/2)+1,fix(Y+PIX/2)]
          UV_IMG[fix(X+PIX/2)-1,fix(Y+PIX/2)-1]=VIS[I,J]*Lw4+UV_IMG[fix(X+PIX/2)-1,fix(Y+PIX/2)-1]
          UV_IMG[fix(X+PIX/2)-1,fix(Y+PIX/2)]=VIS[I,J]*Lw5+UV_IMG[fix(X+PIX/2)-1,fix(Y+PIX/2)]
          UV_IMG[fix(X+PIX/2)-1,fix(Y+PIX/2)+1]=VIS[I,J]*Lw6+UV_IMG[fix(X+PIX/2)-1,fix(Y+PIX/2)+1]
          UV_IMG[fix(X+PIX/2)-1,fix(Y+PIX/2)+2]=VIS[I,J]*Lw7+UV_IMG[fix(X+PIX/2)-1,fix(Y+PIX/2)+2]
          UV_IMG[fix(X+PIX/2),fix(Y+PIX/2)+2]=VIS[I,J]*Lw8+UV_IMG[fix(X+PIX/2),fix(Y+PIX/2)+2]
          UV_IMG[fix(X+PIX/2)+1,fix(Y+PIX/2)+2]=VIS[I,J]*Lw9+UV_IMG[fix(X+PIX/2)+1,fix(Y+PIX/2)+2]
          UV_IMG[fix(X+PIX/2)+2,fix(Y+PIX/2)+2]=VIS[I,J]*Lw10+UV_IMG[fix(X+PIX/2)+2,fix(Y+PIX/2)+2]
          UV_IMG[fix(X+PIX/2)+2,fix(Y+PIX/2)+1]=VIS[I,J]*Lw11+UV_IMG[fix(X+PIX/2)+2,fix(Y+PIX/2)+1]
          UV_IMG[fix(X+PIX/2)+2,fix(Y+PIX/2)]=VIS[I,J]*Lw12+UV_IMG[fix(X+PIX/2)+2,fix(Y+PIX/2)]
          UV_IMG[fix(X+PIX/2)+2,fix(Y+PIX/2)-1]=VIS[I,J]*Lw13+UV_IMG[fix(X+PIX/2)+2,fix(Y+PIX/2)-1]
          UV_IMG[fix(X+PIX/2)+1,fix(Y+PIX/2)-1]=VIS[I,J]*Lw14+UV_IMG[fix(X+PIX/2)+1,fix(Y+PIX/2)-1]
          UV_IMG[fix(X+PIX/2),fix(Y+PIX/2)-1]=VIS[I,J]*Lw15+UV_IMG[fix(X+PIX/2),fix(Y+PIX/2)-1]
          
        END
        ;  ;==============================,case 2, Assign to adjacent 36 points.
        3: begin
          m=6
          Lw=fltarr(m,m)
          Sec0=(X+PIX/2) mod 1
          Sec1=(Y+PIX/2) mod 1
          for ii=0,m-1 do begin
            for jj=0,m-1 do begin
              Lin=sqrt((Sec0+2-ii)^2+(Sec1+2-jj)^2)
              Grid,Lin,Lout
              ;print,Lout
              Lw[ii,jj]=Lout
            endfor
          endfor
          UV_IMG[fix(X+PIX/2)-2:fix(X+PIX/2)+3,fix(Y+PIX/2)-2:fix(Y+PIX/2)+3]=VIS[I,J]*Lw+UV_IMG[fix(X+PIX/2)-2:fix(X+PIX/2)+3,fix(Y+PIX/2)-2:fix(Y+PIX/2)+3]
        end
        else: begin
        end
      endcase
      
    ENDFOR
  ENDFOR
  
  UV_IMG1=shift(UV_IMG,pix/2,pix/2)
  IMG=FFT(UV_IMG1,/center)
  IMG=ROTATE(IMG,5)
  
  taper_uv=fltarr(pix,pix)
  for i=0,pix-1 do begin
    for j=0,pix-1 do begin
      Lin=sqrt((i-pix/2)^2+(j-pix/2)^2)
      Grid,Lin,Lout
      taper_uv[i,j]=Lout
    endfor
  endfor
  taper_uv1=shift(taper_uv,pix/2,pix/2)
  taper=fft(taper_uv1,/center)
  img=img/real_part(taper)
END

;=========================================================================
PRO FFT_beam,U,V,BEAM
  pix=512
  nant=40
  BPIX=2*PIX
  BFov=2./180*!pi
  maxBUV=Bpix/BFov
  
  BEAM=FLTARR(BPIX,BPIX)
  UV_BEAM=FLTARR(BPIX,BPIX)
  
  FOR I=0,Nant-1 DO BEGIN
    FOR J=0,Nant-1 DO BEGIN
      BX=U[I,J]/MAXBUV*BPIX
      BY=V[I,J]/MAXBUV*BPIX
      IF round(abs(BX)+BPIX/2) GE BPIX-3 THEN BX=0
      IF round(abs(BY)+BPIX/2) GE BPIX-3 THEN BY=0
      
      option=3
      case(option) of
        0: begin
          UV_BEAM[round(BX+BPIX/2),round(BY+BPIX/2)]=1.+UV_BEAM[round(BX+BPIX/2),round(BY+BPIX/2)]
        end
        1: begin
          ;  ;==============================,case 1, Assign to adjacent 4 points.
          Sec0=(BX+BPIX/2) mod 1
          Sec1=(BY+BPIX/2) mod 1
          Sec2=1-Sec0
          Sec3=1-Sec1
          L0=sqrt(Sec0^2+Sec1^2)
          L1=sqrt(Sec0^2+Sec3^2)
          L2=sqrt(Sec2^2+Sec3^2)
          L3=sqrt(Sec1^2+Sec2^2)
          TOTAL_L=TOTAL(L0+L1+L2+L3)
          L0=L0/TOTAL_L
          L1=L1/TOTAL_L
          L2=L2/TOTAL_L
          L3=L3/TOTAL_L
          UV_BEAM[fix(BX+BPIX/2),fix(BY+BPIX/2)]=1.*(1-L0)+UV_BEAM[fix(BX+BPIX/2),fix(BY+BPIX/2)]
          UV_BEAM[fix(BX+BPIX/2),fix(BY+BPIX/2)+1]=1.*(1-L1)+UV_BEAM[fix(BX+BPIX/2),fix(BY+BPIX/2)+1]
          UV_BEAM[fix(BX+BPIX/2)+1,fix(BY+BPIX/2)+1]=1.*(1-L2)+UV_BEAM[fix(BX+BPIX/2)+1,fix(BY+BPIX/2)+1]
          UV_BEAM[fix(BX+BPIX/2)+1,fix(BY+BPIX/2)]=1.*(1-L3)+UV_BEAM[fix(BX+BPIX/2)+1,fix(BY+BPIX/2)]
        end
        2: begin
          ;  ;==============================,case 2, Assign to adjacent 16 points.
          Sec0=(BX+BPIX/2) mod 1
          Sec1=(BY+BPIX/2) mod 1
          Sec2=1-Sec0
          Sec3=1-Sec1
          L0=sqrt(Sec0^2+Sec1^2)
          L1=sqrt(Sec0^2+Sec3^2)
          L2=sqrt(Sec2^2+Sec3^2)
          L3=sqrt(Sec1^2+Sec2^2)
          L4=sqrt((1+Sec0)^2+(1+Sec1)^2)
          L5=sqrt((1+Sec0)^2+(Sec1)^2)
          L6=sqrt((1+Sec0)^2+(Sec3)^2)
          L7=sqrt((1+Sec0)^2+(1+Sec3)^2)
          L8=sqrt((Sec0)^2+(1+Sec3)^2)
          L9=sqrt((Sec2)^2+(1+Sec3)^2)
          L10=sqrt((1+Sec2)^2+(1+Sec3)^2)
          L11=sqrt((1+Sec2)^2+(Sec3)^2)
          L12=sqrt((1+Sec2)^2+(Sec1)^2)
          L13=sqrt((1+Sec2)^2+(1+Sec1)^2)
          L14=sqrt((Sec2)^2+(1+Sec1)^2)
          L15=sqrt((Sec0)^2+(1+Sec1)^2)
          
          Grid,L0,Lw0
          Grid,L1,Lw1
          Grid,L2,Lw2
          Grid,L3,Lw3
          Grid,L4,Lw4
          Grid,L5,Lw5
          Grid,L6,Lw6
          Grid,L7,Lw7
          Grid,L8,Lw8
          Grid,L9,Lw9
          Grid,L10,Lw10
          Grid,L11,Lw11
          Grid,L12,Lw12
          Grid,L13,Lw13
          Grid,L14,Lw14
          Grid,L15,Lw15
          
          UV_BEAM[fix(BX+BPIX/2),fix(BY+BPIX/2)]=1.0*Lw0+UV_BEAM[fix(BX+BPIX/2),fix(BY+BPIX/2)]
          UV_BEAM[fix(BX+BPIX/2),fix(BY+BPIX/2)+1]=1.0*Lw1+UV_BEAM[fix(BX+BPIX/2),fix(BY+BPIX/2)+1]
          UV_BEAM[fix(BX+BPIX/2)+1,fix(BY+BPIX/2)+1]=1.0*Lw2+UV_BEAM[fix(BX+BPIX/2)+1,fix(BY+BPIX/2)+1]
          UV_BEAM[fix(BX+BPIX/2)+1,fix(BY+BPIX/2)]=1.0*Lw3+UV_BEAM[fix(BX+BPIX/2)+1,fix(BY+BPIX/2)]
          UV_BEAM[fix(BX+BPIX/2)-1,fix(BY+BPIX/2)-1]=1.0*Lw4+UV_BEAM[fix(BX+BPIX/2)-1,fix(BY+BPIX/2)-1]
          UV_BEAM[fix(BX+BPIX/2)-1,fix(BY+BPIX/2)]=1.0*Lw5+UV_BEAM[fix(BX+BPIX/2)-1,fix(BY+BPIX/2)]
          UV_BEAM[fix(BX+BPIX/2)-1,fix(BY+BPIX/2)+1]=1.0*Lw6+UV_BEAM[fix(BX+BPIX/2)-1,fix(BY+BPIX/2)+1]
          UV_BEAM[fix(BX+BPIX/2)-1,fix(BY+BPIX/2)+2]=1.0*Lw7+UV_BEAM[fix(BX+BPIX/2)-1,fix(BY+BPIX/2)+2]
          UV_BEAM[fix(BX+BPIX/2),fix(BY+BPIX/2)+2]=1.0*Lw8+UV_BEAM[fix(BX+BPIX/2),fix(BY+BPIX/2)+2]
          UV_BEAM[fix(BX+BPIX/2)+1,fix(BY+BPIX/2)+2]=1.0*Lw9+UV_BEAM[fix(BX+BPIX/2)+1,fix(BY+BPIX/2)+2]
          UV_BEAM[fix(BX+BPIX/2)+2,fix(BY+BPIX/2)+2]=1.0*Lw10+UV_BEAM[fix(BX+BPIX/2)+2,fix(BY+BPIX/2)+2]
          UV_BEAM[fix(BX+BPIX/2)+2,fix(BY+BPIX/2)+1]=1.0*Lw11+UV_BEAM[fix(BX+BPIX/2)+2,fix(BY+BPIX/2)+1]
          UV_BEAM[fix(BX+BPIX/2)+2,fix(BY+BPIX/2)]=1.0*Lw12+UV_BEAM[fix(BX+BPIX/2)+2,fix(BY+BPIX/2)]
          UV_BEAM[fix(BX+BPIX/2)+2,fix(BY+BPIX/2)-1]=1.0*Lw13+UV_BEAM[fix(BX+BPIX/2)+2,fix(BY+BPIX/2)-1]
          UV_BEAM[fix(BX+BPIX/2)+1,fix(BY+BPIX/2)-1]=1.0*Lw14+UV_BEAM[fix(BX+BPIX/2)+1,fix(BY+BPIX/2)-1]
          UV_BEAM[fix(BX+BPIX/2),fix(BY+BPIX/2)-1]=1.0*Lw15+UV_BEAM[fix(BX+BPIX/2),fix(BY+BPIX/2)-1]
        END
        ;  ;==============================,case 2, Assign to adjacent 36 points.
        3: begin
          m=6
          Lw=fltarr(m,m)
          Sec0=(BX+BPIX/2) mod 1
          Sec1=(BY+BPIX/2) mod 1
          for ii=0,m-1 do begin
            for jj=0,m-1 do begin
              Lin=sqrt((Sec0+2-ii)^2+(Sec1+2-jj)^2)
              Grid,Lin,Lout
              ;print,Lout
              Lw[ii,jj]=Lout
            endfor
          endfor
          UV_BEAM[fix(BX+BPIX/2)-2:fix(BX+BPIX/2)+3,fix(BY+BPIX/2)-2:fix(BY+BPIX/2)+3]=1.*Lw+UV_BEAM[fix(BX+BPIX/2)-2:fix(BX+BPIX/2)+3,fix(BY+BPIX/2)-2:fix(BY+BPIX/2)+3]
        end
        else: begin
        end
      endcase
      
    ENDFOR
  ENDFOR
  
  UV_BEAM1=shift(UV_BEAM,Bpix/2,Bpix/2)
  BEAM=FFT(UV_BEAM1,/CENTER)
  BEAM=ROTATE(BEAM,5)
  
  taper_uv=fltarr(bpix,bpix)
  for i=0d,bpix-1 do begin
    for j=0d,bpix-1 do begin
      Lin=sqrt((i-bpix/2d)^2+(j-bpix/2d)^2)
      Grid,Lin,Lout
      taper_uv[i,j]=Lout
    endfor
  endfor
  taper_uv1=shift(taper_uv,bpix/2,bpix/2)
  taper=fft(taper_uv1,/center)
  beam=beam/real_part(taper)
  
END

;;;==========No Gradding
;PRO FFT_IMAGE,U,V,VIS,BEAM,IMG
;
;  PIX=512
;  BPIX=PIX*2
;  MAXU=512.*180./!PI; MAXU MEAN 1 DEGREE * 1 DEGREE
;  MAXV=512.*180./!PI
;  
;  IMG=FLTARR(PIX,PIX)
;  BEAM=FLTARR(BPIX,BPIX)
;  UV_IMG=COMPLEXARR(PIX,PIX)
;  UV_BEAM=FLTARR(BPIX,BPIX)
;  
;  FOR I=0,39 DO BEGIN
;    FOR J=0,39 DO BEGIN
;      UV_IMG[U[I,J]/MAXU*PIX+PIX/2,V[I,J]/MAXV*PIX+PIX/2]=UV_IMG[U[I,J]/MAXU*PIX+PIX/2,V[I,J]/MAXV*PIX+PIX/2]+VIS[I,J]
;      UV_BEAM[U[I,J]/MAXU*BPIX+BPIX/2,V[I,J]/MAXV*BPIX+BPIX/2]=1.0
;    ENDFOR
;  ENDFOR
;  BEAM=FFT(UV_BEAM)
;  BEAM=SHIFT(BEAM,BPIX/2,BPIX/2)
;  BEAM=ROTATE(BEAM,5)
;  IMG=FFT(UV_IMG)
;  IMG=SHIFT(IMG,PIX/2,PIX/2)
;  IMG=ROTATE(IMG,5)
;END

;=========================================================================
FUNCTION sdo_aia_scale_hdr,image,index, xllp, yllp , nxp, nyp ,wavelnth=wavelnth, rs_ratio=rs_ratio, no_imgrscale = no_imgrscale, noscale = noscale

  if  n_params() lt 2 then  box_message,'IDL> sdo_aia_scale_hdr,image,index[,xll,yll,nx,ny],wavelnth=wavelnth,rs_ratio=rs_ratio'
  
  if ~keyword_set(rs_ratio) then rs_ratio = 1.0d
  
  r_sun = index[0].r_sun/rs_ratio
  crpix1 = (index[0].crpix1-1)/rs_ratio;2048.5000
  crpix2 = (index[0].crpix2-1)/rs_ratio;2048.5000
  sz_image = size(image);2  4096  4096  2 16777216
  
  if rs_ratio gt 1 and ~keyword_set(no_imgrscale) then image = congrid(image,sz_image[1]/rs_ratio,sz_image[2]/rs_ratio,/interp)
  sz_image = size(image)
  if exist(xllp) then crpix1 = (index[0].crpix1-1-xllp)/rs_ratio
  if exist(yllp) then crpix2 = (index[0].crpix2-1-yllp)/rs_ratio
  
  
  if sz_image[0] eq 2 then begin
    indx = rebin(indgen(sz_image[1]),sz_image[1],sz_image[2]);4096 4096
    indy = rebin(indgen(1,sz_image[2]),sz_image[1],sz_image[2]);4096 4096
    rdist = sqrt((indx-crpix1)^2+(indy-crpix2)^2);create a sun model
    rfilter = (rdist > r_sun)/r_sun-1
    ind_disk = where(rdist le r_sun)
    ind_limb = where(rdist gt r_sun)
  endif else begin
    print,'data must be a two dimensions array'
    return,-1
  endelse
  
  dummy = bytarr(sz_image[1:2])
  
  
  case wavelnth of
    '1600':  begin
      dummy[ind_limb] = bytscl((image[ind_limb]), max = 1000)
      dummy[ind_disk] = bytscl((image[ind_disk]), max = 1000)
      return,dummy
    end
    '1700':  begin
      dummy[ind_limb] = bytscl((image[ind_limb]), max = 2500)
      dummy[ind_disk] = bytscl((image[ind_disk]), max = 2500)
      return,dummy
    end
    '4500':  begin
      dummy[ind_limb] = bytscl((image[ind_limb]<32000.))
      dummy[ind_disk] = bytscl((image[ind_disk]<32000.))
      return,dummy
    end
    '94':  begin
      ; image = image*(sqrt(rfilter*10.0)+1)
      image = image*exp(rfilter*4)
      ; stop
      if keyword_set(noscale) then return,image
      dummy = bytscl(alog10((image>0.5<30)))
      ; dummy = bytscl(sqrt((image>0.1<10)))
      return,dummy
    end
    '131':  begin
      image = image*(sqrt(rfilter*5)+1)
      if keyword_set(noscale) then return,image
      dummy = bytscl(alog10((image>2<200)))
      return,dummy
    end
    '171':  begin
      image = image*exp(rfilter*5)
      if keyword_set(noscale) then return,image
      dummy = bytscl(alog10((image>30<2000)))
      ; image = image*exp(rfilter*5)
      ; dummy = bytscl(sqrt((image>5<1200)))
      return,dummy
    end
    '193':  begin
      image = image*exp(rfilter*3)
      if keyword_set(noscale) then return,image
      dummy = bytscl(alog10((image>30<5000)))
      ; dummy[ind_limb] = bytscl(alog10((image[ind_limb]>30<8000)))
      ; dummy[ind_disk] = bytscl(alog10((image[ind_disk]>40<8000)))
      return,dummy
    end
    '211':  begin
      dummy[ind_limb] = bytscl(alog10((image[ind_limb]>4<2000)))
      dummy[ind_disk] = bytscl(alog10((image[ind_disk]>6<2600)))
      return,dummy
    end
    '304':  begin
      ; stop
      ; image = image*((tanh(rfilter/max(rfilter)*1.*30.-2)+1)*3+1)
      ; image = image*((rfilter)^(0.1)*15+1)
      ; dummy[ind_limb] = bytscl(alog10((image[ind_limb]>1.0<200)))
      ; dummy[ind_disk] = bytscl(alog10((image[ind_disk]>1.0<70)))
      ; image = image*(sqrt(rfilter*100)+1)
      image = image*exp(rfilter*5)
      if keyword_set(noscale) then return,image
      dummy = bytscl(alog10((image>0.8<100)))
      return,dummy
    end
    '335':  begin
      dummy[ind_limb] = bytscl(alog10((image[ind_limb]>0.5<80)))
      dummy[ind_disk] = bytscl(alog10((image[ind_disk]>0.5<150)))
      return,dummy
    end
    '6173':  begin
      dummy[ind_limb] = bytscl(0>(image[ind_limb]<65535.))
      dummy[ind_disk] = bytscl(0>(image[ind_disk]<65535.))
      return,dummy
    end
    '1':  begin
      dummy[ind_limb] = bytscl(-128>image[ind_limb]<128.)
      dummy[ind_disk] = bytscl(-128>image[ind_disk]<128.)
      return,dummy
    end
  endcase
  
  ; '1600':return,bytscl((image), max = 1000)
  ; '1700':return,bytscl((image), max = 2500)
  ; '4500':return,bytscl((image<32000.))
  ; '94':  return,bytscl(sqrt((image>0.2<8)))
  ; '131': return,bytscl(alog10((image>1<200)))
  ; '171': return,bytscl(sqrt((image>2<1000)))
  ; ; '171': return,bytscl(sqrt((image>20<1000)))
  ; '193': return,bytscl(alog10((image>40<1500)))
  ; ; '193': return,bytscl(alog10((image>90<2000)))
  ; '211': return,bytscl(alog10((image>6<2600)))
  ; ; '304': return,bytscl((image>3.5<600)^0.6)
  ; '304': return,bytscl(alog10((image>0.5<30)))
  ; '335': return,bytscl(alog10((image>0.5<1500)))
  ; '6173':return,bytscl(0>(image<65535.))
  ; '1':   return,bytscl(-128>image<128.)
  
  return,-1
END

;=========================================================================
PRO aia_data_prep,filenamea,pix,rsun_radio,imgaia
  ; Name :
  ; aia_data_prep
  ;
  ; Purpose :
  ; the data prep for AIA, transform from original size into the same size with dirty image of MUSER.
  ;
  ; Inputs :
  ; filenamea = the fits file path of the AIA data.
  ; pix = the pixel size of MUSER image, 512 for MUSER-I and 1024 for MUSER-II normally.
  ; rsun_radio = 1 or larger, the value set to be 0 beyond rsun_radio*R_SUN.
  ;
  ; Outputs :
  ; imgaia = transform into the [pix,pix] image.
  ;
  ; History :
  ; Xingyao Chen, June 2017.
  ;
  
  read_sdo, filenamea, indexa, dataa, /UNCOMP_DELETE,/noshell,/use_shared_lib,/HIDE
  dataa=dataa/indexa.EXPTIME
  ;dataa=congrid(dataa,512,512,/interp)
  
  xysize=indexa.NAXIS1
  solar_center=[indexa.NAXIS1/2,indexa.NAXIS2/2]
  dataa_temp=make_array(xysize,xysize)
  for i=0,xysize-1 do begin
    for j=0,xysize-1 do begin
    
      r=sqrt((i-solar_center[0])^2+(j-solar_center[1])^2);pixel
      ;r=sqrt((i-xsize/2)^2+(j-ysize/2)^2);pixel
      if r le ((indexa.R_SUN)*rsun_radio) then dataa_temp[i,j]=dataa[i,j]
      
    endfor
  endfor
  dataa=dataa_temp
  
  pixaia=round(indexa.NAXIS1*indexa.CDELT1/3600.*PIX)
  
  wl=indexa.WAVELNTH
  
  wavelength = ['94','171','193','211','304','335','131','1600']
  mag_scal1 = [1.5,1.5,3,8.,3.,2.,1.5,5.]
  wl_index=where(wavelength eq wl)
  
  imageaia = sdo_aia_scale_hdr(dataa/(mag_scal1[wl_index])[0],indexa,wavelnth=wl,rs_ratio=1.,/no_imgrscale)
  imageaia=congrid(imageaia,pixaia,pixaia,/interp)
  
  edge=(pix-pixaia)/2
  
  imgaia=make_array(pix,pix)
  for xx=edge,edge+pixaia-1 do begin
    for yy=edge,edge+pixaia-1 do begin
      imgaia[xx,yy]=imageaia[xx-edge,yy-edge]
    endfor
  endfor
  
END


;=========================================================================
