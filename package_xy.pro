;.compile package_xy.pro
;
; Name:
; package_xy
;
; History :
; Xingyao Chen, June 2017
;
;;-------------------------name of packages-------------------------
;
;
;
;
;
;
;
;FUNCTION mv_img, input_img, xpixel, ypixel
;FUNCTION slice_pixel_index_smp, xs1, xs2, ys1, ys2
;FUNCTION sl_region_index, p1,p2,p3,p4
;FUNCTION slice_pixel_index, xs1, xs2, ys1, ys2
;FUNCTION aia_slice_pixel, xs1, xs2, ys1, ys2, index, data
;FUNCTION corona_model_Newkirk_density, rx
;FUNCTION corona_model_Allen_density, rx
;FUNCTION corona_model_Mann_density, rx1
;FUNCTION drift_rate, fre, A
;FUNCTION myfunct, X, A
;FUNCTION corona_model_Leblanc, fpe
;FUNCTION corona_model_Newkirk, fpe
;FUNCTION corona_model_Saito, fpe
;FUNCTION corona_model_Mann, fpe
;=========================================================================


;=========================================================================


;=========================================================================


;=========================================================================


;=========================================================================



;=========================================================================



;=========================================================================
FUNCTION mv_img, input_img, xpixel, ypixel
; Name :
; mv_img
;
; Purpose :
; move the image
;
; Examples :
; output_img=mv_img(input_img,10,11)
;
; History :
; Xingyao Chen, 31 Oct 2017.
;-------------------------begins-------------------------
xpixel=xpixel
ypixel=ypixel

xs=(size(input_img))[1]
ys=(size(input_img))[2]
output_img=make_array(xs,ys)

if xpixel ge 0 then output_img[xpixel:(xs-1),*]=input_img[0:(xs-1-xpixel),*]  ;;to the right 
if xpixel lt 0 then output_img[0:(xs-1+xpixel),*]=input_img[-xpixel:(xs-1),*]  ;;to the left
temp_img=output_img

if ypixel ge 0 then output_img[*,ypixel:(ys-1)]=temp_img[*,0:(ys-1-ypixel)]  ;;up
if ypixel lt 0 then output_img[*,0:(ys-1+ypixel)]=temp_img[*,-ypixel:(ys-1)]  ;;down

return,output_img

END

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

;=========================================================================
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
  
END

;=========================================================================
FUNCTION slice_pixel_index, xs1, xs2, ys1, ys2
  ;from axis to pixel for straight lines, from y1 to y2 and x1 to x2
  xs1=float(xs1) & xs2=float(xs2) & ys1=float(ys1) & ys1=float(ys1)
  ;case: parallel to y-axis
  if xs2 ne xs1 then begin
    ks = (ys2-ys1)/(xs2-xs1)
  endif else begin
    num=ceil(ys2-ys1+1)
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
    num=ceil(xs2-xs1+1)
    point=make_array(2,num)
    xs=round(xs1)
    for ii=0,num-1 do begin
      point[0,ii]=xs+ii
      point[1,ii]=ys1
    endfor
  endif
  
  ;case: y plays a leading role
  if ks gt 1 then begin
    num=ceil(ys2-ys1+1)
    point=make_array(2,num)
    ys=round(ys1)
    for ii=0,num-1 do begin
      xs=xs1+ii/ks
      point[0,ii]=round(xs)
      point[1,ii]=ys+ii
    endfor
  endif
  ;case: x plays a leading role
  if ks le 1 and ks gt 0 then begin
    num=ceil(xs2-xs1+1)
    point=make_array(2,num)
    xs=round(xs1)
    for ii=0,num-1 do begin
      ys=ys1+ii*ks
      point[0,ii]=xs+ii
      point[1,ii]=round(ys)
    endfor
  endif
  ;case: y plays a leading role
  if ks lt -1 then begin
    num=ceil(-ys2+ys1)
    point=make_array(2,num)
    ys=round(ys1)
    for ii=0,num-1 do begin
      xs=xs1-ii/ks
      point[0,ii]=round(xs)
      point[1,ii]=ys-ii
    endfor
  endif
  ;case: x plays a leading role
  if ks ge -1 and ks lt 0 then begin
    num=ceil(xs2-xs1)
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

;=========================================================================
FUNCTION aia_slice_pixel, xs1, xs2, ys1, ys2, index, data
  ;from axis to pixel for straight lines, from y1 to y2 and x1 to x2
  pix_xs1=index.crpix1+xs1/index.cdelt1-1
  pix_xs2=index.crpix1+xs2/index.cdelt1-1
  pix_ys1=index.crpix1+ys1/index.cdelt1-1
  pix_ys2=index.crpix1+ys2/index.cdelt1-1
  
  ;case: parallel to y-axis
  if xs2 ne xs1 then begin
    ks = (ys2-ys1)/(xs2-xs1)
  endif else begin
    num=ceil(pix_ys2-pix_ys1)
    point=make_array(4,num)
    ys=round(pix_ys1)
    for ii=0,num-1 do begin
      point[0,ii]=xs1
      point[1,ii]=ys+ii
      point[2,ii]=data[xs1,(ys+ii)]
      point[3,ii]=(ii+1)*(index.rsun_ref/index.r_sun/1e6);Megamiles
    endfor
    GOTO,nafik
  endelse
  ;case: parallel to x-axis
  if ks eq 0 then begin
    num=ceil(pix_xs2-pix_xs1)
    point=make_array(4,num)
    xs=round(pix_xs1)
    for ii=0,num-1 do begin
      point[0,ii]=xs+ii
      point[1,ii]=ys1
      point[2,ii]=data[(xs+ii),ys1]
      point[3,ii]=(ii+1)*(index.rsun_ref/index.r_sun/1e6);Megamiles
    endfor
  endif
  
  ;case: y plays a leading role
  if ks gt 1 then begin
    num=ceil(pix_ys2-pix_ys1)
    ;data_sl=make_array(num)
    point=make_array(4,num)
    ys=round(pix_ys1)
    for ii=0,num-1 do begin
      xs=pix_xs1+ii/ks
      ;data_sl[ii]=data[round(xs),(ys+ii)]
      point[0,ii]=round(xs)
      point[1,ii]=ys+ii
      point[2,ii]=data[round(xs),(ys+ii)]
      point[3,ii]=(ii+1)*sqrt(1+1/(ks^2))*(index.rsun_ref/index.r_sun/1e6);Megamiles
    endfor
  endif
  ;case: x plays a leading role
  if ks le 1 and ks gt 0 then begin
    num=ceil(pix_xs2-pix_xs1)
    ;data_sl=make_array(num)
    point=make_array(4,num)
    xs=round(pix_xs1)
    for ii=0,num-1 do begin
      ys=pix_ys1+ii*ks
      ;data_sl[ii]=data[round(xs),(ys+ii)]
      point[0,ii]=xs+ii
      point[1,ii]=round(ys)
      point[2,ii]=data[(xs+ii),round(ys)]
      point[3,ii]=(ii+1)*sqrt(1+1/(ks^2))*(index.rsun_ref/index.r_sun/1e6);transform into Megamiles
    endfor
  endif
  ;case: y plays a leading role
  if ks lt -1 then begin
    num=ceil(-pix_ys2+pix_ys1)
    ;data_sl=make_array(num)
    point=make_array(4,num)
    ys=round(pix_ys1)
    for ii=0,num-1 do begin
      xs=pix_xs1-ii/ks
      ;data_sl[ii]=data[round(xs),(ys+ii)]
      point[0,ii]=round(xs)
      point[1,ii]=ys-ii
      point[2,ii]=data[round(xs),(ys-ii)]
      point[3,ii]=(ii+1)*sqrt(1+1/(ks^2))*(index.rsun_ref/index.r_sun/1e6);Megamiles
    endfor
  endif
  ;case: x plays a leading role
  if ks ge -1 and ks lt 0 then begin
    num=ceil(pix_xs2-pix_xs1)
    ;data_sl=make_array(num)
    point=make_array(4,num)
    xs=round(pix_xs1)
    for ii=0,num-1 do begin
      ys=pix_ys1+ii*ks
      ;data_sl[ii]=data[round(xs),(ys+ii)]
      point[0,ii]=xs+ii
      point[1,ii]=round(ys)
      point[2,ii]=data[(xs+ii),round(ys)]
      point[3,ii]=(ii+1)*sqrt(1+1/(ks^2))*(index.rsun_ref/index.r_sun/1e6);Megamiles
    endfor
  endif
  nafik:
  RETURN, point
  
END

;=========================================================================
;Use the corona density model of Newkirk,1961, input r, output n(r)
FUNCTION corona_model_Newkirk_density, rx

  ;n(r)=n0*10^4.32/r
  miu=0.6
  G=6.67408e-8;cm3/g/s2
  Ms=1.98855e33;g solar mass
  kB=1.3806503e-16;-16 erg/K;-23 J/K
  T=1e6;K
  Rs=6.958e10;cm
  ns=5.14e9;cm-3
  e=3*1.602176462e-10;C in CGS *3*1e^9
  m=9.10938188e-28;g
  mp=1.67262158e-24;g
  
  n0=4.2e4
  
  nr=n0*(10^(4.32d/rx))
  
  RETURN,nr
  
END

;=========================================================================
;Use the corona density model of Allen,1963, input r, output n(r)
FUNCTION corona_model_Allen_density, rx
  ;n(r)=c1*r^-d1+c2*r^-d2
  miu=0.6
  G=6.67408e-8;cm3/g/s2
  Ms=1.98855e33;g solar mass
  kB=1.3806503e-16;-16 erg/K;-23 J/K
  T=1e6;K
  Rs=6.958e10;cm
  ns=5.14e9;cm-3
  e=3*1.602176462e-10;C in CGS *3*1e^9
  m=9.10938188e-28;g
  mp=1.67262158e-24;g
  
  ;fpe=wpe/2/!pi=sqrt(4*!pi*e^2*n(r)/m)/2/!pi
  c1=1.36e+6
  c2=1.68e+8
  d1=2.14
  d2=6.13
  
  nr=c1*rx^(-d1)+c2*rx^(-d2)
  
  RETURN,nr
  
END

;=========================================================================
;the corona density model of Mann,input r, output n(r)
FUNCTION corona_model_Mann_density, rx1
  ;n(r)=ns*exp(A/Rs*(Rs/r-1))
  ;fpe=wpe/2/!pi=sqrt(4*!pi*e^2*n(r)/m)/2/!pi
  miu=0.6
  G=6.67408e-8;cm3/g/s2
  Ms=1.98855e33;g solar mass
  kB=1.3806503e-16;-16 erg/K;-23 J/K
  T=1e6;K
  Rs=6.958e10;cm
  ns=5.14e9;cm-3
  ;r=Rs+0+x
  e=3*1.602176462e-10;C in CGS *3*1e^9
  m=9.10938188e-28;g
  mp=1.67262158e-24;g
  A=miu*mp*G*Ms/kB/T;*mp
  vc=sqrt(kB*T/miu/mp);vc=1.17292e+07
  rc=G*Ms/2/vc^2;rc=4.82351e+11;C=6.3*10^34
  ;print,(2*!pi*fpe)^2*m/(4*!pi*e^2);n(r)
  ;n(r)=ns*exp(A/Rs*(Rs/r-1))
  ;fpe=wpe/2/!pi=sqrt(4*!pi*e^2*n(r)/m)/2/!pi
  rx1=rx1*Rs
  nr=ns*exp(A/Rs*(Rs/rx1-1))
  
  RETURN, nr
  
END

;=========================================================================
FUNCTION drift_rate, fre, A
  ;       input Frequency and return time and coefficient
  ;
  ;       The function being fit is of the following form:
  ;      t = A[0] + fre^(1-A[2])/A[1]/(1-A[2])
  ;
  ;
  ;       dF(x)/dA(0) = 1
  ;       dF(x)/dA(1) = -X^(1-A[2])/(1-A[2])/A[1]^2
  ;       dF(x)/dA(2) = (X^(1-A[2])*(1-A[2])^(-2) - X^(-A[2]))/A[1]
  ;
  ;       dF(x)/dA(2) = (((1/(1-A[2])-alog(X))*X^(1-A[2]))/(1-A[2])/A[1]
  ;       dF(x)/dA(2) =(X^(1-A[2])*(1/(1-A[2])-alog(X)))/(1-A[2])/A[1]
  ;
  ;       return,[[tim],[dt/dA(0)],[dt/dA(1)],[dt/dA(2)]]
  
  tim = A[0] + fre^(1-A[2])/A[1]/(1-A[2])
  
  ;   RETURN, [[F], [replicate(1d,n_elements(X))], [-X^(1-A[2])/(1-A[2])/A[1]^2], [(X^(1-A[2])*(1-A[2])^(-2) - X^(-A[2]))/A[1]]]
  RETURN, [[tim], [replicate(1d,n_elements(fre))], [-fre^(1-A[2])/(1-A[2])/A[1]^2],[(fre^(1-A[2])*(1/(1-A[2])-alog(fre)))/(1-A[2])/A[1]]]
  
END

;=========================================================================
FUNCTION myfunct, X, A
  ; First, define a return function for LMFIT:
  bx = A[0]*EXP(A[1]*X)
  RETURN,[ [bx+A[2]+A[3]*SIN(X)], [EXP(A[1]*X)], [bx*X], $
    [1.0] ,[SIN(X)] ]
END
;f(x)=a[0] * exp(a[1]*x) + a[2] + a[3] * sin(x)

;=========================================================================
;Use the corona density model of Leblanc
FUNCTION corona_model_Leblanc, fpe
  ;n(r)=c1*r^-d1+c2*r^-d2+c3*r^-d3
  ;c1=3.3*1e5 & d1=2 & c2=4.1*1e6 & d2=4 & c3=8.0*1e7 & d3=6
  miu=0.6
  G=6.67408e-8;cm3/g/s2
  Ms=1.98855e33;g solar mass
  kB=1.3806503e-16;-16 erg/K;-23 J/K
  T=1e6;K
  Rs=6.958e10;cm
  ns=5.14e9;cm-3
  e=3*1.602176462e-10;C in CGS *3*1e^9
  m=9.10938188e-28;g
  mp=1.67262158e-24;g
  
  ;fpe=wpe/2/!pi=sqrt(4*!pi*e^2*n(r)/m)/2/!pi
  c1=3.3e+5
  c2=4.1e+6
  d1=2.
  d2=4.
  c3=8.0e7
  d3=6.
  
  nf=!pi*fpe^2*m/e^2
  
  r1=1.01
  r2=2.00
  r12=(r1+r2)/2.
  num=0
  nr12=c1*(r12^(-d1))+c2*(r12^(-d2))+c3*(r12^(-d3))
  dn=abs(nr12-nf)
  nr1=c1*r1^(-d1)+c2*r1^(-d2)+c3*(r1^(-d3))
  while dn gt 100 do begin
  
    nr1=c1*r1^(-d1)+c2*r1^(-d2)+c3*(r1^(-d3))
    nr2=c1*r2^(-d1)+c2*r2^(-d2)+c3*(r2^(-d3))
    nr12=c1*r12^(-d1)+c2*r12^(-d2)+c3*(r12^(-d3))
    if nr12 lt nf then begin
      r2=r12
    endif else begin
      r1=r12
    endelse
    r12=(r1+r2)/2.
    nr12=c1*r12^(-d1)+c2*r12^(-d2)+c3*(r12^(-d3))
    dn=abs(nr12-nf)
    ;print,r12
    num=num+1
    ;print,num
    ;print,nf,nr12
    if num ge 50 then goto,jump
  endwhile
  
  jump:
  r12=r12*Rs/1e8;Mm
  RETURN,r12
  
END

;=========================================================================
;Use the corona density model of Newkirk,1961
FUNCTION corona_model_Newkirk, fpe

  ;n(r)=n0*10^4.32/r
  miu=0.6
  G=6.67408e-8;cm3/g/s2
  Ms=1.98855e33;g solar mass
  kB=1.3806503e-16;-16 erg/K;-23 J/K
  T=1e6;K
  Rs=6.958e10;cm
  ns=5.14e9;cm-3
  e=3*1.602176462e-10;C in CGS *3*1e^9
  m=9.10938188e-28;g
  mp=1.67262158e-24;g
  
  nf=!pi*fpe^2*m/e^2
  
  n0=4.2e4
  r1=1.01
  r2=2.00
  r12=(r1+r2)/2.
  num=0
  nr12=n0*(10^(4.32d/r12))
  dn=abs(nr12-nf)
  nr1=n0*(10^(4.32d/r1))
  while dn gt 100 do begin
  
    nr1=n0*(10^(4.32d/r1))
    nr2=n0*(10^(4.32d/r2))
    nr12=n0*(10^(4.32d/r12))
    if nr12 lt nf then begin
      r2=r12
    endif else begin
      r1=r12
    endelse
    r12=(r1+r2)/2.
    nr12=n0*(10^(4.32d/r12))
    dn=abs(nr12-nf)
    ;print,r12
    num=num+1
    ;print,num
    ;print,nf,nr12
    if num ge 50 then goto,jump
  endwhile
  
  jump:
  r12=r12*Rs/1e8;Mm
  RETURN,r12
END

;=========================================================================
;Use the corona density model of Saito, 1977
FUNCTION corona_model_Saito, fpe
  ;n(r)=c1*r^-d1+c2*r^-d2
  miu=0.6
  G=6.67408e-8;cm3/g/s2
  Ms=1.98855e33;g solar mass
  kB=1.3806503e-16;-16 erg/K;-23 J/K
  T=1e6;K
  Rs=6.958e10;cm
  ns=5.14e9;cm-3
  e=3*1.602176462e-10;C in CGS *3*1e^9
  m=9.10938188e-28;g
  mp=1.67262158e-24;g
  
  ;fpe=wpe/2/!pi=sqrt(4*!pi*e^2*n(r)/m)/2/!pi
  c1=1.36e+6
  c2=1.68e+8
  d1=2.14
  d2=6.13
  
  nf=!pi*fpe^2*m/e^2
  
  r1=1.01
  r2=2.00
  r12=(r1+r2)/2.
  num=0
  nr12=c1*(r12^(-d1))+c2*(r12^(-d2))
  dn=abs(nr12-nf)
  nr1=c1*r1^(-d1)+c2*r1^(-d2)
  while dn gt 100 do begin
  
    nr1=c1*r1^(-d1)+c2*r1^(-d2)
    nr2=c1*r2^(-d1)+c2*r2^(-d2)
    nr12=c1*r12^(-d1)+c2*r12^(-d2)
    if nr12 lt nf then begin
      r2=r12
    endif else begin
      r1=r12
    endelse
    r12=(r1+r2)/2.
    nr12=c1*r12^(-d1)+c2*r12^(-d2)
    dn=abs(nr12-nf)
    ;print,r12
    num=num+1
    ;print,num
    ;print,nf,nr12
    if num ge 50 then goto,jump
  endwhile
  
  jump:
  r12=r12*Rs/1e8  ;Mm
  RETURN,r12
  
END

;=========================================================================
;the corona density model of Mann
FUNCTION corona_model_Mann, fpe
  ;n(r)=ns*exp(A/Rs*(Rs/r-1))
  ;fpe=wpe/2/!pi=sqrt(4*!pi*e^2*n(r)/m)/2/!pi  
  miu=0.6
  G=6.67408e-8;cm3/g/s2
  Ms=1.98855e33;g solar mass
  kB=1.3806503e-16;-16 erg/K;-23 J/K
  T=1e6;K  
  Rs=6.958e10;cm
  ns=5.14e9;cm-3
  ;r=Rs+0+x
  e=3*1.602176462e-10;C in CGS *3*1e^9
  m=9.10938188e-28;g
  mp=1.67262158e-24;g
  A=miu*mp*G*Ms/kB/T;*mp
  vc=sqrt(kB*T/miu/mp);vc=1.17292e+07
  rc=G*Ms/2/vc^2;rc=4.82351e+11;C=6.3*10^34
  ;print,(2*!pi*fpe)^2*m/(4*!pi*e^2);n(r)
  ;n(r)=ns*exp(A/Rs*(Rs/r-1))
  ;fpe=wpe/2/!pi=sqrt(4*!pi*e^2*n(r)/m)/2/!pi
 
  nrr=alog((2*!pi*fpe)^2*m/(4*!pi*e^2*ns))  
  x=1/(1+nrr*Rs/A)
  x=x*Rs/1e8;Mm
  RETURN, x
     
END 
;=========================================================================

