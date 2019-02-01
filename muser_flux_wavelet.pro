pro tv_img_x1y1, xss, yss, poza
  ; a
  pxs = 90./10.
  pys = 50./10.
  xmgL = 20./10.
  ymgL = 20./10.
  xmgR = 12./10.
  ymgU = 12./10.
  xn = 1. & yn = 1.
  xss = xn*pxs + xmgL + xmgR
  yss = yn*pys + ymgL + ymgU
  xcL = 1./xss
  ycL = 1./yss
  ymin1 = (xmgL + pys*0.)*ycL
  ymax1 = (xmgL + pys*1.)*ycL
  xmin1 = (xmgL + pxs*0.)*xcL
  xmax1 = (xmgL + pxs*1.)*xcL
  poza = [xmin1, ymin1, xmax1, ymax1]
end


pro wr_ch, str, x, y, xy = xy, ch_angle = ch_angle, ch_size = ch_size, ch_thick = ch_thick, ch_color = ch_color

ps_xmin = xy(0)*10. &  ps_xmax = xy(2)*10. &  ps_ymin = xy(1)*10. & ps_ymax = xy(3)*10.

if(n_elements(ch_size) eq 0.)  then ch_size = 1.
if(n_elements(ch_thick) eq 0.) then ch_thick = 1.
if(n_elements(ch_angle) eq 0.) then ch_angle = 0.
if(n_elements(ch_color) eq 0.) then ch_color = 0.

xrio = 1./(ps_xmax-ps_xmin)
yrio = 1./(ps_ymax-ps_ymin)
xpos = 0.0 + xrio*(x-ps_xmin)
ypos = 0.0 + yrio*(y-ps_ymin)

xyouts,  xpos, ypos, str, charsize = ch_size, charthick = ch_thick, orientation = ch_angle, color = ch_color, /normal

end


pro plt_line, x1, y1, x2, y2, xy = xy, ar_color = ar_color, line_s = line_s, l_thick = l_thick
  ps_xmin = xy(0)*10. &  ps_xmax = xy(2)*10. &  ps_ymin = xy(1)*10. & ps_ymax = xy(3)*10.
  if(n_elements(ar_size) eq 0.)  then ar_size = 1.
  if(n_elements(ar_thick) eq 0.) then ar_thick = 1.
  if(n_elements(ar_color) eq 0.) then ar_color = 0.
  xrio = 1./(ps_xmax-ps_xmin)
  yrio = 1./(ps_ymax-ps_ymin)
  xp1 = 0.0 + xrio*(x1-ps_xmin)
  yp1 = 0.0 + yrio*(y1-ps_ymin)
  xp2 = 0.0 + xrio*(x2-ps_xmin)
  yp2 = 0.0 + yrio*(y2-ps_ymin)
  plots, [xp1, xp2], [yp1, yp2], color = ar_color, /normal, line = line_s, thick = l_thick
end


;===================================================================function definition
FUNCTION morlet, $ ;*********************************************** MORLET
    k0,scale,k,period,coi,dofmin,Cdelta,psi0

    IF (k0 EQ -1) THEN k0 = 6d
    n = N_ELEMENTS(k)
    expnt = -(scale*k - k0)^2/2d*(k GT 0.)
    dt = 2*!PI/(n*k(1))
    norm = SQRT(2*!PI*scale/dt)*(!PI^(-0.25))   ; total energy=N   [Eqn(7)]
    morlet = norm*EXP(expnt > (-100d))
    morlet = morlet*(expnt GT -100)  ; avoid underflow errors
    morlet = morlet*(k GT 0.)  ; Heaviside step function (Morlet is complex)
    fourier_factor = (4*!PI)/(k0 + SQRT(2+k0^2)) ; Scale-->Fourier [Sec.3h]
    period = scale*fourier_factor
    coi = fourier_factor/SQRT(2)   ; Cone-of-influence [Sec.3g]
    dofmin = 2   ; Degrees of freedom with no smoothing
    Cdelta = -1
    IF (k0 EQ 6) THEN Cdelta = 0.776 ; reconstruction factor
    psi0 = !PI^(-0.25)
;   PRINT,scale,n,SQRT(TOTAL(ABS(morlet)^2,/DOUBLE))
    RETURN,morlet
END

FUNCTION paul, $ ;************************************************** PAUL
    m,scale,k,period,coi,dofmin,Cdelta,psi0

    IF (m EQ -1) THEN m = 4d
    n = N_ELEMENTS(k)
    expnt = -(scale*k)*(k GT 0.)
    dt = 2d*!PI/(n*k(1))
    norm = SQRT(2*!PI*scale/dt)*(2^m/SQRT(m*FACTORIAL(2*m-1)))
    paul = norm*((scale*k)^m)*EXP(expnt > (-100d))*(expnt GT -100)
    paul = paul*(k GT 0.)
    fourier_factor = 4*!PI/(2*m+1)
    period = scale*fourier_factor
    coi = fourier_factor*SQRT(2)
    dofmin = 2   ; Degrees of freedom with no smoothing
    Cdelta = -1
    IF (m EQ 4) THEN Cdelta = 1.132 ; reconstruction factor
    psi0 = 2.^m*FACTORIAL(m)/SQRT(!PI*FACTORIAL(2*m))
;   PRINT,scale,n,norm,SQRT(TOTAL(paul^2,/DOUBLE))*SQRT(n)
    RETURN,paul
END

FUNCTION dog, $ ;*************************************************** DOG
    m,scale,k,period,coi,dofmin,Cdelta,psi0

    IF (m EQ -1) THEN m = 2
    n = N_ELEMENTS(k)
    expnt = -(scale*k)^2/2d
    dt = 2d*!PI/(n*k(1))
    norm = SQRT(2*!PI*scale/dt)*SQRT(1d/GAMMA(m+0.5))
    I = DCOMPLEX(0,1)
    gauss = -norm*(I^m)*(scale*k)^m*EXP(expnt > (-100d))*(expnt GT -100)
    fourier_factor = 2*!PI*SQRT(2./(2*m+1))
    period = scale*fourier_factor
    coi = fourier_factor/SQRT(2)
    dofmin = 1   ; Degrees of freedom with no smoothing
    Cdelta = -1
    psi0 = -1
    IF (m EQ 2) THEN BEGIN
       Cdelta = 3.541 ; reconstruction factor
       psi0 = 0.867325
    ENDIF
    IF (m EQ 6) THEN BEGIN
       Cdelta = 1.966 ; reconstruction factor
       psi0 = 0.88406
    ENDIF
;   PRINT,scale,n,norm,SQRT(TOTAL(ABS(gauss)^2,/DOUBLE))*SQRT(n)
    RETURN,gauss
END


;****************************************************************** WAVELET
FUNCTION wavelet,y1,dt, $   ;*** required inputs
    S0=s0,DJ=dj,J=j, $   ;*** optional inputs
    PAD=pad,MOTHER=mother,PARAM=param, $
    VERBOSE=verbose,NO_WAVE=no_wave,RECON=recon, $
    LAG1=lag1,SIGLVL=siglvl,DOF=dof,GLOBAL=global, $   ;*** optional inputs
    SCALE=scale,PERIOD=period,YPAD=ypad, $  ;*** optional outputs
    DAUGHTER=daughter,COI=coi, $
    SIGNIF=signif,FFT_THEOR=fft_theor, $
    OCT=oct,VOICE=voice   ;*** defunct inputs

    ON_ERROR,2
    r = CHECK_MATH(0,1)
    n = N_ELEMENTS(y1)
    n1 = n
    base2 = FIX(ALOG(n)/ALOG(2) + 0.4999)   ; power of 2 nearest to N

;....check keywords & optional inputs
    IF (N_ELEMENTS(s0) LT 1) THEN s0 = 2.0*dt
    IF (N_ELEMENTS(voice) EQ 1) THEN dj = 1./voice
    IF (N_ELEMENTS(dj) LT 1) THEN dj = 1./8
    IF (N_ELEMENTS(oct) EQ 1) THEN J = FLOAT(oct)/dj
    IF (N_ELEMENTS(J) LT 1) THEN J=FIX((ALOG(FLOAT(n)*dt/s0)/ALOG(2))/dj)  ;[Eqn(10)]
    IF (N_ELEMENTS(mother) LT 1) THEN mother = 'MORLET'
    IF (N_ELEMENTS(param) LT 1) THEN param = -1
    IF (N_ELEMENTS(siglvl) LT 1) THEN siglvl = 0.95
    IF (N_ELEMENTS(lag1) LT 1) THEN lag1 = 0.0
    lag1 = lag1(0)
    verbose = KEYWORD_SET(verbose)
    do_daughter = KEYWORD_SET(daughter)
    do_wave = NOT KEYWORD_SET(no_wave)
    recon = KEYWORD_SET(recon)
    IF KEYWORD_SET(global) THEN MESSAGE, $
       'Please use WAVE_SIGNIF for global significance tests'

;....construct time series to analyze, pad if necessary
    ypad = y1 - TOTAL(y1)/n    ; remove mean
    IF KEYWORD_SET(pad) THEN BEGIN   ; pad with extra zeroes, up to power of 2
       ypad = [ypad,FLTARR(2L^(base2 + 1) - n)]
       n = N_ELEMENTS(ypad)
    ENDIF

;....construct SCALE array & empty PERIOD & WAVE arrays
    na = J + 1                  ; # of scales
    scale = DINDGEN(na)*dj      ; array of j-values
    scale = 2d0^(scale)*s0      ; array of scales  2^j   [Eqn(9)]
    period = FLTARR(na,/NOZERO) ; empty period array (filled in below)
    wave = COMPLEXARR(n,na,/NOZERO)  ; empty wavelet array
    IF (do_daughter) THEN daughter = wave   ; empty daughter array

;....construct wavenumber array used in transform [Eqn(5)]
    k = (DINDGEN(n/2) + 1)*(2*!PI)/(DOUBLE(n)*dt)
    k = [0d,k,-REVERSE(k(0:(n-1)/2 - 1))]

;....compute FFT of the (padded) time series
    yfft = FFT(ypad,-1,/DOUBLE)  ; [Eqn(3)]

    IF (verbose) THEN BEGIN  ;verbose
       PRINT
       PRINT,mother
       PRINT,'#points=',n1,'   s0=',s0,'   dj=',dj,'   J=',FIX(J)
       IF (n1 NE n) THEN PRINT,'(padded with ',n-n1,' zeroes)'
       PRINT,['j','scale','period','variance','mathflag'], $
         FORMAT='(/,A3,3A11,A10)'
    ENDIF  ;verbose
    IF (N_ELEMENTS(fft_theor) EQ n) THEN fft_theor_k = fft_theor ELSE $
       fft_theor_k = (1-lag1^2)/(1-2*lag1*COS(k*dt)+lag1^2)  ; [Eqn(16)]
    fft_theor = FLTARR(na)

;....loop thru each SCALE
    FOR a1=0,na-1 DO BEGIN  ;scale
       psi_fft=CALL_FUNCTION(mother, $
         param,scale(a1),k,period1,coi,dofmin,Cdelta,psi0)
       IF (do_wave) THEN $
         wave(*,a1) = FFT(yfft*psi_fft,1,/DOUBLE)  ;wavelet transform[Eqn(4)]
       period(a1) = period1   ; save period
       fft_theor(a1) = TOTAL((ABS(psi_fft)^2)*fft_theor_k)/n
       IF (do_daughter) THEN $
         daughter(*,a1) = FFT(psi_fft,1,/DOUBLE)   ; save daughter
       IF (verbose) THEN PRINT,a1,scale(a1),period(a1), $
          TOTAL(ABS(wave(*,a1))^2),CHECK_MATH(0), $
          FORMAT='(I3,3F11.3,I6)'
    ENDFOR  ;scale

    coi = coi*[FINDGEN((n1+1)/2),REVERSE(FINDGEN(n1/2))]*dt   ; COI [Sec.3g]

    IF (do_daughter) THEN $   ; shift so DAUGHTERs are in middle of array
       daughter = [daughter(n-n1/2:*,*),daughter(0:n1/2-1,*)]

;....significance levels [Sec.4]
    sdev = (MOMENT(y1))(1)
    fft_theor = sdev*fft_theor  ; include time-series variance
    dof = dofmin
    signif = fft_theor*CHISQR_CVF(1. - siglvl,dof)/dof   ; [Eqn(18)]

    IF (recon) THEN BEGIN  ; Reconstruction [Eqn(11)]
       IF (Cdelta EQ -1) THEN BEGIN
         y1 = -1
         MESSAGE,/INFO, $
          'Cdelta undefined, cannot reconstruct with this wavelet'
       ENDIF ELSE BEGIN
         y1=dj*SQRT(dt)/(Cdelta*psi0)*(FLOAT(wave) # (1./SQRT(scale)))
         y1 = y1[0:n1-1]
       ENDELSE
    ENDIF

    RETURN,wave(0:n1-1,*)    ; get rid of padding before returning

END


;------------AIA--------------------------------------------------------------

;;.r muser_flux_wavelet.pro
;;.compile package_xy.pro

filename_='./sav/flux_raw_20141217_043302_827.sav',/ver
restore,filename_,/ver

of1=of1
of2=of1
of4=of1
of9=of1
of17=of1
of35=of1

timerange=['17-Dec-14 04:26:00.000','17-Dec-14 04:44:00.000']
timerange=[atime(time[0]),atime(max(time))]
xtickname=[strmid(timerange[0],11,7),' ',' ','04:35:00',' ',' ',strmid(timerange[1],11,7)]

dt=1.

yrange=[10,800]

level=0.95
levels='95%'

smo=501;

tim=time

of1=of1-mean(of1[0:5])
of2=of2-mean(of2[0:5])
of4=of4-mean(of4[0:5])
of9=of9-mean(of9[0:5])
of17=of17-mean(of17[0:5])
of35=of35-mean(of35[0:5])

sf1=smooth(of1,smo)
sf2=smooth(of2,smo)
sf4=smooth(of4,smo)
sf9=smooth(of9,smo)
sf17=smooth(of17,smo)
sf35=smooth(of35,smo)

df1=of1-smooth(of1,smo)
df2=of2-smooth(of2,smo)
df4=of4-smooth(of4,smo)
df9=of9-smooth(of9,smo)
df17=of17-smooth(of17,smo)
df35=of35-smooth(of35,smo)

xy=0 & f1=of1 & f2=of2 & f4=of4 & f9=of9 & f17=of17 & f35=of35
xy=1 & f1=sf1 & f2=sf2 & f4=sf4 & f9=sf9 & f17=sf17 & f35=sf35
xy=2 & f1=df1 & f2=df2 & f4=df4 & f9=df9 & f17=df17 & f35=df35

;----------------------------------------para-----------------------------
time0=(str2arr(anytim(tim[0],/yohkoh,/sec,/time),'.'))[0]

xticks=6
xtitle=strmid(timerange[0],0,10)

hh=float(strmid(time0,0,2))
mm=float(strmid(time0,3,2))
ss=float(strmid(time0,6,2))

bhh=float(strmid(timerange(0),11,2))
bmm=float(strmid(timerange(0),14,2))

ehh=float(strmid(timerange(1),11,2))
emm=float(strmid(timerange(1),14,2))

bn=0L
en=0L

if bhh gt 15 then begin
bn=(bhh-hh-1)*60*60+(bmm+60-mm)*60-ss
en=(ehh-hh-1)*60*60+(emm+60-mm)*60-ss
endif

if bhh lt 15 then begin
bn=(bhh+24-hh-1)*60*60+(bmm+60-mm)*60-ss
en=(ehh+24-hh-1)*60*60+(emm+60-mm)*60-ss
endif

;----------------------------------------plot 01-----------------------------
set_plot,'PS'
tv_img_x1y1, xss, yss, poza
!P.FONT = 0

filenamep='./imagetest/muser_wavelet/'+strmid(filename_,38,26)

device, filename = filenamep+'_sm'+strcompress(smo,/remove_all)+'_flux'+'.eps',xsize=20,ysize=10,BITS_PER_PIXEL=8,$
  SET_FONT='Helvetica',/ENCAPSULATED,scale_factor=2,/color
!p.multi=[0,3,2]

if xy ne 2 then begin
  plot,f1,xstyle=1,ystyle=1
  plot,f2,xstyle=1,ystyle=1
  plot,f4,xstyle=1,ystyle=1
  plot,f9,xstyle=1,ystyle=1,xticks=6,xminor=3,xtickname=xtickname,xtitle=xtitle
  plot,f17,xstyle=1,ystyle=1
  plot,f35,xstyle=1,ystyle=1
endif

if xy eq 2 then begin
  plot,of1,xstyle=1,ystyle=1
  ;oplot,of1,color=cgcolor('red')
  oplot,sf1,color=cgcolor('blue')
  plot,of2,xstyle=1,ystyle=1
  oplot,sf2,color=cgcolor('blue')
  plot,of4,xstyle=1,ystyle=1
  oplot,sf4,color=cgcolor('blue')
  plot,of9,xstyle=1,ystyle=1,xticks=6,xminor=3,xtickname=xtickname,xtitle=xtitle
  oplot,sf9,color=cgcolor('blue')
  plot,of17,xstyle=1,ystyle=1
  oplot,sf17,color=cgcolor('blue')
  plot,of35,xstyle=1,ystyle=1
  oplot,sf35,color=cgcolor('blue')
endif

;----------------------------------------plot 02-----------------------------
device, filename = filenamep+'_sm'+strcompress(smo,/remove_all)+'_wave'+'.eps',xsize=20,ysize=10,BITS_PER_PIXEL=8,$
  SET_FONT='Helvetica',/ENCAPSULATED,scale_factor=2,/color


;----------------------------------------para-----------------------------
S=n_elements(f1)
Yr=findgen(S)

;111===============================================================
R=f1
data=R(0:S-1)
resultR=findgen(S)
resultR(*)=data(0:S-1)
ntime=S
dt=dt ; time resolution 1 year
time=findgen(ntime)*dt
waveR=wavelet(resultR, dt, PERIOD=periodR, COI=coi, /PAD, SIGNIF=signifR)
nscaleR=N_ELEMENTS(periodR)
help,waver
print,signifr
m=max(alog((abs(waveR)^2)))
print,periodR
maxt=max(alog(abs(waveR)^2))

loadct,39
contour, alog(abs(waveR)^2), Yr(*), periodR, position=[0.08, 0.54, 0.38, 0.98], xsty=1, ysty=1, xrange=[0, S-1], $
yrange=yrange, /ylog, nlevels=30, /fill,title=' ',$
ytitle='Period (s)',xticks=6,xminor=3,$ levels=(0.68+indgen(30)/100.)*maxt,c_colors=150+indgen(30)*4,$
xtickname=[replicate(' ',7)],charsize=1.5,charthick=2.5,thick=2.5,/normal
signifR=rebin(transpose(signifR),ntime,nscaleR)
contour,abs(waveR)^2/signifR, Yr(*), periodR, /overplot, level=level, C_ANNOT=levels, xrange=[0, S-1]
PLOTS,Yr(*), coi, noclip=0

;222===============================================================
R=f2
data=R(0:S-1)
resultR=findgen(S)
resultR(*)=data(0:S-1)
ntime=S
dt=dt ; time resolution 1 year
time=findgen(ntime)*dt
waveR=wavelet(resultR, dt, PERIOD=periodR, COI=coi, /PAD, SIGNIF=signifR)
nscaleR=N_ELEMENTS(periodR)
help,waver
print,signifr
m=max(alog((abs(waveR)^2)))
print,periodR
maxt=max(alog(abs(waveR)^2))

loadct,39
contour, alog(abs(waveR)^2), Yr(*), periodR, position=[0.38, 0.54, 0.68, 0.98], xsty=1, ysty=1, xrange=[0, S-1], $
yrange=yrange, /ylog, nlevels=30, /fill,title=' ',$
ytitle=' ',ytickname=[' ',' ',' '],xticks=6,xminor=3,$ levels=(0.68+indgen(30)/100.)*maxt,c_colors=150+indgen(30)*4,$
xtickname=[replicate(' ',7)],charsize=1.5,charthick=2.5,thick=2.5,/normal
signifR=rebin(transpose(signifR),ntime,nscaleR)
contour,abs(waveR)^2/signifR, Yr(*), periodR, /overplot, level=level, C_ANNOT=levels, xrange=[0, S-1]
PLOTS,Yr(*), coi, noclip=0

;333===============================================================
R=f4
data=R(0:S-1)
resultR=findgen(S)
resultR(*)=data(0:S-1)
ntime=S
dt=dt ; time resolution 1 year
time=findgen(ntime)*dt
waveR=wavelet(resultR, dt, PERIOD=periodR, COI=coi, /PAD, SIGNIF=signifR)
nscaleR=N_ELEMENTS(periodR)
help,waver
print,signifr
m=max(alog((abs(waveR)^2)))
print,periodR
maxt=max(alog(abs(waveR)^2))

loadct,39
contour, alog(abs(waveR)^2), Yr(*), periodR, position=[0.68, 0.54, 0.98, 0.98], xsty=1, ysty=1, xrange=[0, S-1], $
yrange=yrange, /ylog, nlevels=30, /fill,title=' ',$
ytitle=' ',ytickname=[' ',' ',' '],xticks=6,xminor=3,$ levels=(0.68+indgen(30)/100.)*maxt,c_colors=150+indgen(30)*4,$
xtickname=[replicate(' ',7)],charsize=1.5,charthick=2.5,thick=2.5,/normal
signifR=rebin(transpose(signifR),ntime,nscaleR)
contour,abs(waveR)^2/signifR, Yr(*), periodR, /overplot, level=level, C_ANNOT=levels, xrange=[0, S-1]
PLOTS,Yr(*), coi, noclip=0

;444===============================================================
R=f9
data=R(0:S-1)
resultR=findgen(S)
resultR(*)=data(0:S-1)
ntime=S
dt=dt ; time resolution 1 year
time=findgen(ntime)*dt
waveR=wavelet(resultR, dt, PERIOD=periodR, COI=coi, /PAD, SIGNIF=signifR)
nscaleR=N_ELEMENTS(periodR)
help,waver
print,signifr
m=max(alog((abs(waveR)^2)))
print,periodR
maxt=max(alog(abs(waveR)^2))

loadct,39
contour, alog(abs(waveR)^2), Yr(*), periodR, position=[[0.08, 0.1, 0.38, 0.54]], xsty=1, ysty=1, xrange=[0, S-1], $
  yrange=yrange, /ylog, nlevels=30, /fill,title=' ',$
  ytitle='Period (s)',xticks=6,xminor=3,$ levels=(0.68+indgen(30)/100.)*maxt,c_colors=150+indgen(30)*4,$
  xtickname=xtickname,charsize=1.5,charthick=2.5,thick=2.5,/normal,xtitle=xtitle
signifR=rebin(transpose(signifR),ntime,nscaleR)
contour,abs(waveR)^2/signifR, Yr(*), periodR, /overplot, level=level, C_ANNOT=levels, xrange=[0, S-1]
PLOTS,Yr(*), coi, noclip=0

;555===============================================================
R=f17
data=R(0:S-1)
resultR=findgen(S)
resultR(*)=data(0:S-1)
ntime=S
dt=dt ; time resolution 1 year
time=findgen(ntime)*dt
waveR=wavelet(resultR, dt, PERIOD=periodR, COI=coi, /PAD, SIGNIF=signifR)
nscaleR=N_ELEMENTS(periodR)
help,waver
print,signifr
m=max(alog((abs(waveR)^2)))
print,periodR
maxt=max(alog(abs(waveR)^2))

loadct,39
contour, alog(abs(waveR)^2), Yr(*), periodR, position=[0.38, 0.1, 0.68, 0.54], xsty=1, ysty=1, xrange=[0, S-1], $
  yrange=yrange, /ylog, nlevels=30, /fill,title=' ',$
  ytitle=' ',ytickname=[' ',' ',' '],xticks=6,xminor=3,$ levels=(0.68+indgen(30)/100.)*maxt,c_colors=150+indgen(30)*4,$
  xtickname=[replicate(' ',7)],charsize=1.5,charthick=2.5,thick=2.5,/normal
signifR=rebin(transpose(signifR),ntime,nscaleR)
contour,abs(waveR)^2/signifR, Yr(*), periodR, /overplot, level=level, C_ANNOT=levels, xrange=[0, S-1]
PLOTS,Yr(*), coi, noclip=0

;666===============================================================
R=f35
data=R(0:S-1)
resultR=findgen(S)
resultR(*)=data(0:S-1)
ntime=S
dt=dt ; time resolution 1 year
time=findgen(ntime)*dt
waveR=wavelet(resultR, dt, PERIOD=periodR, COI=coi, /PAD, SIGNIF=signifR)
nscaleR=N_ELEMENTS(periodR)
help,waver
print,signifr
m=max(alog((abs(waveR)^2)))
print,periodR
maxt=max(alog(abs(waveR)^2))

loadct,39
contour, alog(abs(waveR)^2), Yr(*), periodR, position=[0.68, 0.1, 0.98, 0.54], xsty=1, ysty=1, xrange=[0, S-1], $
  yrange=yrange, /ylog, nlevels=30, /fill,title=' ',$
  ytitle=' ',ytickname=[' ',' ',' '],xticks=6,xminor=3,$ levels=(0.68+indgen(30)/100.)*maxt,c_colors=150+indgen(30)*4,$
  xtickname=[replicate(' ',7)],charsize=1.5,charthick=2.5,thick=2.5,/normal
signifR=rebin(transpose(signifR),ntime,nscaleR)
contour,abs(waveR)^2/signifR, Yr(*), periodR, /overplot, level=level, C_ANNOT=levels, xrange=[0, S-1]
PLOTS,Yr(*), coi, noclip=0

;;====================================
;R=df9
;data=R(0:S-1)
;resultR=findgen(S)
;resultR(*)=data(0:S-1)
;
;ntime=S
;
;
;dt=1 ; time resolution 1 year
;time=findgen(ntime)*dt
;
;waveR=wavelet(resultR, dt, PERIOD=periodR, COI=coi, /PAD, SIGNIF=signifR)
;nscaleR=N_ELEMENTS(periodR)
;help,waver
;print,signifr
;
;m=max(alog((abs(waveR)^2)))
;print,periodR
;
;
;loadct,39
;maxt=max(alog(abs(waveR)^2))
;
;contour, alog(abs(waveR)^2), Yr(*), periodR, position=[0.08, 0.1, 0.38, 0.54], xsty=1, ysty=1, xrange=[0, S-1], $
;yrange=[50,800], /ylog, nlevels=30, /fill,title=' ',$
;ytitle='Period (s)',xticks=6,xminor=3,$ levels=(0.68+indgen(30)/100.)*maxt,c_colors=150+indgen(30)*4,$
;xtickname=['05:00','05:15','05:30','05:45','06:00','06:15',' '],charsize=1.5,charthick=2.5,thick=2.5,/normal
;signifR=rebin(transpose(signifR),ntime,nscaleR)
;contour,abs(waveR)^2/signifR, Yr(*), periodR, /overplot, level=0.95, C_ANNOT='95%', xrange=[0, S-1]
;PLOTS,Yr(*), coi, noclip=0
;

;;====================================
;R=fltarr(S)
;R(0:1079)=lc(2,*)
;data=R(0:S-1)
;resultR=findgen(S)
;resultR(*)=data(0:S-1)
;Yr=findgen(S)
;ntime=S
;
;
;dt=1 ; time resolution 1 year
;time=findgen(ntime)*dt
;
;waveR=wavelet(resultR, dt, PERIOD=periodR, COI=coi, /PAD, SIGNIF=signifR)
;nscaleR=N_ELEMENTS(periodR)
;help,waver
;print,signifr
;
;m=max(alog((abs(waveR)^2)))
;print,periodR
;
;
;loadct,39
;maxt=max(alog(abs(waveR)^2))
;
;contour, alog(abs(waveR)^2), Yr(*), periodR, position=[0.38, 0.1, 0.68, 0.54], xsty=1, ysty=1, xrange=[0, S-1], $
;yrange=[50,800], /ylog, nlevels=30, /fill,title=' ',$
;ytitle=' ',ytickname=[' ',' ',' '],xticks=6,xminor=3,$ levels=(0.68+indgen(30)/100.)*maxt,c_colors=150+indgen(30)*4,$
;xtickname=['05:00','05:15','05:30','05:45','06:00','06:15',' '],charsize=1.5,charthick=2.5,thick=2.5,/normal
;signifR=rebin(transpose(signifR),ntime,nscaleR)
;contour,abs(waveR)^2/signifR, Yr(*), periodR, /overplot, level=0.95, C_ANNOT='95%', xrange=[0, S-1]
;PLOTS,Yr(*), coi, noclip=0
;
;;==================================
;R=fltarr(S)
;R(0:1079)=lc(2,*)
;
;resultR=findgen(S)
;resultR(*)=data(0:S-1)
;Yr=findgen(S)
;ntime=S
;
;
;dt=1 ; time resolution 1 year
;time=findgen(ntime)*dt
;
;waveR=wavelet(resultR, dt, PERIOD=periodR, COI=coi, /PAD, SIGNIF=signifR)
;nscaleR=N_ELEMENTS(periodR)
;help,waver
;print,signifr
;
;m=max(alog((abs(waveR)^2)))
;print,periodR
;
;
;loadct,39
;maxt=max(alog(abs(waveR)^2))
;
;contour, alog(abs(waveR)^2), Yr(*), periodR, position=[0.68, 0.1, 0.98, 0.54], xsty=1, ysty=1, xrange=[0, S-1], $
;yrange=[50,800], /ylog, nlevels=30, /fill,title=' ',$
;ytitle=' ',ytickname=[' ',' ',' '],xticks=6,xminor=3,$ levels=(0.68+indgen(30)/100.)*maxt,c_colors=150+indgen(30)*4,$
;xtickname=['05:00','05:15','05:30','05:45','06:00','06:15','06:30'],charsize=1.5,charthick=2.5,thick=2.5,/normal
;signifR=rebin(transpose(signifR),ntime,nscaleR)
;contour,abs(waveR)^2/signifR, Yr(*), periodR, /overplot, level=0.95, C_ANNOT='95%', xrange=[0, S-1]
;PLOTS,Yr(*), coi, noclip=0
;

;xyouts,0.085,0.94,'(a) 1 GHz',/normal,charthick=2.5
;xyouts,0.385,0.94,'(b) 2 GHz',/normal,charthick=2.5
;xyouts,0.685,0.94,'(c) 4 GHz',/normal,charthick=2.5
;xyouts,0.085,0.5,'(d) 9 GHz',/normal,charthick=2.5
;xyouts,0.385,0.5,'(e) 17 GHz',/normal,charthick=2.5
;xyouts,0.685,0.5,'(f) 35 GHz',/normal,charthick=2.5

device, /close
set_plot,'X'
!p.multi=[0,1,1]
end