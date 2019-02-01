PRO FFT_beam,U,V,BEAM
pix=1024
nant=60
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
;BEAM=ROTATE(BEAM,5)

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