PRO FFT_IMAGE,U,V,VIS,IMG

PIX=1024.
nant=60
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
;IMG=ROTATE(IMG,5)

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