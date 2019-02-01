PRO FFT_IMAGE,U,V,VIS,BEAM,IMG

FREQ=4.1875D  ;GHZ
LAMDA=0.3/FREQ

PIX=1000.
nant=60
BPIX=PIX*2.
Fov=1.0/180*!pi
maxUV=Pix/Fov
maxBUV=Bpix/Fov

IMG=FLTARR(PIX,PIX)
BEAM=FLTARR(BPIX,BPIX)
UV_IMG=COMPLEXARR(PIX,PIX)
UV_BEAM=FLTARR(BPIX,BPIX)

;FOR I=0,Nant-1 DO BEGIN
; FOR J=0,Nant-1 DO BEGIN
;  X=U[I,J]/MAXUV*PIX
;  Y=V[I,J]/MAXUV*PIX
;  BX=U[I,J]/MAXBUV*BPIX
;  BY=V[I,J]/MAXBUV*BPIX
;  IF X GE PIX/2 THEN X=0
;  IF X LE -PIX/2 THEN X=0
;  IF Y GE PIX/2 THEN Y=0
;  IF Y LE -PIX/2 THEN Y=0
;  IF BX GE BPIX/2 THEN BX=0
;  IF BX LE -BPIX/2 THEN BX=0
;  IF BY GE BPIX/2 THEN BY=0
;  IF BY LE -BPIX/2 THEN BY=0
;  UV_IMG[X+PIX/2,Y+PIX/2]=UV_IMG[X+PIX/2,Y+PIX/2]+VIS[I,J]
;  UV_BEAM[BX+BPIX/2,BY+BPIX/2]=1.0
; ENDFOR
;ENDFOR

m=6
WGT_RG=1
big_matrix=complexarr(pix,pix)
FOR I=0,Nant-1 DO BEGIN
 FOR J=0,Nant-1 DO BEGIN
   if i ne j then begin 
     X=U[I,J]/MAXuv*PIX
     Y=V[I,J]/MAXuv*PIX
     BX=U[I,J]/MAXuv*BPIX
     BY=V[I,J]/MAXuv*BPIX
      IF X GE PIX/2-M/2 THEN X=0
      IF X LE -PIX/2+M/2 THEN X=0
      IF Y GE PIX/2-M/2 THEN Y=0
      IF Y LE -PIX/2+M/2 THEN Y=0
      IF BX GE BPIX/2-M/2 THEN BX=0
      IF BX LE -BPIX/2+M/2 THEN BX=0
      IF BY GE BPIX/2-M/2 THEN BY=0
      IF BY LE -BPIX/2+M/2 THEN BY=0
      if x ne 0 and y ne 0 then begin
       gridding,x,y,mini_matrix
       ;print,mini_matrix
       ;big_matrix=big_matrix-big_matrix
       big_matrix[FIX(X)-m/2+pix/2:FIX(X)+m/2+pix/2-1,FIX(y)-m/2+pix/2:FIX(y)+m/2+pix/2-1]=$
         big_matrix[FIX(X)-m/2+pix/2:FIX(X)+m/2+pix/2-1,FIX(y)-m/2+pix/2:FIX(y)+m/2+pix/2-1]+mini_matrix*vis[i,j]
       ;WGT_FACT=TOTAL(UV_IMG[FIX(X)-WGT_RG+pix/2:FIX(X)+WGT_RG+pix/2,FIX(y)-WGT_RG+pix/2:FIX(y)+WGT_RG+pix/2])
       UV_IMG=UV_IMG+big_matrix;/WGT_FACT
      endif
   ENDIF
 ENDFOR
ENDFOR

BEAM=FFT(UV_BEAM)
BEAM=SHIFT(BEAM,BPIX/2,BPIX/2)
BEAM=ROTATE(BEAM,5)
IMG=FFT(UV_IMG)
IMG=SHIFT(IMG,PIX/2,PIX/2)
IMG=ROTATE(IMG,5)
END