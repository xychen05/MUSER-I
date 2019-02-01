pro weighting,u,v,vis
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
end