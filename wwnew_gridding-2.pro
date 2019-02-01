pro gridding,x,y,cfunc

nsize=size(uv_s)
pix=nsize[1]
m=6
locu=dblarr(m,m)
locv=dblarr(m,m)
for i=0,m-1 do begin
  for j=0,m-1 do begin
      locu[i,j]=(x mod 1)+(i-m/2)
      locv[i,j]=(y mod 1)+(j-m/2)
  endfor
endfor

select_func=3
case(select_func) of
  0: begin ;pillbox
    cfunc=0
  
  end
  1: begin ; exponential
    w=1.
    alfa=2.
    cfunc=exp(-(abs(locu)/w)^alfa-(abs(locv)/w)^alfa)
  end
  2: begin ; sinc
    w=1.
    cfunc=sin(locu/w)/locu/w*sin(locv/w)/locv/w
  end
  3: begin ; exponential times sinc
    w1=2.52
    w2=1.55
    alfa=2.
    cfunc=exp(-abs(locu/w1)^alfa-abs(locv/w1)^alfa)*sin(locu/w2)/locu/w2*sin(locv/w2)/locv/w2
  end
  else: begin
  end
endcase

cfunc=cfunc/total(cfunc)
;print,'qq'
end