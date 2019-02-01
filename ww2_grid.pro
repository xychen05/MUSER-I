pro grid,Lin,Lout

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

end