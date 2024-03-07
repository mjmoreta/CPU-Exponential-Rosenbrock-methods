function valz=funhtt(tt,xx)

% This function is used to calculate h_tt(t,x)

valz=sin(xx+tt)-cos(tt+xx)+2*cos(xx+tt)^2-2*sin(xx+tt)^2;

