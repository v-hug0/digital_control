function [omega_0,dmp]=omega_dmp(tm_req,M_req);
%function [omega_0,dmp]=omega_dmp(tm_req,M_req);
%computes natural frequency omega_0 and damping dmp of a continues 
%2nd order system with required properties: maximum overshoot M_req
%and rise time tm_req. The continues system has a standard form:
%                       omega_0
%   G2(s)= ---------------------------------
%          s^2 + 2.dmp.omega_0.s + omega_0^2

%inputs:
%tm_req ... required rise time in seconds
%M_req ... required maximum overshoot
%outputs:
%omega_0 ... natural frequency of the continues system
%dmp ... damping of the continues system
%
%written by:  H. Prochazka, I.D. Landau
%7th june 2002

%damping computing 
precision=0.1;%precision of 0.1%
omega_0=6;%initial value 6 rad/s to compute damping for overshoot
dmp_min=0;
dmp_max=1;
error=precision*10;%initial setting
while error>precision,
	dmp_act=dmp_min+ (dmp_max-dmp_min)/2;
	sys=tf([omega_0^2],[1 2*dmp_act*omega_0 omega_0^2]);
	[resp,tim]=step(sys);
	M_act=(max(resp)-1)*100;
	if M_act<M_req, 
	   dmp_max=dmp_act;
	else
	   dmp_min=dmp_act;
   end;
   error=abs(M_act-M_req);
end;
dmp=dmp_act;

%natural frequency computing
k=1;while resp(k)<0.1, k=k+1;end;
t01=tim(k-1);
k=1;while resp(k)<0.9, k=k+1;end;
t09=tim(k);
tm_act=t09-t01;
omega_0=tm_act*omega_0/tm_req;

