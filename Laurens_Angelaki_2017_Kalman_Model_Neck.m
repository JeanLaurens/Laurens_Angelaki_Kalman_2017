function [Out] = Laurens_Angelaki_2017_Kalman_Model_Neck(t, OmegaN,OmegaB, OmegaN_u,OmegaB_u, dt, param)

% Inputs
% t: time (nx1 vector)
% OmegaN: real head on body angular velocity (nx1 vector)
% OmegaN: real body in space angular velocity (nx1 vector)
% OmegaN_u: head on body angular velocity motor command (nx1 vector)
% OmegaN_u: body in space angular velocity motor command (nx1 vector)
% dt: time step
% Param: matrix of parameters (see below)

% Parameters
% dt=0.01;
% param.sN = 200*pi/180;
% param.sB = 40*pi/180;
% param.sv=10*pi/180;
% param.sP=10*dt*pi/180;
% param.tau=4;

N = cumsum(OmegaN)*dt;
k = 1-dt/param.tau ;
k1 = param.tau/(param.tau+dt) ;
k2 = dt/(param.tau+dt) ;

A = [0 0 0 0;...
     0 0 0 0;...
     0 0 1 0;...
     0 0 0 k1] ;

B = [1 0;...
     0 1;...
     dt 0;...
     k2 k2];

eN = [1 0 dt k2]*param.sN;
eB = [0 1 0 k2]*param.sB;

Q = eN'*eN+eB'*eB ;

H = [1 1 0 -1;0 0 1 0] ;

R = [param.sv*param.sv 0;0 param.sP*param.sP] ;

x = zeros(4,1) ; 
v = 0 ; % 
ht_1 = 0 ;

P = zeros(4,4) ;
[m,n] = size(t) ; n = max([n m]) ;

for i = 1:500
    Pm = A*P*(A') + Q ;
    K = Pm*(H')*((H*Pm*(H')+R)^(-1)) ;
    P = (eye(size(K*H)) - K*H)*Pm ;
end


for i = 1:n
    u = [OmegaN_u(i); OmegaB_u(i)] ;
    a = (OmegaN(i)+OmegaB(i)-ht_1)/dt ; ht_1 = OmegaN(i)+OmegaB(i) ; v = k*v + a*dt ;
    
    z = [v;N(i)] ;
    xm = A*x + B*u ;
    Pm = A*P*(A') + Q ;
    K = Pm*(H')*((H*Pm*(H')+R)^(-1)) ;
    P = (eye(size(K*H)) - K*H)*Pm ;
    x = xm + K*(z-H*xm) ;

    Out(i).time = t(i) ;
    Out(i).Xf = x ;
    Out(i).Xm = xm ;
    Out(i).K = K ;
    Out(i).Zp = H*xm;
    Out(i).Z = z ;
end

