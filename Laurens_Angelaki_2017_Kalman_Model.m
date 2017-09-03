function [Out] = Laurens_Angelaki_2017_Kalman_Model(t, Omega, F, Omega_u, Acc_u, dt, param, Vis)

% Inputs
% t: time (nx1 vector)
% Omega: real head angular velocity (nx1 vector)
% F: real GIA (nx1 vector). Fill with NaN to simulate rotations around an
% earth-vertical axis.
% Omega_u: angular head velocity motor command (nx1 vector)
% Acc_u: linear head acceleration motor command (nx1 vector)
% dt: time step
% Param: matrix of parameters (see below)
% Vis: visual head velocity signals (optional) (nx1 vector)

% Parameters
% dt=0.01;
% param.so = 40*pi/180;
% param.sa=0.3;
% param.sv=10*pi/180;
% param.tau=4;param.Rtau=4;
% param.sf=0.002;
% param.svis = 7*pi/180;
% param.sensory_noise = 0 ;

if ~exist('Vis','var'), Vis = [] ; end

k = 1-dt/param.tau ;
k1 = param.tau/(param.tau+dt) ;
k2 = dt/(param.tau+dt) ;

if isempty(Vis)
    H = [1 -1 0 0;0 0 1 1] ;
    R = [param.sv^2 0;0 param.sf^2] ;
else
    H = [1 -1 0 0;0 0 1 1;1 0 0 0] ;
    R = [param.sv^2 0 0;0 param.sf^2 0;0 0 param.svis^2] ;
end

x = zeros(4,1) ; 
v = 0 ; 
ht_1 = 0 ;

P = zeros(4,4) ;
[m,n] = size(t) ; n = max([n m]) ;
if isnan(F(1)), sw = 0;F(1)=0;else sw = 1;end

A = [0 0 0 0;...
    0 k1 0 0;...
    0 0 1 0;...
    0 0 0 0] ;

B = [1     0;...
    k2     0;...
    sw*dt  0;...
    0      1];

eo = [1 k2 sw*dt 0]*param.so;
ea = [0 0 0 param.sa];
Q = eo'*eo+ea'*ea ;

for i = 1:500
    Pm = A*P*(A') + Q ;
    K = Pm*(H')*((H*Pm*(H')+R)^(-1)) ;
    P = (eye(size(K*H)) - K*H)*Pm ;
end


for i = 1:n
    
    if isnan(F(i)), sw = 0;F(i)=0;else sw = 1;end
    A = [0 0 0 0;...
        0 k1 0 0;...
        0 0 1 0;...
        0 0 0 0] ;
    
    B = [1     0;...
        k2     0;...
        sw*dt  0;...
        0      1];

    u = [Omega_u(i); Acc_u(i)] ;
    a = (Omega(i)-ht_1)/dt ; ht_1 = Omega(i) ; v = (1-dt/param.Rtau)*v + (param.Rtau/(param.Rtau+dt))*a*dt ;

    if isempty(Vis)
        z = [v;F(i)] ;
    else
        if isnan(Vis(i))
            z = [v;F(i);0] ;
            R = [param.sv^2 0 0;0 param.sf^2 0;0 0 10000^2] ;
        else
            z = [v;F(i);Vis(i)] ;
            R = [param.sv^2 0 0;0 param.sf^2 0;0 0 param.svis^2] ;
        end
    end
    
    if param.sensory_noise
        z(1)=z(1)+randn()*param.sv ;
        z(2)=z(2)+randn()*param.sf ;
        if size(z,1)>2
            z(3)=z(3)+randn()*param.svis ;
        end
    end
    
    eo = [1 k2 sw*dt 0]*param.so;
    ea = [0 0 0 param.sa];
    Q = eo'*eo+ea'*ea ;
    
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
