%% EXPLANATION
% These are example of simulation performed with the Kalman filter model from
%     Laurens and Angelaki, eLife 2017
%
% To run simulations:
% Execute Section 1 (Model Parameters) 
% Execute one of the Section 2 (Model inputs)
% Execute Section 3 to perform the simulation and draw a figure
%% Section 1: Model Parameters
dt = 0.01 ;

param.so = 40*pi/180;
param.sa=0.3;
param.sv=10*pi/180;
param.tau=4;
param.sf=0.002;
param.svis = 7*pi/180;
param.sensory_noise = 0 ;

%% Section 2: Short or long duration Yaw rotations around a earth-vertical axis, active or passive

total_duration = 10; % Total duration of the simulation
motion_period = [1 5]; % Interval during which the head rotates
is_active = false ; % Set to true or false to simulate active/passive rotation
in_light = true ; % Set to true or false to simulate rotation in light or darkness

time = (0:dt:total_duration)' ;
Omega = time*0 ; Omega(time>=motion_period(1) & time<=motion_period(2)) = 1 ;
A = time*NaN ; % Real acceleration, set to NaN to tell the model to simulate a rotation around a earth-vertical axis
G = time*NaN ; % Real tilt, set to NaN to tell the model to simulate a rotation around a earth-vertical axis
F = A+G ;
A_u = time*0 ;
if is_active
    Omega_u = Omega ;
else
    Omega_u = Omega*0 ;
end
if in_light
    Omega_vision = Omega ;
else
    Omega_vision = Omega*NaN ;
end

%% Section 2: Short or long duration translation

total_duration = 10; % Total duration of the simulation
motion_period = [1 5]; % Interval during which the head translates
is_active = false ; % Set to true or false to simulate active/passive translation
in_light = false ; % Set to true or false to simulate rotation in light or darkness - note that the model includes visual signals that encode rotations only

time = (0:dt:total_duration)' ;
Omega = time*0 ; 
A = time*0 ; A(time>=motion_period(1) & time<=motion_period(2)) = 1 ;% Real acceleration
G = time*0 ; % Real tilt
F = A+G ;
Omega_u = time*0 ;
if is_active
    A_u = A ;
else
    A_u = A*0 ;
end
if in_light
    Omega_vision = Omega ;
else
    Omega_vision = Omega*NaN ;
end

%% Section 2: Short or long duration tilt

total_duration = 2; % Total duration of the simulation
motion_period = [0.5 1.5]; % Interval during which the head tilts (at a constant velocity)
is_active = false ; % Set to true or false to simulate active/passive tilt
in_light = false ; % Set to true or false to simulate rotation in light or darkness - note that the model includes visual signals that encode rotations only

time = (0:dt:total_duration)' ;
Omega = time*0 ; Omega(time>=motion_period(1) & time<=motion_period(2)) = 1 ;
A = time*0 ; % Real acceleration
G = cumsum(Omega)*dt ; % Real tilt
F = A+G ;

Omega_u = time*0 ;
if is_active
    Omega_u = Omega ;
else
    Omega_u = Omega*0 ;
end
if in_light
    Omega_vision = Omega ;
else
    Omega_vision = Omega*NaN ;
end

%% Section 3: Perform the simulation and draw a figure
[Result] = Laurens_Angelaki_2017_Kalman_Model(time, Omega, F, Omega_u, A_u, dt, param, Omega_vision) ;
clf
X = [Result.X]' ;
Xm = [Result.Xm]' ;
Vestibular_Feedback = [Result.Vestibular_Feedback]' ;
Otolith_Feedback = [Result.Otolith_Feedback]' ;
Visual_Feedback = [Result.Visual_Feedback]' ;
Xf = [Result.Xf]' ;

Z = [Result.Z]' ;
Zp = [Result.Zp]' ;
delta_S = [Result.delta_S]' ;

variable_names = {'\Omega','C','G','A'} ;
variable_colors = [0 0 1;0 1 1;0 0.5 0;1 0 0] ;

sensory_names = {'V','F','Vision'} ;
sensory_colors = [1 0 1;0.5 0.5 0.5;0 1 0] ;

for i = 1:4
   subplot(7,4,(i-1)*4+1) ;
   plot(time,X(:,i),'LineWidth',2,'Color',variable_colors(i,:),'Clipping','off') ;
   if i == 1, title('Real Value') ; end ; if i == 4, xlabel('time (s)') ; end 
   ylabel(variable_names{i});
   
   subplot(7,4,(i-1)*4+2) ;
   plot(time,Xm(:,i),'LineWidth',2,'Color',variable_colors(i,:),'Clipping','off') ;
   if i == 1, title('Predicted Value') ; end ; if i == 4, xlabel('time (s)') ; end 
   
   subplot(7,4,(i-1)*4+3) ;hold on
   plot(time,Vestibular_Feedback(:,i),'LineWidth',2,'Color',sensory_colors(1,:),'Clipping','off') ;
   plot(time,Otolith_Feedback(:,i),'LineWidth',2,'Color',sensory_colors(2,:),'Clipping','off') ;
   plot(time,Visual_Feedback(:,i),'LineWidth',2,'Color',sensory_colors(3,:),'Clipping','off') ;
    if i == 1, title('Kalman Feedback') ; end ; if i == 4, xlabel('time (s)') ; end
    box on

   subplot(7,4,(i-1)*4+4) ;
   plot(time,Xf(:,i),'LineWidth',2,'Color',variable_colors(i,:),'Clipping','off') ;
   if i == 1, title('Final Value') ; end ; if i == 4, xlabel('time (s)') ; end   
end

for i = 1:size(Z,2)
   subplot(7,4,(i+3)*4+1) ;
   plot(time,Z(:,i),'LineWidth',2,'Color',sensory_colors(i,:),'Clipping','off') ;
   if i == 1, title('Real Value') ; end ; if i == size(Z,2), xlabel('time (s)') ; end 
   ylabel(sensory_names{i});
      
   subplot(7,4,(i+3)*4+2) ;
   plot(time,Zp(:,i),'LineWidth',2,'Color',sensory_colors(i,:),'Clipping','off') ;
   if i == 1, title('Predicted Value') ; end ; if i == size(Z,2), xlabel('time (s)') ; end 
   
   subplot(7,4,(i+3)*4+3) ;
   plot(time,delta_S(:,i),'LineWidth',2,'Color',sensory_colors(i,:),'Clipping','off') ;
   if i == 1, title('Sensory Prediction Error') ; end ; if i == size(Z,2), xlabel('time (s)') ; end    
  
end

subplot_handle=[] ;
for i = 1:4, for j = 1:4, subplot_handle(i,j) = subplot(7,4,(i-1)*4+j); end; end
for i = 4+(1:size(Z,2)), for j = 1:3, subplot_handle(i,j) = subplot(7,4,(i-1)*4+j); end; end

linkaxes(subplot_handle(subplot_handle~=0),'x') ;
for i = 1:length(subplot_handle(:))
    if subplot_handle(i)~=0
        y = get(subplot_handle(i),'YLim') ; y(1)=min([y(1) -1]);y(2)=max([y(2) 1]); set(subplot_handle(i),'YLim',y)
    end
end
