% exINS_MEKF is compatible with MATLAB and GNU Octave (www.octave.org).
% exINS_MEKF error-state (indirect) feedback Kalman filter for INS aided by
% position measurements. The attitude is estimated using the MEKF formalism 
% where attitude is parametrized using unit quaternions/Gibbs vector.
%
% The GNSS position measurement frequency f_gnss can be chosen smaller or
% equal to the  sampling frequency f_s, which is equal to the IMU
% measurement frequency. The ratio between the frequencies must be an integer:
%
%     Integer:  Z = f_s/f_gnss >= 1 
%
% The main loop calls ins_mekf(...) for magnetometer measurements and
% ins_mekf_psi(...) for compass measurements:
%
% [x_ins, P_prd] = ins_mekf( ...
%    x_ins, P_prd, mu, h, Qd, Rd, f_imu, w_imu, m_imu, m_ref, y_pos, y_vel)
%
% [x_ins, P_prd] = ins_mekf_psi( ...
%     x_ins, P_prd, mu, h, Qd, Rd, f_imu, y_psi, y_pos, y_vel)
%
% each time a GNSS position y_pos is received at the slow frequency f_gnss.
% For samples without new measurements, the arguments y_pos and y_vel are
% omitted when calling the function ins_ahrs when there are no new GNSS 
% measurements.
%
% Author:    Thor I. Fossen
% Date:      26 Mar 2020
% Revisions: 28 Mar 2020, modified to use an INS signal generator
%            29 Nov 2020, new Kalman filter weights

%% USER INPUTS
T_final = 100;	  % Final simulation time (s)
f_s    = 100;     % Sampling frequency (Hz)
f_gnss = 1;       % GNSS measurement frequency (Hz)

% Flags
mag = 0;          % 0 = compass, 1 = magnetometer 
vel = 1;          % 0 = no velocity meaurement, 1 = velocity aiding

% Parameters
Z = f_s/f_gnss;   % Ratio betwween sampling/IMU frequencies
h  = 1/f_s; 	  % Sampling time: h  = 1/f_s (s) 
h_gnss = 1/f_gnss;  

% Biases
b_acc = [0.1 0.3 -0.1]';
b_ars = [0.05 0.1 -0.05]';
   
% Initial values for x for signal generator
x = [zeros(1,6) b_acc' zeros(1,3) b_ars']';	        

% Initialization of Kalman filter
P_prd = eye(15);

% Process noise weights: v, acc_bias, w, ars_bias
Qd = diag([0.01 0.01 0.01 0.01 0.01 0.01 0.1 0.1 0.1 0.001 0.001 0.001]);
   
if (vel == 0 && mag == 1)  % Position aiding + magnetometer
     
   Rd = diag([1 1 1  1 1 1 0.01 0.01 0.01]);  % p, acc, mag  
   
elseif  (vel == 1 && mag == 1)  % Position/velocity aiding + magnetometer 
     
   Rd = diag([1 1 1  0.1 0.1 0.1  1 1 1  0.01 0.01 0.01]);  % p, v, acc, mag
   
elseif (vel == 0 && mag == 0)  % Position aiding + compass 
    
   Rd = diag([1 1 1  1 1 1  0.01]);  % p, acc, psi
   
else % Position/velocity aiding + compass 
      
   Rd = diag([1 1 1  1 1 1  1 1 1  0.01]);  % p, vel, acc, psi

end

% Initialization of INS
p_ins = [0 0 0]'; 
v_ins = [0 0 0]';
b_acc_ins = [0 0 0]';
q_ins = euler2q(0, 0, 0);
b_ars_ins = [0 0 0]';
x_ins = [p_ins; v_ins; b_acc_ins; q_ins; b_ars_ins];

% Time vector initialization
t = 0:h:T_final;                % Time vector from 0 to T_final          
nTimeSteps = length(t);         % Number of time steps

% WGS-84 gravity model
mu = 63.4305 * pi / 180;    % Lattitude  
g = gravity(mu);  

%% Display
disp('----------------------------------------------------------');
disp('MSS toolbox: Error-state (indirect) feedback Kalman filter');
disp('Attitude parametrization: 2 x Gibbs vector (MEKF)');
if (vel == 0)
   disp(['INS aided by position at ',num2str(f_gnss), ' Hz']);
else
   disp(['INS aided by position and velocity at ',num2str(f_gnss),' Hz']);  
end
   disp(['IMU measurements (specific force and ARS) at ',num2str(f_s),' Hz']);
if (mag == 1)
   disp(['MAGNETOMETER measurements at ',num2str(f_s), ' Hz']);
else
   disp(['COMPASS measurements at ',num2str(f_s), ' Hz']);
end
disp('----------------------------------------------------------');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MAIN LOOP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
simdata = zeros(nTimeSteps,31); % Pre-allocate table of simulation data
ydata = [0 x(1:3)'];            % Pre-allocatetable of position measurements

for i=1:nTimeSteps
    
    % INS signal generator  
    [x, f_imu, w_imu, m_imu, m_ref] = insSignal(x, mu, h, t(i));
    y_psi = x(12);
    
    % GNSS measurements are Z times slower than the sampling time
    if mod( t(i), h_gnss ) == 0
        
        y_pos = x(1:3) + 0.05 * randn(3,1);   % Position measurements
        y_vel = x(4:6) + 0.01 * randn(3,1);   % Optionally velocity meas.
        ydata = [ydata; t(i), y_pos'];        % Store position measurements                  
              
        if (vel == 0 && mag == 1)  % Position aiding + magnetometer
            
            [x_ins,P_prd] = ins_mekf(...
                x_ins,P_prd,mu,h,Qd,Rd,f_imu,w_imu,m_imu,m_ref,y_pos);
            
        elseif (vel == 1 && mag == 1)  % Position/velocity aiding + magnetometer
            
             [x_ins,P_prd] = ins_mekf(...
                 x_ins,P_prd,mu,h,Qd,Rd,f_imu,w_imu,m_imu,m_ref,y_pos,y_vel);
        
        elseif (vel == 0 && mag == 0)  % Position aiding + compass
            
              [x_ins,P_prd] = ins_mekf_psi(...
                 x_ins,P_prd,mu,h,Qd,Rd,f_imu,w_imu,y_psi,y_pos);
             
        else  % Position/velocity aiding + compass
            
              [x_ins,P_prd] = ins_mekf_psi(...
                 x_ins,P_prd,mu,h,Qd,Rd,f_imu,w_imu,y_psi,y_pos,y_vel);
             
        end
        
    else  % No aiding
        
        if (mag == 1)  % Magnetometer
        
           [x_ins,P_prd] = ins_mekf(...
               x_ins,P_prd,mu,h,Qd,Rd,f_imu,w_imu,m_imu,m_ref);
             
        else           % Compass
            
           [x_ins,P_prd] = ins_mekf_psi(...
               x_ins,P_prd,mu,h,Qd,Rd,f_imu,w_imu,y_psi);          
           
        end
        
    end
       
    % Store simulation data in a table (for testing)
    simdata(i,:) = [x' x_ins']; 
    
end

%% PLOTS          
x     = simdata(:,1:15); 
x_hat = simdata(:,16:31); 

phi = zeros(length(t),1);theta = zeros(length(t),1);psi = zeros(length(t),1);
for i = 1:length(t)
    [phi(i), theta(i), psi(i)] = q2euler(x_hat(i,10:13));
end
Theta = [phi theta psi];

t_m = ydata(:,1);              % Slow GNSS measurements
y_m = ydata(:,2:4);

figure(1); figure(gcf)

subplot(311)
h1 = plot(t_m,y_m,'xr'); hold on;
h2 = plot(t,x_hat(:,1:3),'b'); hold off;
xlabel('time (s)'),title('Position [m]'),grid
legend([h1(1),h2(1)],['Measurement at ', num2str(f_gnss), ' Hz'],...
    ['Estimate at ', num2str(f_s), ' Hz'] );

subplot(312)
h1 = plot(t,x(:,4:6),'r'); hold on;
h2 = plot(t,x_hat(:,4:6),'b'); hold off;
xlabel('time (s)'),title('Velocity [m/s]'),grid
legend([h1(1),h2(1)],['True velocity at ', num2str(f_s), ' Hz'],...
    ['Estimate at ', num2str(f_s), ' Hz'] );

subplot(313)
h1 = plot(t,x(:,7:9),'r'); hold on;
h2 = plot(t,x_hat(:,7:9),'b'); hold off;
xlabel('time (s)'),title('Acc bias'),grid
legend([h1(1),h2(1)],['True acc bias at ', num2str(f_s), ' Hz'],...
    ['Estimate at ', num2str(f_s), ' Hz'] );

set(findall(gcf,'type','line'),'linewidth',2)
set(findall(gcf,'type','text'),'FontSize',14)
set(findall(gcf,'type','legend'),'FontSize',14)

figure(2); figure(gcf)

subplot(211)
h1 = plot(t,rad2deg( x(:,10:12) ),'r'); hold on;
h2 = plot(t,rad2deg( Theta ),'b'); hold off;
xlabel('time (s)'),title('Angle [deg]'),grid
legend([h1(1),h2(1)],['Measurement at ', num2str(f_s), ' Hz'],...
    ['Estimate at ', num2str(f_s), ' Hz'] );

subplot(212)
h1 = plot(t,x(:,13:15),'r'); hold on;
h2 = plot(t,x_hat(:,14:16),'b'); hold off;
xlabel('time (s)'),title('ARS bias'),grid
legend([h1(1),h2(1)],['True ARS bias at ', num2str(f_s), ' Hz'],...
    ['Estimate at ', num2str(f_s), ' Hz'] );

set(findall(gcf,'type','line'),'linewidth',2)
set(findall(gcf,'type','text'),'FontSize',14)
set(findall(gcf,'type','legend'),'FontSize',14)
