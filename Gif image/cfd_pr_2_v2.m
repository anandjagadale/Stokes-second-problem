clear all
close all
clc

%% 
ADIABATIC_LOWER_PLATE = false; 
filename = 'OscillatingTemp_velo.gif';

%% PHYSICAL CONSTANTS 
freq = 1000;             
omega = 2 * pi * freq;  
Period = 1 / freq;      

U_inf = 200;            
h_inch = 0.01;          
h = h_inch / 12;        
T0 = 519;               

% Fluid Properties
rho = 0.00237;          
mu = 3.737e-7;          
J = 778.17; g_c = 32.174; cp_btu = 0.24;          
cp = cp_btu * g_c * J; Pr = 0.71; k = (mu * cp) / Pr;     

%% SIMUALTION SETUP 
n = 61;                 
dy = h / (n - 1);       
y = linspace(0, h, n);  

% Stability Time Step
alpha = k / (rho * cp); 
dt_limit_v = 0.5 * (rho * dy^2) / mu;  
dt_limit_t = 0.5 * dy^2 / alpha;       
dt = min(dt_limit_v, dt_limit_t) * 0.9; 

%% MAKING VIDEO
steps_per_cycle = Period / dt;
plot_interval = floor(steps_per_cycle / 60); % set fps here 

fprintf('Total steps per cycle: %.0f\n', steps_per_cycle);

% Creating the figure window
h_fig = figure('Color','black');
axis tight manual 

%% IC
u = zeros(1, n);        
T = ones(1, n) * T0;    
u(n) = 0;          % MAKING A CHANGE HERE   
if ~ADIABATIC_LOWER_PLATE, T(1) = T0; end

%% Main Loop
t_now = 0;
t_max = 2.0 * Period; % three cycles of oscillation 
iter = 0;

u_new = u;
T_new = T;

while (t_now < t_max)
    
    % Momentum eqn
    for i = 2:(n-1)
        diffusion = (mu/rho) * (u(i+1) - 2*u(i) + u(i-1)) / dy^2;
        u_new(i) = u(i) + dt * diffusion;
    end
    
    % Energy eqn
    for i = 2:(n-1)
        diff_T = (k / (rho*cp)) * (T(i+1) - 2*T(i) + T(i-1)) / dy^2;
        du_dy = (u(i+1) - u(i-1)) / (2*dy); 
        visc_heat = (mu / (rho*cp)) * (du_dy)^2;
        T_new(i) = T(i) + dt * (diff_T + visc_heat);
    end
    
    % BC
    u_new(1) = 0; 
    u_new(n) = U_inf * sin(omega * (t_now + dt)); 
    T_new(n) = T0;
    
    if ADIABATIC_LOWER_PLATE
        T_new(1) = T_new(2); 
    else
        T_new(1) = T0;       
    end
    
    % GIF creation
    if mod(iter, plot_interval) == 0
        % Velocity plot
        subplot(1,2,1);
        plot(u_new, y*12, 'b-', 'LineWidth', 2);
        xlim([-U_inf-20, U_inf+20]); ylim([0, h_inch]);
        title('Velocity'); grid on;
        
        % Temperature Plot
        subplot(1,2,2);
        plot(T_new, y*12, 'r-', 'LineWidth', 2);
        xlim([T0 - 0.05, 523.5 + 0.1]); ylim([0, h_inch]);
        title('Temperature'); grid on;
        
        sgtitle(sprintf('t = %.5f s', t_now));
        
        drawnow;
        frame = getframe(gcf); 
        im = frame2im(frame); 
        [imind, cm] = rgb2ind(im, 256); 
        
        % making the video file 
        if iter == 0 
            % First frame
            imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', 0.1); 
        else 
            % Subsequent frames
            imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.1); 
        end 
    end  
    
    u = u_new;
    T = T_new;
    t_now = t_now + dt;
    iter = iter + 1;
end

fprintf(' GIF saved as %s\n', filename);
