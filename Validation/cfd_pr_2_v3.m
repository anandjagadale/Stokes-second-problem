clear all
close all
clc

%% PARAMETERS
freq = 1000;             
omega = 2 * pi * freq;   
Period = 1 / freq;       
U_inf = 200;             % Max Amplitude 
h_inch = 0.01;           
h = h_inch / 12;         % Gap in FEET 
T0 = 519;                % Wall Temp  

rho = 0.00237;          
mu = 3.737e-7;          
J = 778.17;             
g_c = 32.174;           
cp_btu = 0.24;          
cp = cp_btu * g_c * J;  
Pr = 0.71;              
k = (mu * cp) / Pr;     
gamma = 1.4;            
R = 1716;               
ADIABATIC_LOWER_PLATE = false; 

%% GRID INDEPENDENCE STUDY 
Grid_Sizes = [11, 21, 51, 71, 101, 201]; 
colors = {'r', 'b', 'y','g','m','c'};           
line_styles = {'-', '-', '-','-','-','-'};     


SF = 0.5; 


figure(1); clf; hold on; 
title('Velocity Profile at Peak Velocity (t = 1.25 Periods)'); 
xlabel('Velocity (ft/s)'); ylabel('Gap Height (in)');

figure(2); clf; hold on; 
title('Temperature Profile at Peak Velocity (t = 1.25 Periods)'); 
xlabel('Temperature (R)'); ylabel('Gap Height (in)');


%% MAIN LOOP 
for g_idx = 1:length(Grid_Sizes)
    
    n = Grid_Sizes(g_idx);
    
    
    dy = h / (n - 1);       
    y = linspace(0, h_inch, n);  
    
    
    alpha = k / (rho * cp); 
    dt_limit_v = 0.5 * (rho * dy^2) / mu;  
    dt_limit_t = 0.5 * dy^2 / alpha;       
    dt = min(dt_limit_v, dt_limit_t) * SF; 
    
    

    u = zeros(1, n);        
    T = ones(1, n) * T0;    
    u(1) = 0; u(n) = 0; T(n) = T0;
    if ~ADIABATIC_LOWER_PLATE, T(1) = T0; end
    
    % Loop 
    t_now = 0;
    t_max = 2.0 * Period;
    

    captured_peak = false; 
    u_peak = u; 
    T_peak = T; 
    
    % Simulation Loop
    while (t_now < t_max)
        u_old = u;
        T_old = T;
        
        % Momentum
        for i = 2:(n-1)
            diffusion = (mu/rho) * (u_old(i+1) - 2*u_old(i) + u_old(i-1)) / dy^2;
            u(i) = u_old(i) + dt * diffusion;
        end
        
        % Energy
        for i = 2:(n-1)
            diff_T = (k / (rho*cp)) * (T_old(i+1) - 2*T_old(i) + T_old(i-1)) / dy^2;
            du_dy = (u_old(i+1) - u_old(i-1)) / (2*dy); 
            visc_heat = (mu / (rho*cp)) * (du_dy)^2;
            T(i) = T_old(i) + dt * (diff_T + visc_heat);
        end
        
        % Boundary Conditions
        u(1) = 0; 
        u(n) = U_inf * sin(omega * (t_now + dt)); 
        T(n) = T0;
        if ADIABATIC_LOWER_PLATE
            T(1) = T(2); 
        else
            T(1) = T0;       
        end
        
        % at 1.25 Periods 
        target_time = 1.25 * Period;
        if (t_now >= target_time) && (captured_peak == false)
            u_peak = u;
            T_peak = T;
            captured_peak = true; 
        end
        
        t_now = t_now + dt;
    end
    fprintf('Complete.\n');
    
    figure(1);
    plot(u_peak, y, 'Color', colors{g_idx}, 'LineStyle', line_styles{g_idx}, 'LineWidth', 1.5, ...
        'DisplayName', sprintf('Grid N = %d', n));
    
    figure(2);
    plot(T_peak, y, 'Color', colors{g_idx}, 'LineStyle', line_styles{g_idx}, 'LineWidth', 1.5, ...
        'DisplayName', sprintf('Grid N = %d', n));
end

%% POST-PROCESSING

figure(1); legend('show'); grid on;
figure(2); legend('show', 'Location', 'Best'); grid on;
%% VALIDATION: ANALYTICAL SOLUTION vs CFD

y_feet = linspace(0, h, n); 


y_dist = h - y_feet; 


nu = mu / rho;
k_stokes = sqrt(omega / (2 * nu));


t_capture = 1.25 * Period;
u_analytical = U_inf .* exp(-k_stokes .* y_dist) .* ...
               sin(omega * t_capture - k_stokes .* y_dist);

% 5. Plot Comparison
figure('Name', 'Validation', 'Color', 'black');
hold on;

% Plot CFD Result 
plot(u_peak, y_feet * 12, 'r-', 'LineWidth', 2, 'DisplayName', 'CFD Simulation'); 

% Plot Analytical Result
plot(u_analytical, y_feet * 12, 'y--', 'LineWidth', 1.5, 'MarkerSize', 6, ...
    'DisplayName', 'Analytical Solution');

title(' CFD vs. Analytical Solution');
xlabel('Velocity (ft/s)');
ylabel('Gap Height (in)'); 
legend('Location', 'Best');
grid on;

% 6. Calculate Error
error_val = u_peak - u_analytical;
rmse = sqrt(mean(error_val.^2));
fprintf('Root Mean Square Error (RMSE): %.4f ft/s\n', rmse);

%% ERROR DISTRIBUTION ANALYSIS

error_dist = u_peak - u_analytical;
figure('Name', 'Error Distribution', 'Color', 'black');
plot(error_dist, y_feet * 12, 'y-', 'LineWidth', 1.5, 'MarkerSize', 5, ...
    'MarkerFaceColor', 'r');
xline(0, 'c-', 'Zero Error Line', 'LabelVerticalAlignment', 'bottom');

title('Velocity Error Distribution across the Gap');
xlabel('Velocity Error [ft/s]');
ylabel('Gap Height (in)');
grid on;
max_err = max(abs(error_dist));
fprintf('max error: %.4f ft/s\n', max_err);
%% SHEAR STRESS ANALYSIS t = 1.25 Periods

du_dy_num = gradient(u_peak, y_feet); 
tau_cfd = mu .* du_dy_num;



theta = omega * 1.25 * Period - k_stokes .* y_dist;
du_dy_exact = k_stokes * U_inf .* exp(-k_stokes .* y_dist) .* (sin(theta) + cos(theta));

tau_analytical = mu .* du_dy_exact;


figure('Name', 'Shear Stress Distribution', 'Color', 'black');
ax = axes;
set(ax, 'Color', 'black', 'XColor', 'white', 'YColor', 'white');
hold on;

% Plot CFD Shear
plot(tau_cfd, y_feet * 12, 'r-', 'LineWidth', 2, 'DisplayName', 'CFD Shear'); 

% Plot Analytical Shear
plot(tau_analytical, y_feet * 12, 'w--', 'LineWidth', 2, 'DisplayName', 'Analytical Shear');

title('Shear Stress Distribution (\tau_w)', 'Color', 'white');
xlabel('Shear Stress (lb/ft^2)', 'Color', 'white');
ylabel('Channel Height (inches)', 'Color', 'white');
legend('Location', 'Best', 'TextColor', 'white', 'EdgeColor', 'white', 'Color', 'black');
grid on;
ax.GridColor = [0.3, 0.3, 0.3]; 
ax.GridAlpha = 0.5;


tau_wall_cfd = tau_cfd(end);          
tau_wall_exact = tau_analytical(end); 

error_tau = abs(tau_wall_cfd - tau_wall_exact);
pct_error_tau = (error_tau / abs(tau_wall_exact)) * 100;

fprintf('SHEAR STRESS ANALYSIS:\n');
fprintf('Wall Shear (CFD):       %.4f lb/ft^2\n', tau_wall_cfd);
fprintf('Wall Shear (Analytical): %.4f lb/ft^2\n', tau_wall_exact);
fprintf('Error at Wall:           %.2f%%\n', pct_error_tau);
