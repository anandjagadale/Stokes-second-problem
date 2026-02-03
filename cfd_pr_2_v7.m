clear all
close all
clc

%% PARAMETERS
freq = 4000;            
omega = 2 * pi * freq;   
Period = 1 / freq;       
U_inf = 200;             
h_inch = 0.01;           
h = h_inch / 12;        

rho = 0.00237;          
mu = 3.737e-7;          
nu = mu / rho;          

%% 
n = 51;                 
dy = h / (n - 1);       
y = linspace(0, h_inch, n); 

dt = 0.5 * (rho * dy^2) / mu * 0.9; 

%% SETUP
n_cycles = 4.2;
t_max = n_cycles * Period;
num_steps = ceil(t_max / dt);
u_spacetime = zeros(n, num_steps);
t_vector = zeros(1, num_steps);

% Initial Conditions
u = zeros(1, n);        
u(n) = 0; 
t_now = 0;

fprintf('Simulating Spacetime Evolution (%d cycles)...\n', n_cycles);

%% MAIN LOOP
for k = 1:num_steps
    u_old = u;
    
    % Diffusion Equation
    for i = 2:(n-1)
        diffusion = (nu) * (u_old(i+1) - 2*u_old(i) + u_old(i-1)) / dy^2;
        u(i) = u_old(i) + dt * diffusion;
    end
    
    % Boundary Conditions
    u(1) = 0; 
    u(n) = U_inf * sin(omega * (t_now + dt)); 
    
    u_spacetime(:, k) = u'; 
    t_vector(k) = t_now;
    
    t_now = t_now + dt;
end


%% PLOTTING: SPACETIME CONTOUR 
figure('Name', 'Spacetime Diagram', 'Color', 'black');
ax = axes;
set(ax, 'Color', 'black', 'XColor', 'white', 'YColor', 'white');
hold on;
cmap = jet(256);
mid_idx = 128; 
cmap(mid_idx-5:mid_idx+5, :) = 0; 
[C, h_contour] = contourf(t_vector/Period, y, u_spacetime, 100, 'LineColor', 'none');
colormap(jet); 
clim([-U_inf, U_inf]); 

c = colorbar;
c.Color = 'white';
c.Label.String = 'Velocity (ft/s)';
c.Label.Color = 'white';


title('Spacetime Diagram: Velocity Contour', 'Color', 'white');
xlabel('Time (Periods)', 'Color', 'white');
ylabel('Channel Height (inches)', 'Color', 'white');
ylim([0, h_inch]);
xlim([0, n_cycles]);
ax.Layer = 'top'; 
grid on;
ax.GridColor = [1, 1, 1]; 
ax.GridAlpha = 0.3;