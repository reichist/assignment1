clear
clc
% Variables    
m0 = 9.11e-31;
mn = 0.26*m0;
width = 200e-9;
height = 100e-9;
k = 1.381e-23;

% 1 - Electron Modelling
% 1.1 - What is the thermal velocity v_th? Assume T = 300K. 
T = 300;
v_th = sqrt((k * T )/ mn);
fprintf('The thermal velocity is %3.3d\n', v_th);
    
% 1.2 - If the mean time between collisions is tau_mn = 0.2ps, what is the
% mean free path?
tau_mn = 0.2e-12;
lambda = v_th * tau_mn;
fprintf('The original mean time between collisions is %3.3d\n', tau_mn);
fprintf('The original mean free path is %3.3d\n', lambda);

% 1.3 - Write a program that will model the random motion of the electrons.
% Variables    
m0 = 9.11e-31;
mn = 0.26*m0;
width = 200e-9;
height = 100e-9;
k = 1.381e-23;
n = 10000;
% n = 250;
T = 300;
time_interval = 1e-14;
% steps = 1000;
steps = 200;

% Electron Simulation
% Vector Setup
electrons_x = rand(1, n)*width;
electrons_y = rand(1, n)*height;
electrons_vx = (v_th/sqrt(2)).*randn(1, n);
electrons_vy = (v_th/sqrt(2)).*randn(1, n);

new_electrons_x = zeros(1, n);
new_electrons_y = zeros(1, n);
old_z = 0;
old_temperature = 300;
total_temperature = 0;
old_total_temperature = 0;

% Scattering Setup
p_scat = 1-exp(-time_interval/tau_mn);

% Plotting
figure(1)
clf

for z = 1:steps
    % Check for random scattering
    a=rand(1, n);
    electrons_vx(a<p_scat) = (v_th/sqrt(2))*randn(1, length(electrons_vx(a<p_scat))); 
    electrons_vy(a<p_scat) = (v_th/sqrt(2))*randn(1, length(electrons_vx(a<p_scat))); 
        
    % New X&Y position calculations
    new_electrons_x = electrons_x + time_interval*electrons_vx;
    new_electrons_y = electrons_y + time_interval*electrons_vy;
    
    % Check for BCs
    index = new_electrons_x>width;
    new_electrons_x(index) = new_electrons_x(index) - width;
    electrons_x(index) = electrons_x(index) - width;
    
    index = new_electrons_x<0;
    new_electrons_x(index) = new_electrons_x(index) + width;
    electrons_x(index) = electrons_x(index) + width;
    
    index = new_electrons_y>height;
    electrons_vy(index) = -electrons_vy(index);
    
    index = new_electrons_y<0;
    electrons_vy(index) = -electrons_vy(index);
    
    % Drift Velocity in both directions
    V(1, :) = sqrt(electrons_vx(1, :).^2 + electrons_vy(1, :).^2);
    V_mean = mean(V.^2);
    temperature = V_mean*mn/k;
    total_temperature = total_temperature + temperature;
    
    subplot(2, 1, 2)
    plot([old_z z], [old_temperature temperature], 'r');
    hold on; 
    plot([old_z z], [(old_total_temperature/(z-1)) (total_temperature/z)], 'b');
    legend('Temperature/Step', strcat('Average Overall Temperature =', num2str(round(total_temperature/z)), 'K'));
    title('Temperature');
    xlabel('Time Step'); ylabel('Temperature (K)');
    xlim([1 steps]);
    hold on;
    
    % Plot all electrons
    subplot(2, 1, 1)
    plot([electrons_x(1) new_electrons_x(1)], [electrons_y(1) new_electrons_y(1)], 'b');
    hold on;
    plot([electrons_x(2) new_electrons_x(2)], [electrons_y(2) new_electrons_y(2)], 'g');
    hold on;
    plot([electrons_x(3) new_electrons_x(3)], [electrons_y(3) new_electrons_y(3)], 'r');
    hold on;
    plot([electrons_x(4) new_electrons_x(4)], [electrons_y(4) new_electrons_y(4)], 'c');
    hold on;
    plot([electrons_x(5) new_electrons_x(5)], [electrons_y(5) new_electrons_y(5)], 'm');
    hold on;
    plot([electrons_x(6) new_electrons_x(6)], [electrons_y(6) new_electrons_y(6)], 'k');
    hold on;
    title('Electron Modelling');
    xlabel('width'); ylabel('height');
    grid on;
    xlim([0 200e-9]);
    ylim([0 100e-9]);
       
    pause(0.1);
    
    % Update electron coordinates
    electrons_x = new_electrons_x;
    electrons_y = new_electrons_y;
    old_z = z;
    old_temperature = temperature;
    old_total_temperature = total_temperature;
end

% Plot velocity histogram
figure(2)
histogram(V);
hold on;
y = get(gca, 'ylim');
x = v_th;
plot([v_th v_th], y);
hold off;
title('Maxwell-Boltzmann Distribution of Thermal Velocities');
legend('Probability', 'Thermal Velocity');
xlabel('Velocity (m/s)'); ylabel('Probability');
grid on;

fprintf('The calculated mean time between collisions is %3.3d\n', tau_mn);
fprintf('The calculated mean free path is %3.3d\n', lambda);