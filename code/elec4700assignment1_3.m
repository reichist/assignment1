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
fprintf('The mean free path is %3.3d\n', lambda);

% 1.3 - Write a program that will model the random motion of the electrons.
% Variables    
m0 = 9.11e-31;
mn = 0.26*m0;
width = 200e-9;
height = 100e-9;
k = 1.381e-23;
% n = 10;
% n = 100;
n = 10000;
T = 300;
time_interval = 5e-15;
% steps = 1000;
steps = 200;
% steps = 5;
grid_separation = 20e-9;

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

% Box Setup
 % cell 1 is left side of box
 % cell 2 is right side of box
 % cell 3 is bottom side of box
 % cell 4 is top side of box
box1coords = [80e-9 120e-9 0 40e-9];
box2coords = [80e-9 120e-9 60e-9 100e-9];

% Check for Box spawns
 % box 1
index1 = electrons_x > box1coords(1);
index2 = electrons_x < box1coords(2);
index3 = bitand(index1, index2);
index4 = electrons_y < box1coords(4);
index5 = bitand(index3, index4);
electrons_x(index5) = electrons_x(index5) - box1coords(1);

 % box 2
index1 = electrons_x > box2coords(1);
index2 = electrons_x < box2coords(2);
index3 = bitand(index1, index2);
index4 = electrons_y > box1coords(3);
index5 = bitand(index3, index4);
electrons_x(index5) = electrons_x(index5) + box2coords(1);

% Temperature Distribution Setup
temperature_matrix = zeros(20, 10);

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
    
    % Check for box BCs
     % Box 1
      % left side
    index1 = electrons_y < box1coords(4);         % check for electron below the top of box 1
    index2 = electrons_y > box1coords(3);         % check for electron above the bottom of box 1
    index3 = electrons_x < box1coords(1);         % check for old electron location on the left side of the left side of box 1
    index4 = new_electrons_x > box1coords(1);     % check for new electron location on the right side of the left side of box 1
    index5 = bitand(index1, index2);              % compare for y coordinate -> box1coords(3) < electrons_y < box1coords(4)
    index6 = bitand(index3, index4);              % compare for x coordinate -> electrons_x < box1coords(1) < new_electrons_x
    index7 = bitand(index5, index6);              % compare for qualifying electrons 
    electrons_vx(index7) = -electrons_vx(index7);
    
      % right side
    index1 = electrons_y < box1coords(4);         % check for electron below the top of box 1
    index2 = electrons_y > box1coords(3);         % check for electron above the bottom of box 1
    index3 = electrons_x > box1coords(2);         % check for old electron location on the right side of the right side of box 1
    index4 = new_electrons_x < box1coords(2);     % check for new electron location on the left side of the right side of box 1
    index5 = bitand(index1, index2);              % compare for y coordinate -> box1coords(3) < electrons_y < box1coords(4)
    index6 = bitand(index3, index4);              % compare for x coordinate -> electrons_x > box1coords(2) > new_electrons_x
    index7 = bitand(index5, index6);              % compare for qualifying electrons 
    electrons_vx(index7) = -electrons_vx(index7);
    
      % top side
    index1 = electrons_x > box1coords(1);         % check for electron to the right of the left side of box 1
    index2 = electrons_x < box1coords(2);         % check for electron to the left of the right side of box 1
    index3 = electrons_y > box1coords(4);         % check for old electron location above the top side of box 1
    index4 = new_electrons_y < box1coords(4);     % check for new electron location below the top side of box 1
    index5 = bitand(index1, index2);              % compare for y coordinate -> box1coords(1) < electrons_x < box1coords(2)
    index6 = bitand(index3, index4);              % compare for x coordinate -> electrons_y > box1coords(4) > new_electrons_y
    index7 = bitand(index5, index6);              % compare for qualifying electrons 
    electrons_vy(index7) = -electrons_vy(index7);

     % Box 2
      % left side
    index1 = electrons_y < box2coords(4);         % check for electron below the top of box 2
    index2 = electrons_y > box2coords(3);         % check for electron above the bottom of box 2
    index3 = electrons_x < box2coords(1);         % check for old electron location on the left side of the left side of box 2
    index4 = new_electrons_x > box2coords(1);     % check for new electron location on the right side of the left side of box 2
    index5 = bitand(index1, index2);              % compare for y coordinate -> box2coords(3) < electrons_y < box2coords(4)
    index6 = bitand(index3, index4);              % compare for x coordinate -> electrons_x < box2coords(1) < new_electrons_x
    index7 = bitand(index5, index6);              % compare for qualifying electrons 
    electrons_vx(index7) = -electrons_vx(index7);
    
      % right side
    index1 = electrons_y < box2coords(4);         % check for electron below the top of box 2
    index2 = electrons_y > box2coords(3);         % check for electron above the bottom of box 2
    index3 = electrons_x > box2coords(2);         % check for old electron location on the right side of the right side of box 2
    index4 = new_electrons_x < box2coords(2);     % check for new electron location on the left side of the right side of box 2
    index5 = bitand(index1, index2);              % compare for y coordinate -> box2coords(3) < electrons_y < box2coords(4)
    index6 = bitand(index3, index4);              % compare for x coordinate -> electrons_x > box2coords(2) > new_electrons_x
    index7 = bitand(index5, index6);              % compare for qualifying electrons 
    electrons_vx(index7) = -electrons_vx(index7);
    
      % bottom side
    index1 = electrons_x > box2coords(1);         % check for electron to the right of the left side of box 2
    index2 = electrons_x < box2coords(2);         % check for electron to the left of the right side of box 2
    index3 = electrons_y < box2coords(3);         % check for old electron location below the bottom side of box 2
    index4 = new_electrons_y > box2coords(3);     % check for new electron location above the bottom side of box 2
    index5 = bitand(index1, index2);              % compare for y coordinate -> box2coords(1) < electrons_x < box2coords(2)
    index6 = bitand(index3, index4);              % compare for x coordinate -> electrons_y > box2coords(4) > new_electrons_y
    index7 = bitand(index5, index6);              % compare for qualifying electrons 
    electrons_vy(index7) = -electrons_vy(index7);   
    
    % Drift Velocity to Temperature
    V(1, :) = sqrt(electrons_vx(1, :).^2 + electrons_vy(1, :).^2);
    V_mean = mean(V.^2);
    temperature = V_mean*mn/k;
    total_temperature = total_temperature + temperature;
    
%     % Plot velocity histogram
%     subplot(2, 3, 3)
%     histogram(V);
%     hold on;
%     y = get(gca, 'ylim');
%     x = v_th;
%     plot([v_th v_th], y);
%     hold off;
%     title('Maxwell-Boltzmann Distribution of Thermal Velocities');
%     legend('Probability', 'Thermal Velocity');
%     xlabel('Velocity (m/s)'); ylabel('Probability');
%     grid on;
    
    % Plot Temperature Lines
%     subplot(2, 3, 2)
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
%     subplot(2, 3, 1)
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
    
    % Boxes Setup
    box1 = rectangle('Position', [80e-9 0     40e-9 40e-9]); % Creates Box 1 from (80nm, 0nm) to (120nm, 40nm)
    box2 = rectangle('Position', [80e-9 60e-9 40e-9 40e-9]); % Creates Box 2 from (80nm, 60nm) to (120nm, 100nm)
    
    title('Electron Modelling');
    xlabel('Width'); ylabel('Height');
    grid on;
    xlim([0 200e-9]);
    ylim([0 100e-9]);
    
%     % Electron Density Map
%     subplot(2, 3, 4)
%     hist3([electrons_x', electrons_y'], 'CDataMode', 'auto', 'FaceColor', 'interp')
%     colormap('default');
%     colorbar;
%     xlabel('Width'); ylabel('Height');
%     title('Electron Density');
%     view(2);
    
%     % Temperature Distribution Map
%     for y_axis = 1:10
%         y_upper = y_axis*10;
%         y_lower = y_upper-10;
%         for x_axis = 1:20
%             x_upper = x_axis*10;
%             x_lower = x_upper-10;
%             index1 = electrons_x > (x_lower*1e-9);
%             index2 = electrons_x < (x_upper*1e-9);
%             index3 = electrons_y > (y_lower*1e-9);
%             index4 = electrons_y < (y_upper*1e-9);
%             index5 = bitand(index1, index2);
%             index6 = bitand(index3, index4);
%             index7 = bitand(index5, index6);
%             velocity = (electrons_vx(index7).^2) + (electrons_vy(index7).^2);
%             v_mean = mean(velocity);
%             temperature_value = (((v_mean)*(mn))/(k));
%             temperature_matrix(x_axis, y_axis) = temperature_value;
%             temperature_matrix(8:12, 1:4) = NaN;
%             temperature_matrix(8:12, 7:10) = NaN;
%         end       
%     end
% 
%      % Plotting Temperature Distribution
%     subplot(2, 3, 5)
%     surf(transpose(temperature_matrix));
%     xlabel('Width (x*10nm)'); ylabel('Height(x*10nm)');
%     title('Temperature Distribution');
%     colormap('default');
%     colorbar;
%     view(2);
    
    pause(0.1);
    
    % Update electron coordinates
    electrons_x = new_electrons_x;
    electrons_y = new_electrons_y;
    old_z = z;
    old_temperature = temperature;
    old_total_temperature = total_temperature;
end

figure(2)
% Electron Density Map
hist3([electrons_x', electrons_y'], 'CDataMode', 'auto', 'FaceColor', 'interp')
colormap('default');
colorbar;
xlabel('Width'); ylabel('Height');
title('Electron Density');
view(2);

figure(3)
 % Plotting Temperature Distribution
surf(transpose(temperature_matrix));
xlabel('Width (x*10nm)'); ylabel('Height(x*10nm)');
title('Temperature Distribution');
colormap('default');
colorbar;
view(2);