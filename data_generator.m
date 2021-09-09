%% Define fuse shape [m]
len = 0.06; % 6cm
hei = 0.03; % 3cm
r = 0.01; % 1cm   len/2-r>=0 and hei/2-r>=0

%% Define True deviation
L = 0.1; % end-effector length 10cm
rx = -0.0015; % true x % |x|<= len/2 [1.5mm]
ry = -0.0025 - L; % true y % |y+L|<= hei/2 [2.5mm]
phi = 0.63/180*pi; % true pose deviation [2deg]

%% Param Setting
S = 100; % contact data points
% Force
Fy_mean = 5; % [N]
Fy_sigma = 0.1; % control precision [N]
mu = 0.5;
Fy = Fy_mean*ones(1,S) + Fy_sigma*randn(1,S); % true contact D_F_y
Fx = mu/3*randn(1,S).*Fy; % true contact D_F_x  i.e. static friction
% Grip pose
phi0 = 5;  % assume phi in [-phi0,phi0] [5deg]
alpha = linspace(phi0,90-phi0,S)/180*pi; % in [phi0, 90deg-phi0]
% Sensor Precision
% Noise Free Resolution (100Hz) 0.15 N, 0.15 N, 0.005 Nm, 0.002 Nm
sigma_F = 0.15/6; % standard deviation for F sensor noise [N]
sigma_M = 0.005/6; % standard deviation for M sensor noise [Nm]
% Mesh
N = 200; % mesh points for each corner. 4N in total % Mesh 越精细，误差越小
%% Generate Data (only left-bottom contact)

% contact point in fuse frame
B_r_BP = [-(len/2-r)-r*sin(phi+alpha);
          -(hei/2-r)-r*cos(phi+alpha)];
% Sensor Force
true_F = [Fy.*sin(alpha) + Fx.*cos(alpha);
          Fy.*cos(alpha) - Fx.*sin(alpha)];     
noise_F = sigma_F*randn(2,S);
sensor_F = true_F + noise_F;
% Sensor Moment
true_M = ones(1,S);
A_r_AP = ones(2,S);
D_r_AP = ones(2,S);
% A_M = A_r_AP × A_F
% A_r_AP = A_r_AB + C_AB * B_r_BP
A_r_AB = [rx;ry];
for i = 1:S
    C_AB = [cos(phi), -sin(phi);
            sin(phi), cos(phi)];
    A_r_AP(:,i) = A_r_AB + C_AB*B_r_BP(:,i);
    C_DA = [cos(alpha(i)), -sin(alpha(i));
            sin(alpha(i)), cos(alpha(i))];
    D_r_AP(:,i) = C_DA*A_r_AP(:,i);
    A_M = cross([A_r_AP(:,i);0],[true_F(:,i);0]);
    true_M(i) = A_M(3);
end
noise_M = sigma_M*randn(1,S);
sensor_M = true_M + noise_M;
% Gripper position D_r_DA
D_r_DA = zeros(2,S);
D_r_DA(2,:)=-D_r_AP(2,:);

% given data: sensor_F sensor_M C_DA(alpha) D_r_DA

%% fuse Mesh

theta = linspace(0,1,N)*pi/2; 
mesh_p = [(len/2-r)+r*cos(theta);
          (hei/2-r)+r*sin(theta)];
Mesh = [mesh_p(1,:),-mesh_p(1,:),-mesh_p(1,:),mesh_p(1,:);
        mesh_p(2,:),mesh_p(2,:),-mesh_p(2,:),-mesh_p(2,:)];
scatter(Mesh(1,:),Mesh(2,:)); % fuse mesh viz
axis equal
