% given data: sensor_F sensor_M C_DA(alpha) D_r_DA Mesh

%% Param Config
sigma_x = 0.0020; % [m] 2mm
sigma_y = 0.0020; % [m] 2mm
mean_phi = 0;       % deg
sigma_phi = 1;    % deg
Np = 6000; % particles
% contact point standard deviation
s1 = 0.000008; %[m]
s2 = 0.000008;

%% init particles
Rx = sigma_x*randn(1,Np);
Ry = sigma_y*randn(1,Np);
Phi = (mean_phi + sigma_phi*randn(1,Np))/180*pi;
if max(Phi)>phi0/180*pi || min(Phi)<-phi0/180*pi
    disp('Warning, contact assumption break')
end

subplot(1,2,1)
scatter(Rx,Ry,15,'filled')
axis equal
xlim([-0.01 0.01])
ylim([-0.01 0.01])
hold on

subplot(1,2,2)
histogram(Phi/pi*180,200);
hold on

%% Measurement Update
for i = 1:S
    C_da = [cos(alpha(i)), -sin(alpha(i));
            sin(alpha(i)), cos(alpha(i))];
    d_r_da = D_r_DA(:,i);
    fx = sensor_F(1,i);
    fy = sensor_F(2,i);
    moment = sensor_M(i);
    
    b = D_r_DA(2,i);
    x0 = (b*fy*sin(alpha(i))-b*fx*cos(alpha(i))+moment)/(fy*cos(alpha(i))+fx*sin(alpha(i)));
    p_zx = zeros(1,Np);
    
    for p = 1:Np
        a_r_ab = [Rx(p);-L+Ry(p)];
        % find contact point
        C_ab = [cos(Phi(p)), -sin(Phi(p));
                sin(Phi(p)), cos(Phi(p))];
        d_r_dp = (d_r_da + C_da*a_r_ab)*ones(1,4*N) + C_da*C_ab*Mesh;
        %X = [x0*ones(1,4*N); zeros(1,4*N)];
        %dist = sqrt(sum((d_r_dp'-X').^2,2));
        [min_y,min_ind] = min(d_r_dp(2,:));
        %[min_d,min_ind] = min(dist);

        p_zx(p) = mvnpdf([x0;0],d_r_dp(:,min_ind),diag([s1,s2]));
    end
    
    SUM = sum(p_zx);
    beta = p_zx/SUM;
    cumSUM = cumsum(beta);
    ind = zeros(1,Np);
    for p = 1:Np
        r = rand();
        indexs = find(cumSUM>=r);
        ind(p) = indexs(1);
    end
    
    Rx = Rx(ind);
    Ry = Ry(ind);
    Phi = Phi(ind);
    Ux = unique(Rx);
    numbers = length(Ux);
    
% %    roughening
%     if i==S
%         continue
%     end
%  
%     if numbers<100 && i<S-10
%         Kx = 0.001;
%         Ky = 0.0001;
%         Kphi = 0.00001;
%         Ex = max(Rx)-min(Rx);
%         Ey = max(Ry)-min(Ry);
%         Ep = max(Phi)-min(Phi);
%         sigma_rx = Kx * Ex * Np^(-1/3);
%         sigma_ry = Ky * Ey * Np^(-1/3);
%         sigma_phi = Kphi * Ep * Np^(-1/3);
%         Rx = Rx + sigma_x*randn(1,Np);
%         Ry = Ry + sigma_y*randn(1,Np);
%         Phi = Phi +sigma_phi*randn(1,Np);
%     end
  
    disp(numbers)
end

% result analysis
% position
result_x = mean(Rx);
result_y = mean(Ry);
subplot(1,2,1)
scatter(Rx,Ry,10,'r','filled')
scatter(rx,ry+L,50,'d','filled','g');
scatter(result_x,result_y,50,'filled','m');
disp(['estimated pos deviation: [',num2str(result_x),' , ',num2str(result_y),']']);
disp(['position absolute error: [',num2str((result_x-rx)*1000),'mm , ',num2str((result_y-ry-L)*1000),'mm]'])
disp(['position relative error: [',num2str((result_x-rx)/rx*100),'% ,',num2str((result_y-ry-L)/(ry+L)*100),'%]'])

% pose
result_phi = mean(Phi);

subplot(1,2,2)
ylim([0 500])
histogram(Phi/pi*180,400);
scatter(Phi/pi*180,zeros(1,Np),50,'r','filled');
scatter(phi/pi*180,0,50,'d','filled','g');
disp(['estimated phi deviation: ',num2str(result_phi/pi*180),' deg']);
disp(['phi absolute error: ',num2str((result_phi-phi)/pi*180),' deg'])
disp(['phi relative error: ',num2str((result_phi-phi)/phi*100),'%'])


