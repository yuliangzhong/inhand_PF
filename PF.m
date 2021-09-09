% given data: sensor_F sensor_M C_DA(alpha) D_r_DA Mesh

%% Param Config
sigma_x = 0.0020; % [m] 2mm
sigma_y = 0.0020; % [m] 2mm
mean_phi = 0.5;       % deg
sigma_phi = 0.25;    % deg
Np = 6000; % particles
% suppose moment measurement noise 
sigma_w = 0.005/6;
% contact threshold
threshold = 0.00015; % [m] estimated by F/K


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
    p_zx = zeros(1,Np);
    for p = 1:Np
        a_r_ab = [Rx(p);-L+Ry(p)];
        % find contact point
        C_ab = [cos(Phi(p)), -sin(Phi(p));
                sin(Phi(p)), cos(Phi(p))];
        d_r_dp = (d_r_da + C_da*a_r_ab)*ones(1,4*N) + C_da*C_ab*Mesh;
        [min_y,min_ind] = min(d_r_dp(2,:));
        d_r_ap = d_r_dp(:,min_ind)-d_r_da;

        if abs(min_y)>threshold % hard contact
            Mo = zeros(3,1);
        else
            % A_M = a_r_ap ¡Á a_F
            a_r_ap = C_da\d_r_ap;
            Mo = cross([a_r_ap;0],[fx;fy;0]);
        end

        p_zx(p) = normpdf(moment,Mo(3),sigma_w);
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
    
    % roughening
    
    if i==S
        continue
    end

    Kx = 0.1;
    Ky = 0.2;
    Kphi = 0.18;
    Ex = max(Rx)-min(Rx);
    Ey = max(Ry)-min(Ry);
    Ep = max(Phi)-min(Phi);
    sigma_rx = Kx * Ex * Np^(-1/3);
    sigma_ry = Ky * Ey * Np^(-1/3);
    sigma_phi = Kphi * Ep * Np^(-1/3);
    Rx = Rx + sigma_x*randn(1,Np);
    Ry = Ry + sigma_y*randn(1,Np);
    Phi = Phi +sigma_phi*randn(1,Np);
  
    disp([num2str(i/S*100),'%'])
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
% Multimodal distribution!!
Uphi = unique(Phi);
count = hist(Phi,Uphi);
[num ,ind] = max(count);
result_phi = Uphi(ind);

subplot(1,2,2)
ylim([0 500])
histogram(Phi/pi*180,400);
scatter(Phi/pi*180,zeros(1,Np),50,'r','filled');
scatter(phi/pi*180,0,50,'d','filled','g');
disp(['estimated phi deviation: ',num2str(result_phi/pi*180),' deg']);
disp(['phi absolute error: ',num2str((result_phi-phi)/pi*180),' deg'])
disp(['phi relative error: ',num2str((result_phi-phi)/phi*100),'%'])


