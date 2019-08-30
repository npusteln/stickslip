%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXAMPLE 2 :
% nonlinear vs. linear denoising + t_start and t_stop detection 
% for an inertial regime [k = 168 N/m; V = 4300 microns/s].
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

addpath('include/');
addpath('data_raw/');

%% Download data
kstiff      = 168; % 1002, 2254
vitesse = 4300; % 1100, 4300
tit     = ['k',int2str(kstiff),'Nm_v',int2str(vitesse),'.mat'];
load(tit);
data    = Fnorm';
t = t;

%% Parameters 
m = 30.7e-3;  % mass of the slider [kg]
fa = 2e3; % acquisition frequency [Hz}
g = 9.81; % gravitational acceleration [m/s^2]
th_noise = 40; % noise threshold

%% Screen display
msg1='EXAMPLE 2';
msg2='nonlinear vs. linear denoising + t_start and t_stop detection'; 
msg3='for an INERTIAL REGIME [k = 168 N/m; V = 4300 microns/s].';

disp('-------------------------------------------------------------------')
disp(msg1)
disp(msg2)
disp(msg3)
disp('-------------------------------------------------------------------')
disp('Computing nonlinear & linear filtering...');

%% Algorithm L2L1: proposed method, non-linear filtering
param.filter      = 'laplacian'; % Type of filter in the regularization
param.computation = 'direct';    % Convolution type
param.lambda      = 0.8;           % Regularization parameter
param.iter      = 100000;        % Iteration number
param.normL     = 1;
param.mu        = 0;
op.direct       = @(x)opL_1D(x,param);
op.adjoint      = @(x)opLadj_1D(x,param);

param.mu=1;
prox.fidelity = @(y,data,tau) prox_L2(y-data,tau)+data;
objective.fidelity =  @(y,data) 1/2*sum(abs(y-data).^2);
prox.regularization = @(y,tau) prox_L1(y,tau);  
objective.regularization =  @(y,tau) tau*sum(abs(y)); 

[xl2l1,~]= PD_ChambollePock(data, param, op, prox, objective);
%titsave = ['k',int2str(kstiff),'v',int2str(vitesse),'lambda_',num2str(param.lambda),'l2l1.mat'];       
%save(titsave,'xl2l1','param','data');

%% Algorithm L2L2: linear filtering
param.lambda=10000;
H = psf2otf([1 -2 1],[1,length(data)]);
xl2l2 = real(ifft( fft(data)./(1+param.lambda*conj(H).*H)));
%titsave = ['k',int2str(kstiff),'v',int2str(vitesse),'lambda_',num2str(param.lambda),'l2l2.mat'];       
%save(titsave,'xl2l2','param','data');


%% Display t_start and t_stop / Extract tau_m and tau_s for L2L1
Fnorm = xl2l1; % normalized force applied on the slider
vp = -(fa*m*g*1e6/kstiff*(diff(Fnorm))-vitesse); % slider velocity
t = linspace(0,length(Fnorm)-1,length(Fnorm))./fa;
[start,stop,taus,taum] = detect_tstartstop(Fnorm,g,vitesse,kstiff,m,fa, th_noise);

% Display figure
figure(1); ax1 = subplot(121);
hold on; grid on; box on
h11=plot(t,data,'Color',[0.0417 0 0],'DisplayName','data');
h12=plot(t,Fnorm,'Color',[1 0.45 0],'linewidth',1.5,'DisplayName','optimal \lambda');
h13=plot(t(start),Fnorm(start),'ok','MarkerSize',7,'MarkerFaceColor',[1 1 1],'DisplayName','t_{start}');
h14=plot(t(stop),Fnorm(stop),'ok','MarkerSize',7,'MarkerFaceColor',[.7 .7 .7],'DisplayName','t_{stop}');
title 'Nonlinear denoising: L2L1'
xlabel('$t$ [s]','interpreter','Latex','Fontname','Times','FontSize', 16)
ylabel('$F/mg$ [-]','interpreter','Latex','Fontname','Times','FontSize', 16)
hh=legend([h11,h12,h13,h14],{'data','$\lambda=0.8$','$t_{start}$','$t_{stop}$'},'Location','NorthEast');
set(hh,'interpreter','latex','Fontname','Times','Fontsize',11)  % pour l'interpréter en Latex

clear Fnorm vp t start stop taus taum 

%% Display t_start and t_stop / Extract tau_m and tau_s for L2L2
Fnorm = xl2l2; % normalized force applied on the slider
vp = -(fa*m*g*1e6/kstiff*(diff(Fnorm))-vitesse); % slider velocity
t = linspace(0,length(Fnorm)-1,length(Fnorm))./fa;
[start,stop,taus,taum] = detect_tstartstop(Fnorm,g,vitesse,kstiff,m,fa, th_noise);

% Display figure
figure(1); ax2 = subplot(122);
hold on; grid on; box on
h21=plot(t,data,'Color',[0.0417 0 0],'DisplayName','data');
h22=plot(t,Fnorm,'Color',[1 0.45 0],'linewidth',1.5,'DisplayName','optimal \lambda');
h23=plot(t(start),Fnorm(start),'ok','MarkerSize',7,'MarkerFaceColor',[1 1 1],'DisplayName','t_{start}');
h24=plot(t(stop),Fnorm(stop),'ok','MarkerSize',7,'MarkerFaceColor',[.7 .7 .7],'DisplayName','t_{stop}');
title 'Linear denoising: L2L2'
xlabel('$t$ [s]','interpreter','Latex','Fontname','Times','FontSize', 16)
ylabel('$F/mg$ [-]','interpreter','Latex','Fontname','Times','FontSize', 16)
hh=legend([h21,h22,h23,h24],{'data','$\lambda=10^4$','$t_{start}$','$t_{stop}$'},'Location','NorthEast');
set(hh,'interpreter','latex','Fontname','Times','Fontsize',11)  % pour l'interpréter en Latex
%
linkaxes([ax1,ax2],'xy')

% ------------------------------------
% Figure size & properties
% ------------------------------------
set(gcf,'unit','centimeter','Position',[1 7 30 12]);
set(ax1,'Fontname','Times','FontSize', 14);
set(ax2,'Fontname','Times','FontSize', 14);

