%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT SIGNAL EXAMPLE, nonlinear filtering and event detection
% as displayed in Fig.8
% (choose from the command window (k,V)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

addpath('include/');
addpath('data_Fig8/');

%% Screen display
msg1='Plot signal example, nonlinear filtering and event detection';

disp('-------------------------------------------------------------------')
disp(msg1)
disp('-------------------------------------------------------------------')

%% Parameters
m = 30.7e-3;  % et pas 32.7e-3, confirmé par Cristobal le 18/07/2018 !
g = 9.81;
fa = 2e3;
th_noise = 40;
th_event = th_noise*10;

%% Download data
kstiff      = input('Enter the value of k [N/m] (168/1002/2254) --> ') ;
vitesse  = input('Enter the value of V [microns/s] (42/1100/4300) --> ') ;
tit     = ['k',int2str(kstiff),'v',int2str(vitesse),'_lambda0.8.mat'];
load(tit);
%

%% Display t_start and t_stop / Extract tau_m and tau_s for L2L1
Fnorm=x1; % normalized force applied on the slider
vp = -(fa*m*g*1e6/kstiff*(diff(Fnorm))-vitesse); % slider velocity
t = linspace(0,length(Fnorm)-1,length(Fnorm))./fa; % time
%
[start,stop,taus,taum] = detect_tstartstop(Fnorm,g,vitesse,kstiff,m,fa,th_noise);

% ------------------------------------
% Figure size & properties
% ------------------------------------
ha=figure(1); 
set(gcf,'unit','centimeter','Position',[1 7 20 12]);
set(gca,'Fontname','Times','FontSize', 14);
%
hold on; grid on; box on
h11=plot(t,data,'Color',[0.0417 0 0],'DisplayName','data');
h12=plot(t,Fnorm,'Color',[1 0.45 0],'linewidth',1.5,'DisplayName','optimal \lambda');
% 
if kstiff==2254 && vitesse==42
    % No tstart & tstop for continuous sliding
else
    h13=plot(t(start),Fnorm(start),'ok','MarkerSize',7,'MarkerFaceColor',[1 1 1],'DisplayName','t_{start}');
    h14=plot(t(stop),Fnorm(stop),'ok','MarkerSize',7,'MarkerFaceColor',[.7 .7 .7],'DisplayName','t_{stop}');
end
TIT=strcat('$k$=',num2str(kstiff),' N/m - $V$=',num2str(vitesse),' $\mu$m/s');
title(TIT,'interpreter','Latex','Fontname','Times','FontSize',14)
xlabel('$t$ [s]','interpreter','Latex','Fontname','Times','FontSize', 16)
ylabel('$F/mg$ [-]','interpreter','Latex','Fontname','Times','FontSize', 16)
%
if kstiff ==2254 && vitesse==42
    hh=legend([h11,h12],{'data','NL filtering ($\lambda=0.8$)'},'Location','NorthEast');
    set(hh,'interpreter','latex','Fontname','Times','Fontsize',11)  % pour l'interpréter en Latex
else
    hh=legend([h11,h12,h13,h14],{'data','NL filtering ($\lambda=0.8$)','$t_{start}$','$t_{stop}$'},'Location','NorthEast');
    set(hh,'interpreter','latex','Fontname','Times','Fontsize',11)  % pour l'interpréter en Latex
end
%
V=vitesse;
% SET Xlim, Ylim of each signal
if kstiff==168 && V==42
    xlim([70,80]); ylim([.2,.3])
elseif kstiff==1002 && V==42
    xlim([70,80]); ylim([.15,.4])
elseif kstiff==2254 && V==42
    xlim([70,80]); ylim([.1,.6])
elseif kstiff==168 && V==1100
    xlim([5,5.6]); ylim([0.2,.45])
elseif kstiff==1002 && V==1100
    xlim([5,5.6]); ylim([.27,.42])
elseif kstiff==2254 && V==1100    
    xlim([5,5.6]); ylim([0.15,.45])
elseif kstiff==168 && V==4300    
    xlim([0,0.5]); ylim([0.4,.6])
elseif kstiff==1002 && V==4300    
    xlim([0,0.5]); ylim([0.15,.45])
elseif kstiff==2254 && V==4300    
    xlim([0,0.5]); ylim([0.15,.5])
end


