%%LAX Friedrichs
% Alex Davies
% Dr. Powers
clear all;
clc;
close all;
%Scaling = 1; %%% probably we don't need, we need more time discretization
%% Import Density Data
density_data = readmatrix("NGSIM_US101_Density_Data.csv");
density_data = (density_data);
%% Initializing Parameters
t=0;
% only used as a fill in for plot command,
% alternates between 1 and 2 depending on i
dx=20; % 20 feet
% 0.001 works well for 1/1000,
% 0.005 works well for most larger values
dt=5; % 5 sec  we need dt to be smaller
ttot=2695; % Total time 2700 sec

pi=3.14159265;
N=ttot/dt; %Total number of time grid points

x_left= 0; % Left position
x_right= 2060; % Right position
xtot = x_right - x_left;
K= xtot/dx;

x=linspace(x_left,x_right,K+1); %Defining a uniform x
ti=linspace(0,ttot,N+1);
%%%%%%% maximum velocity and maximum density
vm = 80;
um = 0.12;


%% Creating more data points across time
time_scale = 50;
dt = dt/time_scale;
N_b = (N)*time_scale+1; %% new number of data points
u1 = um*ones(K+1,N_b);

%%%%%%%%%% initial condition
u1(:,1) = density_data(:,1);

%% left boundary condition
for i=1:N_b-1
    k1= rem(i-1,time_scale);
    k2= (i-1-k1)/time_scale;
    if k1==0
        u1(1,i) = density_data(1,k2+1);
    else
        u1(1,i) = (density_data(1,k2+1)*(time_scale-k1) + (k1)*density_data(1,k2+2))/time_scale;
    end
    
end
u1(1,N_b) = density_data(1,k2+2);

%% right boundary condition
for i=1:N_b-1
    k1= rem(i-1,time_scale);
    k2= (i-1-k1)/time_scale;
    if k1==0
        u1(K+1,i) = density_data(K+1,k2+1);
    else
        u1(K+1,i) = (density_data(K+1,k2+1)*(time_scale-k1) + (k1)*density_data(K+1,k2+2))/time_scale;
    end
   
end
u1(K+1,N_b) = density_data(K+1,k2+2);








for i=1:N_b-1 %Time loop

    %%%%%% finding flux
    u_d = u1(:,i);
f = vm*u_d.*(1 - u_d./um);



for j=2:K % Space loop

u1(j,i+1) = (1/2)*(u1(j-1,i)+u1(j+1,i))- (dt/(2*dx))*(f(j+1) - f(j-1));


end
%25
end

index = 1:1:540;
scaled_up_index = (index-1)*time_scale+1;
compressed_u = ones(K+1,N+1);
for j=1:K+1
   compressed_u(j,:) = u1(j,scaled_up_index);
end

%%% error computation
U = compressed_u(1:K+1,1:N+1);
D = density_data(1:K+1, 1:N+1);
S =0;
S_e = 0;
for i=1:N+1
    for j=1:K+1
       S = S+ (U(j,i) - D(j,i))^2;
       S_e = S_e + (D(j,i))^2;
    end
end
error = 100*S/S_e

%% Ploting Density
figure;
contourf(ti,x,compressed_u,64,'LineColor', 'none');
%caxis([0 2]);
colorbar;
set(gca,'FontSize',20,'FontName','CMU Serif');
set(gca, 'TickLabelInterpreter','latex');
% Create xlabel
ylabel({'x (feet)'});    
% Create ylabel
xlabel({'time (s)'});

%% Ploting Density
figure;
contourf(ti,x,density_data,64,'LineColor', 'none');
%caxis([0 2]);
colorbar;
set(gca,'FontSize',20,'FontName','CMU Serif');
set(gca, 'TickLabelInterpreter','latex');
% Create xlabel
ylabel({'x (feet)'});    
% Create ylabel
xlabel({'time (s)'});

csvwrite('re_local_updated.csv', compressed_u)
csvwrite('den_101_updated.csv', density_data)