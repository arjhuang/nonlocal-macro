%%LAX Friedrichs
% Alex Davies
% Dr. Powers
clear all;
clc;
close all;
%Scaling = 20;
%% Import Density Data
density_data = readmatrix("NGSIM_US101_Density_Data.csv");

%% Initializing Parameters
t=1;
% only used as a fill in for plot command,
% alternates between 1 and 2 depending on i
dx=20; % 20 feet
% 0.001 works well for 1/1000,
% 0.005 works well for most larger values
dt=5; % 5 sec
ttot=2695; % Total time 2700 sec

pi=3.14159265;
N=ttot/dt; %Total number of time grid points


x_left= 0; % Left position
x_right= 2060; % Right position
xtot = x_right - x_left;
K= xtot/dx;

x=linspace(x_left,x_right,K+1); %Defining a uniform x
ti=linspace(0,ttot,N+1);

%%%%%%%%%%%%% defining the kernel
right_end = 100;  %%% it is the d such that eta is supported on (0,d).
etatot = right_end/dx;
eta_n(1:etatot+1) = 1/right_end; %%  constant kernel

% eta_n = zeros(1,etatot+1);
% for eta_index = 1: etatot+1
%    % eta_n(eta_index) = 1.5*((right_end)^2 - ((eta_index-1)*dx)^2)/((right_end)^3);  %%% quad kernel
%    eta_n(eta_index) = 2*(right_end - (eta_index-1)*dx)/(right_end^2);  %%%%%%%%%%%%% linear kernel
% end

x_collar= x_right - right_end;
xtot_coll = x_collar - x_left;
xk=xtot_coll/dx; %Total number of space grid points of the kernel

%%%%%%% maximum velocity and maximum density
vm = 80;
um = 0.12;




%% Creating more data points across time
time_scale = 50;
dt = dt/time_scale;
N_b = (N)*time_scale+1;
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
for j=xk+1:K+1
for i=1:N_b-1
    k1= rem(i-1,time_scale);
    k2= (i-1-k1)/time_scale;
    if k1==0
        u1(j,i) = density_data(j,k2+1);
    else
        u1(j,i) = (density_data(j,k2+1)*(time_scale-k1) + (k1)*density_data(j,k2+2))/time_scale;
    end
   
end
u1(j,N_b) = density_data(j,k2+2);
end






for i=1:N_b-1 %Time loop

%%%%% finding the u_d
u_d = u1(:,i);
for ud_index=1:xk+1
  vec_1 = u1(ud_index:ud_index+etatot-1,i)'; %%%% vec_1 = u(x) over the interval [x, x+d]
  vec_2 = eta_n(1:etatot);   %%%% eta function 
 u_d(ud_index) = sum(vec_1.*vec_2)*dx;  %%% convolution
end

%%%%%% finding flux
f = vm*u1(:,i).*(1 - u_d./um);



for j=2:xk % Space loop

u1(j,i+1) = (1/2)*(u1(j-1,i)+u1(j+1,i))- (dt/(2*dx))*(f(j+1) - f(j-1));


end
25
end

index = 1:1:540;
scaled_up_index = (index-1)*time_scale+1;
compressed_u = ones(K+1,N+1);
for j=1:K+1
   compressed_u(j,:) = u1(j,scaled_up_index);
end

%%% error computation,
U = compressed_u(1:xk+1,1:N+1);
D = density_data(1:xk+1, 1:N+1);
S =0;
S_e = 0;
for i=1:N+1
    for j=1:xk+1
       S = S+ (U(j,i) - D(j,i))^2;
       S_e = S_e + (D(j,i))^2;
    end
end
error = 100*S/S_e

% %% Ploting Density
figure;
contourf(ti,x,compressed_u,64,'LineColor', 'none');
%caxis([0 2]);
colorbar;
set(gca,'FontSize',20,'FontName','CMU Serif');
set(gca, 'TickLabelInterpreter','latex');
% Create xlabel
xlabel({'x (feet)'});    
% Create ylabel
ylabel({'time (s)'});

%% Ploting Real Density
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


% csvwrite('re_non_line_100_updated.csv', compressed_u)