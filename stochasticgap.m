% 2D cable for gap model

% initial values
dt = 0.01;
dx = 0.05;
num = 30/dx;
L = 1.6;
u = zeros(num,1); % membrane potential
u(1:num) = -0.2824;
u(1:50) = 1;

tmp = zeros(num,1);
w = zeros(num,1); % diffusion factor
w(1:5) = 0.004;
s = ones(num,1); % change to zero in the gap region 

gap = L/dx;
s(num/2-gap/2:num/2+gap/2) = 0;

% constants
I = 0.08;
g_Ca = 1.0;
E_Ca = 1.0;
g_k = 2.0;
E_K = -0.7;
g_L = 0.5;
E_L = -0.5;

v1 = -0.01;
v2 = 0.15;
v3 = 0.1;
v4 = 0.145;
a = 3.0;

% constants for PDE
D = 1;
nch = 100; % change nch
tmax = 100; % ms
tnmax = round(tmax/dt);
resultu = [];
resultw = [];
% w = zeros(num,1);
% temp = zeros(num,1);
w = zeros(num,nch);
resultwavg = [];

for tn = 0:tnmax
    tn/tnmax
    
    m_ss = 0.5.*(1+tanh((u-v1)./v2));
    w_ss = 0.5.*(1+tanh((u-v3)./v4));
    tau_w = a.*1./cosh((u-v3)./(2*v4));
    
    % add stochastic to w gating
    alpha_w = w_ss./tau_w;
    beta_w = (1-w_ss)./tau_w;
    
 
%     dw_dt = 1./tau_w.*(w_ss - w);
%     w = w + dw_dt.*dt;
             r=rand(num,nch);  
             w(w==0 & r < alpha_w*dt)=1;
             
             r=rand(num,nch);
             w(w==1 & r < beta_w*dt)=0;
                
             w_avg = mean(w');
             w_avg = w_avg';
            
                           
    du_dt = I-g_Ca*m_ss.*(u-E_Ca) - g_k*w_avg.*(u-E_K) - g_L*(u-E_L);
%     du_dt = I-g_Ca*m_ss.*(u-E_Ca) - g_k*w.*(u-E_K) - g_L*(u-E_L);
    u = u + du_dt*dt.*s;
 
   % use implict method     
    for t = 1:100
    % solve for diffusion equation
        u(1) = u(3);
        u(num) = u(num-2);

        for c = 2:num-1
            tmp(c) = u(c) + (u(c-1) + u(c+1) - 2*u(c))*D*dt/100./(dx*dx);
        end

        u = tmp;
    end
    
    if (mod(tn,10) == 0)
        resultu = [resultu u];
        resultw = [resultw w_avg];
%         resultw = [resultw w];
    end
end

% space time plot
figure(1)
surf(resultu);
shading interp
% axis([0 inf 0 inf]);
axis off
view(2)
figure(2)
surf(resultw);
shading interp
% axis([0 inf 0 inf]);
axis off
view(2)
    
    

    
    
    
    