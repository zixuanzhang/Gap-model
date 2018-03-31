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
nch = 10; % change nch
tmax = 100; % ms
tnmax = round(tmax/dt);
resultw = [];
w = zeros(num,nch);
resultwavg = [];

resultpro = [];
resultdis = [];

%% main loop
itr = 1000;

for loop = 1:itr
    u = zeros(num,1); % membrane potential
    u(1:num) = -0.2824;
    u(1:50) = 1;
    resultu = [];  
    loop
    for tn = 0:tnmax
        m_ss = 0.5.*(1+tanh((u-v1)./v2));
        w_ss = 0.5.*(1+tanh((u-v3)./v4));
        tau_w = a.*1./cosh((u-v3)./(2*v4));

        % add stochastic to w gating
        alpha_w = w_ss./tau_w;
        beta_w = (1-w_ss)./tau_w;

         r=rand(num,nch);  
         w(w==0 & r < alpha_w*dt)=1;

         r=rand(num,nch);
         w(w==1 & r < beta_w*dt)=0;

         w_avg = mean(w');
         w_avg = w_avg';


        du_dt = I-g_Ca*m_ss.*(u-E_Ca) - g_k*w_avg.*(u-E_K) - g_L*(u-E_L);
        u = u + du_dt*dt.*s;
    
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
        end
    end

%     figure(1)
%     surf(resultu);
%     shading interp
%     % axis([0 inf 0 inf]);
%     axis off
    
    countpro = 0;
    countdis = 0;

        for ti = 1:800
            if resultu(5,ti) > 0.2 & resultu(5,ti+1) < 0.2
                countpro = countpro + 1;
            end

            if resultu(595,ti) > 0.2 & resultu(595,ti+1) < 0.2
                countdis = countdis + 1;
            end

        end
    
    resultpro = [resultpro countpro];
    resultdis = [resultdis countdis];
     
end
        
wavecount = [resultpro;resultdis]
% T = array2table(wavecount','VariableNames',{'pro','dis'})

% reflection pattern
count10 = 0;
count11 = 0;
count21 = 0;
count22 = 0;
count32 = 0;
count33 = 0;

% count waves
for n = 1:itr
    count10(wavecount(1,n)==1 & wavecount(2,n) == 0) = count10+1;
    count11(wavecount(1,n)==1 & wavecount(2,n) == 1) = count11+1;
    count21(wavecount(1,n)==2 & wavecount(2,n) == 1) = count21+1;
    count22(wavecount(1,n)==2 & wavecount(2,n) == 2) = count22+1;
    count32(wavecount(1,n)==3 & wavecount(2,n) == 2) = count32+1;
    count33(wavecount(1,n)==3 & wavecount(2,n) == 3) = count33+1;
end

showset= [count10,count11,count21,count22,count32,count33]
percentage = showset/itr

         
% view(2)
% figure(2)
% surf(resultw);
% shading interp
% % axis([0 inf 0 inf]);
% axis off
% view(2)
    
 

    
    
    
    