% Morris Lecar model
% simple model of excitable media that exhibits true threshold behavior
% Single cell

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

tmax = 50; %s
% dt = 0.0001; 
dt = 0.01;
durn = round(1/dt);
tnmax = round(tmax/dt);

resultt = [];
resultu = [];
resultw = [];
resulttau = [];

u = -0.283;
% w = -0.03;
nch = 1000;
w = zeros(1,nch);

iappl=zeros(1,tnmax);
iappl(1100:1200)=0.065;

for tn = 1:tnmax
%     tn/tnmax
    m_ss = 0.5*(1+tanh((u-v1)/v2));
    w_ss = 0.5*(1+tanh((u-v3)/v4));
    tau_w = a*1./cosh((u-v3)/(2*v4));
    
%     dw_dt = 1/tau_w*(w_ss - w);
%     w = w + dw_dt*dt;
    
     alpha_w = w_ss/tau_w;
     beta_w = (1-w_ss)/tau_w;
     
         r=rand(1,nch);  
         w(w==0 & r < alpha_w*dt)=1;

         r=rand(1,nch);
         w(w==1 & r < beta_w*dt)=0;

         w_avg = mean(w');
         w_avg = w_avg';


        du_dt = I-g_Ca*m_ss.*(u-E_Ca) - g_k*w_avg.*(u-E_K) - g_L*(u-E_L)+iappl(tn);    
%     du_dt = I-g_Ca*m_ss*(u-E_Ca) - g_k*w*(u-E_K) - g_L*(u-E_L);
    u = u + du_dt * dt;
    
    if (mod(tn,10) == 0)
        resultt = [resultt tn*dt];
        resultu = [resultu u];
%         resultw = [resultw w];
        resultw = [resultw w_avg];
    end   
end

% u_t = -1 : 0.01: 0.2;
% w_ssnull = 0.5*(1+tanh((u_t-v2)/v4));
% m_ssnull = 0.5*(1+tanh((u_t-v1)/v2));
% u_ssnull = (I - g_Ca*m_ssnull.*(u_t-E_Ca) - g_L*(u_t-E_L))./(g_k*(u_t-E_K));

figure(1);
plot(resultu,resultw, 'r')
% xlabel('u(t)');
% ylabel('w(t)');
hold on;
% plot(u_t,w_ssnull,'b--')
% plot(u_t,u_ssnull,'k--')
axis([-0.4 0.5 -0.05 0.45]);
% title('w(t) vs u(t) phase plane');
%hold off;

% figure(2);
% plot(resultt,resultu,'b')
% xlabel('time');
% ylabel('u(t)');
% 
% hold on;
% plot(resultt,resultw,'g')
% xlabel('time');
% ylabel('w(t)& u(t)');
% legend('u(t)','w(t)','Location','Northeast');
%hold off
% 
figure
subplot(2,1,1)
plot(resultt,resultu,'b')
% axis([0 50 -0.3 0.5])
xlabel('time')
ylabel('u(t)')

iappl=zeros(1,tnmax/10);
iappl(110:120)=0.065;

subplot(2,1,2)
plot(resultt,iappl,'r')
axis([0,50,0,0.07])
xlabel('time');
ylabel('w(t)');
    
    