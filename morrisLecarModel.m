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

tmax = 200; %s
dt = 0.001;
durn = round(1/dt);
tnmax = round(tmax/dt);

resultt = [];
resultu = [];
resultw = [];
resulttau = [];

resultdw = [];
resultdu = [];
u = -0.25;
w = -0.03;

for tn = 1:tnmax
    
    m_ss = 0.5*(1+tanh((u-v1)/v2));
    w_ss = 0.5*(1+tanh((u-v3)/v4));
    tau_w = a*1./cosh((u-v3)/(2*v4));
    
    dw_dt = 1/tau_w*(w_ss - w);
    w = w + dw_dt*dt;
    
    du_dt = I-g_Ca*m_ss*(u-E_Ca) - g_k*w*(u-E_K) - g_L*(u-E_L);
    u = u + du_dt * dt;
    
    if (mod(tn,100) == 0)
        resultt = [resultt tn*dt];
        resultu = [resultu u];
        resultw = [resultw w];
%         resultdw = [resultdw dw_dt];
%         resultdu = [resultdu du_dt];
    end   
end

u_t = -1 : 0.01: 0.2;
w_ssnull = 0.5*(1+tanh((u_t-v2)/v4));
m_ssnull = 0.5*(1+tanh((u_t-v1)/v2));
u_ssnull = (I - g_Ca*m_ssnull.*(u_t-E_Ca) - g_L*(u_t-E_L))./(g_k*(u_t-E_K));


figure(1);
grid on;
plot(u_t,w_ssnull,'b--')
hold on;
plot(u_t,u_ssnull,'k--')
% legend('w nullcline','u nullcline','action potential')
axis([-0.4 0.4 -0.05 0.45]);
plot(resultu,resultw, 'r','linewidth',1)
legend('w nullcline','u nullcline','action potential')
xlabel('u(t)');
ylabel('w(t)');
title('w(t) vs u(t) phase plane');
hold off;

figure(2);
plot(resultt,resultu,'b','linewidth',2)
xlabel('time');
ylabel('u(t)');

hold on;
plot(resultt,resultw,'k','linewidth',2)
xlabel('time');
ylabel('w(t)& u(t)');
legend('u(t)','w(t)','Location','Northeast');
title('u & w vs time');
axis([0 80 -0.4 0.5]);
hold off



    
    