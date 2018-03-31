% 2D cable for gap model

% initial values
dt = 0.01;
dx = 0.05;
num = 30/dx; % cell number
L = 1.6;
u = zeros(num,1); % membrane potential
u(1:num) = -0.28236;
u(1:20) = 1;

tmp = zeros(num,1);
w = ones(num,1); % diffusion factor
s = ones(num,1); % change to zero in the gap region 

% gap: row 284 to 316

gap = round(L/dx);
s(num/2-gap/2:num/2+gap/2) = 0;

% build coefficient matrix As(num-1,1);
vs = ones(num-1,1);
vd = -2*ones(num,1);
vd(1) = -1; vd(num) = -1;
A = diag(vd);
A = diag(vd) + diag(vs,1) + diag(vs,-1);

D = 1;
k = D*dt/(dx*dx);
I = eye(num);
Q = I - k*A/2;
Q = inv(Q);
B = I + k*A/2;
QB=Q*B;

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
tmax = 100; % ms
tnmax = round(tmax/dt);
resultu = [];
resultw = [];
resultt = [];
w = zeros(num,1);

tcount = 0;
for tn = 0:tnmax
%     tn/tnmax
    
    m_ss = 0.5.*(1+tanh((u-v1)./v2));
    w_ss = 0.5.*(1+tanh((u-v3)./v4));
    tau_w = a.*1./cosh((u-v3)./(2*v4));
     
    dw_dt = 1./tau_w.*(w_ss - w);
                           
    du_dt = I-g_Ca*m_ss.*(u-E_Ca) - g_k*w.*(u-E_K) - g_L*(u-E_L);
    w = w + dw_dt.*dt;            
    u = QB*u+du_dt*dt.*s;
      
    if (mod(tn,10) == 0)
        tn/tnmax
        resultt = [resultt tn/100];
        resultu = [resultu u];
        resultw = [resultw w];
    end
%      if length(u(u>-0.2)) == 0;
%          break
%      end
end

% calculate delay time
edgecell = resultu(600,:);
[m,n] = size(edgecell);
% transtime = 0;
for cell = 1:n
    if edgecell(cell) > 0.2 & edgecell(a+1) < 0.2
        transtime = cell;
        break
    end
end
transtime

proximal = resultu(283,:);
distal317 = resultu(317,:);
proximalw = resultw(283,:);
distalw = resultw(317,:);
% space time plot
figure(1)
surf(resultu);
shading interp
axis off
view(2)
figure(2)
surf(resultw);
shading interp
axis off
view(2)

figure(3) 
plot(resultt,proximal,'r','linewidth',1.5)
hold on;
plot(resultt, distal317,'b','linewidth',1.5)
axis([0 50 -0.35 0.2]);
title('L = 1.6 proximal&distal');
legend('proximal','distal','location','northeast')
xlabel('t');
ylabel('u');
hold off;

% figure(4)
% plot(resultt, proximalw, 'r','linewidth',1.5)
% hold on;
% plot(resultt, distalw, 'b','linewidth',1.5)
% axis([0 50 0 0.25]);
% hold off
   

    
    
    
    
    
    

    
    
    
    