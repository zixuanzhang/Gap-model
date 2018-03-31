% 2D cable for gap model

% initial values
dt = 0.01;
dx = 0.01;
num = 30/dx; % cell number
L = 1.52;
u = zeros(num,1); % membrane potential
u(1:num) = -0.2824;
u(1:20) = 1;

tmp = zeros(num,1);
w = zeros(num,1); % diffusion factor
w(1:5) = 0.004;
s = ones(num,1); % change to zero in the gap region 

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
w = zeros(num,1);

for tn = 0:tnmax
%     tn/tnmax
    
    m_ss = 0.5.*(1+tanh((u-v1)./v2));
    w_ss = 0.5.*(1+tanh((u-v3)./v4));
    tau_w = a.*1./cosh((u-v3)./(2*v4));
     
    dw_dt = 1./tau_w.*(w_ss - w);
    w = w + dw_dt.*dt;        
                           
    du_dt = I-g_Ca*m_ss.*(u-E_Ca) - g_k*w.*(u-E_K) - g_L*(u-E_L);
    u = u + du_dt*dt.*s;      
    u = QB*u;
%     end
      
    if (mod(tn,100) == 0)
        tn/tnmax
        resultu = [resultu u];
        resultw = [resultw w];
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
    
    

    
    
    
    
    
    

    
    
    
    