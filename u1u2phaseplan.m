% u1 u2 phase plane

tmax = 100;
dt = 0.01;
tnmax = round(tmax/dt);
resultu10 = [];
resultu20 = [];

alpha = 0.2;
cg = 0.042;

A = 0.259;
A1 = 0.197;

u10 = 0.35;
u20 = -0.1;

for tn = 0:tnmax
    du10_dt = u10*(1-u10)*(u10-alpha) + cg*(u20-u10);
    u10 = u10 + du10_dt*dt;
    
    du20_dt = u20*(1-u20)*(u20-alpha) + cg*(u10-u20);
    u20 = u20 + du20_dt*dt;
    
    if (mod(tn,10)) == 0
        resultu10 = [resultu10 u10];
        resultu20 = [resultu20 u20];
    end
end

u11 = -0.1;
u21 = A1;
resultu11 = [];
resultu21 = [];

for tn = 0:tnmax
    du11_dt = u11*(1-u11)*(u11-alpha) + cg*(u21-u11);
    u11 = u11 + du11_dt*dt;
    
    du21_dt = u21*(1-u21)*(u21-alpha) + cg*(u11-u21);
    u21 = u21 + du21_dt*dt;
    
    if (mod(tn,10)) == 0
        resultu11 = [resultu11 u11];
        resultu21 = [resultu21 u21];
    end
end

U1 = -0.05:0.01:0.3;
U2 = -0.05:0.01:0.3;
t = -0.2:0.01:0.9;

% null cline
u1_null = U2 - U2.*(1-U2).*(U2-alpha)/cg;
u2_null = U1 - U1.*(1-U1).*(U1-alpha)/cg;

figure(1);
plot(U1, u2_null, 'r--')
hold on;
plot(u1_null, U2, 'b--')
plot(t,t,'k--')
plot(resultu10,resultu20,'m')
plot(resultu11,resultu21,'m')
grid on
axis([-0.1,0.9,-0.1,0.9]);


