% integrate and relax model

tmax = 50;
dt = 0.01;
tnmax = round(tmax/dt);
resultu = [];
t = (0:0.1:tmax);

alpha = 0.2;
u = 0.21;
usub = 0.18;
resultusub = [];
% superthreshold
for tn = 0:tnmax
    
    du_dt = u*(1-u)*(u-alpha);
    u = u + du_dt*dt;   
 
    if u > 0.95
        u = -0.1;
    end
    
    if (mod(tn,10) == 0)
        resultu = [resultu u];
    end
    
end

% subthreshold
for tn = 0:tnmax
    
    dusub_dt = usub*(1-usub)*(usub-alpha);
    usub = usub + dusub_dt*dt;   
 
%     if u > 0.95
%         u = -0.1;
%     end
    
    if (mod(tn,10) == 0)
        resultusub = [resultusub usub];
    end
end

figure(1)
plot(t,resultu,'b');
hold on
plot(t,resultusub,'k--')
xlabel('t');
ylabel('u');
title('integrate and relax model');



    