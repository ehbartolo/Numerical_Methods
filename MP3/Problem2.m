%ES 204 Problem # 2
N = 100/(0.005);
% [y,t] = euler2(0,100,2,0,N);
[y2,t] = heun2(0,100,2,0,N);
plot(t,y2);

legend('implicit heun with h = 0.005')
title('Solution of ODE using Implicit Heun');
xlabel('t');
ylabel('ý')

function [y, t] = euler2(a,b,ya,va,n)
h = (b - a) / n;
y(1) = ya;
v(1) = va;
t(1) = a;

for i = 1 : n
    y(i+1) = y(i) + h * yp(y(i),v(i));
%   v(i+1) = v(i) + h* vp(y(i+1),v(i));                         % Normal euler method
    v(i+1) = (v(i) + h *(-1)*y(i+1))/(1-h*25*(1-(y(i+1))^2));   % modified euler method, implicit
    t(i+1) = t(i) + h;
end

end

function [y, t] = heun2(a,b,ya,va,n)
tol  = 0.1;
h = (b - a) / n;
halfh = h / 2;
y(1) = ya;
v(1) = va;  
t(1) = a;
err = 1000;

for i = 1 : n
    
     g1 = yp(0,v(i));           % f(yi,vi) for y
     g2 = vp(y(i),v(i));        % f(yi,vi) for y'
     
     %% Predictor for y'
     if i == 1
        % initial condition
        z2 = v(i) + h * g2;             % predictor of y'
        z2u = z2;
     else
        % central difference method
        z2 = v(i-1) + 2*h * g2;         % predictor of y', O(h3)
        z2u = z2;
     end
     %% Predictor Modifier
     if i > 1
         z2 = z2 + (4/5)*(v(i)-z1u);
     end
     %% Corrector for y and y'
% for k = 1:1000
         y(i+1) = y(i) + halfh * ( g1 + yp(0,z2) );     % corrector for y
         v(i+1) = v(i) + halfh * ( g2 + vp(y(i+1),z2)); % corrector for y'
         
         err = abs(v(i+1)-z2);
         fprintf('i: %d,err: %f\n',i,err);
%        if err < tol
%           break;
%        else
%           z2 = v(i+1);
%        end   
%     end
     %% Corrector Error Estimate
     v(i+1) = v(i+1) - (1/5)*(v(i+1)-z2);
     %% Update time
     t(i+1) = t(i) + h;
     z1u = z2u;

end
end

function res = yp(A,B)

    res = B;

end

function res = vp(A,B)
    
    res = 25*(1-A^2)*B - A;

end