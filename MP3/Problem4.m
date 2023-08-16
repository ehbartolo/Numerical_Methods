% Problem # 4 Source Code:

xa = -0.5;
xb = 3;
ta = 0;
tb = 2;
e = 0.05;
N = 99;

dt = 0.001;
dx = (xb-xa)/(N+1);
x = linspace(xa+dx,xb-dx, N)';               % elements at x-direction
t = [ta:dt:tb]';                            % elements at t-direction

ut0 = exp(-4*( (x-ones( size(x) )).^2 ));   % initial condition at t = 0;

%% Plot initial condition
plot (x,ut0,'-')
title ('Initial condition for t distribution')
xlabel ('x')
ylabel ('u')
%% Perform ODE-IVP using RK45
u = zeros(size(x,1),size(t,1));             % u(x,t)
xmax = size(x,1);
for i = 1:xmax
    %% To inlclude solver for system of nonlinear ODEs......
       
end

%% RK45 Code:
function [T,Y] = RK45(yah)
  hmin = 0.0001;
  hmax = 10/2;
  tol = 0.00001;
  t0 = 0;
  y0 = 1;
  tf = 10;
  i = 0; %step number
  
  h = hmax;
  yi = y0;
  ti = t0;
    
  T = t0;
  Y = y0;
  %fprintf('Step %d: t = %6.4f, y = %18.15f\n', 0, ti, yi);
  while ti<tf
    %h = min(h, 2-ti);
    
    %Coefficients for the computation of k's
    [k1,k2,k3,k4,k5,k6] = CompKs(h,ti,yi);
    yii = yi + 25*k1/216 + 1408*k3/2565 + 2197*k4/4104 - k5/5;
    zii = yi + 16*k1/135 + 6656*k3/12825 + 28561*k4/56430 - 9*k5/50 + 2*k6/55;
    R = abs(zii-yii)/h;
    s = 0.84*(tol/R)^(1/4);
    if (R<=tol) 
      ti = ti+h;
      yi = yii;
      i = i+1;
      if ti> tf
        ti = tf;
      end
      %fprintf('Step %d: t = %6.4f, y = %18.15f\n', i, ti, yi);
      T = [T;ti];
      Y = [Y;yi];
      h = s*h;

    else
      h  = s*h;
    end
  end
end

%First Derivative
function f = yp(x,y)
    e = 0.05;
%     f = (e/(h^2)) * (uii-2*ui+ui0) - ui*(uii-ui0)/(2*h);
    f = -2*y + exp(-2*((x-6)^2));
  end

function [k1,k2,k3,k4,k5,k6] = CompKs(h,tk,yk)
    c20 = 0.25;
    c21 = 0.25;
    c30 = 3/8;
    c31 = 3/32;
    c32 = 9/32;
    c40 = 12/13;
    c41 = 1932/2197;
    c42 = -7200/2197;
    c43 = 7296/2197;
    c50 = 1;
    c51 = 439/216;
    c52 = -8;
    c53 = 3680/513;
    c54 = -845/4104;
    c60 = 0.5;
    c61 = -8/27;
    c62 = 2;
    c63 = -3544/2565;
    c64 = 1859/4104;
    c65 = -11/40;
    
    k1 = h*yp(tk, yk);
    k2 = h*yp(tk+c20*h, yk+c21*k1);
    k3 = h*yp(tk+c30*h, yk+c31*k1+c32*k2);
    k4 = h*yp(tk+c40*h, yk+c41*k1+c42*k2+c43*k3);
    k5 = h*yp(tk+c50*h, yk+c51*k1+c52*k2+c53*k3+c54*k4);
    k6 = h*yp(tk+c60*h, yk+c61*k1+c62*k2+c63*k3+c64*k4+c65*k5);
end
