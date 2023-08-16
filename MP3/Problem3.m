% Problem #3.
% Disclaimer:   The code is not orignally made by the graduate student.
%               However, considerable work was done to edit the B.C.s and
%               other parameters to come up with the solution.
%               Please see reference in the documentation.

% Author Remarks:
% poisson2.m  -- solve the Poisson problem u_{xx} + u_{yy} = f(x,y)
% on [a,b] x [a,b].  
% 
% The 5-point Laplacian is used at interior grid points.
% This system of equations is then solved using backslash.
% 
% From  http://www.amath.washington.edu/~rjl/fdmbook/chapter3  (2007)


a = 0; 
b = 1;
h = 0.01;
m = (b-a)/h - 1;
% m = 99;
% h = (b-a)/(m+1);
x = linspace(a,b,m+2);   % grid points x including boundaries
y = linspace(a,b,m+2);   % grid points y including boundaries

[X,Y] = meshgrid(x,y);      % 2d arrays of x,y values
X = X';                     % transpose so that X(i,j),Y(i,j) are
Y = Y';                     % coordinates of (i,j) point

Iint = 2:m+1;              % indices of interior points in x
Jint = 2:m+1;              % indices of interior points in y
Xint = X(Iint,Jint);       % interior points
Yint = Y(Iint,Jint);

%f = @(x,y) 1.25*exp(x+y/2);         % f(x,y) function

%rhs = f(Xint,Yint);       % evaluate f at interior points for right hand side
                           % rhs is modified below for boundary conditions.
rhs  = ones(size(Xint))*(-100000/16);

utrue = exp(0*X+0*Y);        % true solution for test problem

% set boundary conditions around edges of usoln array:

usoln = utrue;              % use true solution for this test problem
                            % This sets full array, but only boundary values
                            % are used below.  For a problem where utrue
                            % is not known, would have to set each edge of
                            % usoln to the desired Dirichlet boundary values.

% adjust the rhs to include boundary terms:

alpha = (75*h)/(5*h+1);
beta = 1/(5*h+1);
rhs(:,1) = rhs(:,1) - alpha*ones(size(rhs(:,1)))/h^2;
rhs(:,m) = rhs(:,m) - alpha*ones(size(rhs(:,m)))/h^2;
rhs(1,:) = rhs(1,:) - alpha*ones(size(rhs(1,:)))/h^2;
rhs(m,:) = rhs(m,:) - alpha*ones(size(rhs(m,:)))/h^2;

% convert the 2d grid function rhs into a column vector for rhs of system:
F = reshape(rhs,m*m,1);

% form matrix A:
I = speye(m);                        % mxm identitiy matrix
e = ones(m,1);
T = spdiags([e -4*e e],[-1 0 1],m,m);    
S = spdiags([e e],[-1 1],m,m);

A = (kron(I,T) + kron(S,I)) / h^2;


% Solve the linear system:
uvec = A\F;  

% reshape vector solution uvec as a grid function and 
% insert this interior solution into usoln for plotting purposes:
% (recall boundary conditions in usoln are already set) 

usoln(Iint,Jint) = reshape(uvec,m,m);

% % assuming true solution is known and stored in utrue:
% err = max(max(abs(usoln-utrue)));   
% fprintf('Error relative to true solution of PDE = %10.3e \n',err)

% plot results:

clf
hold on

% plot grid:
% plot(X,Y,'g');  plot(X',Y','g')

% plot solution:
contour(X,Y,usoln,30,'k')

axis([a b a b])
daspect([1 1 1])
title('Contour plot of computed solution')
hold off
