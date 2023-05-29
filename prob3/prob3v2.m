% 240032 ExFinal Q1 2022-23
% Problema 3

% Given the BVP
%
%   (a1(x)u')' + a0(x) u = f(x), x \in (a,b)
%     u'(a) + alpha u(a) = beta,
%                  u'(b) = gamma,
%
% with,
%   a = 0, b = 1,
%   a1(x) = cos(x),
%   a0(x) = a0, some constant value            
%   f(x) = f0 = -50.0,
%   alpha = 0, in parts (a) and (b)
%   alpha = 1, in part (c)
%   beta = 0.5,
%   gamma = 2.5;
%
% and we consider its FEM solution using n=N linear elements
clearvars
close all

% Data:
n = 100;                  % number of divisions = number of linear elements
a0 =  6;                  % EDO's linear part's (constant) coefficient 

a = 0; b = 1;
a1 = @(x) cos(x);
A1 = @(x) sin(x);          % A1(x) = Int(a1(x))
f0 = -50.0;                % f(x) = f0 x, with f0 = 50.0

xp = pi/6;                 % Point at which we the solution is approximated 
                           % by interpolation (in part (b))

beta = -0.5; gamma = 2.5;  % rhs of the (Robin) B.C.
alpha = 1.0;               % coefficient at the lhs of the (Robin) B.C.
                           % for part (c) (for parts (a) and (b) 
                           % alpha = 0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n1 = n+1;                  %number of nodes
h = (b-a)/n;               %length of the elements

nodes = linspace(a,b,n1)'; %position of the nodes
elem = [(1:n)',(2:n1)'];   %connectivity matrix

numNodes = size(nodes,1);
numElem = size(elem,1);

K1 = [1, -1; -1, 1]/h^2; K0 = a0*h*[2, 1; 1, 2]/6; Fe = f0*h*[1;1]/2;
K = zeros(n1);
F = zeros(n1,1);
Q = zeros(n1,1);

for e=1:numElem
    rows = [elem(e,1), elem(e,2)];
    cols = rows;
    x1 = nodes(rows(1,1)); x2 = nodes(rows(1,2));
    Ke = (A1(x2)-A1(x1))*K1+K0;
    if e == 1
        Ke1=Ke(1,1);      %save the stifness matrix of the 1st element 
    end
    K(rows,cols) = K(rows,cols) + Ke;
    F(rows)=F(rows)+Fe;
end

clc
fprintf('\tPROBLEM 3\n')
% Part (A)
%BC
%fixed nodes: none
%Natural BC
Q(1) = -beta*a1(a);
Q(end) = gamma*a1(b);
%Essential BC: none
%Solve the system:
Qm = Q + F; 
u = K\Qm;
[minU,nodMin] = min(u);
fprintf(['Part (a)\n',...
         '*** min u = %.4e\n',...
         '*** Hint. K¹(1,1) = %.4e and u(10) = %.4e\n'],...
         minU,Ke1(1,1),u(10));

% Part (B)
interpU = interp1(nodes,u,xp);
fprintf(['Part (b)\n',...
         '*** The interpolated value at pi/6 is: %.4e\n'],interpU)

% Post-process: plot
% Approximate u(x) by a linear spline using the values at the nodes 
xx = linspace(a,b,1000*n+1);
u1 = interp1(nodes,u,xx);
plot(xx,u1,'-','color','blue','lineWidth',2)
hold on
% Show the minimum value of u
plot(nodes(nodMin),u(nodMin),'o','MarkerFaceColor',...
    'red','MarkerEdgeColor','black','MarkerSize',6)
% Show the interpolated value of u at x = pi/6
plot(xp,interpU,'o','MarkerFaceColor',...
    'green','MarkerEdgeColor','black','MarkerSize',6)
title('Plot of the solution for $\alpha = 0$',...
    'Interpreter','latex','FontSize',14)
xlabel('$x$','Interpreter','latex','FontSize',14)
ylabel('$u$','Interpreter','latex','FontSize',14)
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Part (C)
% for the sake of clarity, ve compute again K, F and Q
clear K F Q u;
K = zeros(n1);
F = zeros(n1,1);
Q = zeros(n1,1);

for e=1:numElem
    rows = [elem(e,1), elem(e,2)];
    cols = rows;
    x1 = nodes(rows(1,1)); x2 = nodes(rows(1,2));
    Ke = (A1(x2)-A1(x1))*K1+K0;
    if e == 1
        Ke1=Ke(1,1);      %save the stifness matrix of the 1st element 
    end
    K(rows,cols) = K(rows,cols) + Ke;
    F(rows)=F(rows)+Fe;
end

Q(1) = -beta*a1(a);
Q(end) = gamma*a1(b);
%Essential BC: none
%Solve the system:
Qm = Q + F; 
K(1,1) = K(1,1)-alpha*a1(a);         %add the "spring-like" term to the lhs
u = K\Qm;
avU = sum(u)/n1;
fprintf(['Part (c)\n',...
         '*** For alpha = %.2f, <u> = %.4e\n',...
         '*** Hint. For this value of alpha, u(10) = %.4e\n'],...
         alpha,avU,u(10))






