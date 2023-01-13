% 240032 ExFinal Q1 2022-23
% Problema 4

clearvars
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kc = 5.0;
f = 22.09;
tempULS = 24.10; %Fixed temperature on the upper left side of the rombus
                 % nodes = 1,2,37
% Convection on the circle boundary
tempInf = 21.0;
beta = 5.0;

% Other data
numElemPartA = 50;
numElemPartC = 31;
bcPartC = [0.25,0.25,0.5];
numNodHintPartC = 50;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%eval('meshcirclediamont');
%save meshcirclediamont.mat nodes elem -mat
load meshcirclediamont.mat;

numNod = size(nodes,1);
numElem  = size(elem,1);

numbering = 1;                        %Set to 1 to identify the nodes on 
                                      %the upper side of the rombus 
plotElements(nodes, elem, numbering); %In this case is better to use
                                      %plotElements rather than 
                                      %plotElementsOld
hold on
[indNodsBd, ~, ~, ~] = boundaryNodes(nodes, elem);
indNodsCirc = find(sqrt((nodes(:,1)-0.8).^2 + (nodes(:,2)-0.5).^2) > 0.43);
indNodsIntBd = setdiff(indNodsBd,indNodsCirc);
indNodsULS = [1;2;37]; %these nodes have been found 'by inspection', 
                       %looking at the figure (set numbering = 1)

plot(nodes(indNodsCirc,1),nodes(indNodsCirc,2),...
    'o','MarkerFaceColor','red','MarkerSize',10)
plot(nodes(indNodsIntBd,1),nodes(indNodsIntBd,2),...
    'o','MarkerEdgeColor','blue','MarkerSize',10)
plot(nodes(indNodsULS,1),nodes(indNodsULS,2),...
    'o','MarkerEdgeColor','green','MarkerSize',10)
hold off

a11=kc;
a12=0;
a21=a12;
a22=a11;
a00=0;
%f=; %already defined
coeff=[a11,a12,a21,a22,a00,f];

K=zeros(numNod);   %global Stiff Matrix
F=zeros(numNod,1); %global internal forces vector
Q=zeros(numNod,1); %global secondary variables vector
for e=1:numElem 
  % compute element system  
    [Ke,Fe]=linearTriangElement(coeff,nodes,elem,e);
  %
  % Assemble the elements
  %
    rows=[elem(e,1); elem(e,2); elem(e,3)];
    colums= rows;	
    K(rows,colums)=K(rows,colums)+Ke; %assembly
    if (coeff(6) ~= 0)
        F(rows)=F(rows)+Fe;
    end
end % end for elements
%we save a copy of the initial F array
%for the postprocess step
Kini=K;
Fini=F;
%Boundary Conditions (BC)
%Impose Boundary Conditions for this example.
fixedNodes=indNodsULS'; %Fixed Nodes (global num.)
%Remark: unique command not necessary, since boundaries ere disjoint
freeNodes = setdiff(1:numNod,fixedNodes); %Complementary of fixed nodes
% ------------ Convection BC
convecNodes=indNodsCirc';
[K,Q]=applyConvTriang(convecNodes,beta,tempInf,K,Q,nodes,elem);
% ------------ Essential BC 
u=zeros(numNod,1);     %initialize u vector
u(indNodsULS)=tempULS; %fixed temperature on the upper left side of 
                       %the rombus nodes = 1,2,37
%modify the linear system  
Km=K(freeNodes,freeNodes);
Im=F(freeNodes)+Q(freeNodes)...
    -K(freeNodes,fixedNodes)*u(fixedNodes); %u can be different from zero
                                            %at fixed nodes
%Compute the solution of the reduced system
um=Km\Im;
u(freeNodes)=um;
%PostProcess: Compute secondary variables and plot results
Q=Kini*u-Fini;
titol='Temperature Distribution';
colorScale='jet';
plotContourSolution(nodes,elem,u,titol,colorScale);

clc
fprintf('\tPROBLEM 4\n')
%PART (A)
[Ke,Fe]=linearTriangElement(coeff,nodes,elem,numElemPartA);
fprintf(['Part (a)\n',...
         '*** K^{%d}_{2,3} = %.4e\n',...
         '*** Hint. K^{%d}_{1,1} = %.4e\n'],...
         numElemPartA,Ke(2,3),numElemPartA,Ke(1,1))

%PART (B)
fprintf(['Part (b)\n',...
         '*** The mean value of the temperature in the nodes of the\n', ...
         '*** romboid-shaped hole boundary is <T> = %.4e\n',...
         '*** Hint. The maximum value of the temperature at the\n',...
         '*** nodes is max T = %.4e\n'],...
         sum(u(indNodsIntBd))/length(indNodsIntBd),max(u))

%PART (C)
uNodElem = u(elem(numElemPartC,:));
tempP = bcPartC*uNodElem;
fprintf(['Part (c)\n',...
         '*** The temperature for the point p of the element %d with\n',... 
         '*** barycentric coordinates (%f, %f, %f)\n',...
         '*** is T = %.4e\n',...
         '*** Hint. The temperature on the node %d is T = %.4e\n'],...
         numElemPartC,bcPartC,tempP,numNodHintPartC,u(numNodHintPartC))