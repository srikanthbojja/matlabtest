close all;
clc;
clear all;
global ne np w_mod N d_coeff dd_coeff A B r index  I r_mod X P P_star mu

tic

% ne  = 10; % ne is the number of elements in the x-direction 
% np  = [20 34 26 28 29 22 23 25 26 27]; % np is the number of internal collocation points in one element 
% w   = [.15 .81 .84 .91 0.92 0.93 0.95 0.97 0.99]; % (ne-1) values
% P_star = 3.5;
% P=40;
% mu=0.0449;

% ne  = 7; % ne is the number of elements in the x-direction 
% np  = [20 34 26 28 29 22 23]; % np is the number of internal collocation points in one element 
% w   = [.15 .81 .84 .91 0.93 0.99]; % (ne-1) values
% P_star = 3.5;
% P=1;
% mu=0.0449;

% ne  = 3; % ne is the number of elements in the x-direction 
% np  = [10 10 10]; % np is the number of internal collocation points in one element 
% w   = [.1 .21]; % (ne-1) values
% P_star = 3.5;
% P=40;
% mu=0.0449;

ne  = 2; % ne is the number of elements in the x-direction 
np  = [10 10]; % np is the number of internal collocation points in one element 
w   = [0.5]; % (ne-1) values
P_star = 3.5;
P=40;
mu=0.0449;

% ne  = 1; % ne is the number of elements in the x-direction 
% np  = [100]; % np is the number of internal collocation points in one element 
% w   = []; % (ne-1) values
% P_star = 3.5;
% P=40;
% mu=0.0449;

w_mod = [0,w,1]; 
index = [1];
X     = [];

for i = 1:ne
    [r{i},A{i},B{i},q{i}] = colloc(np(i),1,1);
    I{i} = eye(np(i)+2,np(i)+2);
    d_coeff(i) =1 /(w_mod(i+1)-w_mod(i)); %First-derivative coefficient(constant) after Linear transformation
    dd_coeff(i) = d_coeff(i)^2; %Second-derivative coefficient(constant) after Linear transformation
    index = [index (sum(np(1:i))+i+1)]; % Cumulative Index
    r_mod{i} = w_mod(i)+r{i}*(w_mod(i+1)-w_mod(i));
    X=[X;r_mod{1,i}];
end

X=unique(X);

N = index(end); %n_tot_int + (ne-1) + 2; %Total collocation points i.e. (Total Internal collocation points + common points + boundary points)

% Initial Guess

Y0=ones(2*N,1);


% Mass Matrix implementation for ode15s & ode23s
M=[];Mq=[];
for i = 1:ne
M = [M,0,ones(1,np(i))];
Mq = [Mq,1,ones(1,np(i))];
end
M=[M,0];
Mq=[Mq,1];
M = diag(M);
Mq = diag(Mq);
M = [M,zeros(size(M));zeros(size(M)),Mq];

options=odeset('Mass',M,'MassSingular','yes','RelTol',1e-3,'AbsTol',1e-6,'Vectorized','off','Stats','on','MaxOrder',5);
[t,y]=ode15s(@modeleqns_problem_1_a,0:0.5:5,Y0,options);

% options=odeset('Mass',M,'RelTol',1e-3,'AbsTol',1e-6);
% [t,y]=ode23s(@righthand,[0 0.1],Y0,options);


figure1=figure(1);
axes1 = axes('Parent',figure1);
hold all 
% labels = cellstr( num2str(t'));
for i=1:size(y)
plot(X,y(i,1:N),'displayname',strcat('t=',num2str(t(i))));
end
% text(X(30,1), y(2,30), labels, 'VerticalAlignment','bottom','HorizontalAlignment','right')
% legend 'show'

figure(2)
plot(t,y(:,1),'r-*');
hold all

toc