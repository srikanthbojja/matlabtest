function dfdt=modeleqns_problem_1_a(t,y)

global ne np N d_coeff dd_coeff A B index P P_star mu

dcdt=zeros(N,1);
% pause

%Left Bounday condition
dcdt(1,1)=y(1)-(1/P)*d_coeff(1)*A{1}(1,:)*y(index(1):index(2),1);
%Equation 1
for i=1:ne
    
    dcdt(index(i)+1 : index(i+1)-1,1)=(1/P) * dd_coeff(i)*B{i}(2:np(i)+1,:)*y(index(i):index(i+1),1) ...
                                       - d_coeff(i)*A{i}(2:np(i)+1,:)*y(index(i):index(i+1),1) ...
                                       - mu*P_star*y(index(i)+1:index(i+1)-1,1).*(1-y(N+index(i)+1:N+index(i+1)-1,1)) - y(N+index(i)+1:N+index(i+1)-1,1);
    
    if i<ne        
        dcdt(index(i+1),1) = d_coeff(i)*A{i}(np(i)+2,:)*y(index(i):index(i+1),1) - d_coeff(i+1)*A{i+1}(1,:)*y(index(i+1):index(i+2),1);
    end
        
end
%Right Bounday condition
dcdt(end,1)=d_coeff(ne)*A{ne}(end,:)*y(index(end-1):index(end),1)
% save dcct
% pause
%Equation 2
dqdt=P_star*y(1:N,1).*(1-y(N+1:2*N,1)) - y(N+1:2*N,1);

dfdt=[dcdt;dqdt];

end

