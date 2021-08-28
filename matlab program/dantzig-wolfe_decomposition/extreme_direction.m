function [d]= extreme_direction(c,A,b)

% created by chongnanli
% generate extreme direction for unbounded subproblem SPk:
%
%   minimize (ck)' * (xk)
% subject to (Ak)  * (xk) <= (bk) 
%            (xk) >= 0
% 
% step 1: convert the subproblem to standard form: 
%    minimize (ck)' * (xk) + 0' * s  
%  subject to (Ak) * (xk) + s = (bk)
%             (xk) >= 0 , s >= 0
% 
% where s is the vector of slack variables  
%
% step 2: generate phase-one problem to generate initial solution
%         and note that initial basic solution might include slack
%         variables
%
%   minimize 0'(xk) + 0's +1'u
% subject to (Ak)(xk) + s + Iu = (bk) (bk>=0)
%            (xk) >= 0 , s >= 0, u >= 0


%% Input: please make sure bk>=0

% c=[-2;-3];
% A=[-1,1;
%    -1,2];
% b=[2;8];
% 
% [m,n]=size(A); 


[m,n]=size(A); 


d_len=n; % n denotes the length of extreme direction

%% convert the problem to standard form LP
A_stan=[A,eye(m,m)]; %standard form
c_stan=[c;zeros(m,1)];

[~,n_stan]=size(A_stan);

%% construct phase-one problem
A1=[A_stan,eye(m,m)];
c1=[zeros(n_stan,1);ones(m,1)];
v1=[n_stan+1:n_stan+m,1:n_stan];

[~,~,v1,~]=m_simplex(c1,A1,b,v1);

%display('******** ¡ü phase I problem ¡ü *****');

% if z_opt > 0 
%     display(sprintf('primal lp has no solution!'));
%     return; 
% end

%display('-----------------------------');

%display('******** ¡ý primal lp  ¡ý *****');

v=zeros(1,n_stan);
k=1;
for i=1:n_stan+m
   if v1(1,i)<=n_stan
       v(k)=v1(1,i);
       k=k+1;
   end
end

[~,~,~,d]=m_simplex(c_stan,A_stan,b,v);


temp=d;
d=zeros(d_len,1);
for i=1:n_stan
    index=v(i);
    if index<=d_len
        d(index)=temp(i);
    end
end
%d


end



