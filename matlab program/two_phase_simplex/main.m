% Created by Chongnan Li
% Two phase method to generate initial basic feasible solution
% and use simplex method to solve the standard form lp:
%
%   minimize c'x
% subject to Ax = b 
%             x>= 0
%
% Phase I problem is:
%
%   minimize 0'x+1'u
% subject to Ax + Iu = b (b>=0)
%                 x >= 0 

clc;

%Input: please make sure b>=0
c=[1;-2;0;0];
A=[1,1,1,0;2,1,0,1];
b=[40;60];

[m,n]=size(A);
%construct phase one problem
A1=[A,eye(m,m)];
c1=[zeros(n,1);ones(m,1)];
v1=[n+1:n+m,1:n];

[xB,z_opt,v1]=simplex(c1,A1,b,v1);

display('******** ¡ü phase I problem ¡ü *****');

if z_opt > 0 
    display(sprintf('primal lp has no solution!'));
    return; 
end
display('-----------------------------');

display('******** ¡ý primal lp  ¡ý *****');

v=zeros(1,n);
k=1;
for i=1:n+m
   if v1(1,i)<=n
       v(k)=v1(1,i);
       k=k+1;
   end
end


[xB,z_opt,v]=simplex(c,A,b,v);



