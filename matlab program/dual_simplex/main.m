% created by chongnan li
% using dual simplex method to solve the standard form LP:
%
% (P)minimize c'x
% subject to  Ax = b
%             x>= 0
%
% constuct the following model to generate initial solution 
% that satisfies dual feasibility:
%   minimize c'x
% subject to Ax = Be 
%             x>= 0
%
% B(nonsingular) is a basis of (P)

clc;

c=[-2;-1;0;0];
A=[1,1,1,0;1,0,0,1];
b=[2;1];


[m,n]=size(A);
B=A(1:m); %If B is singular, choose other m columns of A to form the basis.
b1=B*ones(m,1);
v1=1:n; % If B is singular, v1 should change accordingly!

[xB,z_opt,v]=simplex(c,A,b1,v1);

display('******** ¡ü to generate initial dual feasible solution ¡ü *****');


if z_opt == -inf
   display('primal lp is unbounded below.'); 
   return;
end

display('-----------------------------');

display('******** ¡ý primal lp  ¡ý *****');

[xB,z_opt,v]=dual_simplex(c,A,b,v);