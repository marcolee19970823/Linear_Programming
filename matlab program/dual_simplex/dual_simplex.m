%Created by Chongnan Li

%Using dual simplex method to solve primal LP:             
% (P)   minimize c'x
%     subject to Ax=b
%                x>=0

%First version: assume you have initial solution,
%such that the dual feasibility and complementary slcakness
%conditions are qualified.

%The example problem is
% min  -2*x1-x2
% s.t. x1+x2+x3   =2
%      x1+      x4=1
%      x1,x2,x3,x4>=0
%
%The optimal objective value: -3
%The optimal solution is x=[1,1,0,0]'
%The initial basic variable are x1, x4
%The corresponding reduced costs for x2,x3 are 1,2>=0(dual feasible).

%Input(Standard form lp) 

function [xB,z_opt,v]=dual_simplex(c,A,b,v)
[m,n]=size(A);

B=A(:,v(1:m));
N=A(:,v(m+1:n));
xB=inv(B)*b;

while 1
%Check primal feasible

if xB>=0
   cB=c(v(1:m),1);
   display(sprintf('optimal solution found.'));
   z_opt=cB' * xB
   xB
   v
   return;
end

%Find basic variable to leave the basis
p=find(xB<0,1);

%Find nonbasic variable to enter the basis(minimum ratio test)
t=inv(B);
tildeA=t*N;
w=tildeA(p,:);

if all (w>=0)
   display(sprintf('primal lp is infeasible.')); 
   return;
end

cN=c(v(m+1:n),1);
cB=c(v(1:m),1);
rN2=cN'-cB'*inv(B)*N;
rN=rN2';

q=find(w<0);

lambda=9999;
for i=1:length(q)
   if lambda >= -rN(q(i))/w(q(i))
        lambda = -rN(q(i))/w(q(i));
        s=q(i); % s: the position in nonbasic variables of the selection
   end
end

%Generate the edge direction
d=-t*N(:,s);
%Generate the step length
b_bar=t*b;
alpha=b_bar(p)/w(s);

xB=xB+alpha*d;

%Update
temp=v(p);
v(p)=v(m+s);
v(m+s)=temp;

N=A(:,v(m+1:n));
B=A(:,v(1:m));

xB(p)=alpha;

end

end