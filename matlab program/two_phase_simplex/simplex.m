%Created by Chongnan Li
%Using simplex method to solve standard form lp:
%
%   minimize c'x
% subject to Ax = b
%             x>= 0

% input:
% c : cost efficient vector (column vector)
% A : constraint coefficient matrix
% b : demand vector (column vector)
% v : index set of basic variables and nonbasic variables(row vector)

% output:
% xB : basic variables in optimal solution (column vector)
% z_opt : optimal objective value

% Example:
% min  x1-2*x2
% s.t.  x1+x2+x3    =40
%      2x1+x2  +x4  =60
%       x1,x2,x3,x4>=0

% z_opt=-80
% xB=[40,20] x4, x2 are basic variables
% the initial solution -> v=[3,4,1,2]

% c=[1;-2;0;0];
% A=[1,1,1,0;2,1,0,1];
% b=[40;60];
% v=[3,4,1,2];

function [xB,z_opt,v]=simplex(c,A,b,v)


[m,n]=size(A);

%check optimality

B=A(:,v(1:m));
N=A(:,v(m+1:n));
xB=inv(B)*b;
while 1
    cN=c(v(m+1:n),1);
    cB=c(v(1:m),1);
    rN2=cN'-cB'*inv(B)*N;
    if all(rN2>=0)
        display(sprintf('optimal solution found.'));
        z_opt=cB' * xB
        xB
        v
        if min(rN2)==0
            display(sprintf('lp has infinitely many optimal solutions'));
        end
        return;
    end
    rN=rN2;
    
    %Find nonbasic variable to enter the basis
    %Use Bland's rule to prevent cycling
    q=find(rN<0);
    out=zeros(1,length(q));
    for i = 1:length(q)
        out(i)=v(m+q(i));
    end
    s=min(out); % x_s would be the nonbasic variable to enter the basis
    q=find(v==s);
    q=q-m;      % x_s is qth nonbasic variable
    
    %Generate edge direction
    d=-inv(B)*N(:,q);
    
    if all(d>=0)
        display(sprintf('lp is unbounded below.'));
        z_opt = -inf
        xB
        v
        return;
    end
    
    %Generate step length and the basic variable to exit the basis
    %according to minimum ratio test.
    %Use Bland's rule to prevent cycling.
    lambda=9999;
    for i=1:m
        if lambda > -xB(i)/d(i) && d(i)<0
            lambda =  -xB(i)/d(i);
            p = i; t = v(i);
        elseif lambda == -xB(i)/d(i) && d(i)<0
            p2 = i; t2 = v(i);
            if t2 < t
                p = i; t = v(i);
            end
        end
    end % pth basic variable would leve the basis, which corresponds to x_t
    
    %update
    xB=xB+lambda*d;
    
    xB(p)=lambda;
    
    temp=v(p);
    v(p)=s;
    v(q+m)=temp;
    
    B=A(:,v(1:m));
    N=A(:,v(m+1:n));
end

end


