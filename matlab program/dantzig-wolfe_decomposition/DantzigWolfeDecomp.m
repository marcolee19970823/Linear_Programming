%DantzigWolfeDecom solves a linear programming with a special structure:
%        minimize     c'*x
%      subject to   A*x <= b
%                     x >= 0
%where A is an m*n matrix which can be written as block angular form:
%       __           __        __  __        __  __        __  __
%      | L1 L2 ... LK  |      |  x1  |      |  c1  |      |  b0  |
%      | A1            |      |  x2  |      |  c2  |      |  b1  |
%  A = |    A2         |  x = |  ..  |  c = |  ..  |  b = |  b2  |
%      |       ...     |      |  ..  |      |  ..  |      |  ..  |
%      |__         AK__|      |__xK__|      |__cK__|      |__bK__|
%
%so the LP can be decomposed into a master problem(MP) and K subproblems(SPk)
%we can rewrite the MP as a restricted MP by Resolution Theorem
%
% Inputs:
% mast: a struct includes MP's coefficients, i.e., L1,...,LK, bO
% sub: a struct includes coefficients of SPk's, i.e., c1,...,cK, A1,...,AK
%      b1,...,bK and the initial extreme points v1,...,vK. (c,A,b,v)
% K: the number of subproblems
%
% Outputs:
% x: n*1 vector, optimal solution
% fval: scalar, optimal objective value
% bound: a matrix includes all LBs and UBs for each iteration (each column)
% exit_flag: describes the exit condition of the problem as follows:
%             1 - optimal solution
%            -3 - LP is unbounded below

function [x,fval,bound,exit_flag] = DantzigWolfeDecomp(mast,sub,K)
x=[]; bound=[];
%% Step 0: Initialization
%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate an initial basis B for the master problem.
% Let x_B be the basic variables and B_bar the index set of the basic variables
% Set all other variables to zero to get the restricted master problem.
% Go to STEP 1.
s=mast.b; %slack variables for inequality constraints
x_B=[s;ones(K,1)]; %initial basic variables
%x_Bflag is an index of the basic variables in the restricted master
%problem associated with linking constraints and subproblem k, i.e.,
%slack variables of linking constrints are initially basic and other basic
%variables associated with subproblems are set to 1.
x_Bflag=[zeros(length(s),1);[1:K]'];
f_sub=[];
for k=1:K
    %obtain initial extreme points from struct sub for subproblems there are
    %zero vectors, so initial objective function values will be zero,
    %v(k) is initial extreme point of subproblem k
    v_sub{k}=sub.v{k};
    f_sub=cat(1,f_sub,sub.c{k}'*v_sub{k});
    for a=1:length(s)
        %generating initial extreme point for linking constraints
        %v_L{a,k}= initial extreme point of linking costraint Lk
        %b(k) is the RHS vector in A_k*x^k <= b^k
        v_L{a,k}=zeros(length(sub.b{k}),1);
    end
end
f_s=zeros(length(s),1);
f_B=[f_s;f_sub]; %initial f_B i.e., the objective coefficient of
%the restricted MP
B=eye(length(x_B)); %initial basis for master problem
iter_num=0;         %counter
options=optimset('LargeScale','on'); %choose largescale LP
%solver in linprog
while iter_num>=0
    %% Step 1: Simplex Multiplier Generation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Solve for duals by solving the system B^T*pie=f_B, then Go to STEP2.
    pie=B'\f_B;   %solve B^T*pie=f_B
    pie_sub=pie(end-K+1:end); %duals of kth convexity constraints, pie^2_k
    pie(end-K+1:end)=[];      %duals of linking constraints, pie^1
    %% Step 2: Optimality Check
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % For each k=1,...,K, using the revised simplex method to solve SP_k,
    % i.e.
    %        min  sig_k = [ (c^k)^T - (pie^1)^T*L_k ] * x^k
    %       s.t.  A_k * x^k <= b^k
    %                   x^k >= 0
    % If SP_k is unbounded, then Go to STEP3, else let
    % x^k = (v_i*)^k denote the optimal basic feasible solution,
    % compute (r^k)_* = (sig_k)^* - pie_k^2.
    % If r_min = {(r^k)_*} >= 0, then the current basis B is optimal,
    % else Go to STEP3.
    exitflag=zeros(K,1);
    for k=1:K  
        c_sub=[sub.c{k}'-pie'*mast.L{k}]';  %update the objective coefficient
        [x_sub{iter_num+1, k}, sig{k}, exitflag(k)] = ...
            linprog(c_sub, sub.A{k}, sub.b{k}, [],[], ...
            zeros(length(c_sub),1),[],[],options); %call linprog solver
        if exitflag(k)==-3
            sig{k}=-inf;
        end
    end
    sig2=sig;   %sig is cell; sig2 is vector
    sig2=cell2mat(sig2);
    sig2=sig2';
    sig2=sig2-pie_sub;          %compute (r^k_*)
    r_min=min(sig2);          %minimum ratio test to obtain r_min
    r_minflag=find(sig2==r_min);
    
    if isempty(find(r_min<0)) || abs(r_min) <= 1e-8 % reduced cost >= 0, optimal
        disp('Problem solved!')
        fval=0;
        x_Bflag_s=x_Bflag(1:length(s));
        for k=1:K  %convert to optimal solution for original LP
            x{k,1}=x_B(length(s)+k) * v_sub{k};
            for a=1:length(x_Bflag_s)
                if x_Bflag_s(a)==k
                    x{k}=x{k}+x_B(a)*v_L{a,k};
                end
            end
            fval=fval+sub.c{k}'*x{k}; %convert to optimal objective value for original LP
        end
        x=cell2mat(x);
        exit_flag=1;
        break
    else
        %% Step 3: Column Generation
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % If all subproblems are bounded and r_min = {(r^k)_*}<0,
        % then let t be the index of k in SP_k such that r_min=(r^t)_*.
        % Let a_bar=[q_i*^t e_t]^T=[L_t*v^t_i*  e_t]^T
        % where v^t_i* is the optimal extreme point of SP_t and Go to STEP4.
        % Else there is a subproblem SP_s that is unbounded and so an
        % extreme direction d_j*^s will be generated such that
        % [(c^s)^T-(pie^1)^T*L_s] * d_j*^s < 0 and
        % so let a_bar = [(q_bar)_j*^s 0]^T = [L_s*d_j*^s 0]^T
        % and go to STEP4.
        if length(find(exitflag==1))==K % if all subproblems bounded and r_min<0
            t=r_minflag(1);   %subproblem t such that r_min = r^t_*
            q_t=mast.L{t}*x_sub{iter_num+1,t}; %generate q_i*t
            e=zeros(length(x_B)-length(q_t),1);
            e(t)=1;                   %generate e_t
            a_bar=[q_t;e];        %genrate a_bar
        else
            disp('unbounded subproblem exist')
            unboundflag=find(exitflag==-3);
            t=unboundflag(1);    %subproblem s with extreme direction d_j*^s
            [ex_d]= extreme_direction(sub.c{t},sub.A{t},sub.b{t}); %ex_d is the extreme direction
            x_sub{iter_num+1,t} = ex_d;
            q_t=mast.L{t}*ex_d; %generate (q_bar)_j*^s
            %generate a_bar
            a_bar=[q_t;zeros(K,1)];
        end
    end
    %% Step 4: Descent Direction Generation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Solve for d in the linear system Bd = -a_bar.
    % If d>=0, then the LP is unbounded, STOP, else go to STEP5.
    d=-B\a_bar;        %solve Bd = -a_bar
    d_flag=find(d<0);
    if isempty(d_flag) %if d>=0, unbounded
        disp('LP is unbounded below!')
        fval=-inf;
        exit_flag=-3;
        return
    else %else Go to STEP5
        %% Step 5: Step Length Generation
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Compute the step length alpha by minimum ratio test.
        % Let l* be the index of the basic variables then attains
        % the minimum ratio alpha. Go to STEP6.
        alpha=min(x_B(d_flag)./abs(d(d_flag))); %minimum ratio test
        %% Step 6: Update Basic Variables
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Let x_B := x_B + alpha * d. Go to STEP7.
        x_B=x_B + alpha * d;      %get new basic variables
        delta=1e-30;            %computation error tolerance
        leave=find(abs(x_B)<=delta);   %index of leave variable
        while isempty(leave)
            delta = 0.1* delta;
            leave = find(abs(x_B)<=delta);
        end
        x_B(leave(1))=alpha;
        x_Bflag(leave(1))=t;
        if leave(1) <= length(s)   %update f_S and extreme point
            f_s(leave(1))=sub.c{t}'*x_sub{iter_num+1,t};
            v_L{leave(1),t}=x_sub{iter_num+1,t};
        else
            f_sub(leave(1)-length(s))=sub.c{t}'*x_sub{iter_num+1,t};
            v_sub{leave(1)-length(s)}=x_sub{iter_num+1,t};
        end
        %% Step7: Basis Update
        %%%%%%%%%%%%%%%%%%%%%%
        % Let B_l* be the column associated with the leaving variable x_l*.
        % Update the basis B by removing B_l* and adding the column
        % a_bar, and update B_set.
        B(:,leave(1))=a_bar; %update the basis B
    end 
    iter_num=iter_num+1;
    f_B=[f_s;f_sub];     %update f_B for next iteratin
    bound(:,iter_num)=[f_B'*x_B + sum(sig2); f_B'*x_B]; %new lower/upper bound
end %Go to STEP1