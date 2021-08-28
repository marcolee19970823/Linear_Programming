% created by chongnanli
% dantzig-wolfe decomposition for solving lp in block angular
clc;

% c=[-2;-3;-5;-4];
% A=[1,1,2,0;
%    0,1,1,1;
%    2,1,0,0;
%    1,1,0,0;
%    0,0,1,1;
%    0,0,3,2];
% b=[4;3;4;2;2;5];
% 
% mast.L{1}=A(1:2,1:2);
% mast.L{2}=A(1:2,3:4);
% mast.b=b(1:2);
% 
% sub.c{1}=c(1:2);
% sub.A{1}=A(3:4,1:2);
% sub.b{1}=b(3:4);
% 
% sub.c{2}=c(3:4);
% sub.A{2}=A(5:6,3:4);
% sub.b{2}=b(5:6);
% 
% sub.v{1}=zeros(length(sub.c{1}),1);
% sub.v{2}=zeros(length(sub.c{2}),1);
% 
% K=2;
% 
% 
% [x,fval,bound,exit_flag] = DantzigWolfeDecomp(mast,sub,K)

%%%%%

% c=[-2;-3;-2];
% A=[1,1,1;
%    -1,1,0;
%    -1,2,0;
%    0,0,1];
% b=[12;2;8;1];
% 
% mast.L{1}=A(1,1:2);
% mast.L{2}=A(1,3);
% mast.b=b(1);
% 
% sub.c{1}=c(1:2);
% sub.A{1}=A(2:3,1:2);
% sub.b{1}=b(2:3);
% 
% sub.c{2}=c(3);
% sub.A{2}=A(4,3);
% sub.b{2}=b(4);
% 
% sub.v{1}=zeros(length(sub.c{1}),1);
% sub.v{2}=zeros(length(sub.c{2}),1);
% 
% K=2;
% 
% [x,fval,bound,exit_flag] = DantzigWolfeDecomp(mast,sub,K)

%%%%%

c=[-2;-3;-5;-4];
A=[1,2,3,0;
   0,1,1,1;
   3,2,0,0;
   1,1,0,0;
   0,0,1,1;
   0,0,4,3];
b=[8;7;6;3;3;6];

mast.L{1}=A(1:2,1:2);
mast.L{2}=A(1:2,3:4);
mast.b=b(1:2);

sub.c{1}=c(1:2);
sub.A{1}=A(3:4,1:2);
sub.b{1}=b(3:4);

sub.c{2}=c(3:4);
sub.A{2}=A(5:6,3:4);
sub.b{2}=b(5:6);

sub.v{1}=zeros(length(sub.c{1}),1);
sub.v{2}=zeros(length(sub.c{2}),1);

K=2;

[x,fval,bound,exit_flag] = DantzigWolfeDecomp(mast,sub,K)

%%%%%

%¡ý
% c=[-2;-3;-5;-4];
% A=[1,2,3,0;
%    0,1,1,1;
%    3,2,0,0;
%    1,1,0,0;
%    0,0,1,1;
%    0,0,4,3];
% b=[8;7;6;3;3;6];
% 
% mast.L{1}=A(1:2,1:2);
% mast.L{2}=A(1:2,3:4);
% mast.b=b(1:2);
% 
% sub.c{1}=c(1:2);
% sub.A{1}=A(3:4,1:2);
% sub.b{1}=b(3:4);
% 
% sub.c{2}=c(3:4);
% sub.A{2}=A(5:6,3:4);
% sub.b{2}=b(5:6);
% 
% sub.v{1}=zeros(length(sub.c{1}),1);
% sub.v{2}=zeros(length(sub.c{2}),1);
% 
% K=2;
% 
% [x,fval,bound,exit_flag] = DantzigWolfeDecomp(mast,sub,K)
