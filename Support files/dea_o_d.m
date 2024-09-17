function [effk]=dea_o_d(xk,yk,X,Y,rts,options)

% =======
% dea_o_d
% =======
%
% Computation of effiency for an observation (xk,yk)
% in a DEA dual model for output orientation 
% 
%         Usage
%         [effk]=dea_o_d(xk,yk,X,Y,rts)
%         IN   :
%         ------
%            xk  : Vector of input for observation (xk,yk) (px1)
%            yk  : Vector of input for observation (xk,yk) (qx1)
%            X   : Matrix of input(s)  (n x p)
%            Y   : Matrix of output(s) (n x k)
%            rts : Assumption on returns to scale (text)
%                     'NIRS' = Non increasing rts
%                     'NDRS' = Non decreasing rts
%                     'CRS'  = Constant rts
%                     'VRS'  = Variable rts
%         OUT  :
%         ------
%	effk    : Efficiency scores of observation (xk,yk)
%
% Use : linprog
% Called by : dea
% Proposed by L. Simar
% Adapted by P. Vanden Eeckaut (June 1998), D. Mitrut (October 2001)Kai
% Du(December 2010)
% up-dated and checked by L. Simar (december 2002)
%
% INSTITUT DE STATISTIQUE - UNIVERSITE CATHOLIQUE DE LOUVAIN
%
% Output inefficiency measures of a new DMU xk (px1) and yk(qx1)
% with respect to a fixed other reference set
% determined by a given X (nxp) and Y (nxq)
% Need to input f,A,b
% for solving in x minf'x s.t. Ax<=b
% ==> call x = linprog(f,A,b)
%    -if constraints on x: xLB<= x <= xUB
% ==> call x=linprog(f,A,b,Aeq,Beq,xLB,xUB)
%    -if the first N constraints defined by A and b 
%     are equality constraints, must be writen in first N rows
% ==> call x=linprog(f,A,b,Aeq,Beq,xLB,xUB,[],N)


% -------------------
% Identify dimensions
% -------------------

[n,p] = size(X);       [n,q] = size(Y);
   aq = [];               bq = [];

% -----------------------------
% compute efficiency of new DMU
% -----------------------------

f=[-1.00;zeros(n,1)]; % Objective function

if strcmp('NIRS',rts) % NIRS model
  A=[0,ones(1,n);zeros(p,1),X';yk,-Y'];
  b=[1;xk;zeros(q,1)];
  lbx=zeros(p+q+1,1);
  ubx=ones(p+q+1,1)*inf;
  xeff=linprog(b,-A',f,aq,bq,lbx,ubx,[],options);
  effk=b'*xeff;
  return
elseif strcmp('NDRS',rts) % NDRS model
  A=[0,-ones(1,n);zeros(p,1),X';yk,-Y'];
  b=[-1;xk;zeros(q,1)];
  lbx=zeros(p+q+1,1);
  ubx=ones(p+q+1,1)*inf;
  xeff=linprog(b,-A',f,aq,bq,lbx,ubx,[],options);
  effk=b'*xeff;
  return
elseif strcmp('CRS',rts) % CRS model
  A=[zeros(p,1),X';yk,-Y'];
  b=[xk;zeros(q,1)];
  lbx=zeros(p+q,1);
  ubx=ones(p+q,1)*inf;
  xeff=linprog(b,-A',f,aq,bq,lbx,ubx,[],options);
  effk=b'*xeff;
  return
elseif strcmp('VRS',rts) % VRS model
  A=[0,ones(1,n);zeros(p,1),X';yk,-Y'];
  b=[1;xk;zeros(q,1)];
  lbx=[-inf;zeros(p+q,1)];
  ubx=ones(p+q+1,1)*inf;
  xeff=linprog(b,-A',f,aq,bq,lbx,ubx,[],options);
  effk=b'*xeff;
  return
else
  disp('Incorrect asssumption on RTS')
  return
end
