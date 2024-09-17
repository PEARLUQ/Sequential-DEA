function [effsdea,zksdea] = SequentialDEAls(Xk,Yk,X,Y,ori,rts,prdl,Company,District,Companyk,Districtk)

% ------------------------
% Function for sequential DEA decentralized for lambda^s - code by Z. Wang (2023)
% Prerequisite: Matrices Xk, Yk, X, Y need to be sorted by "District"
% ------------------------
% Following Färe & Zelenyuk (2021, AOR) and based on deaXnewYnew.m, which is acknowleged as below:
%   Proposed by L. Simar
%   Adapted by P. Vanden Eeckaut (June 1998), D. Mitrut (August 2001)
%   Readapted by L. SIMAR (23 november 2002)
%   INSTITUT DE STATISTIQUE - UNIVERSITE CATHOLIQUE DE LOUVAIN
% ------------------------
% Use with lp (with qpsubold) (or linprog), dea_i_p, dea_o_p, dea_i_d, dea_o_d
% ------------------------

% =======================================
% Sequential DEA
% =======================================

options = optimset('Display','off');

% ------------------------
% Prepare intermediate indicators
% ------------------------

uniqueKk = unique(Companyk);
Kk = length(uniqueKk);
Sk = length(unique(Districtk));
Maxc = max(histc(Company(:), unique(Companyk)));
Maxd = max(histc(District(:), unique(Districtk)));

Ksub = [];
for i = 1:Kk
    Ksubtemp =  [Districtk(find (Companyk == uniqueKk(i))) ; zeros((Maxc - length(find (Companyk == uniqueKk(i)))),1)*nan];
    Ksub = [Ksub , Ksubtemp];
end

[~,pp] = size(X);
[~,qq] = size(Y);

DistrictS = District;
uniqueS = unique(DistrictS);
countsS = histc(DistrictS(:), uniqueS);

SS = length(uniqueS);
obsS = [0];
for i = 1:SS
obsS (i+1) = sum (countsS(1:i,:));
end
obsS = obsS';

% ------------------------
% test for orientation
% ------------------------
[n,p] = size(X);
[~,q] = size(Y);

p = [strcmp('I', ori),strcmp('O', ori)];
if sum(p)==0
	disp( ' !!! ori option is not correct use I or O !!!');
	return
end
% test for return to scale
p = [strcmp('CRS', rts),strcmp('NIRS', rts),strcmp('NDRS', rts),strcmp('VRS', rts),strcmp('VRS2', rts)];
if sum(p)==0
	disp( ' !!! rts option is not correct use CRS, NIRS, NDRS or VRS !!!');
	return
end
% test for primal dual
p = [strcmp('P', prdl),strcmp('D', prdl)];
if sum(p)==0
	disp( ' !!! prdl option is not correct use P or D !!!');
	return
end

% ------------------------
% Estimate efficiency score of DMU
% ------------------------

eff=[];
zk=[];

% Currently for Dual and Output-orientation only:
		for i = 1:Kk
            xktemp=[];
            xktemp = Xk (Companyk == uniqueKk(i),:);
            transtemp = [];
            transtempcum = [];
            for s = 1:height(xktemp);
                transtemp = xktemp(s,:)';
                transtempcum = [transtempcum; transtemp];
            end
            xktemp = transtempcum;
            
            yktemp=[];
            yktemp = Yk (Companyk == uniqueKk(i),:);
            transtemp = [];
            transtempcum = [];
            for s = 1:height(yktemp)
                transtemp = yktemp(s,:)';
                transtempcum = [transtempcum; transtemp];
            end
            yktemp = transtempcum;
            
            XStemp = [];
            YStemp = [];
            Kinsub = Ksub(:,i);
            Kinsub(isnan(Kinsub))=[];
            
            for j = 1:length(Kinsub)
                XSt = [X((obsS(Kinsub(j))+1):obsS(Kinsub(j)+1),:) ; zeros(Maxd-countsS(Kinsub(j)),pp)];
                YSt = [Y((obsS(Kinsub(j))+1):obsS(Kinsub(j)+1),:) ; zeros(Maxd-countsS(Kinsub(j)),qq)];
                XStemp = [XStemp, XSt];
                YStemp = [YStemp, YSt];
            end

            % Set up nan for linprog
            XStemp(:, isnan(xktemp)) = [];
            YStemp(:, isnan(yktemp)) = [];
            xktemp = xktemp(~any(isnan(xktemp), 2), :);
            yktemp = yktemp(~any(isnan(yktemp), 2), :);
            XStemp = fillmissing(XStemp,'constant',0);
            YStemp = fillmissing(YStemp,'constant',0);
            
            % Specific A,B,f for Output-oriented VRS:
            [n,p] = size(XStemp);
            [~,q] = size(YStemp);
            aq=[];
            bq=[];
             Skk = length(yktemp);
             f=[repmat(-1,Skk,1);zeros(n,1)];
             
             % Specific matrices for linprog
             % Lambda extended by S but linked with the same zk
             if strcmp('VRS',rts)
             A=[repmat(0,1,Skk),ones(1,n);repmat(0,p,Skk),XStemp';diag(yktemp),-YStemp'];                                                            
             b=[1;xktemp;zeros(q,1)]; 

             elseif strcmp('CRS',rts)
             A=[repmat(0,p,Skk),XStemp';diag(yktemp),-YStemp'];                         
             b=[xktemp;zeros(q,1)]; 
             else
             disp('Unsupported RTS')
             end

             lbx=[repmat(-inf,Skk,1);zeros(n,1)]; 
             ubx=ones(n+Skk,1)*inf;
             xeff=linprog(f,A,b,aq,bq,lbx,ubx);
             effk=xeff(1:Skk); 
             zkt = xeff(Skk+1:end);
         eff=[eff;effk]; 
         zk=[zk;zkt];
        end
effsdea=eff;
zksdea=zk;
return
