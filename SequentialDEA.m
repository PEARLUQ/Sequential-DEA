function effsdea = SequentialDEA(Xk,Yk,X,Y,ori,rts,prdl,Company,District,Companyk,Districtk)

% ------------------------
% Function for sequential DEA estimation - code by Z. Wang (2021)
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

%---------------
% Primal
%---------------
if strcmp('P',prdl)
% Input orientation
	if strcmp('I',ori)
		for i = 1:Kk
            xktemp=[];
            xktemp = Xk (Companyk == uniqueKk(i),:);
            transtemp = [];
            transtempcum = [];
            for s = 1:height(xktemp)
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
	 		xeff=dea_i_p(xktemp,yktemp,XStemp,YStemp,rts);
         eff=[eff;xeff]; 
        end
    end
% output orientation
   if strcmp('O',ori)
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
            for s = 1:height(yktemp);
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
	 		xeff=dea_o_p(xktemp,yktemp,XStemp,YStemp,rts);
         eff=[eff;xeff]; 
        end
       end
end
%---------------
% Dual
%---------------
if strcmp('D',prdl)
% Input orientation
   if strcmp('I',ori)
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
            for s = 1:height(yktemp);
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
	 		xeff=dea_i_d(xktemp,yktemp,XStemp,YStemp,rts);
         eff=[eff;xeff]; 
        end
   end
% Output orientation
   if strcmp('O',ori)
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
                % Set up matrices by DMU, stacking S districts
                % S DEA models solved jointly
                XSt = [X((obsS(Kinsub(j))+1):obsS(Kinsub(j)+1),:) ; zeros(Maxd-countsS(Kinsub(j)),pp)];
                YSt = [Y((obsS(Kinsub(j))+1):obsS(Kinsub(j)+1),:) ; zeros(Maxd-countsS(Kinsub(j)),qq)];
                XStemp = [XStemp, XSt];
                YStemp = [YStemp, YSt];
            end

            % Set up nan for linprog in dea_o_d
            XStemp(:, isnan(xktemp)) = [];
            YStemp(:, isnan(yktemp)) = [];
            xktemp = xktemp(~any(isnan(xktemp), 2), :);
            yktemp = yktemp(~any(isnan(yktemp), 2), :);
            XStemp = fillmissing(XStemp,'constant',0);
            YStemp = fillmissing(YStemp,'constant',0);

	 		xeff=dea_o_d(xktemp,yktemp,XStemp,YStemp,rts,options);
         eff=[eff;xeff]; 
      end
	end
end

effsdea=eff;
return
