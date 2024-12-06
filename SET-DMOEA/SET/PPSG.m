function [priorPS,negatives] = PPSG(Problem,popSize,K,t,P)
%输入：当前问题，种群大小，近邻值，环境变量，当前种群
%输出：先验PS（正样本集），负样本集
d = size(Problem.XLow,1);    %决策变量维度 

if nargin == 4
    RandomP=generateRandomPoints(popSize,Problem);
else
    P = unique(P','rows','stable')';
    RandomP=[generateRandomPoints(popSize,Problem) P];
end

F = [];
for i=1:size(RandomP,2)
    F(:,i) = Problem.FObj(RandomP(:,i)',t);
end
[pRank,~] = fastNonDominatedSort(F');
[~,indxp] = sort(pRank);
ind_NDS = size(find(pRank==1),2);
if(ind_NDS<=popSize/10)
    positives = RandomP(:,indxp(1:popSize/10));
    negatives = RandomP(:,indxp(popSize/10+1:end));
else
    positives = RandomP(:,indxp(1:ind_NDS));
    negatives = RandomP(:,indxp(ind_NDS+1:end));
end

priorPS = positives;
while(size(priorPS,2)<popSize)
    ttyh=pdist2(priorPS',priorPS');
    for i =1:size(priorPS,2)
        %得到k近邻
        [~,Kd] = sort(ttyh(i,:));
        B(i,:) = Kd(2:K+1);
        %计算局部密度
        local_d(i) = mean(pdist2(priorPS(:,i)',priorPS(:,B(i,:))'));
    end
    hreshold_p = mean(local_d);
    SampleP = [];
    for i =1:size(priorPS,2)
        if(local_d(i)>hreshold_p)
            %随机选择一个邻居xp
            xp = B(i,randperm(K,1));
            for mm = 1:d
                xnew(mm)=priorPS(mm,i)+rand*(priorPS(mm,i)-priorPS(mm,xp));
                xnew(mm) = min(max(xnew(mm),Problem.XLow(mm)), Problem.XUpp(mm));
            end
            SampleP = [SampleP xnew'];
        end
    end
    priorPS = [priorPS SampleP];
end
end