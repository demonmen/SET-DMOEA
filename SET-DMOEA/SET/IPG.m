function initPop = IPG(Problem,popSize,priorPS,similarE,Historical,t)
%输入：当前问题，种群大小，先验PS，相似环境索引，历史PS，历史环境的负样本，环境变量
%输出：新环境的初始种群

%构建训练集测试集
trainz = [Historical{similarE}.positives' ones(size(Historical{similarE}.positives,2),1)];
trainf = [Historical{similarE}.negatives' 2*ones(size(Historical{similarE}.negatives,2),1)];
trainPop=[trainz;trainf];
%matlab 存在bug, 需要打乱一点数据才会在迭代的时候生成随机，否则有可能生成一样的随机解
%testSize = randi([24*popSize 30*popSize],1) + randi([3 3*popSize],1)*randi([1 5],1);
testSize = randi([2*popSize 5*popSize],1) + randi([3 3*popSize],1)*randi([1 5],1);
testPop = generateRandomPoints(testSize,Problem)';
%模型训练及分类
[~,y_pred] = EasyTL_kernel(trainPop(:,1:end-1),trainPop(:,end),testPop,2*ones(size(testPop,1),1));
ind_z = find(y_pred==1);
if size(ind_z,2)>popSize
    ind_z = ind_z(1:popSize);
end
initPop = [priorPS testPop(ind_z,:)'];
F=[];
for i=1:size(initPop,2)
    F(:,i) = Problem.FObj(initPop(:,i)',t);
end
[pRank,~] = fastNonDominatedSort(F');
[~,indxp] = sort(pRank);
initPop = initPop(:,indxp(1:popSize));


end