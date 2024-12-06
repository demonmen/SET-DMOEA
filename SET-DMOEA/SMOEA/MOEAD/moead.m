function [PopX,Pareto,POF_iter]=moead( Problem,popSize,MaxIt, t,init_pop)

CostFunction=Problem.FObj;  % Cost Function
nVar=size(Problem.XLow,1);  % Number of Decision Variables
VarSize=[nVar 1];           % Decision Variables Matrix Size 决策变量矩阵大小
VarMin = Problem.XLow;      % Decision Variables Lower Bound
VarMax = Problem.XUpp;      % Decision Variables Upper Bound
nObj=Problem.NObj;

%% MOEA/D Settings


nPop=popSize;               % Population Size (Number of Sub-Problems) 子问题数量
nArchive=popSize;

%T = ceil(popSize/10);
T=max(ceil(0.15*nPop),2);   % Number of Neighbors  邻居数量
T=min(max(T,2),15);

GenOperator_params.gamma=0.8;
GenOperator_params.VarMin=VarMin;
GenOperator_params.VarMax=VarMax;

%% Initialization

% Create Sub-problems
sp=CreateSubProblems(nObj,nPop,T);

% Empty Individual
empty_individual.Position=[];
empty_individual.Cost=[];
empty_individual.g=[];
empty_individual.IsDominated=[];

% Initialize Goal Point
%z=inf(nObj,1);
z=zeros(nObj,1);

% Create Initial Population
pop=repmat(empty_individual,nPop,1);

if nargin == 4
    for i=1:nPop   
        pop(i).Position=unifrnd(VarMin,VarMax,VarSize);
        pop(i).Cost=CostFunction(pop(i).Position',t);
        z=min(z,pop(i).Cost);
    end
elseif nargin == 5
    for i=1:size(init_pop,2)
        pop(i).Position=init_pop(:,i);
        pop(i).Cost=CostFunction(pop(i).Position',t);
        z=min(z,pop(i).Cost);
    end
    for i=size(init_pop,2)+1:nPop
        pop(i).Position=unifrnd(VarMin,VarMax,VarSize);
        pop(i).Cost=CostFunction(pop(i).Position',t);
        z=min(z,pop(i).Cost);
    end
end

% for i=1:nPop
%     pop(i).g=DecomposedCost(pop(i),z,sp(i).lambda);
% end

% Determine Population Domination Status确定人口主导地位
pop=DetermineDomination(pop);

% Initialize Estimated Pareto Front
EP=pop(~[pop.IsDominated]);

%% Main Loop

for it=1:MaxIt%%迭代开始，循环进化种群
    for i=1:nPop
        % Reproduction (Crossover)
        K=randsample(T,2);
        j1=sp(i).Neighbors(K(1));
        p1=pop(j1);
        
        j2=sp(i).Neighbors(K(2));
        p2=pop(j2);
        
        y=empty_individual;
        
        %y.Position=M_Crossover(p1.Position,p2.Position,GenOperator_params);
        y.Position = GeneticOperator(p1.Position,p2.Position,GenOperator_params);%!!!SBX
        y.Cost=CostFunction(y.Position',t);
       
        z=min(z,y.Cost);
        %% 更新邻居
        for j=sp(i).Neighbors
            pop(j).g=DecomposedCost(pop(j),z,sp(j).lambda);
            y.g=DecomposedCost(y,z,sp(j).lambda);
            if y.g<=pop(j).g
                pop(j)=y;
            end
        end
    end
    
    %Determine Population Domination Status
	pop=DetermineDomination(pop);
    ndpop=pop(~[pop.IsDominated]);
    EP=[EP
        ndpop]; %#ok
    %% 重新对POS个体比较，去除支配解
    EP=DetermineDomination(EP);
    EP=EP(~[EP.IsDominated]);

%     % 去除POS中重复元素
%     [~,m,~]=unique([EP.Position]','rows');
%     EP = EP(m);
    %% 若POS数量大于种群数量，随机抽取多余个体，将种群缩小到popsize
    if numel(EP)>nArchive
        Extra=numel(EP)-nArchive;
        ToBeDeleted=randsample(numel(EP),Extra);
        EP(ToBeDeleted)=[];
    end
    %% 保存此次迭代获得的POF
    for arcnum=1: size(EP,1)
        pareto(:,arcnum)=EP(arcnum).Cost;    
    end
    POF_iter{it}=pareto;
end

Pareto.F=[EP.Cost];
Pareto.X=[EP.Position];
PopX = zeros(size(Problem.XLow,1),popSize);
for i=1:popSize
   PopX(:,i) = pop(i).Position;
end

end
