function res=SET_DMOEA(Problem,popSize,MaxIt,T_parameter,group,K)
Historical = [];
for T = 1:T_parameter(group,3)/T_parameter(group,2)
    t= 1/T_parameter(group,1)*(T-1);   
    fprintf(' %d',T);%所处环境变化时刻
    if T==1 
        [priorPS,Historical{T}.negatives] = PPSG(Problem,popSize,K,t);
        initPop = priorPS;
        [PopX,Pareto,POF_iter]= moead(Problem,popSize,MaxIt,t,initPop);
        Historical{T}.positives = Pareto.X;
    else
        %预测方法开始
        % I.生成先验PS
        [priorPS,Historical{T}.negatives] = PPSG(Problem,popSize,K,t,Historical{T-1}.positives);
        % II.相似环境选择
        similarE = SHES(Historical,priorPS);
        % III.生成初始种群
        initPop = IPG(Problem,popSize,priorPS,similarE,Historical,t);
        % MOEA
        [PopX,Pareto,POF_iter]=moead(Problem,popSize,MaxIt,t,initPop); 
        Historical{T}.positives = Pareto.X;
    end
    res{T}.POF_iter=POF_iter;%#ok
    %T时刻预测的帕累托解
    res{T}.POS=Pareto.X;%#ok
    %T时刻对应的真实帕累托前沿
    res{T}.turePOF=getBenchmarkPOF(Problem.Name,group,T,T_parameter);  %#ok
end
