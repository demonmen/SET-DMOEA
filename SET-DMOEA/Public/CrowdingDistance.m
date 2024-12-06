function CrowdDis = CrowdingDistance(PopObj,FrontNo)
if(nargin==1)
    % Calculate the crowding distance of each solution in the same front
    % 计算同一前沿的每个解的拥挤距离
        [N,M]    = size(PopObj); 
        CrowdDis = zeros(1,N);
        Fmax     = max(PopObj,[],1);
        Fmin     = min(PopObj,[],1);
        for i = 1 : M
            [~,rank] = sortrows(PopObj(:,i));
            CrowdDis(rank(1))   = inf;
            CrowdDis(rank(end)) = inf;
            for j = 2 : N-1
                CrowdDis(rank(j)) = CrowdDis(rank(j))+(PopObj(rank(j+1),i)-PopObj(rank(j-1),i))/(Fmax(i)-Fmin(i));
            end
        end
else
    % Calculate the crowding distance of each solution front by front

    [N,M]    = size(PopObj);
    CrowdDis = zeros(1,N);
    Fronts   = setdiff(unique(FrontNo),inf);
    for f = 1 : length(Fronts)
        Front = find(FrontNo==Fronts(f));
        Fmax  = max(PopObj(Front,:),[],1);
        Fmin  = min(PopObj(Front,:),[],1);
        for i = 1 : M
            [~,Rank] = sortrows(PopObj(Front,i));
            CrowdDis(Front(Rank(1)))   = inf;
            CrowdDis(Front(Rank(end))) = inf;
            for j = 2 : length(Front)-1
                CrowdDis(Front(Rank(j))) = CrowdDis(Front(Rank(j)))+(PopObj(Front(Rank(j+1)),i)-PopObj(Front(Rank(j-1)),i))/(Fmax(i)-Fmin(i));
            end
        end
    end
end