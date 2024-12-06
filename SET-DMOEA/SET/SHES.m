function similarE = SHES(Historical,priorPS)
%输入：历史环境PS，先验PS
%输出：相似环境索引
E_mmd=[];
for i=1:length(Historical)-1
    E_mmd(i) = mmd_2(Historical{i}.positives',priorPS','rbf');
end
[~,similarE] = min(E_mmd);
%similarE = T-1;

end