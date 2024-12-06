function POF_Banchmark = getBenchmarkPOF(testfunname,group,T,T_parameter )
if strcmp(testfunname, 'D_PID')
   POF_Banchmark = ones(30,3);
   POF_Banchmark(:,1) = POF_Banchmark(:,1)*10;
   POF_Banchmark(:,2) = POF_Banchmark(:,2)*100;
   POF_Banchmark(:,3) = POF_Banchmark(:,3)*100;
else
%UNTITLED2 Summary of this function goes here
    tempPosition = T;
   % ['./Metrics/pof/measures/pof/' 'POF-nt' num2str(T_parameter(group,1)) '-taut' num2str(T_parameter(group,2)) '-' functions{testfunc} '-' num2str(tempPosition) '.txt']
    POF_Banchmark = importdata(['./Benchmark/pof/' 'POF-nt' num2str(T_parameter(group,1)) '-taut' num2str(T_parameter(group,2)) '-' testfunname '-' num2str(tempPosition) '.txt']);
end
end
