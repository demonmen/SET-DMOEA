function con=configure()
con.T_parameter = [
    10 10 500
    10 5 250
    5 10 500
                    ];%% time parameters   nt tauT tau  
con.TestFunctions = {'DF1'};
%con.TestFunctions = {'FDA1','FDA2','FDA3','FDA4','FDA5','dMOP1','dMOP2','dMOP3','UDF1','UDF2','UDF3','UDF4','UDF5','UDF6','UDF7'};
%con.TestFunctions = {'DF1','DF2','DF3','DF4','DF5','DF6','DF7','DF8','DF9','DF10','DF11','DF12','DF13','DF14'};
con.popSize=100;%种群大小
con.repeat=20;%独立运行次数
con.dec=10; %决策变量维度
end