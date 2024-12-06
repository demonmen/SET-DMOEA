clc
clear
close all
addpath(genpath('../SET-DMOEA'));
warning('off')
con = configure();
%函数加载
functions = con.TestFunctions;
%动态配置加载
T_parameter=con.T_parameter;
%种群大小
popSize=con.popSize;
%独立测试
for rep=1:20
    filename1 = ['result/SET/MIGD-', num2str(rep), '.txt'];
    fid1 = fopen(filename1,'w');
    filename2 = ['result/SET/MHV-', num2str(rep), '.txt'];
    fid2 = fopen(filename2,'w');
    for testFuncNo=1:size(functions,2)
        Problem=TestFunctions(functions{testFuncNo});
        if Problem.NObj==3
            popSize=150;
        end 
        for group=1:size(T_parameter,1) 
            MaxIt=T_parameter(group,2);%taut
            fprintf('\n SET-DMOEA dec:%d runing on: %s, configure: %d, environment:',con.dec,Problem.Name,group);
            reskt=SET_DMOEA(Problem,popSize,MaxIt,T_parameter,group,3); 
            [resIGD,resHV]=computeMetrics(reskt);
            fprintf(fid1,'%f \n',resIGD);
            fprintf(fid2,'%f \n',resHV);
        end
    end
    fclose(fid1);
    fclose(fid2);
end%rep

function [resIGD,resHV]=computeMetrics(resStruct)
     IGD_T = []; 
     HV_T = [];
     for T=1:size(resStruct,2)
        POFIter=resStruct{T}.POF_iter;
        POFbenchmark=resStruct{T}.turePOF;
        pof=POFIter{size(POFIter,2)};
        pof(imag(pof)~=0) = abs(pof(imag(pof)~=0));
        IGD_T(T)=IGD(pof',POFbenchmark);
        HV_T(T)=HV(pof',POFbenchmark);
     end
     resIGD=mean(IGD_T);
     resHV=mean(HV_T);
end