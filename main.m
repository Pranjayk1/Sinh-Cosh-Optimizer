%  Sinh Cosh Optimizer (SCHO)                                                                     
%                                                                                                     
%  Developed in MATLAB R2022a                                                                  
%                                                                                                     
%  programming: Jianfu Bai                                                          
%                                                                                                     
%  e-Mail: Jianfu.Bai@UGent.be, magd.abdelwahab@ugent.be                                                               
%  Soete Laboratory, Department of Electrical Energy, Metals, Mechanical Constructions, and Systems, 
%  Faculty of Engineering and Architecture, Ghent University, Belgium                                                           
%                                                                                                                                                                                                                                                              
%  paper: Jianfu Bai, Yifei Li, Mingpo Zheng, Samir Khatir, Brahim Benaisa, Laith Abualigah, Magd Abdel Wahab, A Sinh Cosh Optimizer, Knowledge-Based Systems (2023).                                                                                      

clear all
clc


maxFunc = 30;
SearchAgents_no = 50;
Max_iteration= 1000; 
runs = 30;

for fn = 3:maxFunc
    
    Function_name=strcat('F',num2str(fn));
    [lb,ub,dim,fobj]=CEC2017(Function_name);    
    Best_score_T = zeros(runs,1);
    AvgConvCurve = zeros(1, Max_iteration);
    Convergence_curve=zeros(runs,Max_iteration);
    for run=1:runs
        [Best_score,Best_pos,cg_curve]=SCHO15(SearchAgents_no,Max_iteration,lb,ub,dim,fobj);
        Best_score_T(run) = Best_score;
        Convergence_curve(run,:)=cg_curve;
    end

    Best_score_Best = min(Best_score_T);
    Best_Score_Mean = mean(Best_score_T);
    Best_Score_std = std(Best_score_T);
    format long
    display([Function_name, ' Best:  ', num2str(Best_score_Best), '     ', 'Mean:  ', num2str(Best_Score_Mean), '     ', 'Std. Deviation:  ', num2str(Best_Score_std)]);
end
