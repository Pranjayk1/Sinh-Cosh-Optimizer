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
%  paper: Jianfu Bai, Yifei Li, Mingpo Zheng, Samir Khatir, Brahim Benaisa, Laith Abualigah, Magd Abdel Wahab, A Sinh Cosh Optimizer, Knowledge-Based Systems(2023).

% This function creates the first random population

function X=initialization(SearchAgents_no,dim,ub,lb)

Boundary_no= size(ub,2); 

% If the boundaries of all variables are equal

if Boundary_no==1
    X=rand(SearchAgents_no,dim).*(ub-lb)+lb;  
end

% If each variable has a different lb and ub
if Boundary_no>1
    for i=1:dim
        ub_i=ub(i);
        lb_i=lb(i);
        X(:,i)=rand(SearchAgents_no,1).*(ub_i-lb_i)+lb_i;
    end
end