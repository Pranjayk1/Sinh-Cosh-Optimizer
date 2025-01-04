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

function X = initialization2(SearchAgents_no, dim, ub, lb)

Boundary_no = size(ub, 2);

% Preallocate the solution matrix
X = zeros(SearchAgents_no, dim);

% If boundaries of all variables are equal (uniform bounds across dimensions)
if Boundary_no == 1
    % Calculate number of intervals based on SearchAgents_no
    intervals = round(SearchAgents_no^(1/dim)); % Estimate intervals per dimension
    linspace_values = linspace(lb, ub, intervals); % Generate evenly spaced values

    % Create a grid of points across the search space
    [grid_points{1:dim}] = ndgrid(linspace_values);
    grid_points = cellfun(@(x) x(:), grid_points, 'UniformOutput', false); % Convert to vector format
    X = [grid_points{:}];

    % If there are more grid points than agents, randomly select a subset
    if size(X, 1) > SearchAgents_no
        selected_indices = randperm(size(X, 1), SearchAgents_no);
        X = X(selected_indices, :);
    end

    % If there are fewer grid points than agents, pad with random samples
    if size(X, 1) < SearchAgents_no
        additional_points = rand(SearchAgents_no - size(X, 1), dim) .* (ub - lb) + lb;
        X = [X; additional_points];
    end
end

% If each variable has different lb and ub (non-uniform bounds)
if Boundary_no > 1
    intervals = round(SearchAgents_no^(1/dim)); % Estimate intervals per dimension
    for i = 1:dim
        ub_i = ub(i);
        lb_i = lb(i);
        linspace_values = linspace(lb_i, ub_i, intervals); % Evenly spaced values per dimension
        
        % Generate grid for each dimension separately
        grid_points = linspace_values(mod((1:SearchAgents_no) - 1, intervals) + 1);
        
        % Assign each agent a grid point, looping back as necessary
        X(:, i) = grid_points(randperm(SearchAgents_no)); % Random permutation for diversity
    end
end
