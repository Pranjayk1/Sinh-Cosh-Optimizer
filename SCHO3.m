function [Destination_fitness, Destination_position, Convergence_curve] = SCHO(N, Max_iteration, lb, ub, dim, fobj)

    % Initialize variables
    Destination_position = zeros(1, dim);
    Destination_fitness = inf;
    Destination_position_second = zeros(1, dim);
    Convergence_curve = zeros(1, Max_iteration);
    Position_sort = zeros(N, dim);
    
    % Initialize SCHO parameters
    u = 0.388;
    m = 0.45;
    n = 0.5;
    p = 10;
    q = 9;
    Alpha = 4.6;
    Beta = 1.55;
    BS = floor(Max_iteration / Beta);
    ct = 3.6;
    T = floor(Max_iteration / ct);
    BSi = 0;
    BSi_temp = 0;
    ub_2 = ub;
    lb_2 = lb;
    
    % Initialize the set of random solutions
    X = initialization(N, dim, ub, lb);
    Objective_values = zeros(1, size(X, 1));
    
    % Calculate the fitness of the first set and find the best one
    for i = 1:size(X, 1)
        Objective_values(1, i) = fobj(X(i, :));
        if Objective_values(1, i) < Destination_fitness
            Destination_position = X(i, :);
            Destination_fitness = Objective_values(1, i);
        end
    end
    Convergence_curve(1) = Destination_fitness;
    t = 2;
    
    % Main loop
    while t <= Max_iteration
        for i = 1:size(X, 1) % in i-th solution
            for j = 1:size(X, 2) % in j-th dimension
                % Update A by using Eq. (17)
                cosh2 = (exp(t / Max_iteration) + exp(-t / Max_iteration)) / 2;
                sinh2 = (exp(t / Max_iteration) - exp(-t / Max_iteration)) / 2;
                r1 = rand();
                A = (p - q * (t / Max_iteration) ^ (cosh2 / (sinh2))) * r1;
                
                % Enhanced bounded search strategy with clustering
                if t == BSi
                    % Perform clustering with `kmeans`
                    cluster_count = min(5, N);  % Use up to 5 clusters
                    [idx, ~] = kmeans(X, cluster_count);  % Cluster the population

                    % Apply bounding strategy to each cluster
                    for cluster = 1:cluster_count
                        cluster_points = X(idx == cluster, :);  % Points in the current cluster
                        
                        if isempty(cluster_points)
                            continue; % Skip empty clusters
                        end
                        
                        % Determine bounding box for each cluster
                        lower_bound = min(cluster_points, [], 1);
                        upper_bound = max(cluster_points, [], 1);
                        
                        % If the cluster is dense, narrow the bounds
                        if size(cluster_points, 1) > N / cluster_count
                            lower_bound = lower_bound + 0.1 * (upper_bound - lower_bound);
                            upper_bound = upper_bound - 0.1 * (upper_bound - lower_bound);
                        end
                        
                        % Reinitialize solutions within these updated bounds
                        X(idx == cluster, :) = lower_bound + ...
                            (upper_bound - lower_bound) .* rand(sum(idx == cluster), dim);
                    end
                    
                    % Reset bounded search parameters
                    BSi_temp = BSi;
                    BSi = 0;
                end
                
                % First phase of exploration and exploitation
                if t <= T
                    r2 = rand();
                    r3 = rand();
                    a1 = 3 * (-1.3 * t / Max_iteration + m);
                    r4 = rand();
                    r5 = rand();
                    if A > 1
                        sinh = (exp(r3) - exp(-r3)) / 2;
                        cosh = (exp(r3) + exp(-r3)) / 2;
                        W1 = r2 * a1 * (cosh + u * sinh - 1);
                        if r5 <= 0.5
                            X(i, j) = Destination_position(j) + r4 * W1 * X(i, j);
                        else
                            X(i, j) = Destination_position(j) - r4 * W1 * X(i, j);
                        end
                    else
                        sinh = (exp(r3) - exp(-r3)) / 2;
                        cosh = (exp(r3) + exp(-r3)) / 2;
                        W3 = r2 * a1 * (cosh + u * sinh);
                        if r5 <= 0.5
                            X(i, j) = Destination_position(j) + r4 * W3 * X(i, j);
                        else
                            X(i, j) = Destination_position(j) - r4 * W3 * X(i, j);
                        end
                    end
                else
                    % Second phase of exploration and exploitation
                    r2 = rand();
                    r3 = rand();
                    a2 = 2 * (-t / Max_iteration + n);
                    W2 = r2 * a2;
                    r4 = rand();
                    r5 = rand();
                    if A < 1
                        sinh = (exp(r3) - exp(-r3)) / 2;
                        cosh = (exp(r3) + exp(-r3)) / 2;
                        X(i, j) = X(i, j) + (r5 * sinh / cosh * abs(W2 * Destination_position(j) - X(i, j)));
                    else
                        if r4 <= 0.5
                            X(i, j) = X(i, j) + (abs(0.003 * W2 * Destination_position(j) - X(i, j)));
                        else
                            X(i, j) = X(i, j) - (abs(0.003 * W2 * Destination_position(j) - X(i, j)));
                        end
                    end
                end
            end
            BSi = BSi_temp;
        end
        
        % Update solutions based on boundaries and objective values
        for i = 1:size(X, 1)
            % Ensure solutions are within bounds
            Flag4ub = X(i, :) > ub_2;
            Flag4lb = X(i, :) < lb_2;
            X(i, :) = (X(i, :) .* ~(Flag4ub + Flag4lb)) + (ub_2 + lb_2) / 2 .* Flag4ub + lb_2 .* Flag4lb;
            
            % Calculate objective values
            Objective_values(1, i) = fobj(X(i, :));
            
            % Update the best solution if found
            if Objective_values(1, i) < Destination_fitness
                Destination_position = X(i, :);
                Destination_fitness = Objective_values(1, i);
            end
        end
        
        % Find the second-best solution periodically
        if t == BS
            BSi = BS + 1;
            BS = BS + floor((Max_iteration - BS) / Alpha);
            Position_sort = sortrows([Objective_values' X], 1); % Sort solutions by fitness
            Destination_position_second = Position_sort(2, 2:end); % Second-best solution
        end
        
        % Update convergence curve and iteration counter
        Convergence_curve(t) = Destination_fitness;
        t = t + 1;
    end
end

% Helper function for initialization
function X = initialization(N, dim, ub, lb)
    X = lb + (ub - lb) .* rand(N, dim);
end
