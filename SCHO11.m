function [Destination_fitness, Destination_position, Convergence_curve] = SCHO11(N, Max_iteration, lb, ub, dim, fobj, function_number)

    % Set the target optimum value based on function number
    optimum_value = function_number * 100;

    % Initialize variables
    Destination_position = zeros(1, dim);
    Destination_fitness = inf;
    Convergence_curve = zeros(1, Max_iteration);

    % Fixed SCHO parameters for initial setup
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

    % Initialize the set of random solutions
    X = initialization(N, dim, ub, lb);
    Objective_values = zeros(1, N);

    % Calculate the fitness of the first set and find the best one
    for i = 1:N
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
        for i = 1:N  % in i-th solution
            for j = 1:dim  % in j-th dimension
                % Update A using fixed parameters
                cosh2 = (exp(t / Max_iteration) + exp(-t / Max_iteration)) / 2;
                sinh2 = (exp(t / Max_iteration) - exp(-t / Max_iteration)) / 2;
                r1 = rand();
                A = (p - q * (t / Max_iteration) ^ (cosh2 / sinh2)) * r1;

                % Distance to optimum adjustment factor
                distance_to_optimum = abs(Destination_fitness - optimum_value) / optimum_value;
                adjustment_factor = 1 - distance_to_optimum;  % Guides movement closer to optimum

                % Exploration and exploitation phases with guidance toward optimum
                if t <= T
                    % First phase of exploration and exploitation
                    r2 = rand();
                    r3 = rand();
                    a1 = 3 * (-1.3 * t / Max_iteration + m);
                    r4 = rand();
                    r5 = rand();
                    if A > 1
                        sinh = (exp(r3) - exp(-r3)) / 2;
                        cosh = (exp(r3) + exp(-r3)) / 2;
                        W1 = r2 * a1 * (cosh + u * sinh - 1) * adjustment_factor;
                        if r5 <= 0.5
                            X(i, j) = Destination_position(j) + r4 * W1 * X(i, j);
                        else
                            X(i, j) = Destination_position(j) - r4 * W1 * X(i, j);
                        end
                    else
                        sinh = (exp(r3) - exp(-r3)) / 2;
                        cosh = (exp(r3) + exp(-r3)) / 2;
                        W3 = r2 * a1 * (cosh + u * sinh) * adjustment_factor;
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
                    W2 = r2 * a2 * adjustment_factor;
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
        end

        % Update solutions based on boundaries and objective values
        for i = 1:N
            % Ensure solutions are within bounds
            Flag4ub = X(i, :) > ub;
            Flag4lb = X(i, :) < lb;
            X(i, :) = (X(i, :) .* ~(Flag4ub + Flag4lb)) + ub .* Flag4ub + lb .* Flag4lb;
            
            % Calculate objective values
            Objective_values(1, i) = fobj(X(i, :));
            
            % Update the best solution if found
            if Objective_values(1, i) < Destination_fitness
                Destination_position = X(i, :);
                Destination_fitness = Objective_values(1, i);
            end
        end

        % Update convergence curve and iteration counter
        Convergence_curve(t) = Destination_fitness;
        t = t + 1;
    end
end

% Helper function for initialization
function X = initialization(SearchAgents_no, dim, ub, lb)

    Boundary_no = size(ub, 2); 

    % If the boundaries of all variables are equal
    if Boundary_no == 1
        X = rand(SearchAgents_no, dim) .* (ub - lb) + lb;  
    end

    % If each variable has a different lb and ub
    if Boundary_no > 1
        for i = 1:dim
            ub_i = ub(i);
            lb_i = lb(i);
            X(:, i) = rand(SearchAgents_no, 1) .* (ub_i - lb_i) + lb_i;
        end
    end
end