function [Destination_fitness, Destination_position, Convergence_curve] = SCHO8(N, Max_iteration, lb, ub, dim, fobj)

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

    % Initialize leaderboard (archive) for top solutions
    num_top_solutions = 5;  % Define the number of top solutions to store
    Top_solutions = Inf(num_top_solutions, dim);
    Top_fitness = Inf(1, num_top_solutions);

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
        % Update the leaderboard with the current population
        for i = 1:N
            fitness = fobj(X(i, :));
            [~, worst_idx] = max(Top_fitness);  % Find the worst solution in the leaderboard
            if fitness < Top_fitness(worst_idx)
                % Replace the worst solution with the current solution if it's better
                Top_solutions(worst_idx, :) = X(i, :);
                Top_fitness(worst_idx) = fitness;
            end
        end

        % Update candidate positions, influenced by leaderboard solutions
        for i = 1:N  % in i-th solution
            for j = 1:dim  % in j-th dimension
                % Update A using fixed parameters
                cosh2 = (exp(t / Max_iteration) + exp(-t / Max_iteration)) / 2;
                sinh2 = (exp(t / Max_iteration) - exp(-t / Max_iteration)) / 2;
                r1 = rand();
                A = (p - q * (t / Max_iteration) ^ (cosh2 / sinh2)) * r1;

                % Influence from leaderboard solutions
                influence = mean(Top_solutions(:, j)) - X(i, j);  % Calculate mean influence of top solutions
                influence_factor = 0.1;  % Scale influence factor as needed
                
                % Add influence from leaderboard solutions to candidate position update
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
                        W1 = r2 * a1 * (cosh + u * sinh - 1) + influence_factor * influence;
                        if r5 <= 0.5
                            X(i, j) = Destination_position(j) + r4 * W1 * X(i, j);
                        else
                            X(i, j) = Destination_position(j) - r4 * W1 * X(i, j);
                        end
                    else
                        sinh = (exp(r3) - exp(-r3)) / 2;
                        cosh = (exp(r3) + exp(-r3)) / 2;
                        W3 = r2 * a1 * (cosh + u * sinh) + influence_factor * influence;
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
                    W2 = r2 * a2 + influence_factor * influence;
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
function X = initialization(N, dim, ub, lb)
    X = lb + (ub - lb) .* rand(N, dim);
end
