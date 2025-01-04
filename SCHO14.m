function [Destination_fitness, Destination_position, Convergence_curve] = SCHO14(N, Max_iteration, lb, ub, dim, fobj)
    % Initialize variables
    Destination_position = zeros(1, dim);
    Destination_fitness = inf;
    Convergence_curve = zeros(1, Max_iteration);
    
    % Fixed SCHO parameters
    u = 0.388;
    m = 0.45;
    n = 0.5;
    p = 10;
    q = 9;
    Beta = 1.55;
    
    % Dynamic ct and Alpha adjustment parameters
    initial_ct = 3.6;
    final_ct = 1.2;
    initial_Alpha = 4.6;
    final_Alpha = 1.2;

    % Opposition-Based Initialization
    X = initialization(N, dim, ub, lb);
    opposite_X = ub + lb - X;

    Objective_values = zeros(1, N);
    for i = 1:N
        fitness_original = fobj(X(i, :));
        fitness_opposite = fobj(opposite_X(i, :));
        
        if fitness_opposite < fitness_original
            X(i, :) = opposite_X(i, :);
            Objective_values(1, i) = fitness_opposite;
        else
            Objective_values(1, i) = fitness_original;
        end

        if Objective_values(1, i) < Destination_fitness
            Destination_position = X(i, :);
            Destination_fitness = Objective_values(1, i);
        end
    end
    Convergence_curve(1) = Destination_fitness;
    t = 2;

    % Main loop
    while t <= Max_iteration
        % Dynamic adjustment of ct and Alpha
        ct = initial_ct + (final_ct - initial_ct) * (1 - t / Max_iteration);
        Alpha = initial_Alpha * (1 - t / Max_iteration) + final_Alpha * (t / Max_iteration);
        
        % Adaptive boundary adjustment
        shrink_factor = 1 - t / Max_iteration;
        adaptive_lb = lb + (1 - shrink_factor) * (Destination_position - lb);
        adaptive_ub = ub - (1 - shrink_factor) * (ub - Destination_position);

        for i = 1:N
            for j = 1:dim
                cosh2 = (exp(t / Max_iteration) + exp(-t / Max_iteration)) / 2;
                sinh2 = (exp(t / Max_iteration) - exp(-t / Max_iteration)) / 2;
                r1 = rand();
                A = (p - q * (t / Max_iteration) ^ (cosh2 / sinh2)) * r1;

                % Exploration and exploitation phases based on ct and Alpha
                if t <= floor(Max_iteration / ct)  % Exploration phase
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
                else  % Exploitation phase
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
        end

        % Evaluate solutions and update best
        for i = 1:N
            X(i, :) = max(min(X(i, :), adaptive_ub), adaptive_lb);
            Objective_values(1, i) = fobj(X(i, :));
            
            if Objective_values(1, i) < Destination_fitness
                Destination_position = X(i, :);
                Destination_fitness = Objective_values(1, i);
            end
        end

        Convergence_curve(t) = Destination_fitness;
        t = t + 1;
    end
end
