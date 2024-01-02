% Rayhan Mohammed MATLAB Project 2
% Part 1(a): Symbolic solution
% Symbolic part
syms t y(t);
ode = diff(y, t, 2) + 3*diff(y, t) + 2*y == 4*sin(t);
cond1 = y(0) == 0;
cond2 = subs(diff(y), 0) == 1;
sol = dsolve(ode, [cond1, cond2]);
y_exact = @(t) subs(sol, t);

% Part 1(b): Rewrite the second-order differential equation as a system
syms x1(t) x2(t);
eq1 = diff(x1, t) == x2;
eq2 = diff(x2, t) == 4*sin(t) - 3*x2 - 2*x1;

% Part 2: Numerical solution using Euler's method
h = pi/50; % time step
tfinal = 20*pi; % final time for computation
X0 = [0; 1]; % initial conditions for x1(t) and x2(t)

n = round(tfinal/h);
t = 0:h:tfinal;

% Euler's method for the system of first-order ODEs
X_num_euler = zeros(2, n+1);
X_num_euler(:, 1) = X0;

for i = 1:n
    F = [X_num_euler(2, i); 4*sin(t(i)) - 3*X_num_euler(2, i) - 2*X_num_euler(1, i)];
    X_num_euler(:, i+1) = X_num_euler(:, i) + h * F;
end

% Part 2(a): Plot numerical and exact solutions for h = pi/50
figure;
plot(t, X_num_euler(1, :), '-b', t, X_num_euler(2, :), '-g', t, y_exact(t), '--r');
xlabel('t');
ylabel('Value');
title('Numerical and Exact Solutions (Euler''s Method)');
legend('Numerical x1', 'Numerical x2', 'Exact y');

% Part 2(b): Compute and display error values for different step sizes
h_values = [pi/50, pi/100, pi/200];
error_values_euler = zeros(1, length(h_values));

disp('Error Table (Euler''s Method):');
disp('   h            e1 (Euler)');
disp('--------    ----------------------');

for j = 1:length(h_values)
    h = h_values(j);
    
    n = round(tfinal / h);
    t = 0:h:tfinal;
    
    X_num_euler = zeros(2, n+1);
    X_num_euler(:, 1) = X0;

    for i = 1:n
        F = [X_num_euler(2, i); 4*sin(t(i)) - 3*X_num_euler(2, i) - 2*X_num_euler(1, i)];
        X_num_euler(:, i+1) = X_num_euler(:, i) + h * F;
    end

    % Calculate and store the error
    error_values_euler(j) = max(abs(y_exact(t) - X_num_euler(1, :)));
    
    % Display the error values
    fprintf('%8.4f    %17.10f\n', h, error_values_euler(j));
end

% Part 2(c): Plot the trajectory in the (x1, x2)-plane for h = pi/200
h = pi/200; % reset time step
n = round(tfinal / h);
t = 0:h:tfinal;

X_num_euler_h200 = zeros(2, n+1);
X_num_euler_h200(:, 1) = X0;

for i = 1:n
    F = [X_num_euler_h200(2, i); 4*sin(t(i)) - 3*X_num_euler_h200(2, i) - 2*X_num_euler_h200(1, i)];
    X_num_euler_h200(:, i+1) = X_num_euler_h200(:, i) + h * F;
end

% Plot the trajectory in the (x1, x2)-plane for h = pi/200
figure;
plot(X_num_euler_h200(1, :), X_num_euler_h200(2, :), '-b');
xlabel('x1');
ylabel('x2');
title('Trajectory in (x1, x2)-plane (Euler''s Method)');

% Part 3: Numerical solution using 2nd-order Runge-Kutta method
h_values = [pi/50, pi/100, pi/200];
error_values_rk2 = zeros(1, length(h_values));

disp('Error Table (2nd-order Runge-Kutta):');
disp('   h            e1 (Runge-Kutta)');
disp('--------    ----------------------');

for j = 1:length(h_values)
    h = h_values(j);
    
    n = round(tfinal / h);
    t = 0:h:tfinal;
    
    X_num_rk2 = zeros(2, n+1);
    X_num_rk2(:, 1) = X0;

    % 2nd-order Runge-Kutta method for the system of first-order ODEs
    for i = 1:n
        k1 = h * [X_num_rk2(2, i); 4*sin(t(i)) - 3*X_num_rk2(2, i) - 2*X_num_rk2(1, i)];
        k2 = h * [X_num_rk2(2, i) + k1(2); 4*sin(t(i) + h) - 3*(X_num_rk2(2, i) + k1(2)) - 2*(X_num_rk2(1, i) + k1(1))];
        X_num_rk2(:, i+1) = X_num_rk2(:, i) + 0.5 * (k1 + k2);
    end

    % Part 3(a): Plot numerical and exact solutions for h = pi/50
    if h == pi/50
        figure;
        plot(t, X_num_rk2(1, :), '-b', t, X_num_rk2(2, :), '-g', t, y_exact(t), '--r');
        xlabel('t');
        ylabel('Value');
        title('Numerical and Exact Solutions (2nd-order Runge-Kutta)');
        legend('Numerical x1', 'Numerical x2', 'Exact y');
    end

% Plot the trajectory in the (x1, x2)-plane for h = pi/200 (Runge-Kutta)
figure;
plot(X_num_rk2(1, :), X_num_rk2(2, :), '-b');
xlabel('x1');
ylabel('x2');
title('Trajectory in (x1, x2)-plane (2nd-order Runge-Kutta)');

    % Part 3(b): Calculate and store the error
    error_values_rk2(j) = max(abs(y_exact(t) - X_num_rk2(1, :)));
    
    % Display the error values
    fprintf('%8.4f    %17.10f\n', h, error_values_rk2(j));
end
