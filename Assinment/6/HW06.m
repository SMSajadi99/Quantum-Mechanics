% Welcome to the Exciting MATLAB World!
% 🌟 Practice Session: HW - 06 🚀
% 🧠 Your Brain Gym: Tackling MATLAB with Seyed Mohammad Sajadi's Magic
% 🎯 Student Number: 402448040
% Ready, set, code! 💻✨

% Clear the workspace and command window
close all
clc

% Initialize input values (commented out for hardcoded values)
% U0 = input("Enter the U0: ");
% m0 = input("Enter the m0: ");
% L = input("Enter the L: ");
% alpha = input("Enter the alpha: ");
U0 = 0.5 * 1.602 * 10^-19;
m0 = 9.11 * 10^-31;
L = 5 * 10^-9;
alpha = 0.07;

% Define constants
h = 4.13567 * 10^-15 * 1.602 * 10^-19; % Planck's constant in J
H = h / (2 * pi); % Reduced Planck's constant (ħ)
M = alpha * m0; % Scaled mass

% Calculate the radius of the circle
radius = sqrt(2 * M * U0 * L^2 / H^2);

% Display the calculated radius
disp(['The calculated radius of the circle is: ', num2str(radius)]);

% Generate x values and compute y values for tan(x) and cot(x)
x = linspace(0, L*pi / 10^-9, 1000000);
y_tan = x .* tan(x);
y_cot = -x .* cot(x);

% Plot the functions
figure;
plot(x, y_tan);
hold on; % Maintain the current plot for additional plots
plot(x, y_cot);
xlim([0 L*pi / (2*10^-9)]);
ylim([0 10]); % Adjust y-axis limits for better visualization

% Plot the circle using parametric equations
theta = linspace(0, 2*pi, 1000000);
x_circle = radius * cos(theta);
y_circle = radius * sin(theta);
plot(x_circle, y_circle);

% Iterate through x values and find intersection points
intersection_points_tan = [];
intersection_points_cot = [];

tolerance = 0.01; % Tolerance for considering points close

for i = 1:length(x)
    % Check if the point is on the circle
    if abs(x(i)^2 + y_tan(i)^2 - radius^2) < tolerance
        % Check if the point is close to any existing point in the array
        is_close = false;
        for j = 1:size(intersection_points_tan, 1)
            if norm([x(i), y_tan(i)] - intersection_points_tan(j, :)) < tolerance
                is_close = true;
                break;
            end
        end
        
        % If the point is not close to any existing point, add it to the array
        if ~is_close
            intersection_points_tan = [intersection_points_tan; x(i), y_tan(i)];
        end
    end

 
    % Check if the point is on the circle
    if abs(x(i)^2 + y_cot(i)^2 - radius^2) < tolerance
        % Check if the point is close to any existing point in the array
        is_close = false;
        for j = 1:size(intersection_points_cot, 1)
            if norm([x(i), y_cot(i)] - intersection_points_cot(j, :)) < tolerance
                is_close = true;
                break;
            end
        end
        
        % If the point is not close to any existing point, add it to the array
        if ~is_close
            intersection_points_cot = [intersection_points_cot; x(i), y_cot(i)];
        end
    end
end


% Display only positive intersection points with y = x*tan(x)
disp('Intersection points with y = x*tan(x):');
positive_intersection_tan = intersection_points_tan(intersection_points_tan(:, 2) > 0, :);
disp(positive_intersection_tan);

% Display only positive intersection points with y = -x*cot(x)
disp('Intersection points with y = -x*cot(x):');
positive_intersection_cot = intersection_points_cot(intersection_points_cot(:, 2) > 0, :);
disp(positive_intersection_cot);


% Perform operations on positive_intersection_tan
tan_col1 = ((positive_intersection_tan(:, 1)).^ 2 * H ^ 2) / (2 * M * (1.602 * 10 ^ -19));
tan_col2 = (((positive_intersection_tan(:, 2)).^ 2 * H ^ 2) / (2 * M * (1.602 * 10 ^ -19)) - U0) * -1;
Energy_positive_intersection_tan = [tan_col1, tan_col2];

% Perform operations on positive_intersection_cot
cot_col1 = ((positive_intersection_cot(:, 1)).^2 * H^2) / (2 * M * (1.602 * 10 ^ -19));
cot_col2 = (((positive_intersection_cot(:, 2)).^2 * H^2) / (2 * M * (1.602 * 10 ^ -19)) - U0) * -1;
Energy_positive_intersection_cot = [cot_col1, cot_col2];

% Display energy 
disp('Energy Intersection points with y = x*tan(x) by e.v:');
disp(Energy_positive_intersection_tan / (1.602*10^-19));

disp('Energy Intersection points with y = -x*cot(x) by e.v:');
disp(Energy_positive_intersection_cot / (1.602*10^-19));

% Count the number of points for tangent and cotangent
num_points_tan = size(positive_intersection_tan, 1);
num_points_cot = size(positive_intersection_cot, 1);

% Display the counts
disp(['Number of points for tangent: ', num2str(num_points_tan)]);
disp(['Number of points for cotangent: ', num2str(num_points_cot)]);
disp(['Total number of points: ', num2str(num_points_tan + num_points_cot)]);


% Plot intersection points
if ~isempty(intersection_points_tan)
    scatter(intersection_points_tan(:, 1), intersection_points_tan(:, 2), 'r', 'filled');
end

if ~isempty(positive_intersection_cot)
    scatter(positive_intersection_cot(:, 1), positive_intersection_cot(:, 2), 'b', 'filled');
end

% Add labels, title, and legend for clarity
xlabel('x');
ylabel('y');
title('Plot of y = x*tan(x), y = -x*cot(x), and a circle with radius r');
legend('y = x*tan(x)', 'y = -x*cot(x)', 'Circle', 'Intersection Points with tan(x)', 'Positive Intersection Points with cot(x)');
grid on;
hold off;


new_positive_intersection_tan = positive_intersection_tan / L;
new_positive_intersection_cot = positive_intersection_cot / L;

% Define the x values for the entire range
x = linspace(-2 * L, 2 * L, 100000); % Adjust the number of points as needed

% Define the corresponding y values for each function
y = zeros(size(x));

row_positive_intersection_cot = size(positive_intersection_cot, 1);
row_positive_intersection_tan = size(positive_intersection_tan, 1);

num_rows = ceil((row_positive_intersection_tan + row_positive_intersection_cot) / 2);
figure;


for i = 1 : (row_positive_intersection_tan + row_positive_intersection_cot)
    subplot(num_rows, 2, i); % Adjust the layout dynamically
    
    if mod(i, 2) == 0
        % Even index
        K = new_positive_intersection_cot(i / 2, 1);
        k = new_positive_intersection_cot(i / 2, 2);

        for j = 1:length(x)
            % First interval: 0 to L
            if x(j) >= 0 && x(j) <= L
                y(j) = sin(K * x(j));
            end
            
            % Second interval: L to 2*L
            if x(j) > L && x(j) <= 2*L
                if sin(K * L) > 0
                    y(j) = exp(-k * (x(j) - 3 * L / 4));
                else
                    y(j) = -exp(-k * (x(j) - 3 * L / 4));
                end
            end
            
            % Third interval: 0 to -L
            if x(j) >= -L && x(j) < 0
                y(j) = sin(K * x(j));
            end
            
            % Fourth interval: -L to -2*L
            if x(j) >= -2*L && x(j) < -L
                if sin(K * (-L)) > 0
                    y(j) = exp(k * (x(j) + 3 * L / 4));
                else
                    y(j) = -exp(k * (x(j) + 3 * L / 4));
                end
            end
        end
        

    else
        % Odd index
        K = new_positive_intersection_tan((i+1) / 2, 1);
        k = new_positive_intersection_tan((i+1) / 2, 2);

        for j = 1:length(x)
            % First interval: 0 to L
            if x(j) >= 0 && x(j) <= L
                y(j) = cos(K * x(j));
            end
            
            % Second interval: L to 2*L
            if x(j) > L && x(j) <= 2*L
                if cos(K * L) > 0
                    y(j) = exp(-k * (x(j) - 3 * L / 4));
                else
                    y(j) = -exp(-k * (x(j) - 3 * L / 4));
                end
            end
            
            % Third interval: 0 to -L
            if x(j) >= -L && x(j) < 0
                y(j) = cos(K * x(j));
            end
            
            % Fourth interval: -L to -2*L
            if x(j) >= -2*L && x(j) < -L
                if cos(K * (-L)) > 0
                    y(j) = exp(k * (x(j) + 3 * L / 4));
                else
                    y(j) = -exp(k * (x(j) + 3 * L / 4));
                end
            end
        end
              
    end
    

    plot(x, y, 'b-', 'LineWidth', 2);
    
    if mod(i, 2) == 0
        title(['Even Subplot ' num2str(i/2)]);
    else
        title(['Odd Subplot ' num2str((i+1)/2)]);
    end
    
    xlabel('x');
    ylabel('ψ(x)');
end
