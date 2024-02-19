% Welcome to the Exciting MATLAB World!
% 🌟 Practice Session: Project - 01 🚀
% 🧠 Your Brain Gym: Tackling MATLAB with Seyed Mohammad Sajadi's Magic
% 🎯 Student Number: 402448040
% Ready, set, code! 💻✨

% Clear the workspace and command window
clear all
close
clc

% Constants
h = 4.13567e-15 * 1.602e-19; % Planck's constant in J
H = h / (2 * pi); % Reduced Planck's constant (ħ)
U0 = 10 * 1.602e-19;
m0 = 9.1e-31;
L = 0.2e-9;
W = 0.1e-9;

% Ask the user for the number of steps
% numSteps = input('Enter the number of steps: ');
numSteps = 4;

% Energy range
E = linspace(0.002 * 1.602e-19, 30 * 1.602e-19, 15000);
Size_E = size(E);

% Initialize array to store results
P11 = zeros(1, Size_E(2));

counter = 1;

% Loop over energy values
for e = 0.002 * 1.602e-19 : 0.002 * 1.602e-19 : 30 * 1.602e-19
    % Initialize the result matrix
    resultMatrix = eye(2);

    % Loop over steps
    for step = 1 : numSteps
        if mod(step, 2) == 1
            % Odd step
            k_odd_FreeSpace = sqrt(2 * m0 * e / (H ^ 2));
            k_odd_PotentialBarrier = sqrt(2 * m0 * (e - U0) / (H ^ 2));

            % Construct matrices for odd steps
            matrix_odd_FreeSpace = [exp(-1i * k_odd_FreeSpace * L), 0;
                                    0, exp(1i * k_odd_FreeSpace * L)];
            matrix_odd_PotentialBarrier = 0.5 * [1 + k_odd_PotentialBarrier / k_odd_FreeSpace, 1 - k_odd_PotentialBarrier / k_odd_FreeSpace;
                                                  1 - k_odd_PotentialBarrier / k_odd_FreeSpace, 1 + k_odd_PotentialBarrier / k_odd_FreeSpace];

            % Update the result matrix
            resultMatrix = resultMatrix * (matrix_odd_FreeSpace * matrix_odd_PotentialBarrier);

        else
            % Even step
            k_Even_FreeSpace = sqrt(2 * m0 * (e - U0) / (H ^ 2));
            k_Even_PotentialBarrier = sqrt(2 * m0 * e / (H ^ 2));

            % Construct matrices for even steps
            matrix_Even_FreeSpace = [exp(-1i * k_Even_FreeSpace * W), 0;
                                     0, exp(1i * k_Even_FreeSpace * W)];
            matrix_Even_PotentialBarrier = 0.5 * [1 + k_Even_PotentialBarrier / k_Even_FreeSpace, 1 - k_Even_PotentialBarrier / k_Even_FreeSpace;
                                                   1 - k_Even_PotentialBarrier / k_Even_FreeSpace, 1 + k_Even_PotentialBarrier / k_Even_FreeSpace];

            % Update the result matrix
            resultMatrix = resultMatrix * (matrix_Even_FreeSpace * matrix_Even_PotentialBarrier);
        end
    end

    % Calculate and store probability
    P11(counter) = (abs(1 / resultMatrix(1)))^2;
    counter = counter + 1;
end

% Plot the results
figure;
plot(E / (1.602e-19), P11, 'LineWidth', 2);
xlabel('Energy (eV)');
ylabel('Transmission Probability');
title('Transmission Probability vs Energy');
grid on;
