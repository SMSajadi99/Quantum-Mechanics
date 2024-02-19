% Welcome to the Exciting MATLAB World!
% 🌟 Practice Session: HW - 08 🚀
% 🧠 Your Brain Gym: Tackling MATLAB with Seyed Mohammad Sajadi's Magic
% 🎯 Student Number: 402448040
% Ready, set, code! 💻✨

% Clear the workspace and command window
close all
clc

% Get user input for the number of iterations
numIterations = input('Enter the number of Levels: ');

% Create a subplot for the Quantum Wave Functions
figure;

% Loop for the specified number of iterations
for i = 0:numIterations
    
    % Define symbolic variable x
    syms x h w

    % Calculate the exponential function
    y = exp(-x^2);
  
    % Calculate the i-th derivative of the exponential function
    derivative = diff(y,i);
    
    % Hermite polynomial expression
    H = ((-1)^i) * (exp(x^2)) * (derivative);

    % Define constants and initial wave function
    A0 = sqrt(1/sqrt(pi));
    wave0 = A0 * exp(-(x^2)/2);

    % Calculate the Quantum Wave Function for the current level
    Wave = (1/(sqrt((2^i) * factorial(i)))) * H * wave0;
    
    % Print the energy expression for each level
    energyExpression = (h * w * i + h * w /2);
    fprintf('Energy Expression for Level %d: E%d = %s\n', i, i, char(energyExpression));

    % Add the Quantum Wave Function to the subplot
    subplot(numIterations+1, 1, i+1);
    fplot(Wave)
    title(['Wave Function - Level ', num2str(i)])
    xlabel('x')
    ylabel('Wave Function Amplitude')
end