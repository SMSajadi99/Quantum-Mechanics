% Welcome to the Exciting MATLAB World!
% 🌟 Practice Session: HW - 02 🚀
% 🧠 Your Brain Gym: Tackling MATLAB with Seyed Mohammad Sajadi's Magic
% 🎯 Student Number: 402448040
% Ready, set, code! 💻✨

% Clear the workspace and command window
close all
clc

% Get a number
userNumber = input('Enter a number of energy level: ');

% Set initial energy E0
syms x E0;

% Create a figure for subplots
figure;

% Loop through userNumber
for n = 1:userNumber
    % Choose the interval [0, L]
    L = 2*pi;

    % Sine function with changing formula based on userNumber
    y = sqrt(2/L)*sin(n*pi*x/L);
    
    % Calculate and print the expression E_n = E0 * userNumber^2
    En = E0 * n^2;
    disp(['E_' num2str(n) ' = ' char(En)]);

    % Create subplots
    subplot(userNumber, 2, (n-1)*2 + 1);
    fplot(y, [0, L]);
    title(['Wave Function, n = ' num2str(n)]);
    xlabel('x');
    ylabel(['sin(' num2str(n) '\pix/' num2str(L) ')']);
    grid on;
    
    subplot(userNumber, 2, (n-1)*2 + 2);
    fplot(y^2, [0, L]); % Plot y^2 for the product with the conjugate
    title(['Probability of Attendance, n = ' num2str(n)]);
    xlabel('x');
    ylabel(['sin(' num2str(n) '\pix/' num2str(L) ')^2']);
    grid on;
end
