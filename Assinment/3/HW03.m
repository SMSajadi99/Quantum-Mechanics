% Welcome to the Exciting MATLAB World!
% 🌟 Practice Session: HW - 03 🚀
% 🧠 Your Brain Gym: Tackling MATLAB with Seyed Mohammad Sajadi's Magic
% 🎯 Student Number: 402448040
% Ready, set, code! 💻✨

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
    y_odd  = sqrt(2/L)*cos(n*pi*x/L);
    y_even = sqrt(2/L)*sin(n*pi*x/L);
    
    % Calculate and print the expression E_n = E0 * userNumber^2
    En_odd  = E0 * n^2;
    En_even = E0 * n^2;

    if mod(n, 2) == 0
        disp(['En_even_' num2str(n) ' = ' char(En_even)]);

        % Create subplots
        subplot(userNumber, 2, (n-1)*2 + 1);
        fplot(y_even, [-L/2, L/2]);
        title(['Wave Function, n = ' num2str(n)]);
        xlabel('x');
        ylabel(['sin(' num2str(n) '\pix/' num2str(L) ')']);
        grid on;
    
        subplot(userNumber, 2, (n-1)*2 + 2);
        fplot(y_even^2, [-L/2, L/2]); % Plot y^2 for the product with the conjugate
        title(['Probability of Attendance, n = ' num2str(n)]);
        xlabel('x');
        ylabel(['sin(' num2str(n) '\pix/' num2str(L) ')^2']);
        grid on;

    else
        disp(['En_odd_' num2str(n) ' = ' char(En_odd)]);

        % Create subplots
        subplot(userNumber, 2, (n-1)*2 + 1);
        fplot(y_odd, [-L/2, L/2]);
        title(['Wave Function, n = ' num2str(n)]);
        xlabel('x');
        ylabel(['cos(' num2str(n) '\pix/' num2str(L) ')']);
        grid on;
    
        subplot(userNumber, 2, (n-1)*2 + 2);
        fplot(y_odd^2, [-L/2, L/2]); % Plot y^2 for the product with the conjugate
        title(['Probability of Attendance, n = ' num2str(n)]);
        xlabel('x');
        ylabel(['cos(' num2str(n) '\pix/' num2str(L) ')^2']);
        grid on;

    end

end
