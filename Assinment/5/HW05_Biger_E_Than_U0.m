% Welcome to the Exciting MATLAB World!
% 🌟 Practice Session: HW - 05 Part 1🚀
% 🧠 Your Brain Gym: Tackling MATLAB with Seyed Mohammad Sajadi's Magic
% 🎯 Student Number: 402448040
% Ready, set, code! 💻✨

% Get a number
U0 = input('Enter the U0: ');

% Set initial data
E = linspace(U0, 7, 50);
h = 4.13567e-15;
H = h / (2 * pi);

m1 = input('Enter the m1: ');
m2 = input('Enter the m2: ');

k1 = sqrt(2 * m1 * E / H^2);
k2 = sqrt(2 * m2 * (E - U0) / H^2);

% Formula transmission, reflection coefficient
Transmission = ((4 * k1 .* k2) / (m1 .* m2)) ./ ((k1 ./ m1) + (k2 ./ m2)).^2;
Reflection = 1 - Transmission;

T = 4 ./ (1 + k2 .* m1 ./ (k1 .* m2)).^2;
R = (1 - k2 * m1 ./ k1 * m2).^2 ./ (1 + k2 * m1 ./ k1 * m2).^2;


% Plotting
subplot(2, 2, 1); % Two rows, two columns, first plot
plot(E, Transmission);
grid on;
xlabel('Energy (E)');
ylabel('Transmission Coefficient (T)');
title('Transmission Coefficient vs Energy');

subplot(2, 2, 2); % Two rows, two columns, second plot
plot(E, Reflection);
grid on;
xlabel('Energy (E)');
ylabel('Reflection Coefficient (R)');
title('Reflection Coefficient vs Energy');

subplot(2, 2, 3); % Two rows, two columns, third plot
plot(E, T);
grid on;
xlabel('Energy (E)');
ylabel('transmission probability');
title('transmission probability vs Energy');

subplot(2, 2, 4); % Two rows, two columns, fourth plot
plot(E, R);
grid on;
xlabel('Energy (E)');
ylabel('reflection probability');
title('reflection probability vs Energy');
