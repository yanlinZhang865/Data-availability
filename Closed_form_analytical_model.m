% Air parameters
u = 1.825e-5; % Dynamic viscosity，unit: Pa∙s
c = 343.6; % Sound speed，unit: m/s
rho = 1.2041; % Air density，unit: kg/m^3

% Impedance tube
D = 0.029; % Impedance tube diameter，unit: m

% Material properties
E = 148500; % Young's modulus，unit: Pa
rho_s = 1070; % Material density，unit: kg/m^3
eta = 0.2; % Loss factor
v = 0.49; % Possion's ratio

% Geometric dimensions
a = 0.004; % Radius of cylindrical shell，unit：m
l = 0.035; % Length of cylindrical shell，unit：m
t = 0.001; % Thickness of cylindrical shell，unit：m
R = 0.0139; % Radius of outer cavity，unit：m
h = 0.0488; % Height of outer cavity，unit：m

% Frequency 
f = 200:2:750; % frequency range, unit: Hz

% Initialize variables
alpha = zeros(size(f));
x_s = zeros(size(f));
y_s = zeros(size(f));

for i = 1:length(f)
    Omega = 2 * pi * f(i); % Angular frequency
    k_v = sqrt(-1i * Omega * rho / u); % Viscous wave number
    S_e = pi * (2 * a + t) * l; % Mid-surface area of the cylindrical shell
    V = pi * R^2 * h - pi * (a + t)^2 * l; % Cavity volume
    S = pi * D^2 / 4; % Cross-sectional area of impedance tube

    % Stiffness enhancement due to Poisson effect
    E_star = E / (1 - v^2); 

    % Acoustic impedance of air column
    Z_a = (1 / (pi * a^2)) * (- u * k_v^2 * l / (1 - 2 * besselj(1, k_v * a) / (k_v * a * besselj(0, k_v * a))) + 2 * sqrt(2 * u * Omega * rho) + 1i * Omega * rho * 0.85 * a *(2 - 1.25 * a / R));

    % Acoustic reactance of cavity
    Z_c = -1i * rho * c^2 / (Omega * V);

    % Equivalent acoustic impedance of shell
    b = a + t / 2; % Mid-surface radius of cylindrical shell，unit: m
    Z_b = E_star * t / (2 * pi * l * Omega * b^3) * (eta + 1i * (rho_s * b^2 * Omega^2 / E_star - 1));
     
    % Total acoustic impedance and sound absorption
    Z_e = tanh(sqrt(Z_a / Z_b)) * sqrt(Z_a * Z_b) + Z_c; 
    alpha(i) = 1 - abs((Z_e * S - rho * c) / (Z_e * S + rho * c))^2;

    % Normalized acoustic impedance
    x_s(i) = real(Z_e * S / (rho * c)); % Normalized acoustic resistance
    y_s(i) = imag(Z_e * S / (rho * c)); % Normalized acoustic reactance
end

% Plot curves
figure;
subplot(3,1,1);
plot(f, alpha);
title('alpha change with f');
xlabel('Frequency (Hz)');
ylabel('Absorption');

subplot(3,1,2);
plot(f, x_s);
title('x_s change with f');
xlabel('Frequency (Hz)');
ylabel('x_s');

subplot(3,1,3);
plot(f, y_s);
title('y_s change with f');
xlabel('Frequency (Hz)');
ylabel('y_s');

% Save data to Excel file
data = table(f', alpha', x_s', y_s', 'VariableNames', {'Frequency', 'Alpha', 'X_s', 'Y_s'});
writetable(data, 'acoustic_impedance_data.xlsx');


