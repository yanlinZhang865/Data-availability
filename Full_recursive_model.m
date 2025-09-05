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
a = 0.0045; % Radius of cylindrical shell，unit：m
l = 0.035; % Length of cylindrical shell，unit：m
t = 0.001; % Thickness of cylindrical shell，unit：m
R = 0.0139; % Radius of outer cavity，unit：m
h = 0.0488; % Height of outer cavity，unit：m

% Frequency 
f = 10:1:1000; % frequency range, unit: Hz

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

    % Initialize variables  
    x_j = zeros(1, 1000);
    q_x_j = zeros(1, 1000);
    k_x_j = zeros(1, 1000);
    R_x_j = zeros(1, 1000);
    C_x_j = zeros(1, 1000);
    M_x_j = zeros(1, 1000);
    Zm_x_j = zeros(1, 1000);
    S_x_j = zeros(1, 1000);
    Zb_x_j = zeros(1, 1000);
    Za_x_j = zeros(1, 1000);

    b = a + t / 2; % Mid-surface radius of cylindrical shell，unit: m

    beta = nthroot((3 * (1 - v^2)) / (b^2 * t^2), 4); % Bending stiffness coefficient 

    % Define sinh(beta * l) and sin(beta * l)
    sinh_bl = sinh(beta * l);
    sin_bl = sin(beta * l);
    cosh_bl = cosh(beta * l);
    cos_bl = cos(beta * l);

    % Calculation of Za_x_j and Zb_x_j
    for j = 1:1000
        x_j(j) = j / 1000 * l;

        % Define
        sinh_bx = sinh(beta * x_j(j));
        sin_bx = sin(beta * x_j(j));
        cosh_bx = cosh(beta * x_j(j));
        cos_bx = cos(beta * x_j(j));
    
        % The spatial modulation function
        q_x_j(j) = 1 + ((sinh_bl^2 + sin_bl^2) * sin_bx * sinh_bx) / (cosh_bl^2 + cos_bl^2) - ...
            (cosh_bl * sinh_bl + cos_bl * sin_bl) * (sin_bx * cosh_bx - cos_bx * sinh_bx) / (cosh_bl^2 + cos_bl^2) - ...
            cos_bx * cosh_bx;

        % Mechanical impedance
        k_x_j(j) = 2 * pi * l * E_star * t / (1000 * b * q_x_j(j));
        R_x_j(j) = (eta * k_x_j(j)) / Omega;
        M_x_j(j) = (rho_s * pi * t * (t + 2 * a) * l) / 1000;
        Zm_x_j(j) = R_x_j(j) + 1j * (Omega * M_x_j(j) - k_x_j(j) / Omega);

        % Conversion to equivalent acoustic impedance
        S_x_j(j) = S_e / 1000;
        Zb_x_j(j) = Zm_x_j(j) / (S_x_j(j))^2;

        % Acoustic impedance
        Za_x_j(j) = Z_a / 1000;
    end

     % Initialize T
    T = zeros(1, 1000); % Create an array of 1000 elements to store T_1 to T_1000

    T(1000) = Zb_x_j(1000) * Za_x_j(1000) / (Zb_x_j(1000) + Za_x_j(1000)) + Za_x_j(999);

    % Recursive calculation from T_999 to T_2
    for n = 999:-1:2
        T(n) = T(n+1) * Zb_x_j(n) / (T(n+1) + Zb_x_j(n)) + Za_x_j(n-1);
    end

    % Calculation of T(1)
    T(1) = T(2) * Zb_x_j(1) / (T(2) + Zb_x_j(1));

    % Total acoustic impedance and sound absorption
    Z_e = T(1) + Z_c;
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


