close all; clc; clearvars;
format long

% Constants and pixel conversions
rod_diameter = 1.25e-3; % Meters
rod_pixel_width = 2144;
pixel_size = rod_diameter/rod_pixel_width; % Meters per pixel reference size
pollen_diameter = 52*pixel_size;
V_off = 0.06; % Volts
gamma = 1620; % Hz
V_ac = 2000;
omega = 50; % Hz
z_err = (2/3)*pollen_diameter;

% Reading data from the excel file
V_data = readmatrix('Paul part 2.xlsx');
V_mon = rmmissing(V_data(:,2));
Z1 = rmmissing(V_data(:,3)*pixel_size);
Z2 = rmmissing(V_data(:,4)*pixel_size);

Z_eq = (Z1 + Z2)/2;
Z_max = (Z2 - Z1 - pollen_diameter)/2;

V_DC = 82*(V_mon - V_off); % Volt
E_DC = 4850*(V_mon - V_off); % Volt/Meter

% Least squares fit for Zmax
[chi_squared_eq,expected_eq,residuals_eq,alpha,alpha_err] = LS_fit(E_DC,Z_eq,z_err);
[chi_squared_max,expected_max,residuals_max,beta,beta_err] = LS_fit(E_DC,Z_max,z_err);

charge_mass = (omega^2)*(beta^2)/(2*alpha)
Z_eff1 = ((alpha*charge_mass*V_ac^2)/(2*gamma^2))^(1/4)
Z_eff2 = sqrt((-beta*omega*V_ac)/(2*gamma))
chi_squared_eq
chi_squared_max


figure (1)
hold on
grid on
plot(E_DC,Z_max,'o')
plot(E_DC,expected_max,'-')
xlabel('E_D_C (Volt/Meter)')
ylabel('Z_m_a_x (Meter)')
legend('',['\alpha = ',num2str(alpha),'\pm',num2str(alpha_err)])


figure (2)
hold on
grid on
errorbar(E_DC,residuals_max,z_err,'r')
xlabel('E_D_C (Volt/Meter)')
ylabel('Z_m_a_x (Meter)')

figure (3)
hold on
grid on
plot(E_DC,Z_eq,'o')
plot(E_DC,expected_eq,'-')
xlabel('E_D_C (Volt/Meter)')
ylabel('Z_e_q (Meter)')
legend('',['\beta = ',num2str(beta),'\pm',num2str(beta_err)])

figure(4)
hold on
grid on
errorbar(E_DC,residuals_eq,z_err,'r')
xlabel('E_D_C (Volt/Meter)')
ylabel('Z_e_q (Meter)')


% Fit according to Barlow, R. J. (1989). Statistics, a Guide to the Use
% of Statistical Methods in the Physical Science, page 100
function [chi_squared,expected,residuals,a,a_err] = LS_fit(x_data,y_data,err)
    x_mean = mean(x_data);
    y_mean = mean(y_data);
    N = numel(y_data);
    Vy = var(y_data);
    sigma_y = err;
    rho_ab = -(x_mean)/sqrt(mean(x_data.^2));

    a = (mean(x_data.*y_data)  - x_mean*y_mean)/(mean(x_data.^2) - x_mean^2);
    b = y_mean - a*x_mean;
    expected = a*x_data + b;
    residuals = y_data - expected;
    a_err = (err^2)/(N*(mean(x_data.^2) - x_mean^2));
    chi_squared = ((N*Vy/(sigma_y^2)))*(1 - rho_ab^2);
end