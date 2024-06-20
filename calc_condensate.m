function condensation_in_pipe
clc;
%% Define parameters
global cp m k D kf A Ta R Tk M
global xT yT

mh = 17.5;          % Mass flow rate in (kg/h) ~= (liter/h)

d = 0.1;            % Diameter of pipe (m)
D = (log((d+2*M)/d));
A = d^2*pi()/4;     % Cross-sectional area of the pipe (m^2)
kf = 0.1;           % Thermal conductivity of the fluid (W/mK)
m = mh/3600;        % Mass flow rate (kg/s)
cp = 2000;          % Specific heat capacity of the fluid (J/kgK)
R = 287;            % Specific gas constant for air (J/kgK)
Ta = 20+273;
k = 0.05;
M = 0.02;
h_kond = 1830550;

p0 = 7*10^5;    % Inlet pressure (Pa)
T0 = 300+273;   % Inlet temperature (K)

Tk = 226.1+273; % Saturation Temperature in p0 (K)

x_span = [0, 50];
%[xT, yT] = ode23s(@T_ODE, x_span, [T0; 0; 1]);
[xT, yT] = ode23s(@T_ODE, x_span, [T0; 0]);

T = yT(:,1);
T_kond = Tk.*ones(length(xT),1);
%X = [xT,yT(:,3)]

T = T - 273;
T_kond = T_kond - 273;

% plot temperature
    figure;
    plot(xT, T, 'r');
    hold on;
    plot(xT, T_kond, 'b');
    xlabel('Distance along pipe (m)');
    ylabel('Temperature (K)');
    title('Temperature Profile along the Pipe');

% save diagram
matlab2tikz('T(x).tex');


m = (50-3.5) * (2 * pi * k) / D * (Tk-Ta) / h_kond  % kg/s
m_min = m * 60                                      % kg/min
m_h = m_min *60                                     % kg/h

M_1year = m * 60 * 60 * 24 * 365 / 1000;    % t/year
cost = M_1year * 31000                      % Ft/year

end

function dydx = T_ODE(x, y)
    global cp m k D kf A Ta Tk T_val M
    % h=cp * m * y(1);
    % Qvez=kf*A*y(2);

    %h_kond = 1830550;
    Qfal= (2 * pi * k) ./ D * (y(1)-Ta);

    % if y(1,1) > Tk
    %     T_val = Ta;
    % else
    %     T_val = Tk;
    % end
    % 
    % Qkond = (2 * pi * k) / D * (T_val-Ta);

    dydx = [y(2)
            -(cp * m * y(2) + Qfal) / (kf * A)];

    % dydx = [y(2)
    %    -(cp * m * y(2) + Qfal - Qkond) / (kf * A)];
end

