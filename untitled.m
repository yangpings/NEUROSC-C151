%% Master File Take #3
tspan = [0:0.0001:50];  % Time span 1 to 50 sec, 0.1ms stamp.

% Initial Conditions / Variables 
V_0 = -71.847;
m_0 = 0.015;
h_0 = 0.981;
m_NaP_0 = 0.002;
h_NaP_0 = 0.797;
n_0 = 0.158;
m_T_0 = 0.001;
h_T_0 = 0.562;
m_P_0 = 0;
m_N_0 = 0.001;
h_N_0 = 0.649;
z_SK_0 = 0;
m_A_0 = 0.057;
h_A_0 = 0.287;
m_H_0 = 0.182;
Ca_i_0 = 0.0604;

y = zeros(16, 1);

y(1) = V_0;
y(2) = m_0;
y(3) = h_0;
y(4) = m_NaP_0;
y(5) = h_NaP_0;
y(6) = n_0;
y(7) = m_T_0;
y(8) = h_T_0;
y(9) = m_P_0;
y(10) = m_N_0;
y(11) = h_N_0;
y(12) = z_SK_0;
y(13) = m_A_0;
y(14) = h_A_0;
y(15) = m_H_0;
y(16) = Ca_i_0;

disp(y)

% Call the ODE solver
[t, y] = ode45(@ode_function, tspan, y);

figure;
plot(t, y(:, 1));
xlabel('Time');
ylabel('Voltage');
title('Voltage vs. Time');
% Plot the variables with respect to time
figure;

subplot(4, 4, 1);
plot(t, y(:, 1));
xlabel('Time');
ylabel('Voltage');
title('Voltage vs. Time');

subplot(4, 4, 2);
plot(t, y(:, 2));
xlabel('Time');
ylabel('m');
title('m vs. Time');

subplot(4, 4, 3);
plot(t, y(:, 3));
xlabel('Time');
ylabel('h');
title('h vs. Time');

subplot(4, 4, 4);
plot(t, y(:, 4));
xlabel('Time');
ylabel('m_{NaP}');
title('m_{NaP} vs. Time');

subplot(4, 4, 5);
plot(t, y(:, 5));
xlabel('Time');
ylabel('h_{NaP}');
title('h_{NaP} vs. Time');

subplot(4, 4, 6);
plot(t, y(:, 6));
xlabel('Time');
ylabel('n');
title('n vs. Time');

subplot(4, 4, 7);
plot(t, y(:, 7));
xlabel('Time');
ylabel('m_{T}');
title('m_{T} vs. Time');

subplot(4, 4, 8);
plot(t, y(:, 8));
xlabel('Time');
ylabel('h_{T}');
title('h_{T} vs. Time');

subplot(4, 4, 9);
plot(t, y(:, 9));
xlabel('Time');
ylabel('m_{P}');
title('m_{P} vs. Time');

subplot(4, 4, 10);
plot(t, y(:, 10));
xlabel('Time');
ylabel('m_{N}');
title('m_{N} vs. Time');

subplot(4, 4, 11);
plot(t, y(:, 11));
xlabel('Time');
ylabel('h_{N}');
title('h_{N} vs. Time');

subplot(4, 4, 12);
plot(t, y(:, 12));
xlabel('Time');
ylabel('z_{SK}');
title('z_{SK} vs. Time');

subplot(4, 4, 13);
plot(t, y(:, 13));
xlabel('Time');
ylabel('m_{A}');
title('m_{A} vs. Time');

subplot(4, 4, 14);
plot(t, y(:, 14));
xlabel('Time');
ylabel('h_{A}');
title('h_{A} vs. Time');

subplot(4, 4, 15);
plot(t, y(:, 15));
xlabel('Time');
ylabel('m_{H}');
title('m_{H} vs. Time');

subplot(4, 4, 16);
plot(t, y(:, 16));
xlabel('Time');
ylabel('Ca_{i}');
title('Ca_{i} vs. Time');