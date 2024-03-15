% Testing Euler's method
clear
x0 = 2;
t0 = 0;
dt = 0.001;
tmax = 20;
tvec = t0:dt:tmax;
xvec = zeros(size(tvec));
xvec(1) = x0;

% Iterative computing of xvec
for i = 2:length(tvec) % there should be no +1
    tval = tvec(i-1);   % should be i not j        
    dxdt = 2*sin(0.1*tval*tval);
    % Euler's method line
    xvec(i) = xvec(i-1) + dxdt*dt;
end

%plot xvec values as time series
plot(tvec,xvec);
