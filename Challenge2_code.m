close all
clear 
clc
%fix isotherm
T = 0.85;
%fix tolerance and maximum number of iteration
tolerance = 1e-2;
maxiter = 900;
%Thermal equation of state in the form of Van der Waals, using
%dimensionless variables
P = @(v) 8*T./(3*v-1)-3./v.^2;
%area between pressure P and p(v)
Area = @(va, vb, P, T) 8/3*T*log(((3*vb-1)./(3*va-1)))+3*(1./vb-1./va);
%Finding minimum and maximum of P, having fixed T
%The derivative is:
dP = @(v) -4*T*v.^3+9*v.^2-6*v+1;
coeff_dP = [-4*T 9 -6 1];
stationary_points = sort(roots(coeff_dP));
vM = stationary_points(end);
vm = stationary_points(end-1);
pM = P(vM);
pm = P(vm);
if pm < 0
    pm = 0;
end
P0 = 1/2*(pM + pm);
Pj = P0;
for jj = 1:maxiter
    coeff_P = [3*Pj -(Pj+8*T) 9 -3];
    zeros_P = sort(roots(coeff_P));
    v1 = zeros_P(1);
    v2 = zeros_P(2);
    v3 = zeros_P(3);
    Area1 = Pj*(v2 - v1) - Area(v1, v2, Pj, T);
    Area2 = Area(v2, v3, Pj, T) - Pj*(v3-v2);
    if abs(Area1-Area2) < tolerance
        break;
    elseif Area1 < Area2
        P1 = 1/2*(pM + Pj);
        flag = 1;
    else 
        P1 = 1/2*(Pj + pm);
        flag = 0;
    end
    coeff_P = [3*P1 -(P1+8*T) 9 -3];
    zeros_P = sort(roots(coeff_P));
    v1 = zeros_P(1);
    v2 = zeros_P(2);
    v3 = zeros_P(3);
    Area1 = P1*(v2 - v1) - Area(v1, v2, P1, T);
    Area2 = Area(v2, v3, P1, T) - P1*(v3-v2);
    if abs(Area1-Area2) < tolerance
        Pj = P1;
        break;
    elseif (flag)&&(Area1 < Area2)
        Pj = 1/2*(pM + P1);
    elseif ((flag)&&(Area1 > Area2))||(not(flag)&&(Area1 < Area2))
        Pj = 1/2*(Pj+P1);
    else
        Pj = 1/2*(P1+pm);
    end
end
output = [jj, Area1, Area2, v1, v2, v3, Pj];
fprintf('The calculated reduced critical pressure is: %f,', Pj)
fprintf('at the reduced molar volumes: %f, %f, of the coexisting liquid and vapour phases ', v1, v3)
fprintf('we have performed: %f iterations and the absolute error is: %f', jj, abs(Area1-Area2))
% Plot
v_int = 0.5:0.001:4;
iso_T = P(v_int);
plot(v_int, iso_T,'b')
hold on
plot(v_int, Pj*ones(size(v_int)),'--r')
plot(zeros_P, Pj*ones(size(zeros_P)),'*g')
plot(vm, pm, 'ok')
plot(vM, pM, 'om')
title('Pressure at liquid-vapour phase transition using Wan der Waals equation')
xlabel('volume [-]')
ylabel('pressure [-]')
legend('Wan der Waals isothermal', 'Pressure of the critical point', 'v1, v2, v3','minimum pressure','maximum pressure')
grid on
hold off