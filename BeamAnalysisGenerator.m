clc; clear;
function reaction_magnitude = reactionsCalc(lp, rp, loads, udls)
    syms r1 r2
    ml = 0;
    for i = 1:size(loads,1)
        ml = ml + loads(i,1)*(loads(i,2) - lp);
    end
    for i = 1:size(udls,1)
        L_udl = udls(i,3) - udls(i,2);
        ml = ml + udls(i,1)*L_udl*((udls(i,2)+udls(i,3))/2 - lp);
    end
    MomEq = r1*lp + r2*rp == ml;
    LoadEq = r1 + r2 == sum(loads(:,1)) + sum(udls(:,1).* (udls(:,3)-udls(:,2)));
    reaction = solve([MomEq, LoadEq], [r1, r2]);
    reaction_magnitude = zeros(2,2);
    reaction_magnitude(1,1) = double(reaction.r1);
    reaction_magnitude(2,1) = double(reaction.r2);
    reaction_magnitude(1,2) = lp;
    reaction_magnitude(2,2) = rp;
end
%% INPUTS
L = input('Enter span length of beam (m): ');
% Point Loads
Total_pt = input('Number of point loads? ');
loads = zeros(Total_pt,2);
for i = 1:Total_pt
    loads(i,1) = input('Load magnitude (N): ');
    loads(i,2) = input('Distance from left support (m): ');
end
% UDLs
Total_udl = input('Number of UDLs? ');
udls = zeros(Total_udl,3);
for i = 1:Total_udl
    udls(i,1) = input('UDL intensity (N/m): ');
    udls(i,2) = input('Start position (m): ');
    udls(i,3) = input('End position (m): ');
end
%% CALCULATE REACTIONS
lp = 0;   % Left support
rp = L;   % Right support
reaction_magnitude = reactionsCalc(lp, rp, loads, udls);
N1 = reaction_magnitude(1,1);
N2 = reaction_magnitude(2,1);
%% DISCRETIZE BEAM
dx = 0.01;
x = 0:dx:L;
V = zeros(size(x));
M = zeros(size(x));
%% CALCULATE SFD & BMD
for i = 1:length(x)
    V(i) = N1;
    M(i) = N1*x(i);
    for j = 1:Total_pt
        if x(i) >= loads(j,2)
            V(i) = V(i) - loads(j,1);
            M(i) = M(i) - loads(j,1)*(x(i)-loads(j,2));
        end
    end
    for j = 1:Total_udl
        if x(i) >= udls(j,2)
            L_udl = min(x(i), udls(j,3)) - udls(j,2);
            if L_udl > 0
                V(i) = V(i) - udls(j,1)*L_udl;
                M(i) = M(i) - udls(j,1)*L_udl*(x(i) - (udls(j,2)+L_udl/2));
            end
        end
    end
end
%% KEY POINTS
[M_max, idx_maxM] = max(M);
x_maxM = x(idx_maxM);
idx_zeroV = find(V(1:end-1).*V(2:end) < 0);
x_zeroV = x(idx_zeroV);
%% MATERIAL PROPERTIES FOR DEFLECTION
E = input('Enter Young''s Modulus E (Pa): ');
I = input('Enter Moment of Inertia I (m^4): ');
%% CALCULATE DEFLECTION
[x_deflection, deflection] = calculateDeflection(L, loads, udls, N1, N2, E, I);
%% PLOTTING
figure;
subplot(3,1,1)
plot(x,V,'b','LineWidth',1.5); hold on;
xlabel('x (m)'); ylabel('Shear Force V (N)');
title('Shear Force Diagram'); grid on;
plot(x_zeroV,zeros(size(x_zeroV)),'ro','MarkerFaceColor','r');
for i=1:Total_pt
    plot(loads(i,2), N1 - sum(loads(1:i,1)),'ks','MarkerFaceColor','k');
end
for i=1:Total_udl
    plot([udls(i,2) udls(i,3)],[0 0],'g^','MarkerFaceColor','g');
end
subplot(3,1,2)
plot(x,M,'r','LineWidth',1.5); hold on;
xlabel('x (m)'); ylabel('Bending Moment M (Nm)');
title('Bending Moment Diagram'); grid on;
plot(x_maxM,M_max,'mo','MarkerFaceColor','m');
for i=1:Total_pt
    idx = find(x>=loads(i,2),1);
    plot(loads(i,2),M(idx),'ks','MarkerFaceColor','k');
end
for i=1:Total_udl
    idx_s = find(x>=udls(i,2),1);
    idx_e = find(x>=udls(i,3),1);
    plot([udls(i,2) udls(i,3)],[M(idx_s) M(idx_e)],'g^','MarkerFaceColor','g');
end
legend('V or M','Shear=0','Point Load','UDL Start/End','Max Moment');
subplot(3,1,3)
plot(x_deflection, deflection*1000, 'LineWidth', 2);  % deflection in mm
xlabel('Position along beam (m)');
ylabel('Deflection (mm)');
title('Beam Deflection Diagram');
grid on;
set(gca,'YDir','reverse');  % Positive downward deflection
[max_def, max_idx] = max(abs(deflection));
fprintf('Maximum deflection: %.2f mm at x = %.2f m\n', max_def*1000, x_deflection(max_idx));
% Functions for deflection calculation
function [x_deflection, deflection] = calculateDeflection(L, pointLoads, udls, Ra, Rb, E, I)
    n = 1000;
    x = linspace(0, L, n);
    M = zeros(1, n);
    for i = 1:n
        M(i) = calculateMomentAt(x(i), L, pointLoads, udls, Ra, Rb);
    end
    dx = L/(n-1);
    slope = zeros(1, n);
    for i = 2:n
        slope(i) = slope(i-1) + (M(i-1)/(E*I) + M(i)/(E*I))*dx/2;
    end
    deflection = zeros(1, n);
    for i = 2:n
        deflection(i) = deflection(i-1) + (slope(i-1) + slope(i))*dx/2;
    end
    % Adjust deflection to satisfy boundary conditions
    deflection_at_end = deflection(end);
    correction = linspace(0, deflection_at_end, n);
    deflection = deflection - correction;
    x_deflection = x;
end
function M = calculateMomentAt(x, L, pointLoads, udls, Ra, Rb)
    M = Ra * x;
    for i = 1:size(pointLoads,1)
        P = pointLoads(i,1);
        a = pointLoads(i,2);
        if x > a
            M = M - P*(x - a);
        end
    end
    for i = 1:size(udls,1)
        w = udls(i,1);
        start_pos = udls(i,2);
        end_pos = udls(i,3);
        if x > start_pos
            if x <= end_pos
                length_under_load = x - start_pos;
                M = M - w * length_under_load * (length_under_load/2);
            else
                length_of_udl = end_pos - start_pos;
                distance_to_center = x - (start_pos + length_of_udl/2);
                M = M - w * length_of_udl * distance_to_center;
            end
        end
    end
end
