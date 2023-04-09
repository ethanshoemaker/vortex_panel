%This is a program hat uses the vortex panel method to determine lifing
%characteristics on an airfoil given an input file with airfoil boundary
%point coordinates, freestream velociy, air density, and angle of attack

%Author: Ethan Shoemaker
%April 2022

%%%IMPORTANT NOTE: THIS PROGRAM ASSUMES AIRFOIL COORDINATES BEGINNING AT
%%%THE TRAILING EDGE AND TRAVELING CLOCKWISE. FOR COUNTERCLOCKWISE TRAVEL,
%%%UNCOMMENT LINES 32-33

clear
clc
close all

tic

%Read text file inputs
data = csvread('NACA2412.txt');
rho = data(1,1);
v_inf = data(1,2);
alpha = data(2,1);
m = data(2,2);
xb = data(3:size(data,1),1);
yb = data(3:size(data,1),2);
xb(m+1) = xb(1);
yb(m+1) = yb(1);


xb = flip(xb);
yb = flip(yb);

transition = find(xb==0); %finds index of boundary points where transition from bottom of airfoil to top of airfoil occurs

%establish empty vectors
cp_top = [];
cp_bottom = [];
A = [];
B = [];
C = [];
D = [];
E = [];
F = [];
G = [];
P = [];
Q = [];
C_n1 = [];
C_n2 = [];
C_t1 = [];
C_t2 = [];
theta = [];
s = [];
xc = [];
yc = [];
A_n = [];
A_t = [];
RHS = [];

%midpoint formula to find control points
for i = 1:m
    xc_new = (xb(i+1)+xb(i))/2;
    xc = [xc xc_new];
    yc_new = (yb(i+1)+yb(i))/2;
    yc = [yc yc_new];
end

%find theta values
for i = 1:m
    y = yb(i+1)-yb(i);
    x = xb(i+1)-xb(i);
    theta_new = atan2d(y,x);
    if x == 0
        theta_new = 0;
    end
    theta = [theta theta_new];
end

%distance formula to find panel lengths
for i = 1:m
    s_new = sqrt((xb(i+1)-xb(i))^2+(yb(i+1)-yb(i))^2);
    s = [s s_new];
end

%create constant matrices. Formulas defined in Kueth and Chow book
for i = 1:m
    for j = 1:m
        if i == j
            A(i,j) = 0;
            B(i,j) = 0;
            C(i,j) = 0;
            D(i,j) = 0;
            E(i,j) = 0;
            F(i,j) = 0;
            G(i,j) = 0;
            P(i,j) = 0;
            Q(i,j) = 0;
            C_n2(i,j) = 1;
            C_n1(i,j) = -1;
            C_t1(i,j) = pi/2;
            C_t2(i,j) = pi/2;
        else
            A_new = (-(xc(i)-xb(j))*cosd(theta(j)))-((yc(i)-yb(j))*sind(theta(j)));
            A(i,j) = A_new;
            B_new = (xc(i)-xb(j))^2 + (yc(i)-yb(j))^2;
            B(i,j) = B_new;
            C_new = sind(theta(i)-theta(j));
            C(i,j) = C_new;
            D_new = cosd(theta(i)-theta(j));
            D(i,j) = D_new;
            E_new = ((xc(i)-xb(j))*sind(theta(j)))-((yc(i)-yb(j))*cosd(theta(j)));
            E(i,j) = E_new;
            F_new = log(1+((s(j)^2+2*A_new*s(j))/(B_new)));
            F(i,j) = F_new;
            G_new = atan2(E_new*s(j),B_new+A_new*s(j));
            G(i,j) = G_new;
            P_new = ((xc(i)-xb(j))*sind(theta(i)-2*theta(j))) + ((yc(i)-yb(j))*cosd(theta(i)-2*theta(j)));
            P(i,j) = P_new;
            Q_new = ((xc(i)-xb(j))*cosd(theta(i)-2*theta(j))) - ((yc(i)-yb(j))*sind(theta(i)-2*theta(j)));
            Q(i,j) = Q_new;
            C_n2(i,j) = D_new + 0.5*Q_new*F_new/s(j) - (A_new*C_new+D_new*E_new)*G_new/s(j);
            C_n1(i,j) = 0.5*D_new*F_new + C_new*G_new - C_n2(i,j);
            C_t2(i,j) = C_new + 0.5*P_new*F_new/s(j) + (A_new*D_new - C_new*E_new)*G_new/s(j);
            C_t1(i,j) = 0.5*C_new*F_new - D_new*G_new - C_t2(i,j);
        end
    end
end

%create A matrices from C matrices
for i = 1:m+1
    for j =1:m+1
        if i < m+1
            RHS(i) = sind(theta(i)-alpha);
            if j == 1
                A_n(i,j) = C_n1(i,j);
                A_t(i,j) = C_t1(i,j);
            elseif j > 1 && j<m+1
                A_n(i,j) = C_n1(i,j)+C_n2(i,j-1);
                A_t(i,j) = C_t1(i,j) + C_t2(i,j-1);
            elseif j == m+1
                A_n(i,j) = C_n2(i,m);
                A_t(i,j) = C_t2(i,m);
            else
                print('error')
            end
        elseif i == m+1
            A_n(i,1) = 1;
            A_n(i,m+1) = 1;
            RHS(i) = 0;
        else
            print('error')
        end
    end
end


%solve set of linear equations to find gamma
gamma_prime = A_n\RHS';
gamma = gamma_prime*2*pi*v_inf;

sum1 = [];
sum2 = 0;

%solve v_i and cp_i
for i = 1:m
    for j = 1:m+1
        new = A_t(i,j)*gamma_prime(j);
        sum2 = sum2 + new;
    end
    sum1(i) = sum2;
    sum2 = 0;
end

for i = 1:m
    v_i(i) = cosd(theta(i)-alpha) + sum1(i);
    c_pi(i) = 1 - v_i(i)^2;
end

%compute panel strengths given gamma at each boundary point
panel_strength = [];
for i = 1:m
    panel_strength(i) = gamma(i)*s(i) + (((gamma(i+1)-gamma(i))*s(i)^2)/(2*s(i)));
end

%compute circulation, lift per unit span, lift coefficient, c_p for top and
%bottom
circulation = sum(panel_strength);
L_prime = rho*v_inf*circulation;
lift_coefficient = L_prime/(0.5*rho*v_inf^2);
cp_bottom = c_pi(1:transition-1);
cp_top = c_pi(transition:m);

%plots
figure();
plot(xb,yb)
ylim([0 1])
hold on
scatter(xc,yc)
hold off

c_l_data = [-0.4 -0.2 0 0.22 0.45 0.63 0.85 1.08 1.26 1.41 1.55 1.69 1.58 1.18];

figure();
plot(flip(xc(1:transition-1)),flip(cp_bottom),'-o')
hold on
plot(xc(transition:m),cp_top,'-o')
title('Pressure Coefficient vs % Chord')
xlabel('% Chord')
ylabel('Pressure Coefficient')
legend('Bottom of airfoil','Top of airfoil')
set(gca, 'YDir','reverse')
hold off


% alpha = [-6 -4 -2 0 2 4 6 8 10 12 14 16 18 20];
% c_l = L_prime/(0.5*rho*v_inf^2)
% figure();
% plot(alpha,c_l)
% hold on
% scatter(alpha,c_l_data)
% title('Predicted Lift Coefficient and Real-World Data vs Angle of Attack')
% xlabel('Angle of Attack (degrees)')
% ylabel('Lift Coefficient')
% legend('Vortex Panel Prediction','Real-World Data')
% hold off
% 
% 
% 
% figure();
% scatter(flip(xc(1:transition-1)),flip(cp_bottom(3,:)),'filled')
% hold on
% scatter(xc(transition:m),cp_top(3,:),'filled')
% scatter(flip(xc(1:transition-1)),flip(cp_bottom(2,:)),'d','filled')
% scatter(xc(transition:m),cp_top(2,:),'d','filled')
% scatter(flip(xc(1:transition-1)),flip(cp_bottom(1,:)),'s','filled')
% scatter(xc(transition:m),cp_top(1,:),'s','filled')
% set(gca, 'YDir','reverse')
% legend('c_p top (alpha = 10°)','c_p bottom (alpha = 10°)','c_p top (alpha = 5°)','c_p bottom (alpha = 5°)','c_p top (alpha = 0°)','c_p bottom (alpha = 0°)')
% title('Pressure Coefficient vs Angle of Attack')
% xlabel('Percent of Total Chord Length')
% ylabel('Pressure Coefficient')
% hold off


%%%outputs to command line%%%
fprintf('Freestream Velocity: %.2f m/s\n',v_inf)
fprintf('Air Density: %.3f kg/m^3\n',rho)
fprintf('Angle of Attack: %.2f degrees\n',alpha)
fprintf('Number of Panels: %d\n\n',m)
fprintf('Total Circulation: %.3f m^2/s\n',circulation)
fprintf('Lift per unit span: %.3f N/m\n',L_prime)
fprintf('Lift Coefficient: %.3f\n\n\n', lift_coefficient)

fprintf('Geometry Summary\n\n')
fprintf('\t    Control Points\t      Boundary Points       theta\n\n')
fprintf('Point |\t  xc\t    yc\t  xc/chord |\txb\t  yb     |  radians\n')
fprintf('----------------------------------------------------------------------\n')
for i = 1:9
    fprintf(' %d    |  %+.3f   %+.3f   %+.3f  |   %+.3f   %+.3f   |  %+.3f\n',i,xc(i),yc(i),xc(i),xb(i),yb(i),theta(i)*(pi/180))
end
for i = 10:m
    fprintf('%d    |  %+.3f   %+.3f   %+.3f  |   %+.3f   %+.3f   |  %+.3f\n',i,xc(i),yc(i),xc(i),xb(i),yb(i),theta(i)*(pi/180))
end
fprintf('%d    |\t\t\t\t   |   %+.3f   %+.3f   |\n\n\n\n\n',m+1,xb(m+1),yb(m+1))



fprintf('Results Summary\n\n')
fprintf('      |       Control Points\t   |  Vortex   | Velocity  |   Pressure\n')
fprintf('      |                            | Strength  |           |  Coefficient\n')
fprintf('--------------------------------------------------------------------------\n')
fprintf('Panel |\t  xc\t    yc\t  xc/chord |   gamma   |    V_i    |    C_p\n')
fprintf('--------------------------------------------------------------------------\n')
for i = 1:9
    fprintf(' %d    |  %+.3f   %+.3f   %+.3f  |  %+.3f   |  %+.3f   |  %+.3f\n',i,xc(i),yc(i),xc(i),gamma_prime(i),v_i(i),c_pi(i))
end
for i = 10:m
    fprintf('%d    |  %+.3f   %+.3f   %+.3f  |  %+.3f   |  %+.3f   |  %+.3f\n',i,xc(i),yc(i),xc(i),gamma_prime(i),v_i(i),c_pi(i))
end
fprintf('      |\t\t\t\t   |  %+.3f   |           |',gamma_prime(m+1))




%text file output%
fileID1 = fopen('vortexpanel_output.txt','w');
fprintf(fileID1,'Freestream Velocity: %.2f m/s\n',v_inf);
fprintf(fileID1,'Air Density: %.3f kg/m^3\n',rho);
fprintf(fileID1,'Angle of Attack: %.2f degrees\n',alpha);
fprintf(fileID1,'Number of Panels: %d\n\n',m);
fprintf(fileID1,'Total Circulation: %.3f m^2/s\n',circulation);
fprintf(fileID1,'Lift per unit span: %.3f N/m\n',L_prime);
fprintf(fileID1,'Lift Coefficient: %.3f\n\n\n', lift_coefficient);

fprintf(fileID1,'Geometry Summary\n\n');
fprintf(fileID1,'\t    Control Points\t     Boundary Points       theta\n\n');
fprintf(fileID1,'Point |    xc      yc     xc/chord |    xb        yb     |  radians\n');
fprintf(fileID1,'----------------------------------------------------------------------\n');
for i = 1:9
    fprintf(fileID1,' %d    |  %+.3f   %+.3f   %+.3f  |   %+.3f   %+.3f   |  %+.3f\n',i,xc(i),yc(i),xc(i),xb(i),yb(i),theta(i)*(pi/180));
end
for i = 10:m
    fprintf(fileID1,'%d    |  %+.3f   %+.3f   %+.3f  |   %+.3f   %+.3f   |  %+.3f\n',i,xc(i),yc(i),xc(i),xb(i),yb(i),theta(i)*(pi/180));
end
fprintf(fileID1,'%d    |                            |   %+.3f   %+.3f   |\n\n\n\n\n',m+1,xb(m+1),yb(m+1));



fprintf(fileID1,'Results Summary\n\n');
fprintf(fileID1,'      |       Control Points       |  Vortex   | Velocity  |   Pressure\n');
fprintf(fileID1,'      |                            | Strength  |           |  Coefficient\n');
fprintf(fileID1,'--------------------------------------------------------------------------\n');
fprintf(fileID1,'Panel |    xc      yc     xc/chord |   gamma   |    V_i    |    C_p\n');
fprintf(fileID1,'--------------------------------------------------------------------------\n');
for i = 1:9
    fprintf(fileID1,' %d    |  %+.3f   %+.3f   %+.3f  |  %+.3f   |  %+.3f   |  %+.3f\n',i,xc(i),yc(i),xc(i),gamma_prime(i),v_i(i),c_pi(i));
end
for i = 10:m
    fprintf(fileID1,'%d    |  %+.3f   %+.3f   %+.3f  |  %+.3f   |  %+.3f   |  %+.3f\n',i,xc(i),yc(i),xc(i),gamma_prime(i),v_i(i),c_pi(i));
end
fprintf(fileID1,'      |                            |  %+.3f   |           |',gamma_prime(m+1));

fclose(fileID1);

%formatted for word%
fileID2 = fopen('vortexpanel_output_word.txt','w');
fprintf(fileID2,'Freestream Velocity: %.2f m/s\n',v_inf);
fprintf(fileID2,'Air Density: %.3f kg/m^3\n',rho);
fprintf(fileID2,'Angle of Attack: %.2f degrees\n',alpha);
fprintf(fileID2,'Number of Panels: %d\n\n',m);
fprintf(fileID2,'Total Circulation: %.3f m^2/s\n',circulation);
fprintf(fileID2,'Lift per unit span: %.3f N/m\n',L_prime);
fprintf(fileID2,'Lift Coefficient: %.3f\n\n\n', lift_coefficient);

fprintf(fileID2,'Geometry Summary\n\n');
fprintf(fileID2,'\t        Control Points\t          Boundary Points       theta\n\n');
fprintf(fileID2,'Point |\txc\t  yc\t   xc/chord |\txb\t    yb     |  radians\n');
fprintf(fileID2,'----------------------------------------------------------------------\n');
for i = 1:9
    fprintf(fileID2,' %d    |  %+.3f   %+.3f   %+.3f  |   %+.3f   %+.3f   |  %+.3f\n',i,xc(i),yc(i),xc(i),xb(i),yb(i),theta(i)*(pi/180));
end
for i = 10:m
    fprintf(fileID2,'%d    |  %+.3f   %+.3f   %+.3f  |   %+.3f   %+.3f   |  %+.3f\n',i,xc(i),yc(i),xc(i),xb(i),yb(i),theta(i)*(pi/180));
end
fprintf(fileID2,'%d    |                            |   %+.3f   %+.3f   |\n\n\n\n\n',m+1,xb(m+1),yb(m+1));



fprintf(fileID2,'Results Summary\n\n');
fprintf(fileID2,'      |       Control Points       |  Vortex   | Velocity  |  Pressure\n');
fprintf(fileID2,'      |                            | Strength  |           | Coefficient\n');
fprintf(fileID2,'-------------------------------------------------------------------------\n');
fprintf(fileID2,'Panel |   xc\t  yc     xc/chord |   gamma   |    V_i    |   C_p\n');
fprintf(fileID2,'-------------------------------------------------------------------------\n');
for i = 1:9
    fprintf(fileID2,' %d    |  %+.3f   %+.3f   %+.3f  |  %+.3f   |  %+.3f   |  %+.3f\n',i,xc(i),yc(i),xc(i),gamma_prime(i),v_i(i),c_pi(i));
end
for i = 10:m
    fprintf(fileID2,'%d    |  %+.3f   %+.3f   %+.3f  |  %+.3f   |  %+.3f   |  %+.3f\n',i,xc(i),yc(i),xc(i),gamma_prime(i),v_i(i),c_pi(i));
end
fprintf(fileID2,'      |                            |  %+.3f   |           |',gamma_prime(m+1));

fclose(fileID2);