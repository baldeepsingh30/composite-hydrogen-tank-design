clear;
clc;

%Properties for Torayca T700 G Carbon fiber
Ef = 588e9; %Pa
Vf = 0.6;
vf = 0.21;
Gf = Ef/(2*(1+vf));
%Properties for Expoxy used in aerospace composites
Em = 3.5e9; %Pa
Vm = 0.4; % volumetric
vm = 0.34; %Poissons Ratio
Gm = Em/(2*(1+vm));

t1 = 0.000625;

ply_pattern = [0, 90, +45, -45, 0, 90]; %0, -15,0 ,0,15, 0

% Define number of times to repeat this pattern
N = 8;  % This gives 6 * N plies (e.g., 240 plies if N = 40)

P = 70e6; %Pa
r = 0.0762;%m
E_al = 68.9e9;%Pa
tal= 0.003;%m

alpha_1 = -1.1e-6;     % Updated from datasheet
alpha_2 = 30e-6;      % in 1/°C
alpha_12 = 0;
alpha = [alpha_1; alpha_2; alpha_12];  % Material coord

beta_1 = 0;
beta_2 = 0.02;
beta_12= 0;
beta = [beta_1;
    beta_2;
    beta_12];

DeltaT = -60;     % Example: thermal change (°C)
DeltaC = 0.2;    % Example: moisture content change




%Elastic Approach E1 calculation
E1_ROM = Ef * Vf + Em * (1 - Vf);

numerator = 2 * Em * Ef * Vf * (vf - vm)^2 * (1 - Vf);

denominator = Ef * (2 * vm^2 * Vf - vm + Vf * vm - Vf - 1) + ...
              Em * (-1 - 2 * Vf * vf^2 + vf - Vf * vf + 2 * vf^2 + Vf);

E1 = E1_ROM - (numerator/denominator);

%************************Elastic Approach E2 calculation******************

nm = 3 - 4 * vm;
nf = 3 - 4 * vf;

% Compute bulk moduli
Kf = Ef / (2 * (1 + vf) * (1 - 2 * vf));
Km = Em / (2 * (1 + vm) * (1 - 2 * vm));

% Compute effective bulk modulus K*
numerator = Km * (Kf + Gm) * Vm + Kf * (Km + Gm) * Vf;
denominator = (Kf + Gm) * Vm + (Km + Gm) * Vf;
K_star = numerator / denominator;

% Compute m
m = 1 + 4 * K_star * (vf^2) / E1;

% Compute A, B, C (from quadratic for G23)
Gf_Gm = Gf / Gm;

A = 3 * Vf * (1 - Vf)^2 * (Gf_Gm - 1) * (Gf_Gm + nf) + ...
    ((Gf_Gm * nm + nf * nm - Gf_Gm * Vf^3 * (Vf * nm * (Gf_Gm - 1) - Gf_Gm * nm - 1)));

B = -3 * Vf * (1 - Vf)^2 * (Gf_Gm - 1) * (Gf_Gm + nf) + ...
    0.5 * (nm * Gf_Gm + Gf_Gm) * Vf + ...
    (nm - 1) * (Gf_Gm + nf) - ...
    2 * (Gf_Gm * nm - Gf_Gm * nf) * Vf^3 + ...
    (Vf / 2) * (nm + 1) * (Gf_Gm - 1) * ...
    (Gf_Gm + nf + Gf_Gm * (nm - nf)) * Vf^3;

C = 3 * Vf * (1 - Vf)^2 * (Gf_Gm - 1) * (Gf_Gm + nf) + ...
    (nm * Gf_Gm + (Gf_Gm - 1) * Vf + 1) * ...
    (Gf_Gm + nf + Gf_Gm * (nm - nf)) * Vf^3;

% Solve the quadratic A*x^2 + 2B*x + C = 0 for x = G23/Gm
coeffs = [A, 2*B, C];
roots_x = roots(coeffs);

% Choose positive, real root for G23
G23 = Gm * max(real(roots_x(real(roots_x) > 0)));

% Compute Poisson’s ratio ν23
v23 = (K_star - m * G23) / (K_star + m * G23);

% Compute E2
E2 = 2 * (1 + v23) * G23;

% Display results
fprintf('G23 = %.2f GPa\n', G23 / 1e9);
fprintf('ν23 = %.4f\n', v23);
fprintf('Transverse modulus E2 = %.2f GPa\n', E2 / 1e9);

% **********************Elasticity approach calculation of v12***********

numerator = Vf * Vm * (vf - vm) * ...
    (2 * Ef * vm^2 + vm * Ef - Ef + Em - Em * Vf - 2 * Em * vf^2);

denominator = (2 * vm^2 * Vf - vm + vm * Vf - 1 - Vf) * Ef + ...
              (2 * vf^2 - Vf * vf - 2 * Vf * vf^2 + Vf + vf - 1) * Em;

v12 = vf*Vf +vm*Vm + (numerator / denominator);

fprintf('Poisson''s ratio v12 = %.6f\n', v12);

%*********************Elasticity approach G12 calculations***************

G12 = Gm *((Gf*(1+Vf)+Gm*(1-Vf))/(Gf*(1-Vf)+Gm*(1+Vf)));

fprintf('Shear modulus G12= %.2f\n',G12);

v21 = v12 * E2 / E1;  % Ensure reciprocal relation
% Material properties
%E1 = 338e9;      % Pa
%E2 = 103e9;      % Pa
%G12 = 7.2e9;     % Pa
%v12 = 0.27;


% Reduced stiffness matrix [Q]
Q11 = E1 / (1 - v12 * v21);
Q12 = v12 * E2 / (1 - v12 * v21);
Q22 = E2 / (1 - v12 * v21);
Q66 = G12;

% Ply angles in degrees



% Build the full theta_deg vector
theta_deg = repmat(ply_pattern, 1, N);

% Count number of plies
nplies = length(theta_deg);

% Initialize 3D matrix to store Q_bar for each ply
Qbar_all = zeros(3, 3, nplies);

for i = 1:nplies
    theta = deg2rad(theta_deg(i));
    c = cos(theta);
    s = sin(theta);

    % Barred Q components from CLT (as in image)
    Q11_bar = Q11 * c^4 + Q22 * s^4 + 2*(Q12 + 2*Q66) * s^2 * c^2;
    Q12_bar = (Q11 + Q22 - 4*Q66) * s^2 * c^2 + Q12 * (s^4 + c^4);
    Q22_bar = Q11 * s^4 + Q22 * c^4 + 2*(Q12 + 2*Q66) * s^2 * c^2;
    Q16_bar = (Q11 - Q12 - 2*Q66) * c^3 * s - (Q22 - Q12 - 2*Q66) * s^3 * c;
    Q26_bar = (Q11 - Q12 - 2*Q66) * c * s^3 - (Q22 - Q12 - 2*Q66) * c^3 * s;
    Q66_bar = (Q11 + Q22 - 2*Q12 - 2*Q66) * s^2 * c^2 + Q66 * (s^4 + c^4);

    % Store the full Q̄ matrix
    Qbar_all(:,:,i) = [Q11_bar, Q12_bar, Q16_bar;
                       Q12_bar, Q22_bar, Q26_bar;
                       Q16_bar, Q26_bar, Q66_bar];
end

% Display result for the first ply
disp('Q̄ matrix for Ply 1 (Pa):');
for i = 1:nplies
disp(Qbar_all(:,:,i));
end
%****************************hcalculation*************************

n = nplies;                        % Number of plies

t_ply = t1;               % Ply thickness in meters
t = t_ply * ones(1, n);        % Vector of ply thicknesses

% Total laminate thickness
h = sum(t);

% Initialize surface position array
h_surface = zeros(1, n+1);     % h_0 to h_n

% Calculate h_0 (top surface)
h_surface(1) = -h/2;

% Loop through to calculate remaining surfaces
for k = 2:n+1
    h_surface(k) = h_surface(k-1) + t(k-1);   % bottom of ply k-1
end

% Display results
disp(['Total laminate thickness h = ', num2str(h), ' m']);
disp('z-positions of ply surfaces (from top h_0 to bottom h_n):');
disp(h_surface');

%**************************** A, B, D matrix calculation *************************

% Initialize A, B, and D matrices
A = zeros(3,3);
B = zeros(3,3);
D = zeros(3,3);

% Loop over each ply to compute contributions
for k = 1:n
    Qk = Qbar_all(:,:,k);             % Q̄ matrix for ply k
    z_k = h_surface(k+1);             % bottom surface of ply k
    z_k_1 = h_surface(k);             % top surface of ply k
    
    delta_z = z_k - z_k_1;
    
    A = A + Qk * delta_z;
    B = B + 0.5 * Qk * (z_k^2 - z_k_1^2);
    D = D + (1/3) * Qk * (z_k^3 - z_k_1^3);
end

% Display results
disp('A matrix [N/m]:');
disp(A);

disp('B matrix [N]:');
disp(B);

disp('D matrix [Nm]:');
disp(D);

%**************************** Inverses of A, B, D *************************

% Inverses of A, B, and D matrices
A_inv = inv(A);
B_inv = inv(B);
D_inv = inv(D);

% Display inverse matrices
disp('Inverse of A matrix [m/N]:');
disp(A_inv);

disp('Inverse of B matrix:');
disp(B_inv);

disp('Inverse of D matrix [1/Nm]:');
disp(D_inv);

% in plane engineering constant of a laminate

Ex = 1/(h*A_inv(1,1));
Ey = 1/(h*A_inv(2,2));
Gxy = 1/(h*A_inv(3,3));
vxy = - A_inv(1,2)/(A_inv(1,1));
vyx = - A_inv(1,2)/(A_inv(2,2));

%Flexural Engineering contants

Efx = 12/(h^3*D_inv(1,1));
Efy = 12/(h^3*D_inv(2,2));
Gfxy = 12/(h^3*D_inv(3,3));
vfxy = - D_inv(1,2)/(D_inv(1,1));
vfyx = - D_inv(1,2)/(D_inv(2,2));

% using the following equation calculated from hoop formula for thin wall
% pressure cylinder, which is Sigal = P*r/(tal+(E_comp * h/E_al))



%************stress(x) developed in alluminum linner due to pressure=70Mpa in cylinder ************

sigma_al = P*r/(tal+(Ex * h/E_al)); % should be less than sigma_al=(sigma_y)ult/SF which is sigma_al=255e6/2.5=1.02e8

%***************************Allowable stress for alluminum linner************
sigma_all=255e6/2;% allowable stress in x direction in alluminum

if sigma_al <= sigma_all
    fprintf('✅ Design is SAFE. sigma_al = %.2f MPa ≤ sigma_all = %.2f MPa\n', sigma_al/1e6, sigma_all/1e6);
else
    fprintf('❌ Design is NOT SAFE. sigma_al = %.2f MPa > sigma_all = %.2f MPa\n', sigma_al/1e6, sigma_all/1e6);
end





%****************Calculation for stress based on chapter 4 in Cabon fiber/epoxy laminate**************



% Assuming ply_pattern is in degrees (like [0 15 -15 ...])
T_inv = zeros(3, 3, nplies);
alpha_all = zeros(3, nplies);  % Store all α in global coord
alpha_new = zeros(3, nplies);

beta_all = zeros(3, nplies);  % Store all α in global coord
beta_new = zeros(3, nplies);
for i = 1:nplies

    theta = deg2rad(theta_deg(i));  % Convert degrees to radians
    c = cos(theta);
    s = sin(theta);

    % Transformation matrix for thermal expansion
    T_inv(:,:,i) = [c^2, s^2, -2*s*c;
                    s^2, c^2,  2*s*c;
                    s*c, -s*c, c^2 - s^2];

    % Apply transformation to α
    alpha_new(:,i) = T_inv(:,:,i) * alpha;
    alpha_all(:,i) = alpha_new(:,i);
 
    alpha_all(3,i) = alpha_new(3,i)/2;

    % Display
    fprintf('Alpha [x y xy] for Ply %d:\n', i);
    disp(alpha_all(:,i));

    beta_new(:,i) = T_inv(:,:,i) * beta;
    beta_all(:,i) = beta_new(:,i);
 
    beta_all(3,i) = beta_new(3,i)/2;

    % Display
    fprintf('beta [x y xy] for Ply %d:\n', i);
    disp(beta_all(:,i));

end

% Hygrothermal force/moment resultants



N_T = zeros(3,1);  % Thermal force
M_T = zeros(3,1);  % Thermal moment
N_C = zeros(3,1);  % Moisture force
M_C = zeros(3,1);  % Moisture moment

for k = 1:nplies

    Qk = Qbar_all(:,:,k);         % Q̄ for ply k
    alpha_k = alpha_all(:,k);     % Global α for ply k
    beta_k  = beta_all(:,k);      % Global β for ply k

    z_k = h_surface(k+1);         % bottom surface
    z_k_1 = h_surface(k);         % top surface

    delta_z = z_k - z_k_1;
    delta_z2 = z_k^2 - z_k_1^2;
    delta_z3 = z_k^3 - z_k_1^3;

    N_T = N_T + DeltaT * Qk * alpha_k * delta_z;
    M_T = M_T + 0.5 * DeltaT * Qk * alpha_k * delta_z2;

    N_C = N_C + DeltaC * Qk * beta_k * delta_z;
    M_C = M_C + 0.5 * DeltaC * Qk * beta_k * delta_z2;

end

% Display results
fprintf('\n--- Thermal Force Resultant N_T [N/m] ---\n');
disp(N_T);
fprintf('--- Thermal Moment Resultant M_T [Nm/m] ---\n');
disp(M_T);

fprintf('\n--- Moisture Force Resultant N_C [N/m] ---\n');
disp(N_C);
fprintf('--- Moisture Moment Resultant M_C [Nm/m] ---\n');
disp(M_C);

NM_T = [N_T; M_T]; %specificly for calculating hygrothermal properties

% Stack moisture force and moment
NM_C = [N_C; M_C]; %specificly for calculating hygrothermal properties

% Add them together
NMTC_total = NM_T + NM_C; %specificly for calculating hygrothermal properties

% Display the result
disp('Total Hygrothermal Force & Moment Resultant [N/m; Nm/m]:');
disp(NMTC_total);

ABD = [A,B;B,D];
disp(ABD);
ABD_inv = inv(ABD);
disp(ABD_inv);
%strain = zeros(3,1);
%*******************Hygrothermal Effect in a Cabon fiber/Epoxy Laminate **********************
%*************Starins in the laminate e0 and k*********************
strain = ABD_inv * NMTC_total;
disp('Starins in the mid plane fue to bending and curvature');
disp(strain);

e0 = [strain(1,1);
    strain(2,1);
    0]; %specificly for calculating hygrothermal properties


k = [0;
    0;
    strain(6,1)];%specificly for calculating hygrothermal properties


%*******************Hygrothermal Effect in a Cabon fiber/Epoxy Laminate **********************
%*************hygrothermal properties for bottom layer*********************
e_bottom = e0 + (h/2) * k;
disp('strain in x,y and xy for the bottom most layer is')
disp(e_bottom);

et_bottom = DeltaT * alpha_all(:,n);
disp('Thermal strain in x,y and xy for the bottom most layer is')
disp(et_bottom);

ec_bottom = DeltaC * beta_all(:,n);
disp('Moisture strain in x,y and xy for the bottom most layer is')
disp(ec_bottom);

em_bottom = e_bottom-et_bottom-ec_bottom;
disp('Mechanical strain in x,y and xy for the bottom most layer is')
disp(em_bottom);

sigma_bottom = Qbar_all(:,:,n) * em_bottom;
disp('stress in x,y and xy for the bottom most layer is')
disp(sigma_bottom);

%*************hygrothermal properties for top layer**********************
e_top = e0 + (-h/2) * k;
disp('strain in x,y and xy for the top most layer is')
disp(e_top);

et_top = DeltaT * alpha_all(:,1);
disp('Thermal strain in x,y and xy for the top most layer is')
disp(et_top);

ec_top = DeltaC * beta_all(:,1);
disp('Moisture strain in x,y and xy for the top most layer is')
disp(ec_top);

em_top = e_top-et_top-ec_top;
disp('Mechanical strain in x,y and xy for the top most layer is')
disp(em_top);

sigma_top = Qbar_all(:,:,1) * em_top;
disp('stress in x,y and xy for the top most layer is')
disp(sigma_top);

% *********** Strain Distribution Through Laminate Thickness no external load**************

% Use z-positions of ply interfaces (nplies + 1 values)
z_points = h_surface;  % Already defined from -h/2 to +h/2 at ply boundaries

% Preallocate strain array [εx; εy; γxy]
strain_through_laminate = zeros(3, length(z_points));

% Compute strain at each interface (z = h_surface(i))
for i = 1:length(z_points)
    z = z_points(i);
    strain_through_laminate(:,i) = e0 + z * k;
end


% *********** Stress Distribution applied Through Laminate Thickness no external load **************

% Use z-positions of ply interfaces (nplies + 1 values)
z_points = h_surface;

% Preallocate stress array [σx; σy; τxy]
stress_through_laminate = zeros(3, length(z_points));

% Compute stress using precomputed strain at each z_point
for i = 1:length(z_points)
    ply_index = min(i, nplies);  % Safely cap the index
    Qk = Qbar_all(:,:,ply_index);  % Use corresponding Qbar matrix
    
    % Compute stress
    stress = Qk * strain_through_laminate(:, i);
    
    % Store result
    stress_through_laminate(:, i) = stress;
end

%----------------chapter 4 finished----------------------------------------

%****************************Chafter 5*************************************

sigma_x_comp = (P * r) / (2 * (tal + (Ex * h / E_al)));
sigma_y_comp = (P * r) / (2 * (tal + (Ey * h / E_al)));
fprintf('⬆️sigma x composite  = %.2e N/m\n ', sigma_x_comp);
fprintf('⬆️sigma y composite = %.2e N/m\n', sigma_y_comp);

Nx = sigma_x_comp * h;

Ny = sigma_y_comp * h;

    
Nxy = 0;

Mx = 0;

My = 0;

Mxy = 0;

NM = [Nx;Ny;Nxy;Mx;My;Mxy];

%e0k = zeros(3);

e0k = ABD_inv * NM;

epsilon_0 = [e0k(1,1);e0k(2,1);e0k(3,1)];
kappa = [e0k(4,1);e0k(5,1);e0k(6,1)];

z_points = 0.5 * (h_surface(1:end-1) + h_surface(2:end));  % mid-plane of each ply

% Initialize strain arrays
strains_g = zeros(3, 1, nplies);  % global strain
strains_m = zeros(3, 1, nplies);  % mechanical strain
stresses_g = zeros(3, 1, nplies);   % resulting stress
stresses_l = zeros(3, 1, nplies);
strains_l = zeros(3,1, nplies);

T = zeros(3, 3, nplies);

R = [1,0,0;
    0,1,0;
    0,0,2]; % reuters matrix
R_inv = inv(R);

for k = 1:nplies
    z = z_points(k);  % mid-plane depth of kth ply

    % Equation (4.69): global strain at depth z
    strains_g(:,:,k) = epsilon_0 + z * kappa;

    % Only apply thermal and moisture strain for ply 1 (top surface)
   

    % Equation (4.70): mechanical strain in kth ply
    strains_m(:,:,k) = strains_g(:,:,k) ;

    % Equation (4.71): mechanical stress in kth ply
    Qk = Qbar_all(:,:,k);
    stresses_g(:,:,k) = Qk * strains_m(:,:,k);

theta = deg2rad(theta_deg(k));

    c = cos(theta);
    s = sin(theta);

    % Transformation matrix for thermal expansion
    T(:,:,k) = [c^2, s^2, 2*s*c;
                    s^2, c^2,  -2*s*c;
                    -s*c, s*c, c^2 - s^2];

    stresses_l(:,:,k) = T(:,:,k)*stresses_g(:,:,k);

   % strains_l (:,:,k) = R * T(:,:,k)* inv(R) * strains_g(:,:,k);
T_strain = [c^2, s^2, c*s;
            s^2, c^2, -c*s;
            -2*c*s, 2*c*s, c^2 - s^2];

strains_l(:,:,k) = T_strain * strains_g(:,:,k);
end


% ------------------- Tsai-Wu Strength Ratio Calculation -------------------

% Constants
sigma1T_u = 2.01e9;   % Longitudinal tensile
sigma1C_u = 0.79e9;   % Longitudinal compressive
sigma2T_u = 34e6;     % Transverse tensile (from 90° test)
sigma2C_u = 246e6;    % Transverse compressive (reasonable assumption)
shear12_u = 55e6;     % Shear strength
%et1_u = 0.015;       % Longitudinal tensile strain
%et2_u = 0.006;       % Transverse tensile strain
  % Shear strength tau12 [Pa]

% Packing-based transverse compressive strain
%d_s = (2*3^(1/3)*Vf/pi())^(1/2);%hexagonal packing
%ecm_u = sigmam_u / Em;

%ec2_u =((d_s)*(Em/Ef)+(1-(d_s)))*ecm_u;
%ef_ult = sigmaf_u / Ef;

% Ultimate stresses
%sigma1T_u = sigmaf_u * Vf + ef_ult * Em * (1 - Vf);
%sigma2T_u = E2 * et2_u;
%sigma1C_u = E1 * et2_u / v12;
%sigma2C_u = E2 * ec2_u;

% Tsai-Wu coefficients
H1 = 1/sigma1T_u - 1/sigma1C_u;
H2 = 1/sigma2T_u - 1/sigma2C_u;
H6 = 0;
H11 = 1 / (sigma1T_u * sigma1C_u);
H22 = 1 / (sigma2T_u * sigma2C_u);
H66 = 1 / (shear12_u^2);
H12 = -0.5 * sqrt(H11 * H22);

% Strength Ratio storage
SR_all = zeros(1, nplies);

for k = 1:nplies
    sigma1 = stresses_l(1,1,k);   % σ1 in local coords
    sigma2 = stresses_l(2,1,k);   % σ2
    tau12  = stresses_l(3,1,k);   % τ12

    % Coefficients of quadratic in SR
    A = H11 * sigma1^2 + H22 * sigma2^2 + H66 * tau12^2 + 2 * H12 * sigma1 * sigma2;
    B = H1 * sigma1 + H2 * sigma2 + H6 * tau12;

    % Solve for SR from Tsai-Wu failure criterion
    SR = (-B + sqrt(B^2 + 4*A)) / (2*A);
    SR_all(k) = SR;
end

% Display results
fprintf('\n--- Tsai-Wu Strength Ratios for Each Ply ---\n');
disp(SR_all');

% Find the minimum strength ratio and the index of the failing ply
[SR_min, min_index] = min(SR_all);

% Update Nx and Ny based on the minimum strength ratio
Nx = Nx * SR_min;
Ny = Ny *SR_min;

% Get ply angle that caused failure
failing_angle = theta_deg(min_index);

% Display information
fprintf('\n*** First Ply Failure ***\n');
fprintf('Minimum Strength Ratio (SR_min): %.4f\n', SR_min);
fprintf('Failing ply index: %d\n', min_index);
fprintf('Ply angle (deg): %d°\n', failing_angle);
fprintf('Updated Nx = %.2e N/m\n', Nx);
fprintf('Updated Ny = %.2e N/m\n', Ny);

Nx_fpf = Nx;
fprintf('\n⚠️ First Ply Failure Stress sigmax = %.4f\n', Nx_fpf);
ex0_fpf = epsilon_0(1);
fprintf('\n⚠️First Ply Failure Strain ex0 = %.4f\n', ex0_fpf);
Ny_fpf = Ny;
fprintf('\n⚠️ First Ply Failure Stress sigmay = %.4f\n', Ny_fpf);
ey0_fpf = epsilon_0(2);
fprintf('\n⚠️ First Ply Failure Strain ey0 = %.4f\n', ey0_fpf);
strain_x_history = ex0_fpf;
Nx_history = Nx_fpf;

% ------------------- Recalculate SR After Ply Failure -------------------

% Define the angle that has failed completely (e.g., 90°)
failed_angle_deg = failing_angle;

% Loop through all plies and zero Qbar if angle is failed
for k = 1:nplies
    if theta_deg(k) == failed_angle_deg
        Qbar_all(:,:,k) = zeros(3);  % Degrade completely
    end
end

% Recalculate A matrix with updated Qbar_all
A = zeros(3,3);
for k = 1:nplies
    z_k = h_surface(k+1);
    z_k_1 = h_surface(k);
    delta_z = z_k - z_k_1;
    A = A + Qbar_all(:,:,k) * delta_z;
end

B = zeros(3,3);
for k = 1:nplies
    z_k = h_surface(k+1);
    z_k_1 = h_surface(k);
    delta_z = z_k^2 - z_k_1^2;
    B = B + 0.5*Qbar_all(:,:,k) * delta_z;
end

D = zeros(3,3);
for k = 1:nplies
    z_k = h_surface(k+1);
    z_k_1 = h_surface(k);
    delta_z = z_k^3 - z_k_1^3;
    D = D + (1/3)*Qbar_all(:,:,k) * delta_z;
end


A_inv = inv(A);
B_inv = inv(B);
D_inv = inv(D);

% Recalculate mid-plane strain using updated loads
ABD = [A, B; B, D];
ABD_inv = inv(ABD);
NM = [Nx; Ny; 0; 0; 0; 0];

strain_k = ABD_inv * NM;
epsilon_0 = strain_k(1:3);
kappa = strain_k(4:6);

% Mid-plane z-points again
z_points = 0.5 * (h_surface(1:end-1) + h_surface(2:end));

% ------------------- Recalculate SR for Each Ply -------------------
SR_updated = ones(1, nplies);  % Default SR = 1 for failed plies

for k = 1:nplies
    if theta_deg(k) == failed_angle_deg
        continue;  % Skip already-failed ply
    end

    z = z_points(k);
    epsilon_g = epsilon_0 + z * kappa;

    theta = deg2rad(theta_deg(k));
    c = cos(theta);
    s = sin(theta);

    T = [c^2, s^2, 2*s*c;
         s^2, c^2, -2*s*c;
         -s*c, s*c, c^2 - s^2];

    T_strain = [c^2, s^2, c*s;
                s^2, c^2, -c*s;
                -2*c*s, 2*c*s, c^2 - s^2];

    epsilon_l = T_strain * epsilon_g;
    Qk = Qbar_all(:,:,k);
    sigma_g = Qk * epsilon_g;
    sigma_l = T * sigma_g;

    sigma1 = sigma_l(1);
    sigma2 = sigma_l(2);
    tau12 = sigma_l(3);

    % Tsai-Wu coefficients
    H1 = 1/sigma1T_u - 1/sigma1C_u;
    H2 = 1/sigma2T_u - 1/sigma2C_u;
    H6 = 0;
    H11 = 1/(sigma1T_u * sigma1C_u);
    H22 = 1/(sigma2T_u * sigma2C_u);
    H66 = 1/(shear12_u^2);
    H12 = -0.5 * sqrt(H11 * H22);

    A_tw = H11 * sigma1^2 + H22 * sigma2^2 + H66 * tau12^2 + 2 * H12 * sigma1 * sigma2;
    B_tw = H1 * sigma1 + H2 * sigma2;
    SR = (-B_tw + sqrt(B_tw^2 + 4*A_tw)) / (2*A_tw);

    SR_updated(k) = SR;
end

% ------------------- Display 2nd SR Values -------------------

fprintf('\n--- Updated Strength Ratios After %d° Ply Failure ---\n', failed_angle_deg);
for k = 1:nplies
    fprintf('Ply %3d (Angle %3d°):  SR = %.4f\n', k, theta_deg(k), SR_updated(k));
end

% Step 1: Remove plies that already failed (SR == 1 or Qbar == 0)
SR_filtered = SR_updated(SR_updated < 1 & SR_updated > 0);

% Step 2: Check if any plies are left to fail
if isempty(SR_filtered)
    fprintf('🛑 All plies have failed or no further failure possible.\n');

end

% Step 3: Get minimum SR (next ply to fail)
SR_min_next = min(SR_filtered);

% Step 4: Apply SR_min to increase applied load
Nx = Nx*SR_min_next ;
Ny = Ny*SR_min_next;

fprintf('\n⚠️ Next failing ply SR = %.6f\n', SR_min_next);
fprintf('⬆️ Updated Nx = %.2e N/m\n', Nx);
fprintf('⬆️ Updated Ny = %.2e N/m\n', Ny);

Nx_2pf = Nx;
fprintf('\n⚠️ second Ply Failure Stress sigmax = %.4f\n', Nx_2pf);
ex0_2pf = epsilon_0(1);
fprintf('\n⚠️ second Ply Failure Strain ex0 = %.4f\n', ex0_2pf);
Ny_2pf = Ny;
fprintf('\n⚠️ second Ply Failure load = %.4f\n', Ny_2pf);
ey0_2pf = epsilon_0(2);
fprintf('\n⚠️ second Ply Failure load= %.4f\n', ey0_2pf);
strain_x_history(end+1) = ex0_2pf;
Nx_history(end+1) = Nx_2pf;

%***********************3rd iteration for third angle*********************
% --- Third Iteration: Handle Failure of ±45° Plies ---

% Generalize: Mark plies with these angles as failed (degrade stiffness)
failed_angles_deg = [45, -45];

for k = 1:nplies
    if any(theta_deg(k) == failed_angles_deg)
        Qbar_all(:,:,k) = zeros(3);  % Fully degrade this ply
    end
end

% Recalculate A matrix
A = zeros(3,3);
for k = 1:nplies
    z_k = h_surface(k+1);
    z_k_1 = h_surface(k);
    A = A + Qbar_all(:,:,k) * (z_k - z_k_1);
end
A_inv = inv(A);

% Recompute mid-plane strain
ABD = [A, zeros(3); zeros(3), eye(3)];
ABD_inv = [A_inv, zeros(3); zeros(3), zeros(3)];
NM = [Nx; Ny; 0; 0; 0; 0];

strain_k = ABD_inv * NM;
epsilon_0 = strain_k(1:3);
kappa = strain_k(4:6);

% Compute SR values again
SR_third = ones(1, nplies);  % Default to 1 for failed plies
z_points = 0.5 * (h_surface(1:end-1) + h_surface(2:end));

for k = 1:nplies
    % Skip failed plies (Q = 0)
    if all(Qbar_all(:,:,k) == 0)
        continue;
    end

    z = z_points(k);
    epsilon_g = epsilon_0 + z * kappa;

    theta = deg2rad(theta_deg(k));
    c = cos(theta); s = sin(theta);

    % Transformation matrices
    T = [c^2, s^2, 2*s*c;
         s^2, c^2, -2*s*c;
         -s*c, s*c, c^2 - s^2];

    T_strain = [c^2, s^2, c*s;
                s^2, c^2, -c*s;
                -2*c*s, 2*c*s, c^2 - s^2];

    epsilon_l = T_strain * epsilon_g;
    Qk = Qbar_all(:,:,k);
    sigma_g = Qk * epsilon_g;
    sigma_l = T * sigma_g;

    sigma1 = sigma_l(1);
    sigma2 = sigma_l(2);
    tau12 = sigma_l(3);

    A_tw = H11 * sigma1^2 + H22 * sigma2^2 + H66 * tau12^2 + 2 * H12 * sigma1 * sigma2;
    B_tw = H1 * sigma1 + H2 * sigma2;

    if A_tw <= 0
        SR = 1;  % If no stiffness, treat as failed/inactive
    else
        SR = (-B_tw + sqrt(B_tw^2 + 4*A_tw)) / (2*A_tw);
    end

    SR_third(k) = SR;
end

% --- Display SRs ---
fprintf('\n--- Third Iteration: SRs for Remaining Plies ---\n');
for k = 1:nplies
    fprintf('Ply %3d (Angle %3d°):  SR = %.4f\n', k, theta_deg(k), SR_third(k));
end

% --- Filter SRs ---
SR_valid = SR_third(SR_third < 1 & SR_third > 0);

if isempty(SR_valid)
    fprintf('✅ All remaining plies are safe. No further failure detected.\n');
else
    SR_min_third = min(SR_valid);

    % Update loads for next iteration
    Nx = Nx*SR_min_third ;
    Ny = Ny*SR_min_third ;

    fprintf('\n⚠️ Third Iteration: Minimum SR = %.4f\n', SR_min_third);
    fprintf('⬆️ Updated Nx = %.2e N/m\n', Nx);
    fprintf('⬆️ Updated Ny = %.2e N/m\n', Ny);
end

Nx_Lpf = Nx;
fprintf('\n⚠️ Last Ply Failure Stress sigmax = %.4f\n', Nx_Lpf);
ex0_Lpf = epsilon_0(1);
fprintf('\n⚠️ Last Ply Failure Strain ex0 = %.4f\n', ex0_Lpf);
Ny_Lpf = Ny;
fprintf('\n⚠️ Last Ply Failure Load = %.4f\n', Ny_Lpf);
ey0_Lpf = epsilon_0(2);
fprintf('\n⚠️ Last Ply Failure load = %.4f\n', ey0_Lpf);

strain_x_history(end+1) = ex0_Lpf;
Nx_history(end+1) = Nx_Lpf;





%***********************4th iteration for fourth angle zero degrees*********************
% Degrade 0° plies: set Q = 0
final_failed_angle = 0;

for k = 1:nplies
    if theta_deg(k) == final_failed_angle
        Qbar_all(:,:,k) = zeros(3);  % Fully fail
    end
end

% Recompute A matrix with all failed plies
A = zeros(3,3);
for k = 1:nplies
    z_k = h_surface(k+1);
    z_k_1 = h_surface(k);
    A = A + Qbar_all(:,:,k) * (z_k - z_k_1);
end
A_inv = inv(A);

% Recompute laminate strain state
ABD = [A, zeros(3); zeros(3), eye(3)];
ABD_inv = [A_inv, zeros(3); zeros(3), zeros(3)];
NM = [Nx; Ny; 0; 0; 0; 0];

strain_k = ABD_inv * NM;
epsilon_0 = strain_k(1:3);
kappa = strain_k(4:6);

% Recompute SRs (should all be 1 or skipped now)
SR_final = ones(1, nplies);
z_points = 0.5 * (h_surface(1:end-1) + h_surface(2:end));

for k = 1:nplies
    if all(Qbar_all(:,:,k) == 0)
        continue;  % Skip already-failed
    end

    z = z_points(k);
    epsilon_g = epsilon_0 + z * kappa;

    theta = deg2rad(theta_deg(k));
    c = cos(theta);
    s = sin(theta);

    T = [c^2, s^2, 2*s*c;
         s^2, c^2, -2*s*c;
         -s*c, s*c, c^2 - s^2];

    T_strain = [c^2, s^2, c*s;
                s^2, c^2, -c*s;
                -2*c*s, 2*c*s, c^2 - s^2];

    epsilon_l = T_strain * epsilon_g;
    Qk = Qbar_all(:,:,k);
    sigma_g = Qk * epsilon_g;
    sigma_l = T * sigma_g;

    sigma1 = sigma_l(1);
    sigma2 = sigma_l(2);
    tau12 = sigma_l(3);

    A_tw = H11 * sigma1^2 + H22 * sigma2^2 + H66 * tau12^2 + 2 * H12 * sigma1 * sigma2;
    B_tw = H1 * sigma1 + H2 * sigma2;

    if A_tw <= 0
        SR = 1;  % Already failed
    else
        SR = (-B_tw + sqrt(B_tw^2 + 4*A_tw)) / (2*A_tw);
    end

    SR_final(k) = SR;
end

% Display SRs (final status)
fprintf('\n--- Final Iteration: 0° Ply SRs ---\n');
for k = 1:nplies
    fprintf('Ply %3d (Angle %3d°):  SR = %.4f\n', k, theta_deg(k), SR_final(k));
end

% Extract minimum SR from active plies
SR_last = SR_final(SR_final < 1 & SR_final > 0);

if isempty(SR_last)
    fprintf('\n✅ All 0° plies now failed. Design cannot carry further load.\n');
    fprintf('🛑 Last Ply Failure Load (Nx) = %.2e N/m\n', Nx);
    fprintf('🛑 Last Ply Failure Load (Ny) = %.2e N/m\n', Ny);
else
    SR_min_last = min(SR_last);

    % Update final failure load
    Nx = Nx*SR_min_last ;
    Ny = Ny*SR_min_last ;

    fprintf('\n⚠️ Final Ply SR = %.4f\n', SR_min_last);
    fprintf('🛑 Last Ply Failure Load (Nx) = %.2e N/m\n', Nx);
    fprintf('🛑 Last Ply Failure Load (Ny) = %.2e N/m\n', Ny);
end

% === Final Plot: Nx/h vs epsilon_x (Stress-Strain Curve) ===
figure;
plot(strain_x_history, Nx_history / h / 1e6, '-o', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Normal Strain \epsilon_x');
ylabel('Stress N_x/h [MPa]');
title('Stress-Strain in X-direction (Progressive Ply Failure)');
grid on;

% Add markers
text(strain_x_history(1), Nx_history(1)/h/1e6, ' First Ply Fail (90°)', 'VerticalAlignment', 'bottom');
text(strain_x_history(2), Nx_history(2)/h/1e6, ' ±45° Fail', 'VerticalAlignment', 'bottom');
text(strain_x_history(3), Nx_history(3)/h/1e6, ' Last Ply Fail (0°)', 'VerticalAlignment', 'bottom');

figure;
plot(strain_x_history, Nx_history / h / 1e6, '-o', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Normal Strain \epsilon_x');
ylabel('Stress N_x/h [MPa]');
title('Stress-Strain in X-direction (Progressive Ply Failure)');
grid on;

% Annotate each point with coordinates
for i = 1:length(strain_x_history)
    str_x = strain_x_history(i);
    str_y = Nx_history(i) / h / 1e6;
    label = sprintf('(%0.4g, %0.2f)', str_x, str_y);
    text(str_x, str_y, label, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'FontSize', 8);
end

% Optionally, add failure labels if desired
text(strain_x_history(1), Nx_history(1)/h/1e6, ' First Ply Fail (90°)', 'VerticalAlignment', 'top', 'FontWeight', 'bold');
text(strain_x_history(2), Nx_history(2)/h/1e6, ' ±45° Fail', 'VerticalAlignment', 'top', 'FontWeight', 'bold');
text(strain_x_history(3), Nx_history(3)/h/1e6, ' Last Ply Fail (0°)', 'VerticalAlignment', 'top', 'FontWeight', 'bold');
