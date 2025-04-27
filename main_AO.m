clear; close all; clc;

%% Set the Algebra
n = 3; % As we are in 3D
clifford_signature(n,0)

%% Set random seed for reproducibility
rng(100);

% Create an output folder for saving figures
outputFolder = 'output_AO_comparison';
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

%% Generate a 3D model point cloud
ptCloud = pcread('apple.ply');
% ptCloud = pcdownsample(ptCloud, 'random', 0.18);

model_Data = ptCloud.Location';

% Plot using scatter3 (transpose data so points are rows)
figure;
scatter3(model_Data(1,:), model_Data(2,:), model_Data(3,:), 5, 'filled');
xlabel('X'); ylabel('Y'); zlabel('Z');
title('Model Point Cloud');
axis equal;
grid on;

%% Simulation Settings
no_test_cases = 200;  % Number of random transformation tests

% --- Preallocate arrays for runtime and error (Frobenius norm error) ---
runtime_CM = zeros(1, no_test_cases);
runtime_CM_GAALOP = zeros(1, no_test_cases);
runtime_CM_GAALOPCSE = zeros(1, no_test_cases);
error_CM = zeros(1, no_test_cases);
error_opt = zeros(1, no_test_cases);
error_optCSE = zeros(1, no_test_cases);

% --- Preallocate cell arrays to store recovered rotation matrices ---
R_CM_all = cell(1, no_test_cases);
R_CM_GAALOP_all = cell(1, no_test_cases);
R_CM_GAALOPCSE_all = cell(1, no_test_cases);
R_GT_all = cell(1, no_test_cases);

%% Main Simulation Loop: Apply Random Transformations and Register
for k = 1:no_test_cases
    fprintf('Test Case: %d\n', k);

    % --- Generate random rotation angles (in degrees) and translation 
    % vector ---
    x_angle = randi([0, 360]);
    y_angle = randi([0, 360]);
    z_angle = randi([0, 360]);
    translation = [10*rand, 10*rand, 10*rand];

    % --- Obtain transformed (measured) data and ground truth transformation ---
    [measured_Data, ~, R_GT, ~] = transformation(model_Data, x_angle, y_angle, z_angle, translation);
    % measured_Data is 3xNpoints, Rotation_GT is 3x3
    
    % --- All Ground Truth Rotations ---
    R_GT_all{k} = R_GT;

    % --- Compute centroids for cross-correlation matrix ---
    c1 = mean(model_Data, 2)';  % 1x3 (model)
    c2 = mean(measured_Data, 2)'; % 1x3 (measured)

    % --- Form cross-correlation matrix H and auto-correlation matrix H2 ---
    H = (model_Data' - c1)' * (measured_Data' - c2); % H: 3x3
    H2 = (measured_Data' - c2)' * (measured_Data' - c2); % H2: 3x3

    % --- Calculate Reciprocals ---
    [F_clifford_up, G_clifford] = calculateReciprocals(n, double(H), double(H2));

    % --- Registration using the Original CM Method ---
    tic;
    [~, ~, R_CM] = CM_method(F_clifford_up, G_clifford);
    runtime_CM(k) = toc;

    % --- Registration using GAALOP---
    tic;
    [R_CM_GAALOP, ~] = CM_method_GAALOP(F_clifford_up, double(H2));
    runtime_CM_GAALOP(k) = toc;

    % --- Registration using GAALOP CSE---
    tic;
    [R_CM_GAALOPCSE, ~] = CM_method_GAALOP(F_clifford_up, double(H2));
    runtime_CM_GAALOPCSE(k) = toc;

    % --- Store the recovered rotations ---
    R_CM_all{k} = R_CM;
    R_CM_GAALOP_all{k} = R_CM_GAALOP;
    R_CM_GAALOPCSE_all{k} = R_CM_GAALOPCSE;

    % --- Compute error as Frobenius norm between recovered and ground
    % truth rotations ---
    error_CM(k) = norm(R_CM - R_GT, 'fro');
    error_opt(k) = norm(R_CM_GAALOP - R_GT, 'fro');
    error_optCSE(k) = norm(R_CM_GAALOPCSE - R_GT, 'fro');
end

%% Plots

% 1. Runtime and Rotation Error Comparison (Figure with two subplots)
hRuntimeError = figure('Name', 'RuntimeRotationComparison_apple', 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
set(0, 'DefaultAxesFontSize', 12);
set(0, 'DefaultTextFontSize', 12);

subplot(2,1,1);
plot(1:no_test_cases, runtime_CM, 'b--s', 'LineWidth', 1.5);
hold on;
plot(1:no_test_cases, runtime_CM_GAALOP, 'g-d', 'LineWidth', 1.5);
plot(1:no_test_cases, runtime_CM_GAALOPCSE, 'r-.*', 'LineWidth', 1.5);
xlabel('Test Case');
ylabel('Runtime (seconds)');
legend('CM Method', 'CM GAALOP', 'CM GAALOP CSE', 'Location', 'best');
title('Runtime Comparison');
grid on;

subplot(2,1,2);
plot(1:no_test_cases, error_CM, 'b--s', 'LineWidth', 1.5);
hold on;
plot(1:no_test_cases, error_opt, 'g-d', 'LineWidth', 1.5);
plot(1:no_test_cases, error_optCSE, 'r-.*', 'LineWidth', 1.5);
xlabel('Test Case');
ylabel('Frobenius Norm Error');
legend('CM Method', 'CM GAALOP', 'CM GAALOP CSE', 'Location', 'best');
title('Rotation Error Comparison');
grid on;

print(hRuntimeError, fullfile(outputFolder, 'Runtime_Rotation_Comparison_apple'), '-dpng', '-r300');


%% SUBPLOTS
% Set larger font size for readability (adjust as needed)
set(0, 'DefaultAxesFontSize', 12);
set(0, 'DefaultTextFontSize', 12);

% Create combined figure with two subplots
hEssentialPlots = figure('Name', 'Essential_AO_Runtimes_apple', 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]); % Adjusted size slightly

% --- Subplot 1: Log Scale Runtime Comparison ---
subplot(2,1,1); % Creates axes in the 1st position of a 2x1 grid
semilogy(1:no_test_cases, runtime_CM, 'b--s', 'LineWidth', 1.5);
hold on;
semilogy(1:no_test_cases, runtime_CM_GAALOP, 'g-d', 'LineWidth', 1.5);
semilogy(1:no_test_cases, runtime_CM_GAALOPCSE, 'r-.*', 'LineWidth', 1.5);
xlabel('Test Case');
ylabel('Runtime (s, log scale)');
% Calculate mean times 
mean_CM = mean(runtime_CM);
mean_GAALOP = mean(runtime_CM_GAALOP);
mean_GAALOPCSE = mean(runtime_CM_GAALOPCSE);
legend({sprintf('CM Method (mean = %.6f s)', mean_CM), ...
        sprintf('CM GAALOP (mean = %.6f s)', mean_GAALOP), ...
        sprintf('CM GAALOP CSE (mean = %.6f s)', mean_GAALOPCSE)}, 'Location', 'best');
grid on;
title('Runtime Comparison (Log Scale)'); % Add title for subplot clarity

% --- Subplot 2: Scatter Plot of Runtime vs. Error ---
subplot(2,1,2); % Creates axes in the 2nd position of a 2x1 grid
plot(runtime_CM, error_CM, 'b--s', 'LineWidth', 1.5, 'MarkerFaceColor', 'b', 'LineStyle', 'none');
hold on;
plot(runtime_CM_GAALOP, error_opt, 'g-d', 'LineWidth', 1.5, 'MarkerFaceColor', 'g', 'LineStyle', 'none');
plot(runtime_CM_GAALOPCSE, error_optCSE, 'r-.*', 'LineWidth', 1.5, 'MarkerFaceColor', 'r', 'LineStyle', 'none');
xlabel('Runtime (s)');
ylabel('Frobenius Norm Error');
legend('CM Method', 'CM GAALOP', 'CM GAALOP CSE', 'Location', 'best');
grid on;
title('Runtime vs Error Scatter Plot'); % Add title for subplot clarity

% Optional: Add overall title to the figure
% sgtitle('Key Performance Metrics for Absolute Orientation Routine'); 

% Save the combined figure with high resolution
print(hEssentialPlots, fullfile(outputFolder, 'Essential_AO_Runtime_Plots_apple'), '-dpng', '-r300'); % Changed filename



%% ----------------- HELPER FUNCTIONS -----------------
function [transformed_data, transform_matrix, Rotation_GT, Translation_GT] = transformation(data, x, y, z, translation)
% Ensure data is 3xN (if not, transpose)
if size(data,1) ~= 3
    data = data';
end

% Convert angles to radians
x = deg2rad(x);
y = deg2rad(y);
z = deg2rad(z);

% Rotation matrices about X, Y, Z axes
Rx = [1 0 0; 0 cos(x) -sin(x); 0 sin(x) cos(x)];
Ry = [cos(y) 0 sin(y); 0 1 0; -sin(y) 0 cos(y)];
Rz = [cos(z) -sin(z) 0; sin(z) cos(z) 0; 0 0 1];

% Combined rotation
Rotation_GT = Rx * Ry * Rz;

% Homogeneous transformation matrix
transform_matrix = [Rotation_GT, translation(:); 0 0 0 1];
Translation_GT = transform_matrix(1:3,4);

% Transform the data
num_points = size(data,2);
homog_data = [data; ones(1, num_points)];
transformed_homog = transform_matrix * homog_data;
transformed_data = transformed_homog(1:3,:);
end


function [F_clifford_up, G_clifford] = calculateReciprocals(n, F_matrix, G_matrix)
% This modified version expects F_matrix (e.g. H) as the first input and
% G_matrix (e.g. H2) as the second input.
%
% The reciprocal frame is computed for F_matrix and the Clifford conversion
% of G_matrix is also performed.

    % Convert inputs to double
    F_matrix = double(F_matrix);
    G_matrix = double(G_matrix);

    % Section 1: Convert Input Matrices into Clifford Vectors
    % Preallocate cell arrays for the Clifford representations.
    F_clifford = cell(1, n);
    G_clifford = cell(1, n);

    % Create a cell array for the Clifford basis vectors (assumed named e1, e2, ..., en)
    basis = cell(1, n);
    for j = 1:n
        basis{j} = eval(sprintf('e%d', j));
    end

    % Convert each column of F_matrix and G_matrix into a Clifford vector.
    for i = 1:n
        F_clifford{i} = 0;
        G_clifford{i} = 0;
        for j = 1:n
            F_clifford{i} = F_clifford{i} + F_matrix(j, i) * basis{j};
            G_clifford{i} = G_clifford{i} + G_matrix(j, i) * basis{j};
        end
    end

    % Section 2: Compute the Reciprocal Frame for F_matrix
    % 2.1: Compute the pseudoscalar of the F frame by wedging all F vectors.
    F_pseudo = F_clifford{1};
    for j = 2:n
        F_pseudo = wedge(F_pseudo, F_clifford{j});
    end

    % 2.2: Compute the unit pseudoscalar of the Clifford basis.
    I_basis = basis{1};
    for j = 2:n
        I_basis = wedge(I_basis, basis{j});
    end

    % 2.3: Compute each reciprocal vector for F.
    F_clifford_up = cell(1, n);
    for i = 1:n
        wedge_excluding_i = [];  % Holds the wedge product of all F_clifford{j} for j ≠ i.
        for j = 1:n
            if j ~= i
                if isempty(wedge_excluding_i)
                    wedge_excluding_i = F_clifford{j};
                else
                    wedge_excluding_i = wedge(wedge_excluding_i, F_clifford{j});
                end
            end
        end
        % Alternating sign and normalization by the pseudoscalar.
        F_clifford_up{i} = (-1)^i * (1 / abs(F_pseudo)) * wedge_excluding_i * I_basis;
    end
end




function [New_Rotor, New_Rotor_Reverse, Rnew] = CM_method(F_clifford_up, G_clifford)
% CM_method Computes the rotor using the Characteristic Multivector (CM) method.
%
% Inputs:
%   F_clifford_up - cell array containing reciprocal frame vectors of F:
%                   {vector1_Fclif_up, vector2_Fclif_up, vector3_Fclif_up}
%   G_clifford    - cell array containing Clifford vectors for G:
%                   {vector1_Gclif, vector2_Gclif, vector3_Gclif}

    % Unpack the cell arrays into individual vectors
    vector1_Fclif_up = F_clifford_up{1};
    vector2_Fclif_up = F_clifford_up{2};
    vector3_Fclif_up = F_clifford_up{3};

    vector1_Gclif = G_clifford{1};
    vector2_Gclif = G_clifford{2};
    vector3_Gclif = G_clifford{3};

    % Calculate the Invariants based on the reciprocal frame of F and G
    Inv1 = vector1_Fclif_up * vector1_Gclif + ...
           vector2_Fclif_up * vector2_Gclif + ...
           vector3_Fclif_up * vector3_Gclif;
       
    Inv2 = wedge(vector2_Fclif_up, vector1_Fclif_up) * wedge(vector1_Gclif, vector2_Gclif) + ...
           wedge(vector3_Fclif_up, vector1_Fclif_up) * wedge(vector1_Gclif, vector3_Gclif) + ...
           wedge(vector3_Fclif_up, vector2_Fclif_up) * wedge(vector2_Gclif, vector3_Gclif);
       
    Inv3 = wedge(vector3_Fclif_up, vector2_Fclif_up, vector1_Fclif_up) * ...
           wedge(vector1_Gclif, vector2_Gclif, vector3_Gclif);
       
    Inv_all = Inv1 + Inv2 + Inv3;

    % Compute the Rotor from the invariants
    Rotor_reverse = 1 + Inv_all;              % Rotor reverse (unnormalized)
    Rotor_reverse_normalized = unit(Rotor_reverse);  % Normalize the rotor reverse
    New_Rotor = reverse(Rotor_reverse_normalized);   % The rotor (normalized)
    New_Rotor_Reverse = Rotor_reverse_normalized;      % Its reverse

    % Check (should be identity)
    R_Rtilde_3D = New_Rotor * New_Rotor_Reverse; %#ok<NASGU>

    % Recover the Rotation Matrix from the Rotor
    f1_F = New_Rotor * e1 * New_Rotor_Reverse;
    f2_F = New_Rotor * e2 * New_Rotor_Reverse;
    f3_F = New_Rotor * e3 * New_Rotor_Reverse;

    coefficients_f1_F = cell2mat(coefficients(f1_F));
    coefficients_f2_F = cell2mat(coefficients(f2_F));
    coefficients_f3_F = cell2mat(coefficients(f3_F));

    % Pad each coefficient vector to length 3 if needed:
    F1 = padCoeffs(coefficients_f1_F);
    F2 = padCoeffs(coefficients_f2_F);
    F3 = padCoeffs(coefficients_f3_F);

    % Assemble the rotation matrix (each row is one of the three padded vectors)
    Rnew = [F1'; F2'; F3'];
    Rnew = reshape(Rnew, 3, 3)';
end


function [RotationMatrix, R] = CM_method_GAALOP(F_clifford_up, G_matrix)

    % Unpack the cell arrays into individual vectors
    vector1_Fclif_up = F_clifford_up{1};
    vector2_Fclif_up = F_clifford_up{2};
    vector3_Fclif_up = F_clifford_up{3};

    % Extract scalar coefficients for the reciprocal frame vectors:
    cr1_coeffs = cell2mat(coefficients(vector1_Fclif_up));
    cr2_coeffs = cell2mat(coefficients(vector2_Fclif_up));
    cr3_coeffs = cell2mat(coefficients(vector3_Fclif_up));
    cr11 = cr1_coeffs(1); cr12 = cr1_coeffs(2); cr13 = cr1_coeffs(3);
    cr21 = cr2_coeffs(1); cr22 = cr2_coeffs(2); cr23 = cr2_coeffs(3);
    cr31 = cr3_coeffs(1); cr32 = cr3_coeffs(2); cr33 = cr3_coeffs(3);

    % Extract components from G_matrix
    % Assume that the columns of G_matrix represent the multivector components
    d1 = G_matrix(:,1); d2 = G_matrix(:,2); d3 = G_matrix(:,3);
    d11 = d1(1); d21 = d1(2); d31 = d1(3);
    d12 = d2(1); d22 = d2(2); d32 = d2(3);
    d13 = d3(1); d23 = d3(2); d33 = d3(3);

	R(1) = ((((cr11 * cr22 - cr12 * cr21) * cr33 + (cr13 * cr21 - cr11 * cr23) * cr32 + (cr12 * cr23 - cr13 * cr22) * cr31) * d11 + cr22 * cr33 - cr23 * cr32) * d22 + (((cr12 * cr21 - cr11 * cr22) * cr33 + (cr11 * cr23 - cr13 * cr21) * cr32 + (cr13 * cr22 - cr12 * cr23) * cr31) * d12 + cr21 * cr33 - cr23 * cr31) * d21 + (cr12 * cr33 - cr13 * cr32) * d12 + (cr11 * cr33 - cr13 * cr31) * d11 + cr33) * d33 + ((((cr12 * cr21 - cr11 * cr22) * cr33 + (cr11 * cr23 - cr13 * cr21) * cr32 + (cr13 * cr22 - cr12 * cr23) * cr31) * d11 - cr22 * cr33 + cr23 * cr32) * d23 + (((cr11 * cr22 - cr12 * cr21) * cr33 + (cr13 * cr21 - cr11 * cr23) * cr32 + (cr12 * cr23 - cr13 * cr22) * cr31) * d13 + cr21 * cr32 - cr22 * cr31) * d21 + (cr13 * cr32 - cr12 * cr33) * d13 + (cr11 * cr32 - cr12 * cr31) * d11 + cr32) * d32 + ((((cr11 * cr22 - cr12 * cr21) * cr33 + (cr13 * cr21 - cr11 * cr23) * cr32 + (cr12 * cr23 - cr13 * cr22) * cr31) * d12 - cr21 * cr33 + cr23 * cr31) * d23 + (((cr12 * cr21 - cr11 * cr22) * cr33 + (cr11 * cr23 - cr13 * cr21) * cr32 + (cr13 * cr22 - cr12 * cr23) * cr31) * d13 - cr21 * cr32 + cr22 * cr31) * d22 + (cr13 * cr31 - cr11 * cr33) * d13 + (cr12 * cr31 - cr11 * cr32) * d12 + cr31) * d31 + ((cr12 * cr23 - cr13 * cr22) * d12 + (cr11 * cr23 - cr13 * cr21) * d11 + cr23) * d23 + ((cr13 * cr22 - cr12 * cr23) * d13 + (cr11 * cr22 - cr12 * cr21) * d11 + cr22) * d22 + ((cr13 * cr21 - cr11 * cr23) * d13 + (cr12 * cr21 - cr11 * cr22) * d12 + cr21) * d21 + cr13 * d13 + cr12 * d12 + cr11 * d11 + 1.0; % 1.0
	R(5) = ((cr23 * cr31 - cr21 * cr33) * d22 + (cr22 * cr33 - cr23 * cr32) * d21 + (cr13 * cr31 - cr11 * cr33) * d12 + (cr12 * cr33 - cr13 * cr32) * d11) * d33 + ((cr21 * cr33 - cr23 * cr31) * d23 + (cr11 * cr33 - cr13 * cr31) * d13 - cr31) * d32 + ((cr23 * cr32 - cr22 * cr33) * d23 + (cr13 * cr32 - cr12 * cr33) * d13 + cr32) * d31 + ((cr13 * cr21 - cr11 * cr23) * d12 + (cr12 * cr23 - cr13 * cr22) * d11) * d23 + ((cr11 * cr23 - cr13 * cr21) * d13 - cr21) * d22 + ((cr13 * cr22 - cr12 * cr23) * d13 + cr22) * d21 - cr11 * d12 + cr12 * d11; % e1 ^ e2
	R(6) = ((cr21 * cr32 - cr22 * cr31) * d22 + (cr11 * cr32 - cr12 * cr31) * d12 - cr31) * d33 + ((cr22 * cr31 - cr21 * cr32) * d23 + (cr23 * cr32 - cr22 * cr33) * d21 + (cr12 * cr31 - cr11 * cr32) * d13 + (cr13 * cr32 - cr12 * cr33) * d11) * d32 + ((cr22 * cr33 - cr23 * cr32) * d22 + (cr12 * cr33 - cr13 * cr32) * d12 + cr33) * d31 + ((cr11 * cr22 - cr12 * cr21) * d12 - cr21) * d23 + ((cr12 * cr21 - cr11 * cr22) * d13 + (cr13 * cr22 - cr12 * cr23) * d11) * d22 + ((cr12 * cr23 - cr13 * cr22) * d12 + cr23) * d21 - cr11 * d13 + cr13 * d11; % e1 ^ e3
	R(7) = ((cr22 * cr31 - cr21 * cr32) * d21 + (cr12 * cr31 - cr11 * cr32) * d11 - cr32) * d33 + ((cr21 * cr33 - cr23 * cr31) * d21 + (cr11 * cr33 - cr13 * cr31) * d11 + cr33) * d32 + ((cr21 * cr32 - cr22 * cr31) * d23 + (cr23 * cr31 - cr21 * cr33) * d22 + (cr11 * cr32 - cr12 * cr31) * d13 + (cr13 * cr31 - cr11 * cr33) * d12) * d31 + ((cr12 * cr21 - cr11 * cr22) * d11 - cr22) * d23 + ((cr11 * cr23 - cr13 * cr21) * d11 + cr23) * d22 + ((cr11 * cr22 - cr12 * cr21) * d13 + (cr13 * cr21 - cr11 * cr23) * d12) * d21 - cr12 * d13 + cr13 * d12; % e2 ^ e3
	f1_F(2) = (R(7) * R(7) - R(6) * R(6) - R(5) * R(5) + R(1) * R(1)) / (R(7) * R(7) + R(6) * R(6) + R(5) * R(5) + R(1) * R(1)); % e1
	f1_F(3) = (-(2.0 * R(6) * R(7) + 2.0 * R(1) * R(5))) / (R(7) * R(7) + R(6) * R(6) + R(5) * R(5) + R(1) * R(1)); % e2
	f1_F(4) = (2.0 * R(5) * R(7) - 2.0 * R(1) * R(6)) / (R(7) * R(7) + R(6) * R(6) + R(5) * R(5) + R(1) * R(1)); % e3
	f2_F(2) = (-(2.0 * R(6) * R(7) - 2.0 * R(1) * R(5))) / (R(7) * R(7) + R(6) * R(6) + R(5) * R(5) + R(1) * R(1)); % e1
	f2_F(3) = (-(R(7) * R(7) - R(6) * R(6) + R(5) * R(5) - R(1) * R(1))) / (R(7) * R(7) + R(6) * R(6) + R(5) * R(5) + R(1) * R(1)); % e2
	f2_F(4) = (-(2.0 * R(1) * R(7) + 2.0 * R(5) * R(6))) / (R(7) * R(7) + R(6) * R(6) + R(5) * R(5) + R(1) * R(1)); % e3
	f3_F(2) = (2.0 * R(5) * R(7) + 2.0 * R(1) * R(6)) / (R(7) * R(7) + R(6) * R(6) + R(5) * R(5) + R(1) * R(1)); % e1
	f3_F(3) = (2.0 * R(1) * R(7) - 2.0 * R(5) * R(6)) / (R(7) * R(7) + R(6) * R(6) + R(5) * R(5) + R(1) * R(1)); % e2
	f3_F(4) = (-(R(7) * R(7) + R(6) * R(6) - R(5) * R(5) - R(1) * R(1))) / (R(7) * R(7) + R(6) * R(6) + R(5) * R(5) + R(1) * R(1)); % e3

    % Assemble the rotation matrix 
    RotationMatrix = [f1_F(2:4), f2_F(2:4), f3_F(2:4)]';
    RotationMatrix = reshape(RotationMatrix, 3, 3);
end

function [RotationMatrix, R] = CM_method_GAALOPCSE(F_clifford_up, G_matrix)

    % Unpack the cell arrays into individual vectors
    vector1_Fclif_up = F_clifford_up{1};
    vector2_Fclif_up = F_clifford_up{2};
    vector3_Fclif_up = F_clifford_up{3};

    % Extract scalar coefficients for the reciprocal frame vectors:
    cr1_coeffs = cell2mat(coefficients(vector1_Fclif_up));
    cr2_coeffs = cell2mat(coefficients(vector2_Fclif_up));
    cr3_coeffs = cell2mat(coefficients(vector3_Fclif_up));
    cr11 = cr1_coeffs(1); cr12 = cr1_coeffs(2); cr13 = cr1_coeffs(3);
    cr21 = cr2_coeffs(1); cr22 = cr2_coeffs(2); cr23 = cr2_coeffs(3);
    cr31 = cr3_coeffs(1); cr32 = cr3_coeffs(2); cr33 = cr3_coeffs(3);

    % Extract components from G_matrix
    % Assume that the columns of G_matrix represent the multivector components
    d1 = G_matrix(:,1); d2 = G_matrix(:,2); d3 = G_matrix(:,3);
    d11 = d1(1); d21 = d1(2); d31 = d1(3);
    d12 = d2(1); d22 = d2(2); d32 = d2(3);
    d13 = d3(1); d23 = d3(2); d33 = d3(3);

	temp_gcse_1 = cr11 * cr33 - cr13 * cr31;
	temp_gcse_2 = (cr11 * cr22 - cr12 * cr21) * cr33 + (cr13 * cr21 - cr11 * cr23) * cr32;
	temp_gcse_4 = cr21 * cr33;
	temp_gcse_5 = cr21 * cr32;
	temp_gcse_6 = cr23 * cr32;
	temp_gcse_7 = cr23 * cr31;
	temp_gcse_9 = cr11 * cr22 - cr12 * cr21;
	temp_gcse_11 = cr11 * cr32 - cr12 * cr31;
	temp_gcse_12 = (cr11 * cr23 - cr13 * cr21) * d11;
	temp_gcse_14 = (cr13 * cr21 - cr11 * cr23) * d13;
	temp_gcse_18 = temp_gcse_1 * d11;
	temp_gcse_20 = cr13 * cr32;
	temp_gcse_21 = (-temp_gcse_1) * d13;
	temp_gcse_23 = (temp_gcse_2 + (cr12 * cr23 - cr13 * cr22) * cr31) * d11 + cr22 * cr33 - temp_gcse_6;
	temp_gcse_25 = (cr12 * cr23 - cr13 * cr22) * cr31;
	temp_gcse_29 = (-temp_gcse_11) * d12;
	temp_gcse_31 = temp_gcse_11 * d11;
	temp_gcse_33 = temp_gcse_9 * d11;
	temp_gcse_34 = (-temp_gcse_9) * d12;
	temp_gcse_38 = temp_gcse_4 - temp_gcse_7;
	temp_gcse_40 = (temp_gcse_2 + temp_gcse_25) * d13 + temp_gcse_5 - cr22 * cr31;
	temp_gcse_42 = temp_gcse_5 - cr22 * cr31;
	temp_gcse_44 = cr22 * cr33;
	temp_gcse_47 = cr13 * cr22;
	temp_gcse_51 = cr12 * cr33;
	temp_gcse_55 = cr12 * cr23;
	temp_gcse_57 = temp_gcse_51 - temp_gcse_20;
	temp_gcse_59 = ((-temp_gcse_2) + (-temp_gcse_25)) * d12 + temp_gcse_38;
	temp_gcse_61 = (temp_gcse_47 - temp_gcse_55) * d13;
	temp_gcse_62 = (temp_gcse_55 - temp_gcse_47) * d12;
	temp_gcse_63 = cr22 * cr33 - temp_gcse_6;
	temp_gcse_74 = temp_gcse_57 * d12;
	temp_gcse_75 = (-temp_gcse_57) * d13;
	R(1) = (temp_gcse_23 * d22 + temp_gcse_59 * d21 + temp_gcse_74 + temp_gcse_18 + cr33) * d33 + ((-temp_gcse_23) * d23 + temp_gcse_40 * d21 + temp_gcse_75 + temp_gcse_31 + cr32) * d32 + ((-temp_gcse_59) * d23 + (-temp_gcse_40) * d22 + temp_gcse_21 + temp_gcse_29 + cr31) * d31 + (temp_gcse_62 + temp_gcse_12 + cr23) * d23 + (temp_gcse_61 + temp_gcse_33 + cr22) * d22 + (temp_gcse_14 + temp_gcse_34 + cr21) * d21 + cr13 * d13 + cr12 * d12 + cr11 * d11 + 1.0; % 1.0
	temp_gcse_8 = (temp_gcse_7 - temp_gcse_4) * d22;
	temp_gcse_13 = (cr13 * cr21 - cr11 * cr23) * d12;
	temp_gcse_22 = (-temp_gcse_1) * d12;
	temp_gcse_54 = (temp_gcse_44 - temp_gcse_6) * d21;
	temp_gcse_60 = (temp_gcse_55 - temp_gcse_47) * d11;
	temp_gcse_76 = temp_gcse_57 * d11;
	R(5) = (temp_gcse_8 + temp_gcse_54 + temp_gcse_22 + temp_gcse_76) * d33 + (temp_gcse_38 * d23 + (-temp_gcse_21) - cr31) * d32 + ((-temp_gcse_63) * d23 + temp_gcse_75 + cr32) * d31 + (temp_gcse_13 + temp_gcse_60) * d23 + ((-temp_gcse_14) - cr21) * d22 + (temp_gcse_61 + cr22) * d21 - cr11 * d12 + cr12 * d11; % e1 ^ e2
	temp_gcse_3 = (cr22 * cr31 - cr21 * cr32) * d23;
	temp_gcse_30 = (-temp_gcse_11) * d13;
	temp_gcse_35 = (-temp_gcse_9) * d13;
	R(6) = (temp_gcse_42 * d22 + (-temp_gcse_29) - cr31) * d33 + (temp_gcse_3 + (-temp_gcse_54) + temp_gcse_30 + (-temp_gcse_76)) * d32 + (temp_gcse_63 * d22 + temp_gcse_74 + cr33) * d31 + ((-temp_gcse_34) - cr21) * d23 + (temp_gcse_35 + (-temp_gcse_60)) * d22 + (temp_gcse_62 + cr23) * d21 - cr11 * d13 + cr13 * d11; % e1 ^ e3
	R(7) = ((-temp_gcse_42) * d21 + (-temp_gcse_31) - cr32) * d33 + (temp_gcse_38 * d21 + temp_gcse_18 + cr33) * d32 + ((-temp_gcse_3) + temp_gcse_8 + (-temp_gcse_30) + temp_gcse_22) * d31 + ((-temp_gcse_33) - cr22) * d23 + (temp_gcse_12 + cr23) * d22 + ((-temp_gcse_35) + temp_gcse_13) * d21 - cr12 * d13 + cr13 * d12; % e2 ^ e3
	temp_gcse_15 = R(7) * R(7) + R(6) * R(6) + R(5) * R(5) + R(1) * R(1);
	temp_gcse_37 = R(7) * R(7) - R(6) * R(6);
	temp_gcse_41 = R(5) * R(5);
	temp_gcse_43 = R(1) * R(1);
	temp_gcse_58 = R(6) * R(6);
	temp_gcse_66 = R(7) * R(7);
	f1_F(2) = (temp_gcse_37 - temp_gcse_41 + temp_gcse_43) / temp_gcse_15; % e1
	temp_gcse_27 = 2.0 * R(6) * R(7);
	temp_gcse_39 = 2.0 * R(1);
	temp_gcse_65 = temp_gcse_39 * R(5);
	f1_F(3) = (-(temp_gcse_27 + temp_gcse_65)) / temp_gcse_15; % e2
	temp_gcse_17 = 2.0 * R(5);
	temp_gcse_36 = temp_gcse_17 * R(7);
	temp_gcse_56 = temp_gcse_39 * R(6);
	f1_F(4) = (temp_gcse_36 - temp_gcse_56) / temp_gcse_15; % e3
	f2_F(2) = (-(temp_gcse_27 - temp_gcse_65)) / temp_gcse_15; % e1
	f2_F(3) = (-(temp_gcse_37 + temp_gcse_41 - temp_gcse_43)) / temp_gcse_15; % e2
	temp_gcse_28 = temp_gcse_17 * R(6);
	temp_gcse_52 = temp_gcse_39 * R(7);
	f2_F(4) = (-(temp_gcse_52 + temp_gcse_28)) / temp_gcse_15; % e3
	f3_F(2) = (temp_gcse_36 + temp_gcse_56) / temp_gcse_15; % e1
	f3_F(3) = (temp_gcse_52 - temp_gcse_28) / temp_gcse_15; % e2
	f3_F(4) = (-(temp_gcse_66 + temp_gcse_58 - temp_gcse_41 - temp_gcse_43)) / temp_gcse_15; % e3

    % Assemble the rotation matrix 
    RotationMatrix = [f1_F(2:4), f2_F(2:4), f3_F(2:4)]';
    RotationMatrix = reshape(RotationMatrix, 3, 3);
end

function padded = padCoeffs(vec)
    % Pads or truncates vec to a 3×1 column vector.
    vec = vec(:);
    n   = length(vec);
    
    if n < 3
        padded = [vec; zeros(3 - n, 1)];
    else
        padded = vec(1:3);
    end
end


