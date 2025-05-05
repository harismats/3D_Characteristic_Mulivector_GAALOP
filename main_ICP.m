clear; close all; clc;
clifford_signature(3,0);

% Set the random seed for reproducibility.
rng(1);

% Define the different point cloud sizes to test.
pointCloudSizes = 1000:500:8000;
noTestCases = 5;    % Number of test cases per cloud size
noIteration = 100;    % Number of ICP iterations per test
numMethods = 4;      % ICP methods: 
                     %   1: SVD, 2: CM, 3: GAALOP, 4: GAALOP CSE,
         
% Preallocate arrays to store average execution times and errors.
numSizes = length(pointCloudSizes);
times = zeros(numSizes, numMethods);
errors = zeros(numSizes, numMethods);

% Loop over each point cloud size.
for idx = 1:numSizes
    N = pointCloudSizes(idx);
    fprintf('--- Testing Point Cloud with %d points ---\n', N);
    % Generate a synthetic point cloud using a helper function.
    ptCloud = generatePointCloud(N,2);
    
    % Preallocate temporary storage for this size.
    tempTime = zeros(noTestCases, numMethods);
    tempError = zeros(noTestCases, numMethods);
    
    for k = 1:noTestCases
        fprintf('Point cloud size: %d, Test Case: %d\n', N, k);
        
        % Generate random transformation parameters (angles in degrees and translation).
        angles = randi([0,30], 1, 3);
        translation = 0.5 * rand(1,3);
        
        % Apply transformation using the dedicated helper function.
        [movedData, Rotation_GT, Translation_GT] = transformation(ptCloud.Location, angles(1), angles(2), angles(3), translation);
        % Create the transformed point cloud; movedData is returned as an Nx3 array.
        ptCloud2 = pointCloud(movedData);
        
        % Preprocess (center) both source and target data.
        [srcData, tgtData] = preprocessData(ptCloud.Location, ptCloud2.Location);
        
        % Loop over each ICP method.
        for method = 1:numMethods
            [finalRot, elapsed, errVal] = runICP(srcData, tgtData, noIteration, Rotation_GT, Translation_GT, method);
            tempTime(k, method) = elapsed;
            tempError(k, method) = errVal;
        end
    end
    
    % Calculate and store average times and errors over test cases.
    times(idx,:) = mean(tempTime, 1);
    errors(idx,:) = mean(tempError, 1);
end


%% Create an output folder for saving figures.
outputFolder = 'output_ICP_comparison';
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

%% Define color and method name arrays for plotting.
allColors = {'k-o', 'b--s', 'g-d', 'r-.*'};
allMethodNames = {'SVD','CM','CM GAALOP','CM GAALOP CSE'};
% For percent difference plots (only methods 1, 3, and 4 will be used):
pdColors = {'k-', 'g-d', 'r-.*'};  % SVD as solid line (baseline), others as indicated.
pdMethodNames = {'SVD','GAALOP','GAALOP CM CSE'};

%% Plot 1: Average Execution Time vs. Number of Points (Linear Scale)
hRuntime = figure('Name', 'Runtime_Rotation_Comparison', 'Units','Normalized', 'OuterPosition', [0 0 1 1]);
set(0, 'DefaultAxesFontSize', 13);
set(0, 'DefaultTextFontSize', 12);
hold on;
for m = 1:numMethods
    plot(pointCloudSizes, times(:, m), allColors{m}, 'LineWidth', 1.5, 'MarkerSize', 8);
end
hold off;
xlabel('Number of Points', 'FontSize', 12);
ylabel('Average Execution Time (s)', 'FontSize', 12);
title('ICP Registration Time Comparison', 'FontSize', 14);
legend(allMethodNames, 'Location','best');
grid on;
print(hRuntime, fullfile(outputFolder, 'ICP_Execution_Time'), '-dpng', '-r300');

%% Plot 2: Average Rotation Error vs. Number of Points
hRotError = figure('Name', 'Final_Angular_Error_Comparison', 'Units','Normalized', 'OuterPosition', [0 0 1 1]);
set(0, 'DefaultAxesFontSize', 13);
set(0, 'DefaultTextFontSize', 12);
hold on;
for m = 1:numMethods
    plot(pointCloudSizes, errors(:, m), allColors{m}, 'LineWidth', 1.5, 'MarkerSize', 8);
end
hold off;
xlabel('Number of Points', 'FontSize', 12);
ylabel('Average Rotation Error (Frobenius Norm)', 'FontSize', 12);
title('ICP Registration Rotation Error Comparison', 'FontSize', 14);
legend(allMethodNames, 'Location','best');
grid on;
print(hRotError, fullfile(outputFolder, 'ICP_Rotation_Error'), '-dpng', '-r300');

%% Plot 3: Percent Difference in Execution Time Relative to SVD
% Use only methods 1 (SVD), 3 (GAALOP) and 4 (GAALOP CM CSE).
percentDiff_SVD = zeros(numSizes, 1);  % Baseline: SVD always 0%
percentDiff_GAALOP = 100 * (times(:,3) ./ times(:,1) - 1);
percentDiff_GAALOPCSE = 100 * (times(:,4) ./ times(:,1) - 1);

mean_PD_GAALOP = mean(percentDiff_GAALOP);
mean_PD_GAALOPCSE = mean(percentDiff_GAALOPCSE);

hPercTime = figure('Name', 'Percent_Difference_Time', 'Units','Normalized', 'OuterPosition', [0 0 1 1]);
set(0, 'DefaultAxesFontSize', 13);
set(0, 'DefaultTextFontSize', 12);
hold on;
plot(pointCloudSizes, percentDiff_SVD, pdColors{1}, 'LineWidth', 1.5, 'MarkerSize', 8);  % SVD baseline
plot(pointCloudSizes, percentDiff_GAALOP, pdColors{2}, 'LineWidth', 1.5, 'MarkerSize', 8);
plot(pointCloudSizes, percentDiff_GAALOPCSE, pdColors{3}, 'LineWidth', 1.5, 'MarkerSize', 8);
hold off;
xlabel('Number of Points', 'FontSize', 12);
ylabel('Percent Difference Relative to SVD (%)', 'FontSize', 12);
title('Percent Difference in Execution Time Compared to SVD', 'FontSize', 14);
grid on;
legendLabel1 = sprintf('SVD (0%%)');
legendLabel2 = sprintf('GAALOP (%.2f%%)', mean_PD_GAALOP);
legendLabel3 = sprintf('GAALOP CM CSE (%.2f%%)', mean_PD_GAALOPCSE);
legend(legendLabel1, legendLabel2, legendLabel3, 'Location', 'best');
print(hPercTime, fullfile(outputFolder, 'Percent_Difference_Time'), '-dpng', '-r300');

%% Plot 4: Percent Difference in Rotation Error Relative to SVD
% For rotation errors: use only methods 1, 3, and 4.
percentDiff_rot_SVD = zeros(numSizes, 1);  % SVD baseline.
percentDiff_rot_GAALOP = 100 * (errors(:,3) ./ errors(:,1) - 1);
percentDiff_rot_GAALOPCSE = 100 * (errors(:,4) ./ errors(:,1) - 1);

mean_PD_rot_GAALOP = mean(percentDiff_rot_GAALOP);
mean_PD_rot_GAALOPCSE = mean(percentDiff_rot_GAALOPCSE);

hPercRot = figure('Name', 'Percent_Difference_Rotation_Error', 'Units','Normalized', 'OuterPosition', [0 0 1 1]);
hold on;
plot(pointCloudSizes, percentDiff_rot_SVD, pdColors{1}, 'LineWidth', 1.5, 'MarkerSize', 8);
plot(pointCloudSizes, percentDiff_rot_GAALOP, pdColors{2}, 'LineWidth', 1.5, 'MarkerSize', 8);
plot(pointCloudSizes, percentDiff_rot_GAALOPCSE, pdColors{3}, 'LineWidth', 1.5, 'MarkerSize', 8);
hold off;
xlabel('Number of Points', 'FontSize', 12);
ylabel('Percent Difference Relative to SVD (%)', 'FontSize', 12);
title('Percent Difference in Rotation Error Compared to SVD', 'FontSize', 14);
grid on;
legendLabel1 = sprintf('SVD (0%%)');
legendLabel2 = sprintf('GAALOP (%.2f%%)', mean_PD_rot_GAALOP);
legendLabel3 = sprintf('GAALOP CM CSE (%.2f%%)', mean_PD_rot_GAALOPCSE);
legend(legendLabel1, legendLabel2, legendLabel3, 'Location', 'best');
print(hPercRot, fullfile(outputFolder, 'Percent_Difference_Rotation_Error'), '-dpng', '-r300');

%% ----------------- HELPER FUNCTIONS -----------------

function ptCloud = generatePointCloud(N, radius)
    % generatePointCloud
    %   Generates a synthetic 3D point cloud of N points that is guaranteed to 
    %   be well conditioned for registration (i.e. its covariance is isotropic).
    %   i.e. create synthetic point clouds that are guaranteed to have 
    %   well-conditioned sample covariance matrices
    %
    %   The function performs the following steps:
    %     1. Generates N points in a cube (in this case, uniformly in [-1,1]^3).
    %     2. Centers the points (mean zero).
    %     3. Whitens the data so that the covariance becomes exactly the identity.
    %     4. Scales the points so that they fill a sphere of the desired radius.
    %     5. Returns a MATLAB pointCloud object.
    %
    %   Inputs:
    %     N      - Number of points (should be >= 3 to get full rank).
    %     radius - The desired spatial extent (default is 1 if omitted).
    %
    %   Output:
    %     ptCloud - A MATLAB pointCloud object with isotropic covariance.
    
    if nargin < 2 || isempty(radius)
        radius = 1;
    end
    
    %--- Step 1: Initial Point Generation ---%
    % Generate N points uniformly in the cube [-1,1]^3.
    pts = 2 * rand(N, 3) - 1;
    
    %--- Step 2: Center the Points ---%
    pts = pts - mean(pts, 1);
    
    %--- Step 3: Whitening Transformation ---%
    % Compute the sample covariance matrix.
    C = cov(pts);
    
    % Do an eigen-decomposition of C.
    [V, D] = eig(C);
    
    % Construct the inverse square-root of D.
    % (For each eigenvalue, take 1/sqrt(lambda).)
    D_inv_sqrt = diag(1 ./ sqrt(diag(D)));
    
    % Apply the whitening transform: X_white = X * V * D_inv_sqrt * V'
    % This results in a point set with covariance equal to the identity.
    pts_white = pts * V * D_inv_sqrt * V';
    
    %--- Step 4: Scaling to the Desired Extent ---%
    % To make the points fill a sphere of radius 'radius', we can scale
    % the cloud so that the maximum Euclidean norm among the points equals 'radius'.
    norms = sqrt(sum(pts_white.^2, 2));
    max_norm = max(norms);
    pts_scaled = pts_white * (radius / max_norm);
    
    % Optionally, re-center the points (they should be nearly centered already).
    pts_final = pts_scaled - mean(pts_scaled, 1);
    
    %--- Step 5: Create the pointCloud Object ---%
    ptCloud = pointCloud(pts_final);
    
    % (Optional) For diagnostics, one can compute:
    % C_final = cov(ptCloud.Location);  % should be nearly (radius/max_norm)^2 * I
    % fprintf('Final covariance eigenvalues: %e, %e, %e\n', sort(eig(C_final)));
end


function [srcData, tgtData] = preprocessData(source, target)
    % preprocessData
    %   Converts the input 3D point sets to 3×N format and centers them.
    if size(source,1) ~= 3
        source = source';
    end
    if size(target,1) ~= 3
        target = target';
    end
    srcData = single(source - mean(source,2));
    tgtData = single(target - mean(target,2));
end

function [movedData, Rotation_GT, Translation_GT] = transformation(data, angleX, angleY, angleZ, translation)
    % transformation
    %   Applies a rotation (angles in degrees) and translation to input 3D data.
    %   If data is in N×3 format, it is transposed to 3×N for processing.
    if size(data, 1) ~= 3 && size(data, 2) == 3
        data = data';
    elseif size(data, 1) ~= 3
        error('Data must be a 3D point set (3×N or N×3).');
    end
    % Create rotation matrices.
    Rx = [1, 0, 0; 0, cosd(angleX), -sind(angleX); 0, sind(angleX), cosd(angleX)];
    Ry = [cosd(angleY), 0, sind(angleY); 0, 1, 0; -sind(angleY), 0, cosd(angleY)];
    Rz = [cosd(angleZ), -sind(angleZ), 0; sind(angleZ), cosd(angleZ), 0; 0, 0, 1];
    % Combine rotations (ZYX order).
    Rotation_GT = Rz * Ry * Rx;
    % Apply rotation and add translation.
    movedData = (Rotation_GT * data)' + translation;
    Translation_GT = translation;
end


function [finalRot, elapsed, errVal] = runICP(sourceData, targetData, noIteration, Rotation_GT, Translation_GT, method)
    % runICP
    %   Runs ICP registration for a single method, measures execution time and rotation error.
    tStart = tic;
    [finalRot, ~, ~, ~, ~, ~, ~, ~, ~] = ...
        icp_registration(sourceData, targetData, noIteration, Rotation_GT, Translation_GT, method, 0);
    elapsed = toc(tStart);
    errVal = norm(eye(3) - Rotation_GT * finalRot', 'fro');
end


function [finalRotation, finalTranslation, iterationDistanceError, transformationHistory, rotationHistory, ...
    angularConvergenceError, normalizedSSED, rmseHistory, translationConvergenceError] = ...
    icp_registration(sourcePoints, targetPoints, numIterations, groundTruthRotation, groundTruthTranslation, methodChoice, verbose)
    % icp_registration
    %   Performs ICP registration between the 3×N source and target point sets for a given
    %   number of iterations using the selected method.
    %
    %   methodChoice:
    %       1: SVD
    %       2: CM
    %       3: GAALOP CM
    %       4: GAALOP CM CSE
    %
    %   Returns the final rotation and translation, as well as error histories.
    angularConvergenceError = zeros(numIterations, 1);
    translationConvergenceError = zeros(numIterations, 1);
    normalizedSSED = zeros(numIterations, 1);
    rmseHistory = zeros(numIterations, 1);
    iterationDistanceError = zeros(1, numIterations);
    transformationHistory = cell(numIterations, 2);
    rotationHistory = cell(numIterations, 1);
    
    Transform_final = eye(4);
    
    for iter = 1:numIterations
        % Find the nearest neighbors between target and source.
        [minIndices, minDists, ~] = match_kNN_search(sourcePoints, targetPoints);
        correspondingPoints = sourcePoints(:, minIndices);
        iterationDistanceError(iter) = sqrt(mean(minDists));
        
        % Compute centroids.
        centerTarget = mean(targetPoints, 2);
        centerCorresponding = mean(correspondingPoints, 2);
        
        % Compute cross-covariance matrix.
        H = bsxfun(@minus, correspondingPoints, centerCorresponding) * bsxfun(@minus, targetPoints, centerTarget)';
        
        % Process H scaling
        s = evalc('disp(H)');
        sTokens = regexp(s, '[\de+-\.]+', 'match');
        if length(sTokens) == 10
            s2 = str2double(sTokens(1));
            H = H ./ s2;
        else
            if floor(H(1,1)) == 2
                H = H ./ 10;
            elseif floor(H(1,1)) > 2
                H = H ./ 100;
            end
        end
        
        % Extract elements of H
        [f11, f12, f13, f21, f22, f23, f31, f32, f33] = extract_fij(double(H));
        
        % Select update method.
        switch methodChoice
            case 1  % SVD-based ICP
                [U, ~, V] = svd(H);
                R_iter = U * V';
                if det(R_iter) < 0
                    R_iter = U * diag([1, 1, det(U*V')]) * V';
                end
                T_iter = centerCorresponding - R_iter * centerTarget;
            case 2  % CM function
                Rotor3D = char_multivector_3D(f11, f12, f13, f21, f22, f23, f31, f32, f33);
                R_iter = single(rotorToRotationMatrix(Rotor3D));
                T_iter = centerCorresponding - R_iter * centerTarget;
            case 3  % GAALOP CM
                [R_iter, ~] = char_multivector_3D_GAALOP(f11, f12, f13, f21, f22, f23, f31, f32, f33);
                R_iter = single(R_iter);
                T_iter = centerCorresponding - R_iter * centerTarget;
            case 4  % GAALOP CM CSE
                [R_iter, ~] = char_multivector_3D_GAALOP_CSE(f11, f12, f13, f21, f22, f23, f31, f32, f33);
                R_iter = single(R_iter);
                T_iter = centerCorresponding - R_iter * centerTarget;
            otherwise
                error('Method choice not implemented.');
        end
        
        % Update the overall transformation.
        T_t = [R_iter, T_iter; 0 0 0 1];
        Transform_final = T_t * Transform_final;
        
        if verbose == 1
            fprintf('Iteration %d of %d\n', iter, numIterations);
            disp('Current Inverse Transform:');
            disp(inv(Transform_final));
        end
        
        % Update target points.
        targetPoints = R_iter * targetPoints + T_iter * ones(1, size(targetPoints, 2));
        
        % Compute convergence errors.
        T_inv = inv(Transform_final);
        currentRotation = T_inv(1:3, 1:3);
        currentTranslation = T_inv(1:3, 4);
        angularConvergenceError(iter) = norm(eye(3) - groundTruthRotation * currentRotation, 'fro');
        translationConvergenceError(iter) = norm(currentTranslation - groundTruthTranslation);
        normalizedSSED(iter) = sum(sum((correspondingPoints - (R_iter * targetPoints + T_iter * ones(1, size(targetPoints, 2)))).^2));
        
        transformationHistory{iter, 1} = currentRotation;
        transformationHistory{iter, 2} = currentTranslation;
        rotationHistory{iter} = R_iter;
    end
    
    normalizedSSED = log10(normalizedSSED / sum(normalizedSSED));
    
    T_inv = inv(Transform_final);
    finalRotation = T_inv(1:3, 1:3);
    finalTranslation = T_inv(1:3, 4);
end

function Rotor = char_multivector_3D(f11, f12, f13, f21, f22, f23, f31, f32, f33)
    % char_multivector_3D
    %   Computes a rotor from the given 3×3 matrix elements.
    f1 = f11 * e1 + f12 * e2 + f13 * e3;
    f2 = f21 * e1 + f22 * e2 + f23 * e3;
    f3 = f31 * e1 + f32 * e2 + f33 * e3;
    Inv1 = e1 * f1 + e2 * f2 + e3 * f3;
    Inv2 = wedge(e2, e1) * wedge(f1, f2) + wedge(e3, e1) * wedge(f1, f3) + wedge(e3, e2) * wedge(f2, f3);
    Inv3 = wedge(e3, e2, e1) * wedge(f1, f2, f3);
    Inv_all = Inv1 + Inv2 + Inv3;
    Rotor_reverse = 1 + Inv_all;
    Rotor_reverse_normalized = unit(Rotor_reverse);
    Rotor = reverse(Rotor_reverse_normalized);
end

function [RotationMatrix] = rotorToRotationMatrix(Rotor_mv)
    % rotorToRotationMatrix
    %   Recovers a 3×3 rotation matrix from the provided rotor via a sandwich product.
    Reverse_Rotor_mv = reverse(Rotor_mv);
    f1_F = Rotor_mv * e1 * Reverse_Rotor_mv;
    f2_F = Rotor_mv * e2 * Reverse_Rotor_mv;
    f3_F = Rotor_mv * e3 * Reverse_Rotor_mv;
    F1 = padCoeffs(cell2mat(coefficients(f1_F)));
    F2 = padCoeffs(cell2mat(coefficients(f2_F)));
    F3 = padCoeffs(cell2mat(coefficients(f3_F)));
    RotationMatrix = reshape([F1; F2; F3], 3, 3)';
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


function [RotationMatrix, R] = char_multivector_3D_GAALOP(f11, f12, f13, f21, f22, f23, f31, f32, f33)
    % char_multivector_3D_GAALOP
    %   Computes the rotation using the GAALOP multivector formulation.
    R(1) = ((f11 + 1.0) * f22 - f12 * f21 + f11 + 1.0) * f33 + (((-f11) - 1.0) * f23 + f13 * f21) * f32 + (f12 * f23 - f13 * f22 - f13) * f31 + (f11 + 1.0) * f22 - f12 * f21 + f11 + 1.0;
    R(5) = (f21 - f12) * f33 + f13 * f32 - f23 * f31 + f21 - f12;
    R(6) = (-(f21 * f32)) + (f22 + 1.0) * f31 + f12 * f23 - f13 * f22 - f13;
    R(7) = (f11 + 1.0) * f32 - f12 * f31 + ((-f11) - 1.0) * f23 + f13 * f21;
    f1_F(2) = (R(7)^2 - R(6)^2 - R(5)^2 + R(1)^2) / (R(7)^2 + R(6)^2 + R(5)^2 + R(1)^2);
    f1_F(3) = (-(2.0 * R(6) * R(7) + 2.0 * R(1) * R(5))) / (R(7)^2 + R(6)^2 + R(5)^2 + R(1)^2);
    f1_F(4) = (2.0 * R(5) * R(7) - 2.0 * R(1) * R(6)) / (R(7)^2 + R(6)^2 + R(5)^2 + R(1)^2);
    f2_F(2) = (-(2.0 * R(6) * R(7) - 2.0 * R(1) * R(5))) / (R(7)^2 + R(6)^2 + R(5)^2 + R(1)^2);
    f2_F(3) = (-(R(7)^2 - R(6)^2 + R(5)^2 - R(1)^2)) / (R(7)^2 + R(6)^2 + R(5)^2 + R(1)^2);
    f2_F(4) = (-(2.0 * R(1) * R(7) + 2.0 * R(5) * R(6))) / (R(7)^2 + R(6)^2 + R(5)^2 + R(1)^2);
    f3_F(2) = (2.0 * R(5) * R(7) + 2.0 * R(1) * R(6)) / (R(7)^2 + R(6)^2 + R(5)^2 + R(1)^2);
    f3_F(3) = (2.0 * R(1) * R(7) - 2.0 * R(5) * R(6)) / (R(7)^2 + R(6)^2 + R(5)^2 + R(1)^2);
    f3_F(4) = (-(R(7)^2 + R(6)^2 - R(5)^2 - R(1)^2)) / (R(7)^2 + R(6)^2 + R(5)^2 + R(1)^2);
    RotationMatrix = reshape([f1_F(2:4), f2_F(2:4), f3_F(2:4)], 3, 3)';
end

function [RotationMatrix, R] = char_multivector_3D_GAALOP_CSE(f11, f12, f13, f21, f22, f23, f31, f32, f33)
    % char_multivector_3D_GAALOP_CSE
    %   Computes the rotation matrix using GAALOP with Common Subexpression Elimination (CSE).
    temp_gcse_4 = f12 * f23 - f13 * f22;
    temp_gcse_8 = (-f11) - 1.0;
    temp_gcse_15 = (f11 + 1.0) * f22 - f12 * f21;
    temp_gcse_17 = temp_gcse_4 - f13;
    temp_gcse_19 = f13 * f21;
    temp_gcse_26 = f11 + 1.0;
    temp_gcse_27 = temp_gcse_8 * f23;
    R(1) = (temp_gcse_15 + f11 + 1.0) * f33 + (temp_gcse_27 + temp_gcse_19) * f32 + temp_gcse_17 * f31 + temp_gcse_15 + f11 + 1.0;
    temp_gcse_24 = f21 - f12;
    R(5) = temp_gcse_24 * f33 + f13 * f32 - f23 * f31 + temp_gcse_24;
    R(6) = (-(f21 * f32)) + (f22 + 1.0) * f31 + temp_gcse_17;
    R(7) = temp_gcse_26 * f32 - f12 * f31 + temp_gcse_27 + temp_gcse_19;
    temp = (R(7)^2 + R(6)^2 + R(5)^2 + R(1)^2);
    f1_F(2) = (R(7)^2 - R(6)^2 - R(5)^2 + R(1)^2) / temp;
    f1_F(3) = (-(2.0 * R(6) * R(7) + 2.0 * R(1) * R(5))) / temp;
    f1_F(4) = (2.0 * R(5) * R(7) - 2.0 * R(1) * R(6)) / temp;
    f2_F(2) = (-(2.0 * R(6) * R(7) - 2.0 * R(1) * R(5))) / temp;
    f2_F(3) = (-(R(7)^2 - R(6)^2 + R(5)^2 - R(1)^2)) / temp;
    f2_F(4) = (-(2.0 * R(1) * R(7) + 2.0 * R(5) * R(6))) / temp;
    f3_F(2) = (2.0 * R(5) * R(7) + 2.0 * R(1) * R(6)) / temp;
    f3_F(3) = (2.0 * R(1) * R(7) - 2.0 * R(5) * R(6)) / temp;
    f3_F(4) = (-(R(7)^2 + R(6)^2 - R(5)^2 - R(1)^2)) / temp;
    RotationMatrix = reshape([f1_F(2:4), f2_F(2:4), f3_F(2:4)], 3, 3)';
end

function [f11, f12, f13, f21, f22, f23, f31, f32, f33] = extract_fij(M)
    % Extracts each element of the 3×3 matrix M into separate outputs f11…f33
    f11 = M(1,1);
    f12 = M(1,2);
    f13 = M(1,3);
    f21 = M(2,1);
    f22 = M(2,2);
    f23 = M(2,3);
    f31 = M(3,1);
    f32 = M(3,2);
    f33 = M(3,3);
end

function [minIndices, minDists, allDists] = match_kNN_search(source, target)
    % For each point in target, finds its nearest neighbor in source
    sourceT = source';
    targetT = target';
    D = pdist2(targetT, sourceT);
    [minDists, idx] = min(D, [], 2);
    minIndices = idx;
    allDists = D;
end

