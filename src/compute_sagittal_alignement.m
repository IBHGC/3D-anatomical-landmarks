clear
clc
close all

%%
% This scripts shows how to calculate sagittal alignement parameters from
% the 3D spinopelvic landmarks.
% For robustness, a plane approximating the endplates of the two vertebrae
% of interest is calculated. Then, the angle between the normals to the
% two endplates is calculated and projected on the plane of interest.

%%
% This is the data from the online repository "3D anatomic landmarks" v1.0
% DOI: https://doi.org/10.57745/TJJLD6
data_url = 'https://entrepot.recherche.data.gouv.fr/api/access/datafile/723128?format=original&gbrecs=true';

websave('repository_data.csv', data_url);

% Load data in script folder
data = readtable('repository_data.csv');
Nsubjects = size(data, 1);

% Helper function to retrieve coordinates of a specific landmark
retrieveCoords = @(N, landmark) [data.([landmark, '_x'])(N), data.([landmark, '_y'])(N), data.([landmark, '_z'])(N)];

%%
% Sagittal parameters that will be calculated
dataLandmarks.T1T12Kyphosis = zeros(Nsubjects, 1);
dataLandmarks.T4T12Kyphosis = zeros(Nsubjects, 1);
dataLandmarks.L1L5Lordosis  = zeros(Nsubjects, 1);
dataLandmarks.pelvicTilt    = zeros(Nsubjects, 1);

headers = data.Properties.VariableNames;


for subjectN = 1 : Nsubjects
    disp(subjectN)

    %% T1-T12 kyphosis
    T1UpperEndplate = extract_endplate_vector(subjectN, 'T1', 'sup', data);
    
    T12LowerEndplate = extract_endplate_vector(subjectN, 'T12', 'inf', data);

    % Angle between endplates
    dataLandmarks.T1T12Kyphosis(subjectN) = sagittalAngleY(T12LowerEndplate, T1UpperEndplate);


    %% T4-T12
    T4UpperEndplate = extract_endplate_vector(subjectN, 'T4', 'sup', data);

    dataLandmarks.T4T12Kyphosis(subjectN) = sagittalAngleY(T12LowerEndplate, T4UpperEndplate);


    %% L1-L5 lordosis
    L1UpperEndplate = extract_endplate_vector(subjectN, 'L1', 'sup', data);
    L5LowerEndplate = extract_endplate_vector(subjectN, 'L5', 'inf', data);

    
    % Angle between endplates
    dataLandmarks.L1L5Lordosis(subjectN) = sagittalAngleY(L5LowerEndplate, L1UpperEndplate);


    %% Pelvic tilt
    % Center of sacrum
    Sacrum = retrieveCoords(subjectN, 'Pelvis_sacrum_sacrum_center');
    
    % Angle with the vertical
    dataLandmarks.pelvicTilt(subjectN) = sagittalAngleY(Sacrum, [0,0,1]);



end

%%
save([fileparts(mfilename("fullpath")), '\dataLandmarks.mat'], 'dataLandmarks');
disp(['File saved: ' fileparts(mfilename("fullpath")), '\dataLandmarks.mat'])





function EndplateVector = extract_endplate_vector(subjectN, Vertebra, endplate, data)
%   EndplateVector = extract_endplate_vector(subjectN, Vertebra, endplate, data)
%   Retrieves the coordinates of the endplate contour for a given subject and
%   vertebra (and upper or lower plateau), calculates an approximating plane
%   of the endplane and returns a normal vector to the plane.
%
%   Inputs:
%     - subjectN  : Index of the subject in the dataset.
%     - Vertebra  : String specifying the vertebra (e.g., 'T1', 'L5').
%     - endplate  : String specifying the endplate ('sup' for superior, 'inf' for inferior).
%     - data      : Struct containing coordinate fields in the format:
%                  'Spine_Vertebra_<Vertebra>_endplate_<endplate>_<point>_x/y/z'
%
%   Output:
%     - EndplateVector : A 1x3 matrix containing the 3D coordinates of endplate normal vector.
%


    % A helper function to extract x, y and z coordinates
    retrieveCoords = @(N, landmark) [data.([landmark, '_x'])(N), data.([landmark, '_y'])(N), data.([landmark, '_z'])(N)];
    
    % Generate the tags for the 12 points of the thoracolumbar endplate
    tags = arrayfun(@(x) sprintf(['Spine_Vertebra_' Vertebra '_endplate_' endplate '_%d'], x), 1:12, 'UniformOutput', false);

    % Retrieve the coordinates
    EndplateVector = zeros(length(tags),3);
    for k = 1 : length(tags)
        EndplateVector(k,:) = retrieveCoords(subjectN, tags{k});
    end

    % On superior endplate, the 9th point is on the right side. On the lower
    % endplate it's the 5th points
    if strcmpi(endplate, 'sup')
        n = '9';
    else
        n = '5';
    end

    % Calculate normal vector
    EndplateVector = planeNormal(EndplateVector, retrieveCoords(subjectN, ['Spine_Vertebra_' Vertebra '_endplate_' endplate '_' n]));

end



function theta = sagittalAngleY(v1, v2)
% sagittalAngleY computes the angle between vectors v1 and v2 projected 
% onto the sagittal plane (XZ-plane)
% Inputs: v1, v2 - 1x3 vectors representing 3D points
% Output: theta - angle in degrees between the projected vectors

    % This vector can be emtpy if the subject did not have a L5 vertebra
    if any(isnan([v1(:) ; v2(:)])) || any([isempty(v1), isempty(v2)])
        theta = NaN;  
        return;   
    end   

    % Project onto Y plane (i.e., zero out Y component)
    v1(2) = 0;
    v2(2) = 0;

    % Normalize vectors
    v1 = v1 / norm(v1);
    v2 = v2 / norm(v2);

     % Compute angle using atan2 for signed result
    crossProd = cross(v1, v2);
    theta = atan2(norm(crossProd), dot(v1, v2)) / pi * 180;

    % Determine sign using Y component of cross product
    if crossProd(2) < 0
        theta = -theta;
    end


end




function n = planeNormal(points, rightPoint)
% planeNormal computes the normal vector to the best-fit plane of 3D points
% Input: points - Nx3 matrix of 3D coordinates
%        rightDirection - a point right of the endplate. This allows
%        correcting the normal pointing "up"
% Output: n - 1x3 normal vector to the approximating plane

    if any(isnan(points))
        n = [];
        return
    end


    % Center the points
    centroid = mean(points, 1);
    centered = points - centroid;

    % Compute covariance matrix and perform SVD
    [~, ~, V] = svd(centered, 0);

    % The normal is the last singular vector
    n = V(:, end)';

    % Must point up, correct if necessary
    rightPoint = rightPoint - centroid;
    
    % This should point forward
    correction = cross(n, rightPoint);

    % If it does not, correct it
    if correction(1) < 0
        n = -n;
    end


end
