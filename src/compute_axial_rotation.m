clear
clc
close all

%%
% This example script computes axial rotations of key vertebrae (T8, T11,
% L2) from the landmarks available in the obline repository.


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

% Most common apical vertebrae according to:
%    Acaroglu et al., 2001. Does Transverse Apex Coincide With Coronal Apex 
%    Levels (Regional or Global) in Adolescent Idiopathic Scoliosis?
%    Spine 26(10): 1143–1146.
%    https://doi.org/10.1097/00007632-200105150-00010
apex_vertebrae = {'T8', 'T11', 'L2'};


% Preallocate
Rotations    = zeros(Nsubjects, 3);


for subjectN = 1 : Nsubjects
    disp(['Processing subject ' num2str(subjectN)])

    for nvert = 1 : length(apex_vertebrae)
        % Landmarks
        vert_name = ['Vertebra_' apex_vertebrae{nvert}];
    
        M = vertebra_frame(subjectN, vert_name, data);
        [Lat,Sag,Ax] = orientation_frame(M);

        Rotations(subjectN,nvert) = Ax;

    end
end

%%
Rotations = cell2table(num2cell(Rotations));
Rotations.Properties.VariableNames = apex_vertebrae;


save('Rotations.mat', 'Rotations');



function M = vertebra_frame(subjectN, vert_name, data)
% Computes a frame of reference for a vertebra

    % Helper function to retrieve coordinates of a specific landmark
    retrieveCoords = @(N, landmark) [data.([landmark, '_x'])(N), data.([landmark, '_y'])(N), data.([landmark, '_z'])(N)];
    
    UpperEndplate = zeros(12,3);
    LowerEndplate = zeros(12,3);
    % Find center of vertebral body
    for n_node = 1 : 12
        UpperEndplate(n_node, :) = retrieveCoords(subjectN, ['Spine_' vert_name '_endplate_sup_' num2str(n_node)]);
        LowerEndplate(n_node, :) = retrieveCoords(subjectN, ['Spine_' vert_name '_endplate_inf_' num2str(n_node)]);
    end

    % Center of the vertebral body is the middle point between
    % endplates
    vertebral_body = mean([mean(UpperEndplate); mean(LowerEndplate)]);

    % Midpoint between pedicles
    pedicles = mean([retrieveCoords(subjectN, ['Spine_' vert_name '_pedicle_right']); ...
                     retrieveCoords(subjectN, ['Spine_' vert_name '_pedicle_left'])]);

    % Anteroposterior vector
    X = vertebral_body - pedicles;
    % Vertical vector
    Z = mean(UpperEndplate) - mean(LowerEndplate);
    % Left vector
    Y = cross(Z/norm(Z), X/norm(X));

    % Recalculate X
    X = cross(Y, Z/norm(Z));

    % Homogeneous coordinate system
    M = [X'/norm(X), Y'/norm(Y), Z'/norm(Z), vertebral_body'; 0 0 0 1];


end



function [Lateral,Sagittal,Axial] = orientation_frame(M)
    % Calculates the orientation of a frame of reference in the Lateral
    % plane, Sagittal and Axial plane
    
    % Invert M only once
    M = inv(M);
    
    Sagittal = asin(M(3,1));
    Lateral = real(asin(-M(3,2)/cos(Sagittal)));
    Axial = real(asin(-M(2,1)/cos(Sagittal)));
    
    % rad to degrees
    Sagittal = Sagittal / pi * 180;
    Lateral = Lateral / pi * 180;
    Axial = Axial / pi * 180;

end  