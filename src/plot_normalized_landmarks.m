clear
clc
close all

%%
% This example script plots all landmarks by splitting the cohort in
% three groupes according to age, and normalizing the landmaarks
% coordinates so that the vertical distance between the orginin (i.e., the 
% center of the femoral heads) and T1 is "100".

%%
% This is the data from the online repository "3D anatomic landmarks" v1.0
% DOI: https://doi.org/10.57745/TJJLD6
data_url = 'https://entrepot.recherche.data.gouv.fr/api/access/datafile/723128?format=original&gbrecs=true';

websave('repository_data.csv', data_url);

% Load data in script folder
data = readtable('repository_data.csv');
Nsubjects = size(data, 1);

% Split cohort by age
determine_group = @(age) find([age < 18, age >= 18 & age < 50, age >= 50]);

% Helper function to compute mean of 12 vertebral landmarks
mean_coords = @(prefix, direction) mean(cell2mat( ...
    arrayfun(@(k) data.(sprintf('%s_%d_%s', prefix, k, direction)), 1:12, 'uni', 0)), 2);


% Compute upper and lower endplates centers
directions = {'x', 'y', 'z'};
upper_endplates = cell2mat(cellfun(@(a) mean_coords('Spine_Vertebra_T1_endplate_sup', a), directions, 'UniformOutput', false));
lower_endplates = cell2mat(cellfun(@(a) mean_coords('Spine_Vertebra_T1_endplate_inf', a), directions, 'UniformOutput', false));

% Compute T1 body center of all subjects as the middle point between 
% upper and lower endplates
T1_body = mean(cat(3, upper_endplates, lower_endplates), 3);


% Extract coordinates from table by removing variables that do not 
% contain "_x", "_y" or "z"
keep_coordinates = (contains(data.Properties.VariableNames, '_x') | ...
                    contains(data.Properties.VariableNames, '_y') |...
                    contains(data.Properties.VariableNames, '_z'));

% tags for plotting later
tags = prepare_tags_by_landmark_type(data, keep_coordinates);

% Prepare figure
figure
subplot(1,3,1)
plot3(0,0,0, 'rx')  % Plot the origin just to intialize the plot
axis equal; hold on;
title('Age < 18')

subplot(1,3,2)
plot3(0,0,0, 'rx')
axis equal; hold on;
title('Age >= 18 & age < 50')

subplot(1,3,3)
plot3(0,0,0, 'rx')
axis equal; hold on;
title('Age >= 50')


groups = cell(1, 3);
for g = 1:3
    groups{g} = struct('Odontoid', [], 'Spine', [], 'Pelvis', [], 'LowerLimb', []);
end

% Collect all points by group
for k = 1:size(data, 1)
    disp(k)
    % Skip if age not available
    if isempty(determine_group(data.Age(k)))
        continue
    end 
    % Keep only datapoints with coordinates for the nth subject
    subject_points = table2array(data(k, keep_coordinates));
    % Rearrange to have [x, y, z]
    subject_points = reshape(subject_points, 3, length(subject_points)/3)';
    
    % Normalize all coordinates according to the vertical height of T1
    normalization = T1_body(k,3);
    subject_points = subject_points ./ normalization * 100;
    
    % Determine which group this subject belongs to
    group_idx = determine_group(data.Age(k));
    
    % Append points to the appropriate group
    groups{group_idx}.Odontoid  = [groups{group_idx}.Odontoid; subject_points(tags.Odontoid, :)];
    groups{group_idx}.Spine     = [groups{group_idx}.Spine; subject_points(tags.Spine, :)];
    groups{group_idx}.Pelvis    = [groups{group_idx}.Pelvis; subject_points(tags.Pelvis, :)];
    groups{group_idx}.LowerLimb = [groups{group_idx}.LowerLimb; subject_points(tags.LowerLimb, :)];
end



% Plot all groups; this is much faster than plotting each subject
% separately!
for g = 1:3
    subplot(1, 3, g)
    hold on
    
    plot3(groups{g}.Odontoid(:,1), groups{g}.Odontoid(:,2), groups{g}.Odontoid(:,3), 'ko', 'MarkerSize', 6)
    plot3(groups{g}.Spine(:,1), groups{g}.Spine(:,2), groups{g}.Spine(:,3), 'r.')
    plot3(groups{g}.Pelvis(:,1), groups{g}.Pelvis(:,2), groups{g}.Pelvis(:,3), 'g.')
    plot3(groups{g}.LowerLimb(:,1), groups{g}.LowerLimb(:,2), groups{g}.LowerLimb(:,3), 'b.')
end



function tags= prepare_tags_by_landmark_type(data, keep_coordinates)
    % Returns tags by landmark type (Spine, Lowerlimb, etc) to plot with  
    % different colours

    % Tags by structure
    tags_names = data.Properties.VariableNames(keep_coordinates);

    % Rearrange tags and keep only tag of each (x,y,z) triplet
    tags_names = tags_names(1:3:end);  % keep first 

    % tags
    tags.Odontoid = contains(tags_names, 'Odontoid');
    tags.Spine = contains(tags_names, 'Spine');
    tags.Pelvis = contains(tags_names, 'Pelvis');
    tags.LowerLimb = contains(tags_names, 'LowerLimb');
end

