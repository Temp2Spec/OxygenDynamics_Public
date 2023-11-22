%% Script for finding relationships between vasculature and hypoxic pockets in the mouse barrel cortex


% Uses the output of OxygenDynamics_Master analysis script and annotated
% png images for arteries and veins 
% Run the script on the masterfolder that contains all the recordings. 
% In the masterfolder have a csv file with the paths of the recordings'
% folders. 
% In each recording folder you need:
% (1) A subfolder named ManualCurOxySinkData that contains a .mat file from
% the OxygenDynamics_Master analysis script
% (2) A png file with the manually annotated Veins that has 'Veins' on its
% name
% (3) A png file with the manually annotated Arteries that has 'Arteriess' on its
% name

% The output of this script are four coloumn vectors that are added into the ManualCuration .mat file. 
% These vectors contain the distance from the closest vein and artery caclulated from the centroid and any pixel belonging to each sink

% Antonis Asiminas, PhD
% Center for Translational Neuromedicine
% University of Copenhagen, November 2021 - January 2023

%%

InputD=readtable('DataPathsExample.csv','Delimiter',',');

%For some reason sometimes the csv is not being read correctly and it InputD ends up being
%a single columns of all info This looks like it is an Excel thing..
if size(InputD,2)<2
    InputD=readtable('DataPathsExample.csv');
end

Paths=table2cell(InputD(:,{'Paths'}));

%% Get the data first
Masterfolder=pwd;

for datai=1:length(Paths)

cd(Paths{datai});
% Open the png images containing the info for artery and vein morphology
% and from rgb convert to grayscale and subsequently to logical

Arteries_path=findFilesWithWord('./', 'Arteries');
Veins_path=findFilesWithWord('./', 'Veins');
Arteries = logical(rgb2gray(imread(Arteries_path{1}))); % './' means the current folder. I can later modify this. 
Veins = logical(rgb2gray(imread(Veins_path{1})));
%% Open the folder with the oxysinks output and load the mat file
OxySinks_dir=findFilesWithWord('./','ManualCurOxySinksData');
cd(OxySinks_dir{1})

OxySinksdata_dir=findFilesWithWord('./','ManualCuration');

OxySinksdata = load(OxySinksdata_dir{1});

Table_OxygenSinks_Out=OxySinksdata.Table_OxygenSinks_Out;

% %% Filterring pockets 
% %If there are curated data use that
% if isfield(OxySinksdata, 'Included_Pocs') %if I'm going to use currated data I need to get the logical vector to filter the oxygen sink data
%                   
%     Table_OxygenSinks_Out=OxySinksdata.Table_OxygenSinks_Out(Included_Pocs,:);
%        
% else
%             
%     %a disclaimer if curated data were chosen to be analysed but the            
%     %filter from the curator is not present
%             
%     fprintf(['No Oxygen sink currated data were fround in data folder', Paths{datai},'\n'])
%             
%             
%     %then if the .mat file contains a vector indicating potential            
%     %noise sinks (if the master script is from after            
%     %February 2022)
%                         
%     if isfield(OxySinksdata, 'Trace_PotentialNoise') 
%                 
% 
%         %filter the data based on this                
%         Table_OxygenSinks_Out=OxySinksdata.Table_OxygenSinks_Out(OxySinksdata.Trace_PotentialNoise,:);
%       
%     end
%                    
% end

% Getting the map of oxygen recordings
Map=OxySinksdata.OxySink_Map;

%calculate the scalling factor 
rescaleF_Veins=size(Veins,1)/size(Map,1);
rescaleF_Arter=size(Arteries,1)/size(Map,1);
%I assume that ceil(size(Veins,1)/rescaleF) == size(Map,1). It is possible
%that for some cases this will not be the case. I need to include some
%checks here for the script to run smoothly. 

%%
% Finding the all the "islands" of non zero that correspond to the vasculature
Vein_regions=regionprops(Veins,'PixelIdxList');
Arte_regions=regionprops(Arteries,'PixelIdxList');

% Getting the indices of the pockets' centroids
Centroids_Idx=[round(Table_OxygenSinks_Out.MeanCentroid_y), round(Table_OxygenSinks_Out.MeanCentroid_x)];



Veins_Distance_all=nan(height(Table_OxygenSinks_Out),1);
Veins_Distance_centroid=nan(height(Table_OxygenSinks_Out),1);
Arteries_Distance_all=nan(height(Table_OxygenSinks_Out),1);
Arteries_Distance_centroid=nan(height(Table_OxygenSinks_Out),1);
% Go through the different pockets
for i=1:height(Table_OxygenSinks_Out)


    Pixels=Table_OxygenSinks_Out.OxySink_Pxls_all{i};

    [row1, col1] = ind2sub(size(Map), Pixels);
    % the coordinates of all pixels for this pocket
    coords1=[col1, row1]; 
    % the coordinates of mean centroid for this pocket
    coords_c=[round(Table_OxygenSinks_Out.MeanCentroid_y(i)),round(Table_OxygenSinks_Out.MeanCentroid_x(i))];

    mindistveins_all=nan(length(Vein_regions),1);
    mindistveins_centroid=nan(length(Vein_regions),1);
    for j=1:length(Vein_regions) %now go through all the identified vein areas

        [row2, col2] = ind2sub(size(Veins), Vein_regions(j).PixelIdxList);


         coords2=[col2, row2]; % the coordinates of all pixels for this vein area

         %adjust coordinates to match the resolution of the oxygen imaging
         coords2=ceil(coords2/rescaleF_Veins);

         % Calculate the pairwise Euclidean distances between the all the coordinates
         distances_pixels_all = pdist2(coords1, coords2);
         % repeat for centroid
         distances_pixels_c = pdist2(coords_c, coords2);


         mindistveins_all(j)=min(distances_pixels_all(:));
         mindistveins_centroid(j)=min(distances_pixels_c(:));

    end

    Veins_Distance_all(i)=min(mindistveins_all);
    Veins_Distance_centroid(i)=min(mindistveins_centroid);


    mindistarter_all=nan(length(Arte_regions),1);
    mindistarter_centroid=nan(length(Arte_regions),1);
    for j=1:length(Arte_regions) %now go through all the identified vein areas

        [row2, col2] = ind2sub(size(Veins), Arte_regions(j).PixelIdxList);


         coords2=[col2, row2]; % the coordinates of all pixels for this vein area

         %adjust coordinates to match the resolution of the oxygen imaging
         coords2=ceil(coords2/rescaleF_Arter);

         % Calculate the pairwise Euclidean distances between the all the coordinates
         distances_pixels_all = pdist2(coords1, coords2);
         % repeat for centroid
         distances_pixels_c = pdist2(coords_c, coords2);


         mindistarter_all(j)=min(distances_pixels_all(:));
         mindistarter_centroid(j)=min(distances_pixels_c(:));

    end

    Arteries_Distance_all(i)=min(mindistarter_all);
    Arteries_Distance_centroid(i)=min(mindistarter_centroid);

end


clear row1 col1 row2 col2 i j Pixels coords1 distances_pixels_all distances_pixels_c mindistarter_all mindistarter_centroid mindistveins_all mindistveins_centroid


%% Amend the .mat file

% Open the .mat file in read-write mode
matObj = matfile(OxySinksdata_dir{1}, 'Writable', true);

% Add or modify variables within the .mat file
matObj.Veins_Distance_all = Veins_Distance_all;  % Add a new variables
matObj.Veins_Distance_centroid = Veins_Distance_centroid;  % Add a new variables
matObj.Arteries_Distance_all = Arteries_Distance_all;  % Add a new variables
matObj.Arteries_Distance_centroid = Arteries_Distance_centroid;  % Add a new variables

%Also export as a csv
exporttable=table(Arteries_Distance_centroid,Arteries_Distance_all,Veins_Distance_centroid,Veins_Distance_all);

Output_xlsx=sprintf('Vascular_Analysis.xlsx');
writetable(exporttable,Output_xlsx);

%%
cd(Masterfolder) %return to the masterfolder before going to the next data folder

end

clear all

%% FUNCTIONS

%% Find files or folder that contain a word of interest on their name
function matchingFilePaths = findFilesWithWord(folderPath, searchWord)
    % Get the list of all files and folders in the specified directory
    fileList = dir(folderPath);

    % Getting the names of the files/folders that contain the word of interest
    matchingFilePaths={fileList(cellfun(@(x) contains(x, searchWord),{fileList.name})).name};
    
    % Adding the name of the path to create full filepaths instead of file
    % names
    matchingFilePaths=cellfun(@(x) fullfile(fileList(1).folder,x),matchingFilePaths,'UniformOutput',false);
    
end

