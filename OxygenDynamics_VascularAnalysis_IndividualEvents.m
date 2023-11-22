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

InputD=readtable('VasculatureData.csv','Delimiter',',');

%For some reason sometimes the csv is not being read correctly and it InputD ends up being
%a single columns of all info This looks like it is an Excel thing..
if size(InputD,2)<2
    InputD=readtable('VasculatureData.csv');
end

Paths=table2cell(InputD(:,{'Paths'}));
Pixelsizes=table2cell(InputD(:,{'Pixelsize'}));

%%

SinkFold_OldORNew = questdlg('Do you want to use the oldest or the most recent data for oxygen sinks?', ...
	'Which folder to use?', ...
	'Oldest','Recent','Recent');

%% Get the data first
Masterfolder=pwd;

for datai=1:length(Paths)

cd(Paths{datai});
PxlSz = Pixelsizes{datai};
% Open the png images containing the info for artery and vein morphology
% and from rgb convert to grayscale and subsequently to logical
%%
Arteries_path=findFilesWithWord('./', 'Arteries');
Veins_path=findFilesWithWord('./', 'Veins');
Arteries = logical(rgb2gray(imread(Arteries_path{1}))); % './' means the current folder. I can later modify this. 
Veins = logical(rgb2gray(imread(Veins_path{1})));
%% Open the folder with the oxysinks output and load the mat and tif files

 % finding the folders that contain the sink outputs from master analysis script         
    
list=dir;  
list=list([list.isdir]);
list=list(3:end);

wrapper = @(x) contains(x,'OxygenSinks_Output'); 
SinkFolders=list(cell2mat(cellfun(wrapper, {list(:).name}, 'UniformOutput', false))); 

if isempty(SinkFolders)

   continue
end

    
if strcmp(SinkFold_OldORNew,'Oldest')
        
   % choose the older folder with that name          
   [~,I]=min([SinkFolders.datenum]); 
   SinksDataFolder= [SinkFolders(I).folder,'\',SinkFolders(I).name];
 
else  
   % or choose the most recent
   [~,I]=max([SinkFolders.datenum]);
   SinksDataFolder= [SinkFolders(I).folder,'\',SinkFolders(I).name];
end

   
cd(SinksDataFolder) %go in the selected manual curation folder
d= dir('*mat'); %get the name of the mat file in this folder  
Tifffile = dir('*.tif'); %get the name of the logical matric in this folder
%load the table with the info on the sinks
%%
wrapper = @(x) contains(x,'Urefined');
load(d(cellfun(wrapper, {d(:).name})).name,"Table_OxygenSinks_Out");
load(d.name,'OxySinkArea_all');
DIR = fullfile(Tifffile.folder,'\', Tifffile.name);
Pixelsize=cell(height(Table_OxygenSinks_Out),1);Pixelsize(:)={PxlSz};
Table_OxygenSinks_Out = addvars(Table_OxygenSinks_Out,Pixelsize,'Before','RecDuration');

[IM_BW,~,~]=loadtiff(DIR);

%converting to logical to make my life easier with size
IM_BW=logical(IM_BW);

%%

%%%%%%%%%%%%%%%%%%%%%
%%Get the hypoxic even specific metrics
OxygenSinksInfo_rec=cell(size(IM_BW,3),1);

for i=1:size(IM_BW,3)

    OxygenSinksInfo_rec{i}=regionprops(IM_BW(:,:,i),'Area','PixelIdxList');

end

sink_overlap_thres=0.6; %percentage of overlap that is needed between sinks in two different frames



Overall_OxygenSinks_Pxllist_rec=cell(size(OxySinkArea_all));
    
    for i=1:size(OxySinkArea_all,1) %for every pocket

        for ij=1:size(OxySinkArea_all,2) %%for every frame

            
           if ~isempty(OxySinkArea_all{i,ij}) %check if the pocket had an ongoing event in that frame


               % comparing with the size of the identified pockets in the
               % same frame.
               %if there is only pocket with the same size
               if sum([OxygenSinksInfo_rec{ij}(:).Area]==OxySinkArea_all{i,ij})==1
               
               
                   %just capture the pixel indices for this frame
                   Overall_OxygenSinks_Pxllist_rec{i,ij}=OxygenSinksInfo_rec{ij}([OxygenSinksInfo_rec{ij}(:).Area]==OxySinkArea_all{i,ij}).PixelIdxList;


                   % replace the area with nan so it cannot be added to
                   % another pocket
                   OxygenSinksInfo_rec{ij}([OxygenSinksInfo_rec{ij}(:).Area]==OxySinkArea_all{i,ij}).Area=nan;

    
               %if there are more pixel lists than pockets in that frame
               elseif length(OxygenSinksInfo_rec{ij})>sum(~cellfun(@isempty,OxySinkArea_all(:,ij)))
                   
                   %create all possible combos to add
                   allcombos=nchoosek(1:length(OxygenSinksInfo_rec{ij}),2);

                   allcomboAreas=[OxygenSinksInfo_rec{ij}(allcombos(:,1)).Area]+[OxygenSinksInfo_rec{ij}(allcombos(:,2)).Area];

                   if sum(allcomboAreas==OxySinkArea_all{i,ij})==1 % check if there is a combo that gives the area of Sink I'm asfter


                       constituentpockets=allcombos(allcomboAreas==OxySinkArea_all{i,ij},:);

                       for t=constituentpockets
                       
                           %capture the pixel indices from all constituents of
                           %this frame into one
                           Overall_OxygenSinks_Pxllist_rec{i,ij}=[Overall_OxygenSinks_Pxllist_rec{i,ij};[OxygenSinksInfo_rec{ij}(t).PixelIdxList]];

                       end

                       % replace the area of all constituent pockets with nans so it cannot be added to
                       % another pocket
                       for t=constituentpockets
                       
                           OxygenSinksInfo_rec{ij}(t).Area=nan;

                       end


                   end
               
               
               else %if neither a single area or a combination matched the area from Sinks_Area, find the one with the most overlap


                   tempintersect=nan(length(OxygenSinksInfo_rec{ij}),1);
                   
                   for k=1:length(OxygenSinksInfo_rec{ij}) %go through the pixles lists of the pockets of this frame


                       %calculate the intersection between the pixels of
                       %the pockets for that frame and the unique pixels
                       %for that pocket across recording.
                       tempintersect(k)=length(intersect([OxygenSinksInfo_rec{ij}(k).PixelIdxList],Table_OxygenSinks_Out.OxySink_Pxls_all{i}));


                   end



                   if sum(logical(tempintersect))==1


                        %and get the pixels form the pocket with the highest
                        %intersection             
                        Overall_OxygenSinks_Pxllist_rec{i,ij}=OxygenSinksInfo_rec{ij}(tempintersect==max(tempintersect)).PixelIdxList;

                        % replace the area with nan so it cannot be added to
                        % another pocket
                        OxygenSinksInfo_rec{ij}(tempintersect==max(tempintersect)).Area=nan;

                   end


               end
           end

        end

        
    end


%%

    Area_um=nan(sum(Table_OxygenSinks_Out.NumOxySinkEvents(:)),1);
    Area_norm=nan(sum(Table_OxygenSinks_Out.NumOxySinkEvents(:)),1);
    FilledArea_um=nan(sum(Table_OxygenSinks_Out.NumOxySinkEvents(:)),1);
    FilledArea_norm=nan(sum(Table_OxygenSinks_Out.NumOxySinkEvents(:)),1);
    Centroid_x=nan(sum(Table_OxygenSinks_Out.NumOxySinkEvents(:)),1);
    Centroid_y=nan(sum(Table_OxygenSinks_Out.NumOxySinkEvents(:)),1);
    Circularity=nan(sum(Table_OxygenSinks_Out.NumOxySinkEvents(:)),1);
    Perimeter_um=nan(sum(Table_OxygenSinks_Out.NumOxySinkEvents(:)),1);
    Diameter_um=nan(sum(Table_OxygenSinks_Out.NumOxySinkEvents(:)),1);
    OxySink_Pxls_event=cell(sum(Table_OxygenSinks_Out.NumOxySinkEvents(:)),1);

    counter=1;
    for k=1:height(Table_OxygenSinks_Out) %for every pocket of this recording

        %create a cell vector with all the pixels of hypoxic events of this pocket
        temp_pocket=Overall_OxygenSinks_Pxllist_rec(k,:);
        %create a logical vector that has tru values in only the frame that
        %this pocket had an event
        temp_pocket_bool= ~cellfun(@isempty,temp_pocket);

        % finding the individual events
        hpevents_temp=regionprops(temp_pocket_bool,'Area','PixelIdxList');

        if length(hpevents_temp)>Table_OxygenSinks_Out.NumOxySinkEvents(k)

            hpevents_temp=[hpevents_temp(1:Table_OxygenSinks_Out.NumOxySinkEvents(k))];

        end
        
        for j=1:length(hpevents_temp) %for every event produced from that pocket 


            %initiate a struct to capture the parameters of the event
            props_event_temp = struct(...
            'Area', [], ...
            'PixelIdxList', [], ...
            'Centroid', [], ...
            'Circularity', [], ...
            'FilledArea', [], ...
            'EquivDiameter', [], ...
            'Perimeter', []);

            for f=1:length(hpevents_temp(j).PixelIdxList) %go through all the frames belonging to that hypoxic event

                % Create a new empty logical vector to have all pocket-frames (events) that has the original dimensions
                Single_BW=false(size(IM_BW(:,:,1)));


                Single_BW(Overall_OxygenSinks_Pxllist_rec{k,hpevents_temp(j).PixelIdxList(f)})=true;

                props_eventframe_single=regionprops(Single_BW,'Area','PixelIdxList','Centroid','Circularity',"FilledArea","EquivDiameter","Perimeter");


                if isempty([props_event_temp.Area])
                    
                    props_event_temp=[props_eventframe_single(:)];

                else

                    props_event_temp=cat(1, [props_event_temp(:)], [props_eventframe_single(:)]);
                end
            
            end


            centroidxy=[props_event_temp.Centroid];
            

            % capturing the results
            Area_um(counter)=mean([props_event_temp.Area])*Table_OxygenSinks_Out.Pixelsize{k}^2;
            Area_norm(counter)=(mean([props_event_temp.Area])*Table_OxygenSinks_Out.Pixelsize{k}^2)/Table_OxygenSinks_Out.RecAreaSize{k};
            FilledArea_um(counter)=mean([props_event_temp.FilledArea])*Table_OxygenSinks_Out.Pixelsize{k}^2;
            FilledArea_norm(counter)=(mean([props_event_temp.FilledArea])*Table_OxygenSinks_Out.Pixelsize{k}^2)/Table_OxygenSinks_Out.RecAreaSize{k};
            Centroid_x(counter)=mean(centroidxy(1:2:end));
            Centroid_y(counter)=mean(centroidxy(2:2:end));
            Circularity(counter)= mean([props_event_temp.Circularity]);
            Perimeter_um(counter)=mean([props_event_temp.Perimeter])*Table_OxygenSinks_Out.Pixelsize{k};
            Diameter_um(counter)=mean([props_event_temp.EquivDiameter])*Table_OxygenSinks_Out.Pixelsize{k};
            OxySink_Pxls_event{counter}=unique(cell2mat({props_event_temp.PixelIdxList}'));

            counter=counter+1;

        end

    end

    Eventspecificmetrics_rec=table(Area_um,Area_norm,FilledArea_um,FilledArea_norm,Centroid_x,Centroid_y,...
        Circularity,Perimeter_um,Diameter_um,OxySink_Pxls_event);


%%%%%%%%%%%%%%%%%%%%
%%

%calculate the scalling factor 
rescaleF_Veins=size(Veins,1)/size(IM_BW,1);
rescaleF_Arter=size(Arteries,1)/size(IM_BW,1);
%I assume that ceil(size(Veins,1)/rescaleF) == size(Map,1). It is possible
%that for some cases this will not be the case. I need to include some
%checks here for the script to run smoothly. 

%%
% Finding the all the "islands" of non zero that correspond to the vasculature
Vein_regions=regionprops(Veins,'PixelIdxList');
Arte_regions=regionprops(Arteries,'PixelIdxList');

% Getting the indices of the pockets' centroids
Centroids_Idx=[round(Eventspecificmetrics_rec.Centroid_y), round(Eventspecificmetrics_rec.Centroid_x)];



Veins_Distance_all=nan(height(Eventspecificmetrics_rec),1);
Veins_Distance_centroid=nan(height(Eventspecificmetrics_rec),1);
Arteries_Distance_all=nan(height(Eventspecificmetrics_rec),1);
Arteries_Distance_centroid=nan(height(Eventspecificmetrics_rec),1);
% Go through the different pocket events
for i=1:height(Eventspecificmetrics_rec)


    Pixels=Eventspecificmetrics_rec.OxySink_Pxls_event{i};

    [row1, col1] = ind2sub(size(IM_BW(:,:,1)), Pixels);
    % the coordinates of all pixels for this pocket
    coords1=[col1, row1]; 
    % the coordinates of mean centroid for this pocket
    coords_c=[round(Eventspecificmetrics_rec.Centroid_y(i)),round(Eventspecificmetrics_rec.Centroid_x(i))];

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


%% Export as a csv
exporttable=table(Arteries_Distance_centroid,Arteries_Distance_all,Veins_Distance_centroid,Veins_Distance_all);

Output_xlsx=sprintf('Vascular_Analysis_Events.xlsx');
writetable(exporttable,Output_xlsx);

%%
cd(Masterfolder) %return to the masterfolder before going to the next data folder

end

clear all



%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTIONS

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

%% Function to load tiff files

function [oimg,Miu,SD]= loadtiff(path)
% Copyright (c) 2012, YoonOh Tak
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are 
% met:
% 
%     * Redistributions of source code must retain the above copyright 
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright 
%       notice, this list of conditions and the following disclaimer in 
%       the documentation and/or other materials provided with the distribution
%     * Neither the name of the Gwangju Institute of Science and Technology (GIST), Republic of Korea nor the names 
%       of its contributors may be used to endorse or promote products derived 
%       from this software without specific prior written permission.
%       

tStart = tic;
warn_old = warning('off', 'all'); % To ignore unknown TIFF tag.

% Check directory and file existence
path_parent = pwd;
[pathstr, ~, ~] = fileparts(path);
if ~isempty(pathstr) && ~exist(pathstr, 'dir')
    error 'Directory is not exist.';
end
if ~exist(path, 'file')
    error 'File is not exist.';
end

% Open file
file_opening_error_count = 0;
while ~exist('tiff', 'var')
    try
        tiff = Tiff(path, 'r');
    catch
        file_opening_error_count = file_opening_error_count + 1;
        pause(0.1);
        if file_opening_error_count > 5 % automatically retry to open for 5 times.
            reply = input('Failed to open the file. Do you wish to retry? Y/n: ', 's');
            if isempty(reply) || any(upper(reply) == 'Y')
                file_opening_error_count = 0;
            else
                error(['Failed to open the file ''' path '''.']);
            end
        end
    end
end

% Load image information
tfl = 0; % Total frame length
tcl = 1; % Total cell length
while true
    try
        tfl = tfl + 1; % Increase frame count
        iinfo(tfl).w       = tiff.getTag('ImageWidth');
        iinfo(tfl).h       = tiff.getTag('ImageLength');
        iinfo(tfl).spp     = tiff.getTag('SamplesPerPixel');
        iinfo(tfl).bps     = tiff.getTag('BitsPerSample');
        iinfo(tfl).color   = iinfo(tfl).spp > 2; % Grayscale: 1(real number) or 2(complex number), Color: 3(rgb), 4(rgba), 6(rgb, complex number), or 8(rgba, complex number)
        iinfo(tfl).complex = any(iinfo(tfl).spp == [2 6 8]);

        if tiff.getTag('Photometric') == Tiff.Photometric.Palette
            disp('Warning: Palette information has been removed.');
        end
    catch
        error 'Failed to load image';
    end
    
    %% Sample format
    switch tiff.getTag('SampleFormat')
        % Unsupported Matlab data type: char, logical, cell, struct, function_handle, class.
        case Tiff.SampleFormat.UInt
            switch tiff.getTag('BitsPerSample')
                case 8
                    iinfo(tfl).sf = 'uint8';
                case 16
                    iinfo(tfl).sf = 'uint16';
                case 32
                    iinfo(tfl).sf = 'uint32';
            end
        case Tiff.SampleFormat.Int
            switch tiff.getTag('BitsPerSample')
                case 8
                    iinfo(tfl).sf = 'int8';
                case 16
                    iinfo(tfl).sf = 'int16';
                case 32
                    iinfo(tfl).sf = 'int32';
            end
        case Tiff.SampleFormat.IEEEFP
            switch tiff.getTag('BitsPerSample')
                case 32
                    iinfo(tfl).sf = 'single';
                case 64
                    iinfo(tfl).sf = 'double';
            end
        otherwise
            % (Unsupported)Void, ComplexInt, ComplexIEEEFP
            error 'Unsupported Matlab data type. (char, logical, cell, struct, function_handle, class)';
    end

    if tfl > 1
        % If tag information is changed, make a new cell
        if iinfo(tfl-1).w ~= iinfo(tfl).w || ...
            iinfo(tfl-1).h ~= iinfo(tfl).h || ...
            iinfo(tfl-1).spp ~= iinfo(tfl).spp || ...
            any(iinfo(tfl-1).sf ~= iinfo(tfl).sf) || ...
            iinfo(tfl-1).color ~= iinfo(tfl).color || ...
            iinfo(tfl-1).complex ~= iinfo(tfl).complex
            tcl = tcl + 1; % Increase cell count
            iinfo(tfl).fid = 1; % First frame of this cell
        else
            iinfo(tfl).fid = iinfo(tfl-1).fid + 1;
        end
    else
        iinfo(tfl).fid = 1; % Very first frame of this file
    end
    iinfo(tfl).cid = tcl; % Cell number of this frame
    
    if tiff.lastDirectory(), break; end;
    tiff.nextDirectory();
end

% Load image data
if tcl == 1 % simple image (no cell)
    if ~iinfo(tfl).color
       oimg = zeros(iinfo(tfl).h,iinfo(tfl).w,iinfo(tfl).fid, iinfo(tfl).sf); % Grayscale image
    else
       oimg = zeros(iinfo(tfl).h,iinfo(tfl).w,iinfo(tfl).spp,iinfo(tfl).fid, iinfo(tfl).sf); % Color image
    end
    for tfl = 1:tfl
        tiff.setDirectory(tfl);
        temp = tiff.read();
        if iinfo(tfl).complex
            temp = temp(:,:,1:2:end-1,:) + temp(:,:,2:2:end,:)*1i;
        end
        if ~iinfo(tfl).color
            oimg(:,:,iinfo(tfl).fid) = temp; % Grayscale image
        else
            oimg(:,:,:,iinfo(tfl).fid) = temp; % Color image
        end
    end
else % multiple image (multiple cell)
    oimg = cell(tcl, 1);
    for tfl = 1:tfl
        tiff.setDirectory(tfl);
        temp = tiff.read();
        if iinfo(tfl).complex
            temp = temp(:,:,1:2:end-1,:) + temp(:,:,2:2:end,:)*1i;
        end
        if ~iinfo(tfl).color
            oimg{iinfo(tfl).cid}(:,:,iinfo(tfl).fid) = temp; % Grayscale image
        else
            oimg{iinfo(tfl).cid}(:,:,:,iinfo(tfl).fid) = temp; % Color image
        end
    end
end

%% Close file
tiff.close();
cd(path_parent);
warning(warn_old);

Miu = squeeze(mean(mean(double(oimg),1),2))';
SD = squeeze(std(double(oimg),[],[1,2]))';
oimg=single(oimg);

display(sprintf('New file was loaded successfully. Elapsed time : %.3f s.', toc(tStart)));
end

%% Unique for cell array of string
%version 1.1.0.0 (1.05 KB) by Wei-Rong Chen
% This function 'uniqueStrCell' performs 'unique' for cell array of string.

function out=uniqueStrCell(inputStrCell)
% This function 'uniqueStrCell' performs 'UNIQUE' for cell array of string.
% The output cell 'out' will include only string cells and numeric cells converted to strings
% , and exclude NaN and empty cells.
% Example: 
% inputStrCell={'ek','wekf', 29, NaN, [],'we'};
% out = uniqueStrCell(inputStrCell);
% >> out = {'ek'    'we'    'wekf'    '29'}
%
% Acknowledgement: 
% This function greatly benefits from Jan Simon's comments. The
% previous version was errorful. 
% See 'unique' for more information
%Weirong Chen   Apr-14-2014
out=[];
if nargin<1, display('Not enough input argument!'); return;end;
A=cellfun('isclass', inputStrCell, 'char'); 
B=cellfun(@isnumeric, inputStrCell); 
C=cellfun(@isnan, inputStrCell,'UniformOutput',false);
C=logical(cell2mat(cellfun(@sum, C,'UniformOutput',false)));
D=cellfun(@isempty, inputStrCell);
numCell = cellfun(@num2str,inputStrCell(logical(B-C-D)),'UniformOutput',false);
numCell = unique(numCell);
strCell = unique(inputStrCell(A));
out = [strCell numCell];
end

%% RemoveSheet123 - removes the sheets that are automatically added to excel
% file. 

function RemoveSheet123(excelFileName,sheetName)

% When Matlab writes data to a new Excel file, the Excel software
% automatically creates 3 sheets (the names are depended on the user
% languade). This appears even if the user has defined the sheet name to be
% added. 
%
% Usage:
% RemoveSheet123(excelFileName) - remove "sheet1", "sheet2","sheet3" from
% the excel file. excelFileName is a string of the Excel file name.
% RemoveSheet123(excelFileName,sheetName) - enables the user to enter the
% sheet name when the language is other than English.
% sheetName is the default sheet name, without the number.
%
%
%                       Written by Noam Greenboim
%                       www.perigee.co.il
%
%% check input arguments
if nargin < 1 || isempty(excelFileName)
    error('Filename must be specified.');
end
if ~ischar(excelFileName)
    error('Filename must be a string.');
end
try
    excelFileName = validpath(excelFileName);
catch 
    error('File not found.');
end
if nargin < 2
    sheetName = 'Sheet'; % EN: Sheet, DE: Tabelle, HE: âéìéåï , etc. (Lang. dependent)
else
    if ~ischar(sheetName)
        error('Default sheet name must be a string.');
    end
end
%%
% Open Excel file.
objExcel = actxserver('Excel.Application');
objExcel.Workbooks.Open(excelFileName); % Full path is necessary!
% Delete sheets.
try
      % Throws an error if the sheets do not exist.
      objExcel.ActiveWorkbook.Worksheets.Item([sheetName '1']).Delete;
      fprintf('\nsheet #1 - deleted.')
      objExcel.ActiveWorkbook.Worksheets.Item([sheetName '2']).Delete;
      fprintf('\nsheet #2 - deleted.')
      objExcel.ActiveWorkbook.Worksheets.Item([sheetName '3']).Delete;
      fprintf('\nsheet #3 - deleted.\n')
catch
    fprintf('\n')
    O=objExcel.ActiveWorkbook.Worksheets.get;
    if O.Count==1
        error('Can''t delete the last sheet. Excel file must containt at least one sheet.')
    else
      warning('Problem occured. Check excel file.'); 
    end
end
% Save, close and clean up.
objExcel.ActiveWorkbook.Save;
objExcel.ActiveWorkbook.Close;
objExcel.Quit;
objExcel.delete;
end
function filenameOut = validpath(filename)
    % VALIDPATH builds a full path from a partial path specification
    %   FILENAME = VALIDPATH(FILENAME) returns a string vector containing full
    %   path to a file. FILENAME is string vector containing a partial path
    %   ending in a file or directory name. May contain ..\  or ../ or \\. The
    %   current directory (pwd) is prepended to create a full path if
    %   necessary. On UNIX, when the path starts with a tilde, '~', then the
    %   current directory is not prepended.
    %
    %   See also XLSREAD, XLSWRITE, XLSFINFO.
    
    %   Copyright 1984-2012 The MathWorks, Inc.
    
    %First check for wild cards, since that is not supported.
    if strfind(filename, '*') > 0
        error(message('MATLAB:xlsread:Wildcard', filename));
    end
    
    % break partial path in to file path parts.
    [Directory, file, ext] = fileparts(filename);
    if ~isempty(ext)
        filenameOut = getFullName(filename);
    else
        extIn = matlab.io.internal.xlsreadSupportedExtensions;
        for i=1:length(extIn)
            try                                                                %#ok<TRYNC>
                filenameOut = getFullName(fullfile(Directory, [file, extIn{i}]));
                return;
            end
        end
        error(message('MATLAB:xlsread:FileDoesNotExist', filename));    
    end
end
function absolutepath=abspath(partialpath)
    
    % parse partial path into path parts
    [pathname, filename, ext] = fileparts(partialpath);
    % no path qualification is present in partial path; assume parent is pwd, except
    % when path string starts with '~' or is identical to '~'.
    if isempty(pathname) && strncmp('~', partialpath, 1)
        Directory = pwd;
    elseif isempty(regexp(partialpath,'(.:|\\\\)', 'once')) && ...
            ~strncmp('/', partialpath, 1) && ...
            ~strncmp('~', partialpath, 1);
        % path did not start with any of drive name, UNC path or '~'.
        Directory = [pwd,filesep,pathname];
    else
        % path content present in partial path; assume relative to current directory,
        % or absolute.
        Directory = pathname;
    end
    
    % construct absolute filename
    absolutepath = fullfile(Directory,[filename,ext]);
end
function filename = getFullName(filename)
    FileOnPath = which(filename);
    if isempty(FileOnPath)
        % construct full path to source file
        filename = abspath(filename);
        if isempty(dir(filename)) && ~isdir(filename)
            % file does not exist. Terminate importation of file.
            error(message('MATLAB:xlsread:FileDoesNotExist', filename));
        end
    else
        filename = FileOnPath;
    end
end
