%% Script for creating 3D plots of oxygen sinks and surges

% All subsequent scripts from Antonis Asiminas, PhD
% Center for Translational Neuromedicine
% University of Copenhagen, August 2023

%%
close all;
clear;

%% Choose the data and colour
uiwait(msgbox("Select a tif file with the masks for sinks or surges", 'Select tiff file'));
[filename, filepath] = uigetfile('*.*', 'Select a File');

c = uisetcolor([0 0.4470 0.7410]);

answer1 = questdlg('Do you want to save a vector? (takes long time)',...
    'Vector graphics?',...
    'Yes', 'No','No');


%% Load the data

DIR = fullfile(filepath,'\', filename);
[RecordingBW,~,~]=loadtiff(DIR);

% Dimensions of the matrix
[M, N, P] = size(RecordingBW);
%%
Epoch = questdlg('Do you want a different colour for a portion of the recording',...
    'What data',...
    'Yes', 'No','No');


switch Epoch
    case 'Yes'

        c2 = uisetcolor([1 0.4470 0.7410]);

        % Ask for user input for start and end along the z-axis
        prompt = 'Enter the start and end of the period with different colour in seconds (e.g. 6301,600).';
        dlgTitle = 'Epoch start-end';
        numLines = 2;
        defaultInput = {'301,600'};  % Default values in the input fields

        userInput = inputdlg(prompt, dlgTitle, numLines, defaultInput);

        % Check if the user clicked Cancel or entered an empty value
        if isempty(userInput)
            disp('Operation canceled or no value entered.');
        else
            % Convert the user input to numerical values
            Epochs_Limits = str2double(strsplit(userInput{1}, ','));
    
            % Check if the input is valid
            if ~any(isnan(Epochs_Limits))
                disp('You entered the following numerical values:');
                disp(Epochs_Limits);

                % Convert user input to numerical values
                startZ=Epochs_Limits(1:2:end);
                endZ=Epochs_Limits(2:2:end);

                if startZ>endZ
                    error('Invalid input. Start index should be < than end index.');
                elseif startZ<1
                    error('Invalid input. Start index should be > 1.');
                elseif endZ>size(RecordingBW,3)
                    error('Invalid input. End index should be < than length of recording.');
                end

            else
               
                error('Invalid input. Please enter numerical values.');
            end
  
        end

        
        
end



%% Generating and plotting 

newMatrix = single(zeros(M, N, 2*P - 1));
% Fill the new matrix with the original matrix and duplicates
for i = 1:P
    newMatrix(:, :, 2*i - 1) = RecordingBW(:, :, i); % Odd-indexed pages
    if i < P
        newMatrix(:, :, 2*i) = RecordingBW(:, :, i); % Even-indexed pages (except the last one)
    end
end
%%
% Define coordinates for the vertices of the cubes

[X, Y, Z] = meshgrid(1:N, 1:M, 1:0.5:P);
fig = figure('Visible', 'off');
% Define the aspect ratio. We can adapt this to what we want
aspect_ratio = [1 1 1];
% Set the aspect ratio of the axes
daspect(aspect_ratio);
s=isosurface(X, Y, Z, newMatrix,0.1);
p = patch(s);
isonormals(X, Y, Z, newMatrix,p)
% By typing this on comand line while looking at the figure
% [caz,cel] = view
% Get the azimuth (caz) and elevation (cel) angles for this plot.
view(3); % This can be updated to view(caz cel) to get a consistent point of view for the vector graphics

switch Epoch
    case 'No'
        
        set(p,'FaceColor',c);  
        set(p,'EdgeColor','none');
        camlight;
        lighting gouraud;
        zlabel('Time (s)');

    case 'Yes'
        %This code sets a different color for the specified section of the isosurface plot along the z-axis while leaving the rest of the isosurface with the default color.
        % Set the default color for the entire isosurface
        set(p, 'FaceColor', 'flat');
        set(p, 'EdgeColor', 'none');
        camlight;
        lighting gouraud;
        zlabel('Time (s)');

        % Get the isosurface vertices
        vertices = p.Vertices;

        % Get the faces of the isosurface
        faces = p.Faces;

        % Compute the z-values of the vertices
        zValues = vertices(:, 3);

        % Find the vertices within the specified z-range
        indices = find(zValues >= startZ & zValues <= endZ);

        % Create a color matrix for the vertices based on their Z-values
        vertexColors = repmat(c, size(vertices, 1), 1);

        % Set the different color for the specified section
        vertexColors(indices, :) = repmat(c2, length(indices), 1);

        % Update the 'FaceVertexCData' property to apply the specified colors
        p.FaceVertexCData = vertexColors;

        % Set the 'FaceColor' property to 'interp' to use FaceVertexCData
        set(p, 'FaceColor', 'interp');

end





%%
exportgraphics(fig,[filename(1:end-4),'.png'],'Resolution',1200)
save([filename(1:end-4),'.mat'], 'fig');
% You can load the interactive figure and replot it by typing this on the
% command window 
% loadedData = load('3D_sinks_interactive.mat');
% loadedFig = loadedData.fig;
% figure(loadedFig);


switch answer1

    case 'Yes'

        exportgraphics(fig,[filename(1:end-4),'.pdf'],'ContentType','vector')
end


%%clearvars -except RecordingBW
%%%%%%%%%%%%%%%% Functions %%%%%%%%%%%%%%%%%%%%%%%%%

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
    
    if tiff.lastDirectory(), break; end
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

fprintf('New file was loaded successfully. Elapsed time : %.3f s \n', toc(tStart));
end
