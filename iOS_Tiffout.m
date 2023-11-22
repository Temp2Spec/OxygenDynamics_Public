% 26 September 2023
%% Initial parameters

smooth = 10;  % This is used to define how big the convolution frame will be. This value shows the number of pixels around each pixel.

if exist('PiSz','var')
    PixelSize = PiSz; %this is taken from the wrapper. Change to specific number if running the script alone
else
    fprintf('Script was call individually. No info from wrapper available! Pixel size is set to 2.5um \n');
    PixelSize =2.5;
end

if exist('SFs','var')
    fs=SFs;    % your frequency of the sampling [Hz]. This is taken from the wrapper. Change to specific number if running the script alone
else
    fprintf('Script was call individually. No info from wrapper available! Sampling Freq is set to 1Hz \n');
    fs =1;
end


Pixel_frame=smooth*2; %this is a frame that will be clipped around the imaging data to exclude artifacts resulting from motion correction.
ROIsize=20; %this is the size (in microns) of edge for the square bins the recorded tissue will be divided
PercentileDetectionThres=99; %this is the statistical threshold for the detection of oxygen sinks in the detrended/z-transformed/convolved/smoothed data
%% STEP 1 Read files and get some basic stats of the signal
fprintf('Loading data... \n');
Tifffiles = dir('*.tif'); 
IM_Raw=[];

for k=1:height(Tifffiles) %go through the itf files in the data folder
    
    DIR = fullfile(Tifffiles(k).folder,'\', Tifffiles(k).name);
  
    [IM_Raw,Miu,~]=loadtiff(DIR); 
    
end

%this is to be used later when naming output files
DatafileID=Tifffiles(1).folder(max(strfind(Tifffiles(1).folder,'\'))+1:end);

%that is the recording duration 
RecDur=size(IM_Raw,3)*fs;
 
clear k DIR 
%% STEP 2 resize the figure to the same dimensions as the oxygen bioluminescence recordings (512x512pxl)
%% Also, invert the data to get bright spots where high concentration of Haemoglobin exist
IM_resize=imcomplement(imresize(IM_Raw,512/size(IM_Raw,1)));

%% STEP 3 Converting the signal to df/f

[IM_Raw_dFoF_global,IM_Raw_dFoF_frame] = normalize_to_dFoF(IM_resize);

IM_Raw_dFoF_global = uint8((IM_Raw_dFoF_global - min(IM_Raw_dFoF_global(:))) * (255 / (max(IM_Raw_dFoF_global(:)) - min(IM_Raw_dFoF_global(:)))));
IM_Raw_dFoF_frame = uint8((IM_Raw_dFoF_frame - min(IM_Raw_dFoF_frame(:))) * (255 / (max(IM_Raw_dFoF_frame(:)) - min(IM_Raw_dFoF_frame(:)))));

%% STEP 4 Detrending 
fprintf('Detrending... This might take some time! \n');

tic;
%detrend the raw data. That's were the analysis will retun to extract the signal. 
IM_Notrend_dFoF_global=detrend_custom(single(IM_Raw_dFoF_global),3);

IM_Notrend_dFoF_global = single(IM_Notrend_dFoF_global);
toc;
%% STEP 6 Normalize according to the standard deviation, to separate white and black spots
% I produced three different z-scored methods. It looks like the 'IM_Zframetime' is better but I need to investigate
fprintf('Z-scoring data... \n');

tic;

[~, IM_Zframetime] = normalize_to_z_stat(IM_Notrend_dFoF_global);

toc;


%% STEP 7 Multiply the z-scored amplitude of each pixel with the mean amplitude of a 21x21 pixel window around it ('smooth' in four directions from each pixel).

fprintf('Convolving data... \n');

tic;
% defining the convolution window. The '+1'makes sure there is alway a center of the window
ConvWin=ones(2*smooth+1,2*smooth+1)/(2*smooth+1)^2;
% Some more info on why I constructed the convolution kernel like that
% https://se.mathworks.com/matlabcentral/answers/358411-i-need-a-code-that-produce-a-moving-average-matrix-with-a-5-5-window

IM_Zframetime_conv=single(zeros(size(IM_Zframetime)));
for z=1:size(IM_Zframetime,3) 
   
     IM_Zframetime_conv(:,:,z)=single(conv2(IM_Zframetime(:,:,z),ConvWin,'same'));
    
end
clear z
toc;

%% STEP 8 Smooth the signal a bit more 
% This step maybe redundant when working for denoised data
disp('Smoothing data...');
IM_Zframetime_Conv_smoothed=smoothdata(IM_Zframetime_conv,3,'gaussian',smooth/2);

% IM_Zframetime_smoothed = (IM_Zframetime_smoothed - min(IM_Zframetime_smoothed(:))) * (255 / (max(IM_Zframetime_smoothed(:)) - min(IM_Zframetime_smoothed(:))));
% IM_Zframetime_smoothed = single(IM_Zframetime_smoothed)/255*(2^16-1);
IM_Zframetime_Conv_smoothed = uint8((IM_Zframetime_Conv_smoothed - min(IM_Zframetime_Conv_smoothed(:))) * (255 / (max(IM_Zframetime_Conv_smoothed(:)) - min(IM_Zframetime_Conv_smoothed(:)))));


%% Saving outputs
disp('Saving manipulated data matrices');


 if not(isfolder('OxygenSinks_Output')) || strOW=='Y' %if there not a folder with previous oxygen sinks output or we want to overwrite previous data

     OxySinksfolder='OxygenSinks_Output';
     ManualCurOxySinksfolder='ManualCurOxySinksData';    
     OxySurgesfolder='OxygenSurges_Output';
     Images_Proce_folder='Images_Processed';
 
 else %if there is a folder with previous output and we do not want to overwrite.

     OxySinksfolder=['OxygenSinks_Output_',datestr(now,30)]; %adds an identifier in the form of yyyymmddTHHMMSS
     ManualCurOxySinksfolder=['ManualCurOxySinksData_',datestr(now,30)];
     OxySurgesfolder=['OxygenSurges_Output',datestr(now,30)];
     Images_Proce_folder=['Images_Processed',datestr(now,30)];
 end 
 
mkdir(OxySinksfolder)
mkdir(ManualCurOxySinksfolder)
mkdir(OxySurgesfolder)
mkdir(Images_Proce_folder)

clear options;
options.overwrite = true;



%export the dF/F matrix from detrended (based on the global mean)
IM_Notrend_dFoF_global = uint8((IM_Notrend_dFoF_global - min(IM_Notrend_dFoF_global(:))) * (255 / (max(IM_Notrend_dFoF_global(:)) - min(IM_Notrend_dFoF_global(:)))));
saveastiff(uint16(single(IM_Notrend_dFoF_global)/255*(2^16-1)), [Tifffiles(1).folder,'\',Images_Proce_folder,'\', 'IM_Notrend_dFoF', DatafileID, '.tif'],options);

%export the dF/F matrix from raw (based on the global mean)
IM_Raw_dFoF_global = uint8((IM_Raw_dFoF_global - min(IM_Raw_dFoF_global(:))) * (255 / (max(IM_Raw_dFoF_global(:)) - min(IM_Raw_dFoF_global(:)))));
saveastiff(uint16(single(IM_Raw_dFoF_global)/255*(2^16-1)), [Tifffiles(1).folder,'\',Images_Proce_folder,'\', 'IM_Raw_dFoF', DatafileID, '.tif'],options);

%export the conv/smoothed matrix
IM_Zframetime_Conv_smoothed = uint8((IM_Zframetime_Conv_smoothed - min(IM_Zframetime_Conv_smoothed(:))) * (255 / (max(IM_Zframetime_Conv_smoothed(:)) - min(IM_Zframetime_Conv_smoothed(:)))));
saveastiff(uint16(single(IM_Zframetime_Conv_smoothed)/255*(2^16-1)), [Tifffiles(1).folder,'\',Images_Proce_folder,'\','IM_Conv', DatafileID, '.tif'],options);

     
%% FUNCTIONS %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Function to fit and remove trends for each pixel timeseries

% Antonis Asiminas

function IM2=detrend_custom(IM,dimensions)

IM2=nan(size(IM));

switch dimensions

case 3  

    % Detrend the original data on a pixel by pixel basis
    for xx=1:size(IM,1)    
        parfor yy=1:size(IM,2)
        
            S1=(squeeze(IM(xx,yy,:)))';       
            [p,s,mu] = polyfit(1:numel(S1), S1, 3);            
            signaltrend=polyval(p,1:numel(S1),[],mu);         
            IM2(xx,yy,:) = S1- signaltrend;    
        end
    end 
    clear xx yy S1 p

case 2

    for xx=1:size(IM,1) 

        S1=IM(xx,:);
        [p,s,mu] = polyfit(1:numel(S1), S1, 5);
        signaltrend=polyval(p,1:numel(S1),[],mu);
        IM2(xx,:) = S1- signaltrend;

    end

end

end

%% Function to z-score data

function [IM2_zsc_frame, IM2_zsc_frame_time] =normalize_to_z_stat(IM)


Miu_frame = squeeze(mean(mean(double(IM),1),2))';
SD_frame = squeeze(std(double(IM),[],[1,2]))';

%H = kstest(IM(:)); %test for the normality 
% I commented this out because the data are never normaly distributed (H=1) therefore this step was taking unecessary time and power.
%The distribution is now z-scored and spread around 0
IM2_zsc_frame=NaN(size(IM));
%IM2_zsc_time=NaN(size(IM));
IM2_zsc_frame_time=NaN(size(IM));


% if H==0
%     for ii=1:length(Miu_frame)
% %         IM2_zsc(:,:,ii)=IM2(:,:,ii)-(Miu3(ii)- SD3(ii));  %major in Felix2/Felix3 file
%           IM2_zsc_frame(:,:,ii)=(IM(:,:,ii)-Miu_frame(ii))./sqrt(SD_frame(ii));
% %         if wrt==1
% %             imwrite(uint16(IM2_zsc(:,:,ii)),[DIR 'Test_zscored_std2.tif'],'Resolution',[96 96], 'Compression','none','writemode','append')
% %         end
%     end
%     
% end
% %i splitted this to increase performance
% 
% if H~=0
    for ii=1:length(Miu_frame)
%         IM2_zsc2(:,:,ii)=IM2(:,:,ii)-(Miu3(ii)- sqrt(SD3(ii))); %major in Felix2/Felix3 file
          IM2_zsc_frame(:,:,ii)=(IM(:,:,ii)-Miu_frame(ii))/SD_frame(ii);
%         if wrt==1
%             imwrite(uint16(IM2_zsc(:,:,ii)),[DIR 'Test_zscored_sqrt2.tif'],'Resolution',[96 96], 'Compression','none','writemode','append')
%         end
    end
% end; clear ii s;
%hist(IM2_zsc(:),100);


%instead of the frame by frame zscoring, I also produce the zscored data over time

%Miu_time=mean(IM, 3);
%SD_time=std(IM,[], 3);

Miu_Zframe=mean(IM2_zsc_frame, 3);
SD_Zframe=std(IM2_zsc_frame,[], 3);

for ii=1:length(Miu_frame)
    
    %IM2_zsc_time(:,:,ii)=(IM(:,:,ii)-Miu_time)./sqrt(SD_time);
    
    IM2_zsc_frame_time(:,:,ii)=(IM2_zsc_frame(:,:,ii)-Miu_Zframe)./sqrt(SD_Zframe);
end

IM2_zsc_frame_time=single(IM2_zsc_frame_time);

IM2_zsc_frame=single(IM2_zsc_frame);
end



%% Function to df/f data

function [IM_dFoF_global,IM_dFoF_frame] =normalize_to_dFoF(IM)

IM_dFoF_global=NaN(size(IM));
IM_dFoF_frame=NaN(size(IM));
%compute the mean across all pixels/frames
Miu=mean(IM(:));
%demean and scale to this gloabl mean 
IM_dFoF_global(:,:,:)=(IM(:,:,:)-Miu)./Miu;

%compute the mean across all pixels for each frame
Miu_frame=mean(mean(IM, 1), 2);

for ii=1:size(IM,3)
       
    IM_dFoF_frame(:,:,ii)=(IM(:,:,ii)-Miu_frame(ii))/Miu_frame(ii);
end

%Replace the values where IM_dFoF==-inf with the minimum real number.
IM_dFoF_global(IM_dFoF_global==-inf) = min(IM_dFoF_global(isfinite(IM_dFoF_global)));
%Replace the values where IM_dFoF==+inf with the maximum real number.
IM_dFoF_global(IM_dFoF_global==inf)  = max(IM_dFoF_global(isfinite(IM_dFoF_global)));

%Replace the values where IM_dFoF==-inf with the minimum real number.
IM_dFoF_frame(IM_dFoF_frame==-inf) = min(IM_dFoF_frame(isfinite(IM_dFoF_frame)));
%Replace the values where IM_dFoF==+inf with the maximum real number.
IM_dFoF_frame(IM_dFoF_frame==inf)  = max(IM_dFoF_frame(isfinite(IM_dFoF_frame)));

end


%% Function to save tiff files

function res = saveastiff(data, path, options)
% options.color
%   : true or FALSE
%   : If this is true, third dimension should be 3 and the data is saved as a color image.
% options.compress
%   : 'no', 'lzw', 'jpeg' or 'adobe'.
%     Compression type.
%       'no'    : Uncompressed(Default)
%       'lzw'   : lossless LZW
%       'jpeg'  : lossy JPEG (When using JPEG compression, ImageWidth,
%                 ImageLength, and RowsPerStrip must be multiples of 16.)
%       'adobe' : lossless Adobe-style
% options.jpegquality
%   : JPEG compression qualtiy. A value between 1 and 100
% options.message
%   : TRUE or false.
%     If this is false, all messages are skipped. 
% options.append
%   : true or FALSE
%     If path is exist, the data is appended to an existing file.
%     If path is not exist, this options is ignored.
% options.overwrite
%   : true or FALSE
%     Overwrite to an existing file.
% options.big 
%   : true or FALSE, 
%     Use 64 bit addressing and allows for files > 4GB
% 
% Defalut value of 'options' is
%     options.color     = false;
%     options.compress  = 'no';
%     options.message   = true;
%     options.append    = false;
%     options.overwrite = false;
%     options.big       = false;
% 
% res : Return value. It is 0 when the function is finished with no error.
%       If an error is occured in the function, it will have a positive
%       number (error code).
%
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
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
% POSSIBILITY OF SUCH DAMAGE.

tStart = tic;
errcode = 0;
try
%% Init options parameter    
if nargin < 3 % Use default options
    options.color = false;
    options.compress = 'no';
    options.message = true;
    options.append = false;
    options.overwrite = false;
end
if ~isfield(options, 'message'),   options.message   = true; end
if ~isfield(options, 'append'),    options.append    = false; end
if ~isfield(options, 'compress'),  options.compress  = 'no';  end
if ~isfield(options, 'color'),     options.color     = false; end
if ~isfield(options, 'overwrite'), options.overwrite = false; end
if  isfield(options, 'big') == 0,  options.big       = false; end

switch class(data)
    case {'uint8', 'uint16', 'uint32', 'int8', 'int16', 'int32', 'single', 'double', 'uint64', 'int64'}
    otherwise
        errcode = 5; assert(false);
end

if isempty(data), errcode = 1; assert(false); end
if (options.color == false && ndims(data) > 3) || ...
   (options.color == true && ndims(data) > 4)
    % Maximum dimension of a grayscale image is 3 of [height, width, frame]
    % Maximum dimension of a color image is 4 of [height, width, color, frame]
    errcode = 2; assert(false);
end

%% Get image informations
% http://www.awaresystems.be/imaging/tiff/tifftags/photometricinterpretation.html
if ~options.color
    if ndims(data) >= 4, errcode = 2; assert(false); end
    [height, width, depth] = size(data);
    tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
%     tagstruct.Photometric = Tiff.Photometric.MinIsWhite;
%     tagstruct.Photometric = Tiff.Photometric.Mask;
%     tagstruct.Photometric = Tiff.Photometric.Separated;
else
    if ndims(data) >= 5, errcode = 2; assert(false); end
    [height, width, cc, depth] = size(data); % cc: color channels. 3: rgb, 4: rgb with alpha channel
    if cc ~= 3 && cc ~= 4, errcode = 3; assert(false); end
    tagstruct.Photometric = Tiff.Photometric.RGB;
%     tagstruct.Photometric = Tiff.Photometric.CIELab;
%     tagstruct.Photometric = Tiff.Photometric.ICCLab;
%     tagstruct.Photometric = Tiff.Photometric.ITULab;
%     (Unsupported)tagstruct.Photometric = Tiff.Photometric.Palette;
%     (Unsupported)tagstruct.Photometric = Tiff.Photometric.YCbCr;
end
tagstruct.ImageLength = height;
tagstruct.ImageWidth = width;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky; % (RGB RGB,RGB RGB,RGB RGB), http://www.awaresystems.be/imaging/tiff/tifftags/planarconfiguration.html
% (Unsupported)tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Separate; % (RRR RRR, GGG GGG, BBB BBB), 
% http://www.awaresystems.be/imaging/tiff/tifftags/planarconfiguration.html

%% Complex number
% http://www.awaresystems.be/imaging/tiff/tifftags/samplesperpixel.html
if ~options.color && isreal(data) % Grayscale image with real numbers
    tagstruct.SamplesPerPixel = 1;
    data = reshape(data, height, width, 1, depth);
elseif ~options.color && ~isreal(data) % Grayscale image with complex numbers
    tagstruct.SamplesPerPixel = 2;
    data = reshape([real(data) imag(data)], height, width, 2, depth);
elseif options.color && isreal(data) % Color image with real numbers
    tagstruct.SamplesPerPixel = cc;
    if cc == 4
        tagstruct.ExtraSamples = Tiff.ExtraSamples.AssociatedAlpha; % The forth channel is alpha channel
    end
    data = reshape(data, height, width, cc, depth);
elseif options.color && ~isreal(data) % Color image with complex numbers
    tagstruct.SamplesPerPixel = cc * 2;
    if cc == 3
        tagstruct.ExtraSamples = repmat(Tiff.ExtraSamples.Unspecified, 1, 3); % 3(real)+3(imag) = 6 = 3(rgb) + 3(Extra)
    else
        tagstruct.ExtraSamples = repmat(Tiff.ExtraSamples.Unspecified, 1, 5); % 4(real)+4(imag) = 8 = 3(rgb) + 5(Extra)
    end
    data = reshape([real(data) imag(data)], height, width, cc*2, depth);
end

%% Image compression
% http://www.awaresystems.be/imaging/tiff/tifftags/compression.html
switch lower(options.compress)
    case 'no'
        tagstruct.Compression = Tiff.Compression.None;
    case 'lzw'
        tagstruct.Compression = Tiff.Compression.LZW;
    case {'jpeg', 7}
        tagstruct.Compression = Tiff.Compression.JPEG;
        if mod(height, 16) ~= 0 || mod(width, 16) ~= 0
            tagstruct.Compression = Tiff.Compression.AdobeDeflate;
            disp('Warning: Image width and height must be multiples of 16 when using JPEG compression. The compression method has been automatically changed to AdobeDeflate.')
        else
            if isfield(options, 'jpegquality') && 1 <= options.jpegquality && options.jpegquality <= 100
                tagstruct.JPEGQuality = options.jpegquality;
            end
        end
    case 'adobe'
        tagstruct.Compression = Tiff.Compression.AdobeDeflate;
    otherwise
        % Use tag nubmer in http://www.awaresystems.be/imaging/tiff/tifftags/compression.html
        tagstruct.Compression = options.compress;
end

%% Sample format
% http://www.awaresystems.be/imaging/tiff/tifftags/sampleformat.html
switch class(data)
    % Unsupported Matlab data type: char, logical, cell, struct, function_handle, class.
    case {'uint8', 'uint16', 'uint32'}
        tagstruct.SampleFormat = Tiff.SampleFormat.UInt;
    case {'int8', 'int16', 'int32'}
        tagstruct.SampleFormat = Tiff.SampleFormat.Int;
        if options.color
            errcode = 4; assert(false);
        end
    case {'uint64', 'int64'}
        tagstruct.SampleFormat = Tiff.SampleFormat.IEEEFP;
        data = double(data);
    case {'single', 'double'}
        tagstruct.SampleFormat = Tiff.SampleFormat.IEEEFP;
    otherwise
        % (Unsupported)Void, ComplexInt, ComplexIEEEFP
        errcode = 5; assert(false);
end

%% Bits per sample
% http://www.awaresystems.be/imaging/tiff/tifftags/bitspersample.html
switch class(data)
    case {'uint8', 'int8'}
        tagstruct.BitsPerSample = 8;
    case {'uint16', 'int16'}
        tagstruct.BitsPerSample = 16;
    case {'uint32', 'int32'}
        tagstruct.BitsPerSample = 32;
    case {'single'}
        tagstruct.BitsPerSample = 32;
    case {'double', 'uint64', 'int64'}
        tagstruct.BitsPerSample = 64;
    otherwise
        errcode = 5; assert(false);
end

%% Rows per strip
tagstruct.RowsPerStrip = 512; % http://www.awaresystems.be/imaging/tiff/tifftags/rowsperstrip.html

%% Overwrite check
if exist(path, 'file') && ~options.append
    if ~options.overwrite
        errcode = 6; assert(false);
    end
end

%% Save path configuration
path_parent = pwd;
[pathstr, fname, fext] = fileparts(path);
if ~isempty(pathstr)
    if ~exist(pathstr, 'dir')
        mkdir(pathstr);
    end
    cd(pathstr);
end

%% Write image data to a file
file_opening_error_count = 0;
while ~exist('tfile', 'var')
    try
        if ~options.append % Make a new file
            s=whos('data');
            if s.bytes > 2^32-1 || options.big
                tfile = Tiff([fname, fext], 'w8'); % Big Tiff file
            else
                tfile = Tiff([fname, fext], 'w');
            end
        else
            if ~exist([fname, fext], 'file') % Make a new file
                s=whos('data');
                if s.bytes > 2^32-1 || options.big
                    tfile = Tiff([fname, fext], 'w8'); % Big Tiff file
                else
                    tfile = Tiff([fname, fext], 'w');
                end
            else % Append to an existing file
                tfile = Tiff([fname, fext], 'r+');
                while ~tfile.lastDirectory() % Append a new image to the last directory of an exiting file
                    tfile.nextDirectory();
                end
                tfile.writeDirectory();
            end
        end
    catch
        file_opening_error_count = file_opening_error_count + 1;
        pause(0.1);
        if file_opening_error_count > 5 % automatically retry to open for 5 times.
            reply = input('Failed to open the file. Do you wish to retry? Y/n: ', 's');
            if isempty(reply) || any(upper(reply) == 'Y')
                file_opening_error_count = 0;
            else
                errcode = 7;
                assert(false);
            end
        end
    end
end

for d = 1:depth
    tfile.setTag(tagstruct);
    tfile.write(data(:, :, :, d));
    if d ~= depth
       tfile.writeDirectory();
    end
end

tfile.close();
if exist('path_parent', 'var'), cd(path_parent); end

tElapsed = toc(tStart);
if options.message
    sprintf('The file was saved successfully. Elapsed time : %.3f s.', tElapsed);
end

catch exception
%% Exception management
    if exist('tfile', 'var'), tfile.close(); end
    switch errcode
        case 1
            if options.message, error '''data'' is empty.'; end
        case 2
            if options.message, error 'Data dimension is too large.'; end
        case 3
            if options.message, error 'Third dimesion (color depth) should be 3 or 4.'; end
        case 4
            if options.message, error 'Color image cannot have int8, int16 or int32 format.'; end
        case 5
            if options.message, error 'Unsupported Matlab data type. (char, logical, cell, struct, function_handle, class)'; end
        case 6
            if options.message, error 'File already exists.'; end
        case 7
            if options.message, error(['Failed to open the file ''' path '''.']); end
        otherwise
            if exist('fname', 'var') && exist('fext', 'var')
                delete([fname fext]);
            end
            if exist('path_parent', 'var'), cd(path_parent); end
            rethrow(exception);
    end
    if exist('path_parent', 'var'), cd(path_parent); end
end
res = errcode;
end
