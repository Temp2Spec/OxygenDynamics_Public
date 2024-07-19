%% Script for the detection and analysis of oxygen dynamics in the mouse barrel cortex

% Provisional title: Optical investigation of oxygen dynamics in the murine cortex

% Each data folder (e.g 'Recording1') should contain ONLY one tif stack of original motion corrected data and one tif stack of denoised data.
% IMPORTANT!! We use this denoising approach 
% https://royerlab.github.io/aydin/getting_started/install.html denoising
% The denoised data should have the tag 'DENOISED' in the name of the tif stack.

% Overview of script (high level)
% Data are imported (both motion corrected and denoised)
% !Future version! Intrinsic Optical Signal Imaging data are imported along with a way to map them spatially to the oxygen dynamics data
% Preprocessing of data including detrending, z-scoring and convolution
% Identification of oxygen sinks seen as spatially confined hypoxic pockets 
% Identification of oxygen surges
% Calculation of events
% The output contains a series tables with data as well as ploted time traces of signal from identified hypoxic poxic
% The output can be curated further with the use of a MATLAB app and create ground-truth data for the training of detection algorithms

% All needed, custom fuctions are contained within the script. This is not good practice in general but it allows for the script to be portable.

% Initial script and selected functions from Ryszard Gomolka, dr. eng.

% All subsequent scripts from Antonis Asiminas, PhD
% Center for Translational Neuromedicine
% University of Copenhagen, November 2021

% Last updated 12 Sep 2018

%% Initial parameters

smooth = 10;  % This is used to define how big the convolution frame will be.
ThresholdMinsize = 100; % These are additional thresholds for the size of the hypoxic pocket (in pixels). 
ThresholdMaxsize = 6400; % !!I may need to have a size in um since the lense and magnification may change in future recordings!!
CircularityThres=0.3;
ThresholdMinsize_Surges=400; %!!!there is no real basis for this threshold (I need to check with felix)
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

Thesholdtime = 20*fs; % This was set after discussions with Felix. Hypoxic events from the same pocket do not happen in close temporal promimity (>20sec, 1Hz Fs).
ThresholMinddur = 3*fs; % This was set after discussions with Felix. Hypoxic events are sharp decreases in signal that last for several seconds (>3sec, 1Hz Fs)
ThresholMaxddur = 150*fs; % This is to exlude some long-lasting movement-related noise in the signal (currently >2.5min)
ThresholMinddur_Surges=10*fs; %this is a bit random but loosely connected to functional hyperemia temporal dynamics
Pixel_frame=smooth*2; %this is a frame that will be clipped around the imaging data to exclude artifacts due to motion correction.
PercentileDetectionThres=99; %this is the statistical threshold for the detection of hypoxic pockets in the detrended/z-transformed/convolved/smoothed data
%% STEP 1 Read files and get some basic stats of the signal
fprintf('Loading data... \n');
Tifffiles = dir('*.tif'); 
IM_Raw=[];
IM_NoNoise=[];
for k=1:height(Tifffiles) %go through the itf files in the data folder
    
    DIR = fullfile(Tifffiles(k).folder,'\', Tifffiles(k).name);
    if ~isempty(strfind(Tifffiles(k).name, 'DENOISED')) %if the file is the denoised
    
        [IM_NoNoise,Miu_NoNoise,~]=loadtiff(DIR); % load the tiff stack using a function that does not care about the size of the stack
    else
        [IM_Raw,Miu,~]=loadtiff(DIR); 
    end
end

%this is to be used later when naming output files
DatafileID=Tifffiles(1).folder(max(strfind(Tifffiles(1).folder,'\'))+1:end);

%that is the recording duration 
RecDur=size(IM_Raw,3)/fs;
 
clear k DIR 

%% STEP 2 Detrending data
fprintf('Detrending... This might take some time! \n');

tic;
%detrend the raw data. That's were the analysis will retun to extract the signal. 
IM_Raw_Notrend=detrend_custom(IM_Raw,3);
%Check if the denoised data is available 
if ~isempty(IM_NoNoise) %if so
    %Detrend this and all subsequenct analyses will be based on the detrended/denoised data .
    IM_Notrend=detrend_custom(IM_NoNoise,3); 
else %if not
    
    fprintf('Cannot find denoised data... Analysis will continue with motion-corrected-only data! \n')
    IM_Notrend=IM_Raw_Notrend; 
end

toc;

%% STEP 3 Defining the area of tissue that is being recorded.
% This is not as clearcut as one might though simply because of the tissue movement. 

%first calculate the mean signal for each frame across
Miu_Notrend=squeeze(mean(mean(IM_Notrend,1),2))';
AboveMean=NaN(size(IM_Notrend));
for f=1:length(Miu_Notrend)
    % I go through all the frames and ask which pixels are above mean of that frame
    AboveMean(:,:,f)=IM_Notrend(:,:,f)>Miu_Notrend(f);     
end
% convert to single to save some memory
AboveMean=single(AboveMean);
% I then collapse across the frames by getting the average value for each
% pixel. If it is 1 this means it was always above mean if it is closest to
% 0 this means that was never above mean. 
CollapsedArea=squeeze(mean(AboveMean,3));

% Get the pixels that are above the 30th percentile. This is a bit arbitrary and could change. 
% In that way pixels with 
Aboveback=CollapsedArea>prctile(CollapsedArea(:),30);

if size(AboveMean,3)<600 %for short recordings (less than 10minutes) modify the definition. 
    % Get the pixels that are above the 20th percentile. This is a bit arbitrary and could change. 
    Aboveback=CollapsedArea>prctile(CollapsedArea(:),20);

end


% Copying Aboveback variable to use it for filterring out sinks and sources
% that are not part of therec area
RecAreaFilter=Aboveback;

% Count how many pixels that area has and multiply with the squared pixelsize to get the surface in um squared
RecArea=(sum(Aboveback(:)))*PixelSize^2;

%divide the rec area in as many square ROIs as possible
%the edge size of the ROI is 50um.
counter=1;
RecArea_bins=[];
for i=1:size(Aboveback,1) %go pixel by pixel
    for j=1:size(Aboveback,2)
        if i+fix(25/PixelSize)*2-1<=size(Aboveback,1) && j+fix(25/PixelSize)*2-1<=size(Aboveback,2) %if the ROI is within the limits of the matrix
            if length(find(Aboveback(i:i+fix(25/PixelSize)*2-1,j:j+fix(25/PixelSize)*2-1)))>=((fix(25/PixelSize)*2)^2)*0.9 %and if at least 90% of the pixels are true
                RecArea_bins(counter,1)=i; %capturing the upper left corner of the ROI
                RecArea_bins(counter,2)=j;
                %eliminating the pixels of the frame from the pool
                Aboveback(i:i+fix(25/PixelSize)*2-1,j:j+fix(25/PixelSize)*2-1)=false;

                counter=counter+1;
            else
                continue
            end
        else
            continue
        end
    end
end

clear f Miu_Notrend CollapsedArea Aboveback AboveMean counter i j 

IM_Notrend = single(IM_Notrend);
IM_Raw_Notrend = single(IM_Raw_Notrend);

%% STEP 4 Normalize according to the standard deviation, to separate white and black spots
% I produced three different z-scored methods. It looks like the 'IM_Zframetime' is better but I need to investigate
fprintf('Z-scoring data... \n');

tic;

[~, IM_Zframetime] = normalize_to_z_stat(IM_Notrend);

toc;

%% STEP 5 Multiply the z-scored amplitude of each pixel with the mean amplitude of a 21x21 pixel window around it ('smooth' in four directions from each pixel).

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

%% STEP 6 Smooth the signal a bit more 
% This step maybe redundant when working for denoised data
disp('Smoothing data...');
IM_Zframetime_smoothed=smoothdata(IM_Zframetime_conv,3,'gaussian',smooth/2);

% IM_Zframetime_smoothed = (IM_Zframetime_smoothed - min(IM_Zframetime_smoothed(:))) * (255 / (max(IM_Zframetime_smoothed(:)) - min(IM_Zframetime_smoothed(:))));
% IM_Zframetime_smoothed = single(IM_Zframetime_smoothed)/255*(2^16-1);

%% STEP 7  Detect the oxygen sink events/pockets in the filtered, smoothed data             
fprintf('Detecting putative oxygen sinks in every frame... \n'); 
tic;
%The putative hypoxic pockets are identified based on singal intensity.
%To avoid weird results at the borders (because of motion correction) I clip a 20 pixel frame around the field of view.
IMclip=IM_Zframetime_smoothed(Pixel_frame+1:end-Pixel_frame,Pixel_frame+1:end-Pixel_frame,:);
%I also clip the recording are so the linear indices to match when I refine
RecAreaFilterclip=RecAreaFilter(Pixel_frame+1:end-Pixel_frame,Pixel_frame+1:end-Pixel_frame);
%this is a cell vector with legth equal to the duration of the recording that has the info for the 
OxygenSinksInfo_all=cell(size(IMclip,3),1);
for i=1:size(IMclip,3)%for every frame of the recording

    %store the frame in temporary variable and invert the values using the function imcomplement since I am
    %looking for oxygen sinks
    temp=imcomplement(IMclip(:,:,i));
    %apply filter to take only pixels that are above the choosen percentile
    %(this is a bit arbitrary)
    BW=temp>prctile(temp(:),PercentileDetectionThres);  
    
    %Define the different putative pockets in the frame
    SinksInfo=regionprops(BW,'Area','PixelIdxList','Circularity');
       
 %if a pocket is outside the predefined size limits or is            
    % irregular shape              
    %eliminate it from the logical array because it is most likely noise
        
    BW(cat(1, SinksInfo(squeeze(cell2mat({SinksInfo.Area}))<ThresholdMinsize).PixelIdxList))=false;        
    BW(cat(1, SinksInfo(squeeze(cell2mat({SinksInfo.Area}))>ThresholdMaxsize).PixelIdxList))=false;     
    BW(cat(1, SinksInfo(squeeze(cell2mat({SinksInfo.Circularity}))<CircularityThres).PixelIdxList))=false;
    
    
    for j=1:length(SinksInfo)  
        %redundant
%         % removing pixels that are outside the defined recording area   
%         SinksInfo(j).PixelIdxList= SinksInfo(j).PixelIdxList(ismember(SinksInfo(j).PixelIdxList, find(RecAreaFilter)));

        % if more than 50% of the sink does not belong in the recording
        % area, remove it from the pool
        if length(SinksInfo(j).PixelIdxList(ismember(SinksInfo(j).PixelIdxList, find(~RecAreaFilterclip))))/length(SinksInfo(j).PixelIdxList)>0.5
            
            BW(SinksInfo(j).PixelIdxList)=false;

        end
     end   
             
    %update the regionprops after the refinement and save them info 
    OxygenSinksInfo_all{i}=regionprops(BW,'Area','PixelIdxList');

end

clear i j BW temp SinksInfo IMclip RecAreaFilterclip
toc;
%% STEP 8  Tracking oxygen sinks across time
% For the remaining oxygen sinks of each frame track them accross time by checking if they overlap with oxygen sinks from subsequent frames.
fprintf('Tracking putative oxygen sinks across imaging session... \n'); 
tic;
%this will contain putative unique pockets(rows) and the pixels belonging to events for every frame of the recording.
Overall_OxygenSinks_Pxllist=cell(1,length(OxygenSinksInfo_all));
sink_overlap_thres=0.6; %percentage of overlap that is needed between sinks in two different frames
counter=1;
for i=1:length(OxygenSinksInfo_all)-ThresholMinddur %for every video frame      
    %Notice that I do not go all the way to the end. That is because by definition any new putative events starting late in recording will not have enough time to...  
    %reach the minimum duration limit.           
    for ij=1:length(OxygenSinksInfo_all{i,1}) %for every pocket of the frame (aka pocket-frame). As you can see I am starting from the smaller and go to the larger to avoid large pockets hogging the subsequent pockets
        if ~isnan(OxygenSinksInfo_all{i,1}(ij).PixelIdxList) %if this pocket has not been removed from the pool
        
            %create a temporary variable containing the pixels (linear
            %indexing) of this oxygen sink frame
            pixl_temp=OxygenSinksInfo_all{i,1}(ij).PixelIdxList;

            Overall_OxygenSinks_Pxllist{counter,i}=OxygenSinksInfo_all{i,1}(ij).PixelIdxList; %start an entry on that identified pocket-frame
            %and then...
            for f=i+1:length(OxygenSinksInfo_all)%...sweep forward to all the subsequent frames
                        
                for fj=1:length(OxygenSinksInfo_all{f,1}) % and check the different oxygen sinks in those frames. Again going from smaller to larger
            
                    if ~isnan(OxygenSinksInfo_all{f,1}(fj).PixelIdxList) %if a subsequent oxygen sink frame has not been removed from the pool 
                        
                        %get the number of pixels belonging to the start pocket that correspond to the overlap thershold
                        F_overlap=length(pixl_temp)*sink_overlap_thres; 
                        %get the number of pixels belonging to the current pocket that correspond to the overlap thershold
                        NextF_overlap=length(OxygenSinksInfo_all{f,1}(fj).PixelIdxList)*sink_overlap_thres; 
                        
                        if length(intersect(pixl_temp,OxygenSinksInfo_all{f,1}(fj).PixelIdxList))>NextF_overlap || length(intersect(pixl_temp,OxygenSinksInfo_all{f,1}(fj).PixelIdxList))>F_overlap    
                          %and if there is an significant overlap between oxygen sink frames of interest and subsequent oxygen sink frame... 

                          %add the pixel list of that pocket-frame in the list
                          % I use vercat in the eventuality the preceding
                          % sink frame(s) overlaps with more than one sink
                          % area (fj) in the current frame (f)
                          Overall_OxygenSinks_Pxllist{counter,f}=vertcat(Overall_OxygenSinks_Pxllist{counter,f},OxygenSinksInfo_all{f,1}(fj).PixelIdxList);                                   
                                                                            
                          %Now I get the unique pixels of frames that this oxygen sink exists and check the overlap with oxygen sinks in subsequent frames                                                     
                          %The reason I do that is that an oxygen sink in                                                      
                          %frames n and n+3 might not have overlap above                                                       
                          %the minimum but n+2 and n+3 may have. There is a                              
                          %drift in the centroid of the pocket.                    
                          pixl_temp=unique(vertcat(pixl_temp,OxygenSinksInfo_all{f,1}(fj).PixelIdxList));
  
                          OxygenSinksInfo_all{f,1}(fj).PixelIdxList=NaN; %removing that pocket-frame from the pool

                        end
                    end
                    
                end
            end
  
            OxygenSinksInfo_all{i,1}(ij).PixelIdxList=NaN; %remove that pocket-frame from the pool before going to the next
            
            counter=counter+1; %count up for my indexing of the identified pockets
   
        end
    end      
end

clear i ij f fj counter OxygenSinksInfo_all pixl_temp
toc;
%% STEP 9 Refine the identified oxygen sinks based on event duration thresholds

fprintf('Refining putative oxygen sinks based on event duration thresholds... \n'); 
tic;
%this creates a logical array with TRUE value at every frame whne a pocket was detected
Overall_OxySinks_logical=~(cellfun(@isempty, Overall_OxygenSinks_Pxllist));
for i=1:size(Overall_OxySinks_logical,1) %for every putative pocket
    
    %finding the events 
    eventstemp=regionprops(Overall_OxySinks_logical(i,:),'Area','PixelIdxList');    
    
    for q=1:length(eventstemp)       
        
        if eventstemp(q).Area<ThresholMinddur || eventstemp(q).Area>ThresholMaxddur %if the duration is too short or too long, 
            
            %erase the putative event from the logical vector
            Overall_OxySinks_logical(i,eventstemp(q).PixelIdxList)=false;            
            %erase the putative event from the pixel list
            Overall_OxygenSinks_Pxllist(i,eventstemp(q).PixelIdxList)={[]};
            
        end  
    end
    
    %finding the different events after the refinement
    eventstemp=regionprops(Overall_OxySinks_logical(i,:),'Area','PixelIdxList');    
    %just a check if there is only one event there is no reason to go to next step
    if length(eventstemp)==1    
        continue
    end
    
    pairs=nchoosek(1:length(eventstemp),2);    
    %Based on observations pockets are not very close to each other. Therefore, I need to refine.
    for q=1:size(pairs,1) %I go through the events                     
    
        %reality check in case the putative events have been removed in a previous iteration
        if Overall_OxySinks_logical(i,eventstemp(pairs(q,1)).PixelIdxList(1))==true && Overall_OxySinks_logical(i,eventstemp(pairs(q,2)).PixelIdxList(1))==true
 
        if eventstemp(pairs(q,2)).PixelIdxList(1)-eventstemp(pairs(q,1)).PixelIdxList(end)< Thesholdtime %if two events are too close to each other I need to exclude one                                            
            % a few conditionals that will allow me to choose what to exclude.
            % the idea here is that when one of the two events are close to the intuitively imposed duration limits that is the most likely false signal
            % this is of course not a perfect approach but it is a start
            if abs(eventstemp(pairs(q,2)).Area-ThresholMinddur)<abs(eventstemp(pairs(q,1)).Area-ThresholMinddur)... 
                    &&   abs(eventstemp(pairs(q,2)).Area-ThresholMinddur)<abs(eventstemp(pairs(q,1)).Area-ThresholMaxddur)                                             
                %if the following event is too close to the min allowed duration (compared to difference between the leading and min or max)
                Overall_OxySinks_logical(i,eventstemp(pairs(q,2)).PixelIdxList)=false; %exclude it 
                Overall_OxygenSinks_Pxllist(i,eventstemp(pairs(q,2)).PixelIdxList)={[]};

            elseif abs(eventstemp(pairs(q,2)).Area-ThresholMaxddur)<abs(eventstemp(pairs(q,1)).Area-ThresholMinddur)... 
                    &&   abs(eventstemp(pairs(q,2)).Area-ThresholMaxddur)<abs(eventstemp(pairs(q,1)).Area-ThresholMaxddur)
                %if the following event is too close to the max allowed duration (compared to difference between the leading and min or max)
                Overall_OxySinks_logical(i,eventstemp(pairs(q,2)).PixelIdxList)=false; %exclude it  
                Overall_OxygenSinks_Pxllist(i,eventstemp(pairs(q,2)).PixelIdxList)={[]};
                
            elseif abs(eventstemp(pairs(q,1)).Area-ThresholMinddur)<abs(eventstemp(pairs(q,2)).Area-ThresholMinddur)... 
                    &&   abs(eventstemp(pairs(q,1)).Area-ThresholMinddur)<abs(eventstemp(pairs(q,2)).Area-ThresholMaxddur)   
                %if the leading event is too close to the min allowed duration (compared to difference between the following and min or max)
                Overall_OxySinks_logical(i,eventstemp(pairs(q,1)).PixelIdxList)=false; %exclude it              
                Overall_OxygenSinks_Pxllist(i,eventstemp(pairs(q,1)).PixelIdxList)={[]};
            
            elseif abs(eventstemp(pairs(q,1)).Area-ThresholMaxddur)<abs(eventstemp(pairs(q,2)).Area-ThresholMinddur)... 
                    &&   abs(eventstemp(pairs(q,1)).Area-ThresholMaxddur)<abs(eventstemp(pairs(q,2)).Area-ThresholMaxddur)
                %if the leading event is too close to the max allowed duration (compared to difference between the following and min or max)
                Overall_OxySinks_logical(i,eventstemp(pairs(q,1)).PixelIdxList)=false; %exclude it            
                Overall_OxygenSinks_Pxllist(i,eventstemp(pairs(q,1)).PixelIdxList)={[]};
            
            elseif eventstemp(pairs(q,1)).Area > eventstemp(pairs(q,2)).Area  
                % if none of the previous cases is true then go by the size and exclude the smaller event
                Overall_OxySinks_logical(i,eventstemp(pairs(q,2)).PixelIdxList)=false; %exclude the following
                Overall_OxygenSinks_Pxllist(i,eventstemp(pairs(q,2)).PixelIdxList)={[]};                
            else                 
                Overall_OxySinks_logical(i,eventstemp(pairs(q,1)).PixelIdxList)=false; %exclude the leading
                Overall_OxygenSinks_Pxllist(i,eventstemp(pairs(q,1)).PixelIdxList)={[]};
            end      
        end   
        
        end
    end 
 
end

%Finding which putative pockets have any events after the refinement
Anypoc=any(Overall_OxySinks_logical,2);
%Excluding the cell array rows(pockets) that did not have any "real" events
Overall_OxygenSinks_Pxllist=Overall_OxygenSinks_Pxllist(Anypoc,:);
Overall_OxySinks_logical=Overall_OxySinks_logical(Anypoc,:);


clear i q eventstemp pairs Anypoc
toc;

%% STEP 10 Extracting the mean trace for each identified pocket
fprintf('Extracting mean traces from Oxygen Sinks... \n');
% I have to go pocket by pocket

OxySink_Map=false(size(IM_Zframetime_smoothed(:,:,1)));
OxySink_Pxls_all=cell(size(Overall_OxygenSinks_Pxllist,1),1);
Mean_OxySink_TraceZ=NaN(size(Overall_OxygenSinks_Pxllist,1),size(IM_Zframetime,3));
Mean_OxySink_Trace_Convo=NaN(size(Overall_OxygenSinks_Pxllist,1),size(IM_Zframetime_smoothed,3));
for i=1:size(Overall_OxygenSinks_Pxllist,1) %go through the pockets
    
    %Get all the frames that are part of en event for that pocket
    Pxls=Overall_OxygenSinks_Pxllist(i,Overall_OxySinks_logical(i,:));

    %getting the unique pixels (linear index) across all frames of the same sink
    Pxls_Unique=cellfun(@transpose,Pxls,'UniformOutput',false); 
    Pxls_Unique = [Pxls_Unique{:}];
    Pxls_Unique = (unique(Pxls_Unique))';

    %getting the subscript index for these oxygen sink pixels
    %remember that the sinks have been identified in the clipped image so I
    %need to reduce the size of OxySink_Map to the correct size
    [idx_temp(:,1), idx_temp(:,2)]=ind2sub(size(OxySink_Map)-2*Pixel_frame,Pxls_Unique);
    
    %Correcting back to the original image dimensions and storing 
    idx_temp=idx_temp+Pixel_frame;
    
    linIndx=sub2ind(size(OxySink_Map),idx_temp(:,1), idx_temp(:,2)); %change back to linear index
    
    OxySink_Pxls_all{i}=linIndx;
    
    OxySink_Map(linIndx)=true;

    %Now get the singal for all 1200 frames from the IM_Raw_Zframetime
    
    %create a temporary vector as long as the video to store the sum of positive pixels  
    sum_poc=zeros(size(IM_Zframetime,3),1);
    sum_poc_convo=single(zeros(size(IM_Zframetime_smoothed,3),1));

    for k=1:size(idx_temp,1) %going through the pixels of each putative pocket
                    
        sum_poc=sum_poc + squeeze(IM_Zframetime(idx_temp(k,1), idx_temp(k,2),:)); %summing the trace for the pixels
        sum_poc_convo= sum_poc_convo + squeeze(single(IM_Zframetime_smoothed(idx_temp(k,1), idx_temp(k,2),:)));  
    end
    
    sum_poc=sum_poc/size(idx_temp,1); %dividing by the number of pixels in a given pocket to get the mean trace
    sum_poc_convo=sum_poc_convo/size(idx_temp,1);

    %here I can add the check point about the correlation with other
    %pockets and if these have signiticant over and lump together. 
    %I need to be careful with the rest of the script because I will end up
    % conflicting with previsous number of pockets. The script needs to be
    % addapted appropriatelly. 

    Mean_OxySink_TraceZ(i,:)=sum_poc;
    Mean_OxySink_Trace_Convo(i,:)=sum_poc_convo;


    clear idx_temp
end
clear i k Pxls Pxls_Unique idx_temp linIndx sum_poc sum_poc_convo

%% Detrending again the individual traces
% When getting the average pocket singal I end up combining pixels with
% different trends which have not been completely removed when detrending
% pixel by pixel the data matrix

Mean_OxySink_TraceZ=detrend_custom(Mean_OxySink_TraceZ,2);
Mean_OxySink_Trace_Convo=detrend_custom(Mean_OxySink_Trace_Convo,2);

%% Here I calculate whether the trace of a oxygen sink pocket is correlated with other identified pockets.
% I first combine traces with good correlation and spatial overlap creating
% a new average trace.
% This is the first time!!
fprintf('Refining putative oxygen sinks based on the correlation of traces with other sinks... \n'); 
tic;

%flip the traces to use the crosscorrelation function in matlab
FlippedTracestemp=Mean_OxySink_Trace_Convo';
%getting the pearsons correlation of all the different oxygen sink traces
Trace_Convo_CrossCor=corr(FlippedTracestemp);
% for each oxygen sink trace I am calculating the number of oxygen sink traces that they correlate
% well will (>0.8). The minimum is 1 (the autocorrelation)
Trace_Convo_CrossCor_High_sum=sum(Trace_Convo_CrossCor>0.8,2);
% If the total number of traces a give oxygen sink trace correlates well
% with more than the 95th percentile it is likely that it is noise. 
Trace_PotentialNoise=~(Trace_Convo_CrossCor_High_sum>prctile(Trace_Convo_CrossCor_High_sum,90));
%
% Initiate the cell array that will take the pixels of the currated Sinks
OxySink_Pxls_Currated = cell(size(Overall_OxygenSinks_Pxllist));
% Now go through the traces to combine if needed
for i=1:size(Trace_Convo_CrossCor,1)

    if sum(Trace_Convo_CrossCor(i,:)>0.8)>1 && Trace_PotentialNoise(i) % if there is good correlation with more than the autocorrelation and it is not potentially noise
       
        %getting the indices of the well correlated sink traces (columns in
        %the cross correlation matrix) 
        CombinedIndx=find(Trace_Convo_CrossCor(i,:)>0.8)';
               
        % get the actual correlation values and sort them from larger (=1
        % autocorrelation) to smallest  
        Corrs=Trace_Convo_CrossCor(i,Trace_Convo_CrossCor(i,:)>0.8);
        [~,IndxSort]=sort(Corrs,'descend');
        
        % Delete the line from the autocorrelation matrix corresponding to the largest correlation (autocorrelation) 
        Trace_Convo_CrossCor(CombinedIndx(IndxSort(1)),:)=nan;
        % Remember that IndxSort contains the sorted indices of traces that correlate well with trace i.
        % These are effectivelly number of columns in the Trace_Convo_CrossCor correlation matrix. 
        % Given that I am going trace by trace, IndxSort(1) IS i
        % (autocorrelation)

        % Adding the pixels for each frame of the trace in question
        OxySink_Pxls_Currated(i,:)=Overall_OxygenSinks_Pxllist(CombinedIndx(IndxSort(1)),:);

        %empty the cells with the pixels from the seed area so they will
        %not be added in other sinks later on.
        Overall_OxygenSinks_Pxllist(CombinedIndx(IndxSort(1)),:)=({[]});
        %getting the unique pixels of the seed sink across all pixels of the recording              
        Pxls_seed=OxySink_Pxls_Currated(CombinedIndx(IndxSort(1)),Overall_OxySinks_logical(IndxSort(1),:));                               
        Pxls_seed=cellfun(@transpose,Pxls_seed,'UniformOutput',false);                 
        Pxls_seed = [Pxls_seed{:}];           
        Pxls_seed = (unique(Pxls_seed))';

        for j=2:length(IndxSort) %going through the different traces from sink area with the second highest trace correlation 

            %getting the unique pixels for the sink trace in question
            %(across all frames)
            Pxls_test=Overall_OxygenSinks_Pxllist(CombinedIndx(IndxSort(j)),Overall_OxySinks_logical(CombinedIndx(IndxSort(j)),:));                                      
            Pxls_test=cellfun(@transpose,Pxls_test,'UniformOutput',false);                         
            Pxls_test = [Pxls_test{:}];                 
            Pxls_test = (unique(Pxls_test))';
                
            %if the overlap is larger that 70% of (practically) the smaller
            %sink
            if length(intersect(Pxls_seed,Pxls_test))>length(Pxls_test)*0.7 ...
                    || length(intersect(Pxls_seed,Pxls_test))>length(Pxls_seed)*0.7
                 % go frame by frame combine pixels between the seed and
                 % the sink in question. Store the unique ones
                 % if more sink traces correlate well with the seed then
                 % their unique pixels for every pixel will be added
                for k=1:size(OxySink_Pxls_Currated(i,:),2)
                            
                    OxySink_Pxls_Currated(i,k)={unique(vertcat(OxySink_Pxls_Currated{i,k},Overall_OxygenSinks_Pxllist{CombinedIndx(IndxSort(j)),k}))}; %getting the unique pixels for that frame              
                                                     
                    % Delete the line of correlations for the overlapping sink area so we                                                          
                    % will not repeat the process the opposite way                
                    Trace_Convo_CrossCor(CombinedIndx(IndxSort(j)),:)=nan;
                    % Delete the pixels from the list so they will not be
                    % added to other things
                    Overall_OxygenSinks_Pxllist(CombinedIndx(IndxSort(j)),:)=({[]});

                end
                
            end
        end
   
    elseif sum(Trace_Convo_CrossCor(i,:)>0.8)==1 %if the trace correlates well only with it self

        OxySink_Pxls_Currated(i,:)=Overall_OxygenSinks_Pxllist(i,:);      
        Overall_OxygenSinks_Pxllist(i,:)=({[]});

    end
end

%Finding which putative pockets have any events (not empty) after this
%refinement step
Anypoc=any(~(cellfun(@isempty,OxySink_Pxls_Currated)),2);
OxySink_Pxls_Currated=OxySink_Pxls_Currated(Anypoc,:);
% Update the cell array that is going to be used in the 
Overall_OxygenSinks_Pxllist=OxySink_Pxls_Currated;
%this updates the logical array with TRUE value at every frame whne a pocket was detected
Overall_OxySinks_logical=~(cellfun(@isempty, Overall_OxygenSinks_Pxllist));

toc;
clear FlippedTracestemp GoodCorrTraces SeedTrace SeedTracePxls Corrs CombinedPxls IndxSort GoodCorrIndx SeedIndx OxySink_Pxls_Currated

%% Repeat the extraction of mean sink traces for the refined sink areas
% I have to go pocket by pocket

OxySink_Map=false(size(IM_Zframetime_smoothed(:,:,1)));
OxySink_Pxls_all=cell(size(Overall_OxygenSinks_Pxllist,1),1);
Mean_OxySink_TraceZ=NaN(size(Overall_OxygenSinks_Pxllist,1),size(IM_Zframetime,3));
Mean_OxySink_Trace_Convo=NaN(size(Overall_OxygenSinks_Pxllist,1),size(IM_Zframetime_smoothed,3));
for i=1:size(Overall_OxygenSinks_Pxllist,1) %go through the pockets
    
    %Get all the frames that are part of en event for that pocket
    Pxls=Overall_OxygenSinks_Pxllist(i,Overall_OxySinks_logical(i,:));

    %getting the unique pixels (linear index) across all frames of the same sink
    Pxls_Unique=cellfun(@transpose,Pxls,'UniformOutput',false); 
    Pxls_Unique = [Pxls_Unique{:}];
    Pxls_Unique = (unique(Pxls_Unique))';

    %getting the subscript index for these oxygen sink pixels
    %remember that the sinks have been identified in the clipped image so I
    %need to reduce the size of OxySink_Map to the correct size
    [idx_temp(:,1), idx_temp(:,2)]=ind2sub(size(OxySink_Map)-2*Pixel_frame,Pxls_Unique);
    
    %Correcting back to the original image dimensions and storing 
    idx_temp=idx_temp+Pixel_frame;
    
    linIndx=sub2ind(size(OxySink_Map),idx_temp(:,1), idx_temp(:,2)); %change back to linear index
    
    OxySink_Pxls_all{i}=linIndx;
    
    OxySink_Map(linIndx)=true;

    %Now get the singal for all 1200 frames from the IM_Raw_Zframetime
    
    %create a temporary vector as long as the video to store the sum of positive pixels  
    sum_poc=zeros(size(IM_Zframetime,3),1);
    sum_poc_convo=zeros(size(IM_Zframetime_smoothed,3),1);

    for k=1:size(idx_temp,1) %going through the pixels of each putative pocket
                    
        sum_poc=sum_poc + squeeze(IM_Zframetime(idx_temp(k,1), idx_temp(k,2),:)); %summing the trace for the pixels
        sum_poc_convo= sum_poc_convo + squeeze(single(IM_Zframetime_smoothed(idx_temp(k,1), idx_temp(k,2),:)));  
    end
    
    sum_poc=sum_poc/size(idx_temp,1); %dividing by the number of pixels in a given pocket to get the mean trace
    sum_poc_convo=sum_poc_convo/size(idx_temp,1);

    %here I can add the check point about the correlation with other
    %pockets and if these have signiticant over and lump together. 
    %I need to be careful with the rest of the script because I will end up
    % conflicting with previsous number of pockets. The script needs to be
    % addapted appropriatelly. 

    Mean_OxySink_TraceZ(i,:)=sum_poc;
    Mean_OxySink_Trace_Convo(i,:)=sum_poc_convo;


    clear idx_temp
end
clear i k Pxls Pxls_Unique idx_temp linIndx sum_poc sum_poc_convo

%% Detrending again the individual traces
% When getting the average pocket singal I end up combining pixels with
% different trends which have not been completely removed when detrending
% pixel by pixel the data matrix

Mean_OxySink_TraceZ=detrend_custom(Mean_OxySink_TraceZ,2);
Mean_OxySink_Trace_Convo=detrend_custom(Mean_OxySink_Trace_Convo,2);

%% Here I calculate whether the trace of a oxygen sink pocket is correlated with other identified pockets.
% I first combine traces with good correlation and spatial overlap creating
% a new average trace.
% This is the sencond time!!
fprintf('Refining putative oxygen sinks based on the correlation of traces with other sinks for a second time... \n'); 
tic;

%flip the traces to use the crosscorrelation function in matlab
FlippedTracestemp=Mean_OxySink_Trace_Convo';
%getting the pearsons correlation of all the different oxygen sink traces
Trace_Convo_CrossCor=corr(FlippedTracestemp);
% for each oxygen sink trace I am calculating the number of oxygen sink traces that they correlate
% well will (>0.8). The minimum is 1 (the autocorrelation)
Trace_Convo_CrossCor_High_sum=sum(Trace_Convo_CrossCor>0.8,2);
% If the total number of traces a give oxygen sink trace correlates well
% with more than the 95th percentile it is likely that it is noise. 
Trace_PotentialNoise=~(Trace_Convo_CrossCor_High_sum>prctile(Trace_Convo_CrossCor_High_sum,90));
%
% Initiate the cell array that will take the pixels of the currated Sinks
OxySink_Pxls_Currated = cell(size(Overall_OxygenSinks_Pxllist));
% Now go through the traces to combine if needed
for i=1:size(Trace_Convo_CrossCor,1)

    if sum(Trace_Convo_CrossCor(i,:)>0.8)>1 && Trace_PotentialNoise(i) % if there is good correlation with more than the autocorrelation and it is not potentially noise
       
        %getting the indices of the well correlated sink traces (columns in
        %the cross correlation matrix) 
        CombinedIndx=find(Trace_Convo_CrossCor(i,:)>0.8)';
               
        % get the actual correlation values and sort them from larger (=1
        % autocorrelation) to smallest  
        Corrs=Trace_Convo_CrossCor(i,Trace_Convo_CrossCor(i,:)>0.8);
        [~,IndxSort]=sort(Corrs,'descend');
        
        % Delete the line from the autocorrelation matrix corresponding to the largest correlation (autocorrelation) 
        Trace_Convo_CrossCor(CombinedIndx(IndxSort(1)),:)=nan;
        % Remember that IndxSort contains the sorted indices of traces that correlate well with trace i.
        % These are effectivelly number of columns in the Trace_Convo_CrossCor correlation matrix. 
        % Given that I am going trace by trace, IndxSort(1) IS i
        % (autocorrelation)

        % Adding the pixels for each frame of the trace in question
        OxySink_Pxls_Currated(i,:)=Overall_OxygenSinks_Pxllist(CombinedIndx(IndxSort(1)),:);

        %empty the cells with the pixels from the seed area so they will
        %not be added in other sinks later on.
        Overall_OxygenSinks_Pxllist(CombinedIndx(IndxSort(1)),:)=({[]});
        %getting the unique pixels of the seed sink across all pixels of the recording              
        Pxls_seed=OxySink_Pxls_Currated(CombinedIndx(IndxSort(1)),Overall_OxySinks_logical(IndxSort(1),:));                               
        Pxls_seed=cellfun(@transpose,Pxls_seed,'UniformOutput',false);                 
        Pxls_seed = [Pxls_seed{:}];           
        Pxls_seed = (unique(Pxls_seed))';

        for j=2:length(IndxSort) %going through the different traces from sink area with the second highest trace correlation 

            %getting the unique pixels for the sink trace in question
            %(across all frames)
            Pxls_test=Overall_OxygenSinks_Pxllist(CombinedIndx(IndxSort(j)),Overall_OxySinks_logical(CombinedIndx(IndxSort(j)),:));                                      
            Pxls_test=cellfun(@transpose,Pxls_test,'UniformOutput',false);                         
            Pxls_test = [Pxls_test{:}];                 
            Pxls_test = (unique(Pxls_test))';
                
            %if the overlap is larger that 70% of (practically) the smaller
            %sink
            if length(intersect(Pxls_seed,Pxls_test))>length(Pxls_test)*0.7 ...
                    || length(intersect(Pxls_seed,Pxls_test))>length(Pxls_seed)*0.7
                 % go frame by frame combine pixels between the seed and
                 % the sink in question. Store the unique ones
                 % if more sink traces correlate well with the seed then
                 % their unique pixels for every pixel will be added
                for k=1:size(OxySink_Pxls_Currated(i,:),2)
                            
                    OxySink_Pxls_Currated(i,k)={unique(vertcat(OxySink_Pxls_Currated{i,k},Overall_OxygenSinks_Pxllist{CombinedIndx(IndxSort(j)),k}))}; %getting the unique pixels for that frame              
                                                     
                    % Delete the line of correlations for the overlapping sink area so we                                                          
                    % will not repeat the process the opposite way                
                    Trace_Convo_CrossCor(CombinedIndx(IndxSort(j)),:)=nan;
                    % Delete the pixels from the list so they will not be
                    % added to other things
                    Overall_OxygenSinks_Pxllist(CombinedIndx(IndxSort(j)),:)=({[]});

                end
                
            end
        end
   
    elseif sum(Trace_Convo_CrossCor(i,:)>0.8)==1 %if the trace correlates well only with it self

        OxySink_Pxls_Currated(i,:)=Overall_OxygenSinks_Pxllist(i,:);      
        Overall_OxygenSinks_Pxllist(i,:)=({[]});

    end
end

%Finding which putative pockets have any events (not empty) after this
%refinement step
Anypoc=any(~(cellfun(@isempty,OxySink_Pxls_Currated)),2);
OxySink_Pxls_Currated=OxySink_Pxls_Currated(Anypoc,:);
% Update the cell array that is going to be used in the 
Overall_OxygenSinks_Pxllist=OxySink_Pxls_Currated;
%this updates the logical array with TRUE value at every frame whne a pocket was detected
Overall_OxySinks_logical=~(cellfun(@isempty, Overall_OxygenSinks_Pxllist));

toc;
clear FlippedTracestemp GoodCorrTraces SeedTrace SeedTracePxls Corrs CombinedPxls IndxSort GoodCorrIndx SeedIndx OxySink_Pxls_Currated

%% Repeat the extraction of mean sink traces for the refined sink areas
% I have to go pocket by pocket

OxySink_Map=false(size(IM_Zframetime_smoothed(:,:,1)));
OxySink_Pxls_all=cell(size(Overall_OxygenSinks_Pxllist,1),1);
Mean_OxySink_TraceZ=NaN(size(Overall_OxygenSinks_Pxllist,1),size(IM_Zframetime,3));
Mean_OxySink_Trace_Convo=NaN(size(Overall_OxygenSinks_Pxllist,1),size(IM_Zframetime_smoothed,3));
for i=1:size(Overall_OxygenSinks_Pxllist,1) %go through the pockets
    
    %Get all the frames that are part of en event for that pocket
    Pxls=Overall_OxygenSinks_Pxllist(i,Overall_OxySinks_logical(i,:));

    %getting the unique pixels (linear index) across all frames of the same sink
    Pxls_Unique=cellfun(@transpose,Pxls,'UniformOutput',false); 
    Pxls_Unique = [Pxls_Unique{:}];
    Pxls_Unique = (unique(Pxls_Unique))';

    %getting the subscript index for these oxygen sink pixels
    %remember that the sinks have been identified in the clipped image so I
    %need to reduce the size of OxySink_Map to the correct size
    [idx_temp(:,1), idx_temp(:,2)]=ind2sub(size(OxySink_Map)-2*Pixel_frame,Pxls_Unique);
    
    %Correcting back to the original image dimensions and storing 
    idx_temp=idx_temp+Pixel_frame;
    
    linIndx=sub2ind(size(OxySink_Map),idx_temp(:,1), idx_temp(:,2)); %change back to linear index
    
    OxySink_Pxls_all{i}=linIndx;
    
    OxySink_Map(linIndx)=true;

    %Now get the singal for all 1200 frames from the IM_Raw_Zframetime
    
    %create a temporary vector as long as the video to store the sum of positive pixels  
    sum_poc=zeros(size(IM_Zframetime,3),1);
    sum_poc_convo=zeros(size(IM_Zframetime_smoothed,3),1);

    for k=1:size(idx_temp,1) %going through the pixels of each putative pocket
                    
        sum_poc=sum_poc + squeeze(IM_Zframetime(idx_temp(k,1), idx_temp(k,2),:)); %summing the trace for the pixels
        sum_poc_convo= sum_poc_convo + squeeze(single(IM_Zframetime_smoothed(idx_temp(k,1), idx_temp(k,2),:)));  
    end
    
    sum_poc=sum_poc/size(idx_temp,1); %dividing by the number of pixels in a given pocket to get the mean trace
    sum_poc_convo=sum_poc_convo/size(idx_temp,1);

    %here I can add the check point about the correlation with other
    %pockets and if these have signiticant over and lump together. 
    %I need to be careful with the rest of the script because I will end up
    % conflicting with previsous number of pockets. The script needs to be
    % addapted appropriatelly. 

    Mean_OxySink_TraceZ(i,:)=sum_poc;
    Mean_OxySink_Trace_Convo(i,:)=sum_poc_convo;


    clear idx_temp
end
clear i k Pxls Pxls_Unique idx_temp linIndx sum_poc sum_poc_convo

%% Detrending again the individual traces
% When getting the average pocket singal I end up combining pixels with
% different trends which have not been completely removed when detrending
% pixel by pixel the data matrix

Mean_OxySink_TraceZ=detrend_custom(Mean_OxySink_TraceZ,2);
Mean_OxySink_Trace_Convo=detrend_custom(Mean_OxySink_Trace_Convo,2);

%% Here I calculate whether the trace of a oxygen sink pocket is correlated with other identified pockets.
% I first combine traces with good correlation and spatial overlap creating
% a new average trace.
% This is the third time!!
fprintf('Refining putative oxygen sinks based on the correlation of traces with other sinks for a third time... \n'); 
tic;

%flip the traces to use the crosscorrelation function in matlab
FlippedTracestemp=Mean_OxySink_Trace_Convo';
%getting the pearsons correlation of all the different oxygen sink traces
Trace_Convo_CrossCor=corr(FlippedTracestemp);
% for each oxygen sink trace I am calculating the number of oxygen sink traces that they correlate
% well will (>0.8). The minimum is 1 (the autocorrelation)
Trace_Convo_CrossCor_High_sum=sum(Trace_Convo_CrossCor>0.8,2);
% If the total number of traces a give oxygen sink trace correlates well
% with more than the 95th percentile it is likely that it is noise. 
Trace_PotentialNoise=~(Trace_Convo_CrossCor_High_sum>prctile(Trace_Convo_CrossCor_High_sum,90));
%
% Initiate the cell array that will take the pixels of the currated Sinks
OxySink_Pxls_Currated = cell(size(Overall_OxygenSinks_Pxllist));
% Now go through the traces to combine if needed
for i=1:size(Trace_Convo_CrossCor,1)

    if sum(Trace_Convo_CrossCor(i,:)>0.8)>1 && Trace_PotentialNoise(i) % if there is good correlation with more than the autocorrelation and it is not potentially noise
       
        %getting the indices of the well correlated sink traces (columns in
        %the cross correlation matrix) 
        CombinedIndx=find(Trace_Convo_CrossCor(i,:)>0.8)';
               
        % get the actual correlation values and sort them from larger (=1
        % autocorrelation) to smallest  
        Corrs=Trace_Convo_CrossCor(i,Trace_Convo_CrossCor(i,:)>0.8);
        [~,IndxSort]=sort(Corrs,'descend');
        
        % Delete the line from the autocorrelation matrix corresponding to the largest correlation (autocorrelation) 
        Trace_Convo_CrossCor(CombinedIndx(IndxSort(1)),:)=nan;
        % Remember that IndxSort contains the sorted indices of traces that correlate well with trace i.
        % These are effectivelly number of columns in the Trace_Convo_CrossCor correlation matrix. 
        % Given that I am going trace by trace, IndxSort(1) IS i
        % (autocorrelation)

        % Adding the pixels for each frame of the trace in question
        OxySink_Pxls_Currated(i,:)=Overall_OxygenSinks_Pxllist(CombinedIndx(IndxSort(1)),:);

        %empty the cells with the pixels from the seed area so they will
        %not be added in other sinks later on.
        Overall_OxygenSinks_Pxllist(CombinedIndx(IndxSort(1)),:)=({[]});
        %getting the unique pixels of the seed sink across all pixels of the recording              
        Pxls_seed=OxySink_Pxls_Currated(CombinedIndx(IndxSort(1)),Overall_OxySinks_logical(IndxSort(1),:));                               
        Pxls_seed=cellfun(@transpose,Pxls_seed,'UniformOutput',false);                 
        Pxls_seed = [Pxls_seed{:}];           
        Pxls_seed = (unique(Pxls_seed))';

        for j=2:length(IndxSort) %going through the different traces from sink area with the second highest trace correlation 

            %getting the unique pixels for the sink trace in question
            %(across all frames)
            Pxls_test=Overall_OxygenSinks_Pxllist(CombinedIndx(IndxSort(j)),Overall_OxySinks_logical(CombinedIndx(IndxSort(j)),:));                                      
            Pxls_test=cellfun(@transpose,Pxls_test,'UniformOutput',false);                         
            Pxls_test = [Pxls_test{:}];                 
            Pxls_test = (unique(Pxls_test))';
                
            %if the overlap is larger that 70% of (practically) the smaller
            %sink
            if length(intersect(Pxls_seed,Pxls_test))>length(Pxls_test)*0.7 ...
                    || length(intersect(Pxls_seed,Pxls_test))>length(Pxls_seed)*0.7
                 % go frame by frame combine pixels between the seed and
                 % the sink in question. Store the unique ones
                 % if more sink traces correlate well with the seed then
                 % their unique pixels for every pixel will be added
                for k=1:size(OxySink_Pxls_Currated(i,:),2)
                            
                    OxySink_Pxls_Currated(i,k)={unique(vertcat(OxySink_Pxls_Currated{i,k},Overall_OxygenSinks_Pxllist{CombinedIndx(IndxSort(j)),k}))}; %getting the unique pixels for that frame              
                                                     
                    % Delete the line of correlations for the overlapping sink area so we                                                          
                    % will not repeat the process the opposite way                
                    Trace_Convo_CrossCor(CombinedIndx(IndxSort(j)),:)=nan;
                    % Delete the pixels from the list so they will not be
                    % added to other things
                    Overall_OxygenSinks_Pxllist(CombinedIndx(IndxSort(j)),:)=({[]});

                end
                
            end
        end
   
    elseif sum(Trace_Convo_CrossCor(i,:)>0.8)==1 %if the trace correlates well only with it self

        OxySink_Pxls_Currated(i,:)=Overall_OxygenSinks_Pxllist(i,:);      
        Overall_OxygenSinks_Pxllist(i,:)=({[]});

    end
end

%Finding which putative pockets have any events (not empty) after this
%refinement step
Anypoc=any(~(cellfun(@isempty,OxySink_Pxls_Currated)),2);
OxySink_Pxls_Currated=OxySink_Pxls_Currated(Anypoc,:);
% Update the cell array that is going to be used in the 
Overall_OxygenSinks_Pxllist=OxySink_Pxls_Currated;
%this updates the logical array with TRUE value at every frame whne a pocket was detected
Overall_OxySinks_logical=~(cellfun(@isempty, Overall_OxygenSinks_Pxllist));
fprintf([num2str(size(Overall_OxygenSinks_Pxllist,1)),' putative hypoxic pockets identified! \n'])
toc;
clear FlippedTracestemp GoodCorrTraces SeedTrace SeedTracePxls Corrs CombinedPxls IndxSort GoodCorrIndx SeedIndx OxySink_Pxls_Currated

%% Repeat the extraction of mean sink traces for the refined sink areas
% I have to go pocket by pocket

OxySink_Map=false(size(IM_Zframetime_smoothed(:,:,1)));
OxySink_Pxls_all=cell(size(Overall_OxygenSinks_Pxllist,1),1);
Mean_OxySink_TraceZ=NaN(size(Overall_OxygenSinks_Pxllist,1),size(IM_Zframetime,3));
Mean_OxySink_Trace_Convo=NaN(size(Overall_OxygenSinks_Pxllist,1),size(IM_Zframetime_smoothed,3));
for i=1:size(Overall_OxygenSinks_Pxllist,1) %go through the pockets
    
    %Get all the frames that are part of en event for that pocket
    Pxls=Overall_OxygenSinks_Pxllist(i,Overall_OxySinks_logical(i,:));

    %getting the unique pixels (linear index) across all frames of the same sink
    Pxls_Unique=cellfun(@transpose,Pxls,'UniformOutput',false); 
    Pxls_Unique = [Pxls_Unique{:}];
    Pxls_Unique = (unique(Pxls_Unique))';

    %getting the subscript index for these oxygen sink pixels
    %remember that the sinks have been identified in the clipped image so I
    %need to reduce the size of OxySink_Map to the correct size
    [idx_temp(:,1), idx_temp(:,2)]=ind2sub(size(OxySink_Map)-2*Pixel_frame,Pxls_Unique);
    
    %Correcting back to the original image dimensions and storing 
    idx_temp=idx_temp+Pixel_frame;
    
    linIndx=sub2ind(size(OxySink_Map),idx_temp(:,1), idx_temp(:,2)); %change back to linear index
    
    OxySink_Pxls_all{i}=linIndx;
    
    OxySink_Map(linIndx)=true;

    %Now get the singal for all 1200 frames from the IM_Raw_Zframetime
    
    %create a temporary vector as long as the video to store the sum of positive pixels  
    sum_poc=zeros(size(IM_Zframetime,3),1);
    sum_poc_convo=zeros(size(IM_Zframetime_smoothed,3),1);

    for k=1:size(idx_temp,1) %going through the pixels of each putative pocket
                    
        sum_poc=sum_poc + squeeze(IM_Zframetime(idx_temp(k,1), idx_temp(k,2),:)); %summing the trace for the pixels
        sum_poc_convo= sum_poc_convo + squeeze(single(IM_Zframetime_smoothed(idx_temp(k,1), idx_temp(k,2),:)));  
    end
    
    sum_poc=sum_poc/size(idx_temp,1); %dividing by the number of pixels in a given pocket to get the mean trace
    sum_poc_convo=sum_poc_convo/size(idx_temp,1);

    %here I can add the check point about the correlation with other
    %pockets and if these have signiticant over and lump together. 
    %I need to be careful with the rest of the script because I will end up
    % conflicting with previsous number of pockets. The script needs to be
    % addapted appropriatelly. 

    Mean_OxySink_TraceZ(i,:)=sum_poc;
    Mean_OxySink_Trace_Convo(i,:)=sum_poc_convo;


    clear idx_temp
end
clear i k Pxls Pxls_Unique idx_temp linIndx sum_poc sum_poc_convo

%% Detrending again the individual traces
% When getting the average pocket singal I end up combining pixels with
% different trends which have not been completely removed when detrending
% pixel by pixel the data matrix

Mean_OxySink_TraceZ=detrend_custom(Mean_OxySink_TraceZ,2);
Mean_OxySink_Trace_Convo=detrend_custom(Mean_OxySink_Trace_Convo,2);

%% Now I correlated the new traces to see which of them correlate with too many other traces, therefore being potential noise

%flip the traces to use the crosscorrelation function in matlab
FlippedTracestemp=Mean_OxySink_Trace_Convo';
%getting the pearsons correlation of all the different oxygen sink traces
Trace_Convo_CrossCor=corr(FlippedTracestemp);
% for each oxygen sink trace I am calculating the number of oxygen sink traces that they correlate
% well will (>0.8). The minimum is 1 (the autocorrelation)
Trace_Convo_CrossCor_High_sum=sum(Trace_Convo_CrossCor>0.8,2);
% If the total number of traces a give oxygen sink trace correlates well
% with is more than the selected percentile (default 95%) it is likely that it is noise. 
Trace_PotentialNoise=~(Trace_Convo_CrossCor_High_sum>prctile(Trace_Convo_CrossCor_High_sum,90));

clear FlippedTracestemp Trace_Convo_CrossCor Trace_Convo_CrossCor_High_sum 
%% STEP 11  Get the properties of sinks and their events in equal size column cells/vectors to export them as a table
% Gathering oxygen sink pocket parameters
tic;
OxySinkArea_all=cell(size(Overall_OxygenSinks_Pxllist));
OxySinkFilledArea_all=cell(size(Overall_OxygenSinks_Pxllist));
OxySinkDiameter_all=cell(size(Overall_OxygenSinks_Pxllist));
OxySinkPerimeter_all=cell(size(Overall_OxygenSinks_Pxllist));
Circularity_all=cell(size(Overall_OxygenSinks_Pxllist));
Centroid_x_all=cell(size(Overall_OxygenSinks_Pxllist));
Centroid_y_all=cell(size(Overall_OxygenSinks_Pxllist));
OxySinkBoundingBox_all=cell(size(Overall_OxygenSinks_Pxllist));

% Create a new empty logical vector to have all pocket-frames (events) that has the original dimensions
IM_OxySinks_BW=false(size(IM_Zframetime));
for j=1:size(Overall_OxygenSinks_Pxllist,2) %for every frame
    IM_BWtemp=false(size(IM_OxySinks_BW(:,:,j)));
    for i=1:size(Overall_OxygenSinks_Pxllist,1) %for every pocket        
        %this is a sham frame that will only contain a single pocket-frame at a time. To be able to calculate the mean size metrics of that pocket-frame
        IM_BWsingle=false(size(IM_OxySinks_BW(:,:,j))); 
        if ~isempty(Overall_OxygenSinks_Pxllist{i,j}) %if that pocket had an event on that frame

            [idx_temp(:,1), idx_temp(:,2)]=ind2sub(size(IM_OxySinks_BW(:,:,1))-2*Pixel_frame,Overall_OxygenSinks_Pxllist{i,j}); %get the subscript index for the pixels of that pocket-frame
            idx_temp=idx_temp+Pixel_frame; %corrrecting the index to the original size
            linIndx=sub2ind(size(IM_OxySinks_BW(:,:,i)),idx_temp(:,1), idx_temp(:,2)); %change back to linear index

            IM_BWsingle(linIndx)=true;
            
            props_temp=regionprops(IM_BWsingle,'Area','FilledArea','Centroid','Circularity','Perimeter','EquivDiameter','BoundingBox');

            %during refinement I may combin sink traces that their unique
            %pixels (across the entire recording) have overlap and go trace
            %correlation. However these may NOT overlap at all times and
            %there may be a split...
            if length(props_temp)> 1
            
                OxySinkArea_all{i,j}= sum([props_temp.Area]);            
                OxySinkFilledArea_all{i,j}= sum([props_temp.FilledArea]);                               
                OxySinkDiameter_all{i,j}= sum([props_temp.EquivDiameter]);           
                OxySinkPerimeter_all{i,j}= sum([props_temp.Perimeter]);          
                Circularity_all{i,j}= mean([props_temp.Circularity]);;          
                c = [props_temp.Centroid];
                x = 1:2:length(c);
                y = 2:2:length(c);                           
                c = round([mean(c(x)) mean(c(y))]);           
                Centroid_x_all{i,j}= c(1);           
                Centroid_y_all{i,j}= c(2); 
                c = [props_temp.BoundingBox];
                x = 1:4:length(c);
                y = 2:4:length(c);
                z = 3:4:length(c);
                t = 4:4:length(c);
                OxySinkBoundingBox_all{i,j}=[median(c(x)) median(c(y)) median(c(z)) median(c(t))];

            else
                OxySinkArea_all{i,j}= props_temp.Area;            
                OxySinkFilledArea_all{i,j}= props_temp.FilledArea;                               
                OxySinkDiameter_all{i,j}= props_temp.EquivDiameter;           
                OxySinkPerimeter_all{i,j}= props_temp.Perimeter;          
                Circularity_all{i,j}= props_temp.Circularity;          
                c = round(props_temp.Centroid);           
                Centroid_x_all{i,j}= c(1);           
                Centroid_y_all{i,j}= c(2);            
                OxySinkBoundingBox_all{i,j}=props_temp.BoundingBox;
            end
            IM_BWtemp(linIndx)=true; %replace the false with true values in the pixels of that pocket-frame
            
            clear idx_temp
        end
        
        
    end
    IM_OxySinks_BW(:,:,j)=IM_BWtemp;
end
clear i j c x y z t IM_BWtemp idx_temp linIndx IM_BWsingle props_temp

% Mean pocket parameters
MeanOxySinkArea_um=NaN(size(Overall_OxygenSinks_Pxllist,1),1);
MeanOxySinkFilledArea_um=NaN(size(Overall_OxygenSinks_Pxllist,1),1);
MeanOxySinkDiameter_um=NaN(size(Overall_OxygenSinks_Pxllist,1),1);
MeanOxySinkPerimeter_um=NaN(size(Overall_OxygenSinks_Pxllist,1),1);
MeanCircularity=NaN(size(Overall_OxygenSinks_Pxllist,1),1);
MeanCentroid_x=NaN(size(Overall_OxygenSinks_Pxllist,1),1);
MeanCentroid_y=NaN(size(Overall_OxygenSinks_Pxllist,1),1);
MeanOxySinkBoundingBox=cell(size(Overall_OxygenSinks_Pxllist,1),1);
for i=1:size(OxySinkArea_all,1)
      
    % multiplying the area in pixel with (um/pixel)^2 value to get size in um^2
    MeanOxySinkArea_um(i)= mean(cell2mat(OxySinkArea_all(i,Overall_OxySinks_logical(i,:))))*PixelSize^2; 
    MeanOxySinkFilledArea_um(i)= mean(cell2mat(OxySinkFilledArea_all(i,Overall_OxySinks_logical(i,:))))*PixelSize^2;
    % multiplying the distances in pixel with um/pixel value to get size in um
    MeanOxySinkDiameter_um(i)=mean(cell2mat(OxySinkDiameter_all(i,Overall_OxySinks_logical(i,:))))*PixelSize;
    MeanOxySinkPerimeter_um(i)=mean(cell2mat(OxySinkPerimeter_all(i,Overall_OxySinks_logical(i,:))))*PixelSize;
    MeanCircularity(i)=mean(cell2mat(Circularity_all(i,Overall_OxySinks_logical(i,:))));
    MeanCentroid_x(i)=mean(cell2mat(Centroid_x_all(i,Overall_OxySinks_logical(i,:))));
    MeanCentroid_y(i)=mean(cell2mat(Centroid_y_all(i,Overall_OxySinks_logical(i,:))));
    MeanOxySinkBoundingBox{i}=mean(vertcat(OxySinkBoundingBox_all{i,Overall_OxySinks_logical(i,:)}),1);
    %making sure the values make sense after averaging
    MeanOxySinkBoundingBox{i}(1)=floor(MeanOxySinkBoundingBox{i}(1)) + floor((MeanOxySinkBoundingBox{i}(1)-floor(MeanOxySinkBoundingBox{i}(1)))/0.5) * 0.5;
    MeanOxySinkBoundingBox{i}(2)=floor(MeanOxySinkBoundingBox{i}(2)) + floor((MeanOxySinkBoundingBox{i}(2)-floor(MeanOxySinkBoundingBox{i}(2)))/0.5) * 0.5;
    MeanOxySinkBoundingBox{i}(3)=round(MeanOxySinkBoundingBox{i}(3));
    MeanOxySinkBoundingBox{i}(4)=round(MeanOxySinkBoundingBox{i}(4));
end


clear OxySinkFilledArea_all OxySinkDiameter_all OxySinkPerimeter_all Circularity_all Centroid_x_all Centroid_y_all OxySinkBoundingBox_all


%% 
%Now the Overall_Pocket_Pxllist has all the putative oxygen sinks (row)and the list of pixels for each frame an event took place
fprintf('Collating properties of oxygen sinks and events... \n')

Experiment=cell(size(Overall_OxygenSinks_Pxllist,1),1); Experiment(:)={DatafileID};
Mouse=cell(size(Overall_OxygenSinks_Pxllist,1),1);
if exist('Mous','var')
    Mouse(:)={Mous};%that is coming from the wrapper
else
    Mouse(:)={'Unknown'};
    fprintf('Script was call individually. No info from wrapper available! Mouse ID is unknown \n');
end
Condition=cell(size(Overall_OxygenSinks_Pxllist,1),1); 
if exist('Cond','var')
    Condition(:)={Cond};%that is coming from the wrapper
else 
    fprintf('Script was call individually. No info from wrapper available! Condition is unknown \n');
    Condition(:)={'Unknown'};
end

DrugID=cell(size(Overall_OxygenSinks_Pxllist,1),1);  
if exist('Drug','var')
    DrugID(:)={Drug}; %that is coming from the wrapper
else  
    fprintf('Script was call individually. No info from wrapper available! DrugID is unknown \n');
    DrugID(:)={'Unknown'};
end

Genotype=cell(size(Overall_OxygenSinks_Pxllist,1),1); 
if exist('Gen','var')
    Genotype(:)={Gen}; %that is coming from the wrapper
else  
    fprintf('Script was call individually. No info from wrapper available! Genotype is unknown \n');
    Genotype(:)={'Unknown'};
end
Promoter=cell(size(Overall_OxygenSinks_Pxllist,1),1); 
if exist('Promo','var')
    Promoter(:)={Promo}; %that is coming from the wrapper
else   
    fprintf('Script was call individually. No info from wrapper available! Genetic promoter is unknown \n');
    Promoter(:)={'Unknown'};
end
RecDuration=cell(size(Overall_OxygenSinks_Pxllist,1),1); RecDuration(:)={RecDur};
RecAreaSize=cell(size(Overall_OxygenSinks_Pxllist,1),1); RecAreaSize(:)={RecArea}; 
%% Collecting metrics for the events

%Event parameters
NumOxySinkEvents=NaN(size(Overall_OxygenSinks_Pxllist,1),1);
Start=cell(size(Overall_OxygenSinks_Pxllist,1),1);
Duration=cell(size(Overall_OxygenSinks_Pxllist,1),1);
NormOxySinkAmp=cell(size(Overall_OxygenSinks_Pxllist,1),1);
Size_modulation=cell(size(Overall_OxygenSinks_Pxllist,1),1);

for i=1:size(Overall_OxygenSinks_Pxllist,1) %for every oxygen sink area
    

    temptrace=Mean_OxySink_Trace_Convo(i,:);    
    [p,~,mu] = polyfit(1:numel(temptrace), temptrace, 7);        
    tracetrend=polyval(p,1:numel(temptrace),[],mu);

    %get the details of events for that sink area
    eventstemp=regionprops(Overall_OxySinks_logical(i,:),'Area','PixelIdxList');            
    NumOxySinkEvents(i)=length(eventstemp);
    
    Start(i)={NaN(1,length(eventstemp))};
    Duration(i)={NaN(1,length(eventstemp))};
    NormOxySinkAmp(i)={NaN(1,length(eventstemp))};
    Size_modulation(i)={NaN(1,length(eventstemp))};

    for q=1:length(eventstemp) %go through the events
        

        if eventstemp(q).PixelIdxList(1)>1 && eventstemp(q).PixelIdxList(end)<size(temptrace,2) %if the event did not happen at the very beggining or end
        
            StartIndx=find((abs(tracetrend(1:eventstemp(q).PixelIdxList(1))-temptrace(1:eventstemp(q).PixelIdxList(1)))<0.015),1,'last');
            if isempty(StartIndx) %just in case is too close to the start and there are not enough elements to find one that is close enough to the base line.          
                StartIndx=eventstemp(q).PixelIdxList(1);
            end           
            Start{i}(q)=StartIndx;

            EndIndx=eventstemp(q).PixelIdxList(end)+find((abs(tracetrend(eventstemp(q).PixelIdxList(end):end)-temptrace(eventstemp(q).PixelIdxList(end):end))<0.015),1,'first')-1;            
            if isempty(EndIndx) %just in case is too close to the end and there are not enough elements to find one that is close enough to the base line.          
                EndIndx=eventstemp(q).PixelIdxList(end);
            end
              
            Duration{i}(q)=EndIndx-Start{i}(q);            
            % Getting the normalised amplitude of each event as the mmaximum distance from the trendline during the duration of the event.     
            NormOxySinkAmp{i}(q) = abs(abs(min(temptrace(Start{i}(q):EndIndx)))-abs(min(tracetrend(Start{i}(q):EndIndx))));

        elseif eventstemp(q).PixelIdxList(1)==1 %if it is at the start
            
            Start{i}(q)=eventstemp(q).PixelIdxList(1);

            EndIndx=eventstemp(q).PixelIdxList(end)+find((abs(tracetrend(eventstemp(q).PixelIdxList(end):end)-temptrace(eventstemp(q).PixelIdxList(end):end))<0.015),1,'first')-1;
            if isempty(EndIndx)        
                EndIndx=eventstemp(q).PixelIdxList(end);
            end

            Duration{i}(q)=EndIndx-Start{i}(q);            
            % Getting the normalised amplitude of each event as the mmaximum distance from the trendline during the duration of the event.     
            NormOxySinkAmp{i}(q) = abs(abs(min(temptrace(Start{i}(q):EndIndx)))-abs(min(tracetrend(Start{i}(q):EndIndx))));

        else

            StartIndx=find((abs(tracetrend(1:eventstemp(q).PixelIdxList(1))-temptrace(1:eventstemp(q).PixelIdxList(1)))<0.015),1,'last');
            if isempty(StartIndx)        
                StartIndx=eventstemp(q).PixelIdxList(1);
            end           
            Start{i}(q)=StartIndx;
           
            EndIndx=eventstemp(q).PixelIdxList(end);
            
            Duration{i}(q)=EndIndx-Start{i}(q);            
            % Getting the normalised amplitude of each event as the mmaximum distance from the trendline during the duration of the event.     
            NormOxySinkAmp{i}(q) = abs(abs(min(temptrace(Start{i}(q):EndIndx)))-abs(min(tracetrend(Start{i}(q):EndIndx))));

        end
        
        
        %Calculating an index (size of pocket end minuse the end of event)/size at the start of the event
        Size_modulation{i}(q)= (length(Overall_OxygenSinks_Pxllist{i,eventstemp(q).PixelIdxList(end)})-length(Overall_OxygenSinks_Pxllist{i,eventstemp(q).PixelIdxList(1)}))...
            /length(Overall_OxygenSinks_Pxllist{i,eventstemp(q).PixelIdxList(1)});
        
    end

    if ~(any(NormOxySinkAmp{i}>2.5)) %if that pocket did not have events that were stronger than 3 times the baseline is probably noise
    
        Trace_PotentialNoise(i)=false;
    end

end

clear i q eventstemp temptrace p mu tracetrend StartIndx EndIndx

toc;
%% STEP 12 Detecting Oxygen surges
fprintf('Identifying oxygen surges now... \n');

%for the oxygen surges I need to follow a different approach. 
% the spatial distribution is less importnant 
% there is no maximum size
% the percentile threshold will be different. I might need to check the original script from Ryszard
tic;
%The putative oxygen surges are identified based on singal intensity.

%this is a cell vector with legth equal to the duration of the recording that has the info for the surges 
OxygenSurgesInfo_all=cell(size(IM_Zframetime_smoothed,3),1);
for i=1:size(IM_Zframetime_smoothed,3)%for every frame of the recording

    %store the frame in temporary variable 
    temp=mat2gray(IM_Zframetime_smoothed(:,:,i));
    
    %apply filter to take only pixels that are above the choosen percentile (this is a hard coded threshold that I might need to revisit)
    BW=temp>prctile(temp(:),90);  
    
    %Define the different putative pockets in the frame
    SurgesInfo=regionprops(BW,'Area','PixelIdxList','Circularity');
    
    %if a pocket is outside the predefined size limits 
    %eliminate it from the logical array because it is most likely noise
    BW(cat(1, SurgesInfo(squeeze(cell2mat({SurgesInfo.Area}))<ThresholdMinsize_Surges).PixelIdxList))=false;        
    BW(cat(1, SurgesInfo(squeeze(cell2mat({SurgesInfo.Circularity}))<0.1).PixelIdxList))=false;


    % removing surges with more than 50% pixels outside the defined recording area
    for j=1:length(SurgesInfo)

        if length(SurgesInfo(j).PixelIdxList(ismember(SurgesInfo(j).PixelIdxList, find(~RecAreaFilter))))/length(SurgesInfo(j).PixelIdxList)>0.5
            
            BW(SurgesInfo(j).PixelIdxList)=false;

        end

    end
    
    %update the regionprops after the refinement and save them info 
    OxygenSurgesInfo_all{i}=regionprops(BW,'Area','PixelIdxList');

end

clear i j BW temp SurgesInfo 
toc;

%% STEP 13 Tracking surges through frames
fprintf('Tracking putative oxygen surges across the imaging session... \n'); 
tic;
%this will contain putative oxygen surges(rows) and the pixels belonging to events for every frame of the recording.
Overall_OxygenSurges_Pxllist=cell(1,length(OxygenSurgesInfo_all));
counter=1;
for i=1:length(OxygenSurgesInfo_all)-ThresholMinddur_Surges %for every video frame (the minimum duration does not really apply to the surges)      

    for j=1:length(OxygenSurgesInfo_all{i,1}) %for every oxygen surge of the frame (aka signal peak in the frame)
        if ~isnan(OxygenSurgesInfo_all{i,1}(j).PixelIdxList) %if this peak has not been removed from the pool already
        
            Overall_OxygenSurges_Pxllist{counter,i}=OxygenSurgesInfo_all{i,1}(j).PixelIdxList; %start an entry on that identified surge-frame
            %and then...
            for f=i+1:length(OxygenSurgesInfo_all)%...sweep forward to all the subsequent frames
                    
                for fj=1:length(OxygenSurgesInfo_all{f,1}) % and check the different surge pockets in each of those frames
            
                    if ~isnan(OxygenSurgesInfo_all{f,1}(fj).PixelIdxList) %if this surge-frame has not been removed from the pool

                        if length(intersect(OxygenSurgesInfo_all{i,1}(j).PixelIdxList,OxygenSurgesInfo_all{f,1}(fj).PixelIdxList))>ThresholdMinsize_Surges-smooth        
                          
                            %and there is an significant overlap between surge-frame of interest and subsequent surge-frame...                          
                            %The overlap threshold is the minimum size of a
                            %pocket minus the convolution radius

                            Overall_OxygenSurges_Pxllist{counter,f}=OxygenSurgesInfo_all{f,1}(fj).PixelIdxList; %add the pixel list of that surge-frame in the list
                            OxygenSurgesInfo_all{f,1}(fj).PixelIdxList=NaN; %removing that pocket-frame from the pool
                            break 

                        end
                    end  
                end

            end
  
            OxygenSurgesInfo_all{i,1}(j).PixelIdxList=NaN; %remove that pocket-frame from the pool before going to the next
            
            counter=counter+1; %count up for my indexing of the identified pockets
   
        end
    end      
end

clear i j f fj counter OxygenSurgesInfo_all
toc;

%% STEP 14 Refining the identified oxygen surges based on event duration thresholds
fprintf('Refining putative oxygen surges based on event duration threshold... \n'); 
tic;
%this creates a logical array with TRUE value at every frame whne a pocket was detected
Overall_OxySurges_logical=~(cellfun(@isempty, Overall_OxygenSurges_Pxllist));
for i=1:size(Overall_OxySurges_logical,1) %for every putative pocket
    
    %finding the events 
    eventstemp=regionprops(Overall_OxySurges_logical(i,:),'Area','PixelIdxList');    
    for q=1:length(eventstemp)       
        
        if eventstemp(q).Area<ThresholMinddur_Surges %if the duration is too short 
            
            %erase the putative event from the logical vector
            Overall_OxySurges_logical(i,eventstemp(q).PixelIdxList)=false;            
            %erase the putative event from the pixel list
            Overall_OxygenSurges_Pxllist(i,eventstemp(q).PixelIdxList)={[]};
            
        end  
    end
    
end

%Finding which putative pockets have any events after the refinement
Anypoc=any(Overall_OxySurges_logical,2);
%Excluding the cell array rows(pockets) that did not have any "real" events
Overall_OxygenSurges_Pxllist=Overall_OxygenSurges_Pxllist(Anypoc,:);
Overall_OxySurges_logical=Overall_OxySurges_logical(Anypoc,:);
fprintf([num2str(size(Overall_OxygenSurges_Pxllist,1)),' putative oxygen surge loci identified! \n'])

clear i q eventstemp pairs Anypoc
toc;

%% Gathering oxygen surge loci parameters

OxySurgeArea_all=cell(size(Overall_OxygenSurges_Pxllist));
OxySurgeFilledArea_all=cell(size(Overall_OxygenSurges_Pxllist));
OxySurgeDiameter_all=cell(size(Overall_OxygenSurges_Pxllist));
OxySurgePerimeter_all=cell(size(Overall_OxygenSurges_Pxllist));
Circularity_Surge_all=cell(size(Overall_OxygenSurges_Pxllist));
Centroid_Surge_x_all=cell(size(Overall_OxygenSurges_Pxllist));
Centroid_Surge_y_all=cell(size(Overall_OxygenSurges_Pxllist));
OxySurgeBoundingBox_all=cell(size(Overall_OxygenSurges_Pxllist));

% Create a new empty logical vector to have all pocket-frames (events) that has the original dimensions
IM_OxySurges_BW=false(size(IM_Zframetime));
for j=1:size(Overall_OxygenSurges_Pxllist,2) %for every frame
    IM_BWtemp=false(size(IM_OxySurges_BW(:,:,j)));
    for i=1:size(Overall_OxygenSurges_Pxllist,1) %for every surge loci        
        %this is a sham frame that will only contain a single surge-frame at a time. To be able to calculate the mean size metrics of that pocket-frame
        IM_BWsingle=false(size(IM_OxySurges_BW(:,:,j))); 
        if ~isempty(Overall_OxygenSurges_Pxllist{i,j}) %if that pocket had an event on that frame


            IM_BWsingle(Overall_OxygenSurges_Pxllist{i,j})=true;
            %of course props_temp will have length=1 since IM_BWsingle contains a single pocket
            props_temp=regionprops(IM_BWsingle,'Area','FilledArea','Centroid','Circularity','Perimeter','EquivDiameter','BoundingBox');

            OxySurgeArea_all{i,j}= props_temp.Area;
            OxySurgeFilledArea_all{i,j}= props_temp.FilledArea;                    
            OxySurgeDiameter_all{i,j}= props_temp.EquivDiameter;
            OxySurgePerimeter_all{i,j}= props_temp.Perimeter;
            Circularity_Surge_all{i,j}= props_temp.Circularity;
            c = round(props_temp.Centroid);
            Centroid_Surge_x_all{i,j}= c(1);
            Centroid_Surge_y_all{i,j}= c(2);
            OxySurgeBoundingBox_all{i,j}=props_temp.BoundingBox;
            
            IM_BWtemp(Overall_OxygenSurges_Pxllist{i,j})=true; %replace the false with true values in the pixels of that pocket-frame
            
            clear idx_temp
        end
        
        
    end
    IM_OxySurges_BW(:,:,j)=IM_BWtemp;
end
clear i j c IM_BWtemp idx_temp linIndx IM_BWsingle props_temp

% Mean pocket parameters
MeanOxySurgeArea_um=NaN(size(Overall_OxygenSurges_Pxllist,1),1);
MeanOxySurgeFilledArea_um=NaN(size(Overall_OxygenSurges_Pxllist,1),1);
MeanOxySurgeDiameter_um=NaN(size(Overall_OxygenSurges_Pxllist,1),1);
MeanOxySurgePerimeter_um=NaN(size(Overall_OxygenSurges_Pxllist,1),1);
MeanCircularity_Surge=NaN(size(Overall_OxygenSurges_Pxllist,1),1);
MeanCentroid_Surge_x=NaN(size(Overall_OxygenSurges_Pxllist,1),1);
MeanCentroid_Surge_y=NaN(size(Overall_OxygenSurges_Pxllist,1),1);
MeanOxySurgeBoundingBox=cell(size(Overall_OxygenSurges_Pxllist,1),1);
for i=1:size(OxySurgeArea_all,1)
      
    % multiplying the area in pixel with (um/pixel)^2 value to get size in um^2
    MeanOxySurgeArea_um(i)= mean(cell2mat(OxySurgeArea_all(i,Overall_OxySurges_logical(i,:))))*PixelSize^2; 
    MeanOxySurgeFilledArea_um(i)= mean(cell2mat(OxySurgeFilledArea_all(i,Overall_OxySurges_logical(i,:))))*PixelSize^2;
    % multiplying the distances in pixel with um/pixel value to get size in um
    MeanOxySurgeDiameter_um(i)=mean(cell2mat(OxySurgeDiameter_all(i,Overall_OxySurges_logical(i,:))))*PixelSize;
    MeanOxySurgePerimeter_um(i)=mean(cell2mat(OxySurgePerimeter_all(i,Overall_OxySurges_logical(i,:))))*PixelSize;
    MeanCircularity_Surge(i)=mean(cell2mat(Circularity_Surge_all(i,Overall_OxySurges_logical(i,:))));
    MeanCentroid_Surge_x(i)=mean(cell2mat(Centroid_Surge_x_all(i,Overall_OxySurges_logical(i,:))));
    MeanCentroid_Surge_y(i)=mean(cell2mat(Centroid_Surge_y_all(i,Overall_OxySurges_logical(i,:))));
    MeanOxySurgeBoundingBox{i}=mean(vertcat(OxySurgeBoundingBox_all{1,Overall_OxySurges_logical(1,:)}),1);
    %making sure the values make sense after averaging
    MeanOxySurgeBoundingBox{i}(1)=floor(MeanOxySurgeBoundingBox{i}(1)) + floor((MeanOxySurgeBoundingBox{i}(1)-floor(MeanOxySurgeBoundingBox{i}(1)))/0.5) * 0.5;
    MeanOxySurgeBoundingBox{i}(2)=floor(MeanOxySurgeBoundingBox{i}(2)) + floor((MeanOxySurgeBoundingBox{i}(2)-floor(MeanOxySurgeBoundingBox{i}(2)))/0.5) * 0.5;
    MeanOxySurgeBoundingBox{i}(3)=round(MeanOxySurgeBoundingBox{i}(3));
    MeanOxySurgeBoundingBox{i}(4)=round(MeanOxySurgeBoundingBox{i}(4));
end


clear OxySurgeFilledArea_all OxySurgeDiameter_all OxySurgePerimeter_all Circularity_Surge_all Centroid_Surge_x_all Centroid_Surge_y_all OxySurgeBoundingBox_all

%defined the metadata vectors based on the length of oxygen surge events
Experiment_Surge=cell(size(Overall_OxygenSurges_Pxllist,1),1); Experiment_Surge(:)={DatafileID};
Mouse_Surge=cell(size(Overall_OxygenSurges_Pxllist,1),1);
if exist('Mous','var')
    Mouse_Surge(:)={Mous};%that is coming from the wrapper
else
    Mouse_Surge(:)={'Unknown'};
    fprintf('Script was call individually. No info from wrapper available! Mouse ID is unknown \n');
end
Condition_Surge=cell(size(Overall_OxygenSurges_Pxllist,1),1); 
if exist('Cond','var')
    Condition_Surge(:)={Cond};%that is coming from the wrapper
else 
    fprintf('Script was call individually. No info from wrapper available! Condition is unknown \n');
    Condition_Surge(:)={'Unknown'};
end

DrugID_Surge=cell(size(Overall_OxygenSinks_Pxllist,1),1);  
if exist('Drug','var')
    DrugID_Surge(:)={Drug}; %that is coming from the wrapper
else  
    fprintf('Script was call individually. No info from wrapper available! DrugID is unknown \n');
    DrugID_Surge(:)={'Unknown'};
end

Genotype_Surge=cell(size(Overall_OxygenSurges_Pxllist,1),1); 
if exist('Gen','var')
    Genotype_Surge(:)={Gen}; %that is coming from the wrapper
else  
    fprintf('Script was call individually. No info from wrapper available! Genotype is unknown \n');
    Genotype_Surge(:)={'Unknown'};
end
Promoter_Surge=cell(size(Overall_OxygenSurges_Pxllist,1),1); 
if exist('Promo','var')
    Promoter_Surge(:)={Promo}; %that is coming from the wrapper
else   
    fprintf('Script was call individually. No info from wrapper available! Genetic promoter is unknown \n');
    Promoter_Surge(:)={'Unknown'};
end
RecDuration_Surge=cell(size(Overall_OxygenSurges_Pxllist,1),1); RecDuration_Surge(:)={RecDur};
RecAreaSize_Surge=cell(size(Overall_OxygenSurges_Pxllist,1),1); RecAreaSize_Surge(:)={RecArea}; 

%% STEP 14

%Event parameters
fprintf('Extracting mean traces from Oxygen Surges... \n');
% I have to go pocket by pocket

OxySurge_Map=false(size(IM_Zframetime_smoothed(:,:,1)));
OxySurge_Pxls_all=cell(size(Overall_OxygenSurges_Pxllist,1),1);
Mean_OxySurge_TraceZ=NaN(size(Overall_OxygenSurges_Pxllist,1),size(IM_Zframetime,3));
for i=1:size(Overall_OxygenSurges_Pxllist,1) %go through the pockets
    
    %Get all the frames with that are part of en event for that pocket
    Pxls=Overall_OxygenSurges_Pxllist(i,Overall_OxySurges_logical(i,:));

    Pxls_Unique=unique(Pxls{1}); %getting the unique pixels for the first frame. This is to initiate that variable.
    for j=2:length(Pxls)
        
        Pxls_Unique=unique(vertcat(Pxls_Unique,Pxls{j})); %getting the unique with all other event frames of the same pocket
    end
    OxySurge_Pxls_all{i}=Pxls_Unique;
    OxySurge_Map(Pxls_Unique)=true;
    
    %getting the subscript index for the unique pixels
    [idx_temp(:,1), idx_temp(:,2)]=ind2sub(size(OxySurge_Map),Pxls_Unique);

    %Now get the singal for all 1200 frames from the IM_Zframetime
    
    %create a temporary vector as long as the video to store the sum of positive pixels  
    sum_poc=zeros(size(IM_Zframetime,3),1);

    for k=1:size(idx_temp,1) %going through the pixels of each putative pocket
                    
        sum_poc=sum_poc + squeeze(IM_Zframetime(idx_temp(k,1), idx_temp(k,2),:)); %summing the trace for the pixels
            
    end
    
    Mean_OxySurge_TraceZ(i,:)=sum_poc/size(idx_temp,1); %dividing by the number of pixels in a given pocket to get the mean trace
    
    clear idx_temp
end
clear i k Pxls Pxls_Unique idx_temp sum_poc

%% Calculating SurgeEvent parameters
NumOxySurgeEvents=NaN(size(Overall_OxygenSurges_Pxllist,1),1);
Start_Surge=cell(size(Overall_OxygenSurges_Pxllist,1),1);
Duration_Surge=cell(size(Overall_OxygenSurges_Pxllist,1),1);
NormOxySurgeAmp=cell(size(Overall_OxygenSurges_Pxllist,1),1);
Size_Surge_modulation=cell(size(Overall_OxygenSurges_Pxllist,1),1);

for i=1:size(Overall_OxygenSurges_Pxllist,1)
    

    eventstemp=regionprops(Overall_OxySurges_logical(i,:),'Area','PixelIdxList');        
    NumOxySurgeEvents(i)=length(eventstemp);
    Start_Surge(i)={NaN(1,length(eventstemp))};
    Duration_Surge(i)={NaN(1,length(eventstemp))};
    NormOxySurgeAmp(i)={NaN(1,length(eventstemp))};
    Size_Surge_modulation(i)={NaN(1,length(eventstemp))};

    for q=1:length(eventstemp)
        
        Start_Surge{i}(q)=eventstemp(q).PixelIdxList(1);
        Duration_Surge{i}(q)=eventstemp(q).Area;
        
        %getting the normalised amplitude of each event as the mean surge loci signal during the duration of the event divided by...
        %the mean surge loci signal during the minimum inter-event period before(or after) the event         
        if eventstemp(q).PixelIdxList(1)-20<1 && eventstemp(q).PixelIdxList(end)+20<size(Mean_OxySurge_TraceZ(i,:),2)+1 %checking if the event happened at the very beggining and finished not too close to the end
            %if it did take the 20 frames after the event as a baseline
            NormOxySurgeAmp{i}(q)=abs(mean(Mean_OxySurge_TraceZ(i,eventstemp(q).PixelIdxList))/mean(Mean_OxySurge_TraceZ(i,eventstemp(q).PixelIdxList(end)+1:eventstemp(q).PixelIdxList(end)+20)));
        elseif eventstemp(q).PixelIdxList(end)+20>size(Mean_OxySurge_TraceZ(i,:),2) && eventstemp(q).PixelIdxList(1)-20>0 % check if it happend too close in the end and started not too close to the start 
            %if it did take the 20 frames before the start as a baseline
            NormOxySurgeAmp{i}(q)=abs(mean(Mean_OxySurge_TraceZ(i,eventstemp(q).PixelIdxList))/mean(Mean_OxySurge_TraceZ(i,eventstemp(q).PixelIdxList(1)-20:eventstemp(q).PixelIdxList(1)-1)));
        elseif eventstemp(q).PixelIdxList(1)-20>0 && eventstemp(q).PixelIdxList(end)+20<size(Mean_OxySurge_TraceZ(i,:),2)+1 %check if the event was neith close to begining nor the end of the recording
            %if it did take the 20 frames before the start as a baseline
            NormOxySurgeAmp{i}(q)=abs(mean(Mean_OxySurge_TraceZ(i,eventstemp(q).PixelIdxList))/mean(Mean_OxySurge_TraceZ(i,eventstemp(q).PixelIdxList(1)-20:eventstemp(q).PixelIdxList(1)-1)));
        else % if the recording is too short there is a posibility that the surge event can be close to the start and the end 

            %check from which side I have more space to take baseline 
            if eventstemp(q).PixelIdxList(1)-20 > size(Mean_OxySurge_TraceZ(i,:),2)-(eventstemp(q).PixelIdxList(end)+20) % if there is more space at the start
            
                %take baseline from the start of recording to the start of
                %the surge
                NormOxySurgeAmp{i}(q)=abs(mean(Mean_OxySurge_TraceZ(i,eventstemp(q).PixelIdxList))/mean(Mean_OxySurge_TraceZ(i,1:eventstemp(q).PixelIdxList(1)-1)));

            else % if not

                %take baseline from the end of the surge to the end of recording
                NormOxySurgeAmp{i}(q)=abs(mean(Mean_OxySurge_TraceZ(i,eventstemp(q).PixelIdxList))/mean(Mean_OxySurge_TraceZ(i,eventstemp(q).PixelIdxList(end)+1:end)));
            end

        end
        %Calculating an index (size of pocket end minuse the end of event)/size at the start of the event
        Size_Surge_modulation{i}(q)= (length(Overall_OxygenSurges_Pxllist{i,eventstemp(q).PixelIdxList(end)})-length(Overall_OxygenSurges_Pxllist{i,eventstemp(q).PixelIdxList(1)}))...
            /length(Overall_OxygenSurges_Pxllist{i,eventstemp(q).PixelIdxList(1)});
        
    end
 
end

clear i q eventstemp


%% extracting traces for spatial bins
Mean_ROI_TraceZ=NaN(size(RecArea_bins,1),size(IM_Zframetime,3));
for i=1:size(RecArea_bins,1) %go through the different ROIs. Remember that RecArea_bins has the index of the upper left corner of the rectangular ROIs  

    %take the average signal for that ROI for each frame
    Mean_ROI_TraceZ(i,:)=squeeze(mean(mean(IM_Zframetime(RecArea_bins(i,1):RecArea_bins(i,1)+fix(25/PixelSize)*2-1,RecArea_bins(i,2):RecArea_bins(i,2)+fix(25/PixelSize)*2-1,:),1),2));

end
 clear i

%% The output tables 
Table_OxygenSinks_Out = table(Experiment,Mouse,Condition,DrugID,Genotype,Promoter,RecDuration,RecAreaSize,...
    MeanCentroid_x,MeanCentroid_y,MeanOxySinkArea_um,MeanOxySinkFilledArea_um,MeanOxySinkDiameter_um,MeanOxySinkPerimeter_um,MeanCircularity,...
    MeanOxySinkBoundingBox,NumOxySinkEvents,Start,Duration, NormOxySinkAmp,Size_modulation,OxySink_Pxls_all);

Table_OxygenSurges_Out = table(Experiment_Surge,Mouse_Surge,Condition_Surge,Genotype_Surge,Promoter_Surge,RecDuration_Surge,RecAreaSize_Surge,...
    MeanCentroid_Surge_x,MeanCentroid_Surge_y,MeanOxySurgeArea_um,MeanOxySurgeFilledArea_um,MeanOxySurgeDiameter_um,MeanOxySurgePerimeter_um,MeanCircularity_Surge,...
    MeanOxySurgeBoundingBox,NumOxySurgeEvents,Start_Surge,Duration_Surge,NormOxySurgeAmp,Size_Surge_modulation,OxySurge_Pxls_all);
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

%export the detrended matrix
IM_Notrend = uint8((IM_Notrend - min(IM_Notrend(:))) * (255 / (max(IM_Notrend(:)) - min(IM_Notrend(:)))));
IM_Notrend = uint16(single(IM_Notrend)/255*(2^16-1));
saveastiff(uint16(single(IM_Notrend)/255*(2^16-1)), [Tifffiles(1).folder,'\',Images_Proce_folder,'\', 'IM_Notrend', DatafileID, '.tif'],options);

%export the conv/smoothed matrix
IM_Zframetime_smoothed = uint8((IM_Zframetime_smoothed - min(IM_Zframetime_smoothed(:))) * (255 / (max(IM_Zframetime_smoothed(:)) - min(IM_Zframetime_smoothed(:)))));
saveastiff(uint16(single(IM_Zframetime_smoothed)/255*(2^16-1)), [Tifffiles(1).folder,'\',Images_Proce_folder,'\','IM_Conv', DatafileID, '.tif'],options);

%export the OxySinks BW matrix
saveastiff(uint16(IM_OxySinks_BW), [Tifffiles(1).folder,'\',OxySinksfolder,'\', 'IM_OxySinks_BW', DatafileID, '.tif'],options);

%export the OxySurges BW matrix
saveastiff(uint16(IM_OxySurges_BW), [Tifffiles(1).folder,'\',OxySurgesfolder,'\', 'IM_OxySurges_BW', DatafileID, '.tif'],options);
     
%% Save .AVI video to access in the manual curation app
%this takes waaaay too long and the benefits are not that big.
%Also the loading on the manual curation app is not great and I need to
%work on the callback functions to be more meaninful.
% disp('Saving detrended matrix as an AVI');
% 
% v = VideoWriter([Tifffiles(1).folder,'\',ManualCurOxySinksfolder,'\','IM_Notrend','_video_',DatafileID,'.avi']);
% open(v);
% set(gca,'visible','off'); %hide the current axes
% set(get(gca,'children'),'visible','off'); %hide the current axes contents
% for i=1:size(IM_Notrend,3)
%     f = figure('visible', 'off');    
%     imshow(IM_Notrend(:,:,i));
%     frame = getframe(gca);
%     writeVideo(v,frame);    
% end
% close(f);
% close(v);
% 
% clear v f i frame
% Check these related to the manual curation app and how to open the video in that
% https://se.mathworks.com/matlabcentral/answers/499047-callback-functions-for-displaying-videos 
% https://se.mathworks.com/matlabcentral/answers/558596-displaying-video-file-on-app-designer-axes
%% Saving variables

% Saving table as a .mat file to help with the curation
save([Tifffiles(1).folder,'\',OxySinksfolder,'\','OxygenSinks_Urefined', DatafileID, '.mat'], 'Table_OxygenSinks_Out','OxySinkArea_all','Mean_ROI_TraceZ','Mean_OxySink_TraceZ','Mean_OxySink_Trace_Convo','OxySink_Map','Trace_PotentialNoise');

% Saving it also as a .mat file to help with the curation
save([Tifffiles(1).folder,'\',OxySurgesfolder,'\','OxygenSurges', DatafileID, '.mat'], 'Table_OxygenSurges_Out','OxySurgeArea_all','Mean_ROI_TraceZ','Mean_OxySurge_TraceZ','OxySurge_Map');

% Saving info for manual curation GUI and further analysis
save([Tifffiles(1).folder,'\',ManualCurOxySinksfolder,'\','ManualCuration', DatafileID, '.mat'], 'Table_OxygenSinks_Out','OxySinkArea_all', 'Mean_OxySink_TraceZ','Mean_OxySink_Trace_Convo','OxySink_Map','Trace_PotentialNoise');

%% FUNCTIONS %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%% Function to fit and remove trend based on the mean field of view timeseries

% Antonis Asiminas

% function IM2=mean_detrend(Miu,IM)
% 
% [p,s,mu] = polyfit(1:numel(Miu), Miu, 7);
%     
% signaltrend=polyval(p,1:numel(Miu),[],mu);
% IM2=nan(size(IM));
% % Final detrend on the original data
% for xx=1:size(IM,1)
%     parfor yy=1:size(IM,2)
%         S1=(squeeze(IM(xx,yy,:)))';
%         IM2(xx,yy,:) = S1- signaltrend;
%     end
% end; clear xx yy S1 p
% 
% end

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


%% Function that can find peaks in noisy data. I use it only at the begging where I try to find pixels with peaks.


function varargout = peakfinder(x0, sel, thresh, extrema, includeEndpoints, interpolate)
%PEAKFINDER Noise tolerant fast peak finding algorithm
%   INPUTS:
%       x0 - A real vector from the maxima will be found (required)
%       sel - The amount above surrounding data for a peak to be,
%           identified (default = (max(x0)-min(x0))/4). Larger values mean
%           the algorithm is more selective in finding peaks.
%       thresh - A threshold value which peaks must be larger than to be
%           maxima or smaller than to be minima.
%       extrema - 1 if maxima are desired, -1 if minima are desired
%           (default = maxima, 1)
%       includeEndpoints - If true the endpoints will be included as
%           possible extrema otherwise they will not be included
%           (default = true)
%       interpolate - If true quadratic interpolation will be performed
%           around each extrema to estimate the magnitude and the
%           position of the peak in terms of fractional indicies. Note that
%           unlike the rest of this function interpolation assumes the
%           input is equally spaced. To recover the x_values of the input
%           rather than the fractional indicies you can do:
%           peakX = x0 + (peakLoc - 1) * dx
%           where x0 is the first x value and dx is the spacing of the
%           vector. Output peakMag to recover interpolated magnitudes.
%           See example 2 for more information.
%           (default = false)
%
%   OUTPUTS:
%       peakLoc - The indicies of the identified peaks in x0
%       peakMag - The magnitude of the identified peaks
%
%   [peakLoc] = peakfinder(x0) returns the indicies of local maxima that
%       are at least 1/4 the range of the data above surrounding data.
%
%   [peakLoc] = peakfinder(x0,sel) returns the indicies of local maxima
%       that are at least sel above surrounding data.
%
%   [peakLoc] = peakfinder(x0,sel,thresh) returns the indicies of local
%       maxima that are at least sel above surrounding data and larger
%       (smaller) than thresh if you are finding maxima (minima).
%
%   [peakLoc] = peakfinder(x0,sel,thresh,extrema) returns the maxima of the
%       data if extrema > 0 and the minima of the data if extrema < 0
%
%   [peakLoc] = peakfinder(x0,sel,thresh,extrema, includeEndpoints)
%       returns the endpoints as possible extrema if includeEndpoints is
%       considered true in a boolean sense
%
%   [peakLoc, peakMag] = peakfinder(x0,sel,thresh,extrema,interpolate)
%       returns the results of results of quadratic interpolate around each
%       extrema if interpolate is considered to be true in a boolean sense
%
%   [peakLoc, peakMag] = peakfinder(x0,...) returns the indicies of the
%       local maxima as well as the magnitudes of those maxima
%
%   If called with no output the identified maxima will be plotted along
%       with the input data.
%
%   Note: If repeated values are found the first is identified as the peak
%
% Example 1:
% t = 0:.0001:10;
% x = 12*sin(10*2*pi*t)-3*sin(.1*2*pi*t)+randn(1,numel(t));
% x(1250:1255) = max(x);
% peakfinder(x)
%
% Example 2:
% ds = 100;  % Downsample factor
% dt = .001; % Time step
% ds_dt = ds*dt; % Time delta after downsampling
% t0 = 1;
% t = t0:dt:5 + t0;
% x = 0.2-sin(0.01*2*pi*t)+3*cos(7/13*2*pi*t+.1)-2*cos((1+pi/10)*2*pi*t+0.2)-0.2*t;
% x(end) = min(x);
% x_ds = x(1:ds:end); % Downsample to test interpolation
% [minLoc, minMag] = peakfinder(x_ds, .8, 0, -1, false, true);
% minT = t0 + (minLoc - 1) * ds_dt; % Take into account 1 based indexing
% p = plot(t,x,'-',t(1:ds:end),x_ds,'o',minT,minMag,'rv');
% set(p(2:end), 'linewidth', 2); % Show the markers more clearly
% legend('Actual Data', 'Input Data', 'Estimated Peaks');
% Copyright Nathanael C. Yoder 2015 (nyoder@gmail.com)
% Perform error checking and set defaults if not passed in
narginchk(1, 6);
nargoutchk(0, 2);
s = size(x0);
flipData =  s(1) < s(2);
len0 = numel(x0);
if len0 ~= s(1) && len0 ~= s(2)
    error('PEAKFINDER:Input','The input data must be a vector')
elseif isempty(x0)
    varargout = {[],[]};
    return;
end
if ~isreal(x0)
    warning('PEAKFINDER:NotReal','Absolute value of data will be used')
    x0 = abs(x0);
end
if nargin < 2 || isempty(sel)
    sel = (max(x0)-min(x0))/4;
elseif ~isnumeric(sel) || ~isreal(sel)
    sel = (max(x0)-min(x0))/4;
    warning('PEAKFINDER:InvalidSel',...
        'The selectivity must be a real scalar.  A selectivity of %.4g will be used',sel)
elseif numel(sel) > 1
    warning('PEAKFINDER:InvalidSel',...
        'The selectivity must be a scalar.  The first selectivity value in the vector will be used.')
    sel = sel(1);
end
if nargin < 3 || isempty(thresh)
    thresh = [];
elseif ~isnumeric(thresh) || ~isreal(thresh)
    thresh = [];
    warning('PEAKFINDER:InvalidThreshold',...
        'The threshold must be a real scalar. No threshold will be used.')
elseif numel(thresh) > 1
    thresh = thresh(1);
    warning('PEAKFINDER:InvalidThreshold',...
        'The threshold must be a scalar.  The first threshold value in the vector will be used.')
end
if nargin < 4 || isempty(extrema)
    extrema = 1;
else
    extrema = sign(extrema(1)); % Should only be 1 or -1 but make sure
    if extrema == 0
        error('PEAKFINDER:ZeroMaxima','Either 1 (for maxima) or -1 (for minima) must be input for extrema');
    end
end
if nargin < 5 || isempty(includeEndpoints)
    includeEndpoints = true;
end
if nargin < 6 || isempty(interpolate)
    interpolate = false;
end
x0 = extrema*x0(:); % Make it so we are finding maxima regardless
thresh = thresh*extrema; % Adjust threshold according to extrema.
dx0 = diff(x0); % Find derivative
dx0(dx0 == 0) = -eps; % This is so we find the first of repeated values
ind = find(dx0(1:end-1).*dx0(2:end) < 0)+1; % Find where the derivative changes sign
% Include endpoints in potential peaks and valleys as desired
if includeEndpoints
    x = [x0(1);x0(ind);x0(end)];
    ind = [1;ind;len0];
    minMag = min(x);
    leftMin = minMag;
else
    x = x0(ind);
    minMag = min(x);
    leftMin = min(x(1), x0(1));
end
% x only has the peaks, valleys, and possibly endpoints
len = numel(x);
if len > 2 % Function with peaks and valleys
    % Set initial parameters for loop
    tempMag = minMag;
    foundPeak = false;
    if includeEndpoints
        % Deal with first point a little differently since tacked it on
        % Calculate the sign of the derivative since we tacked the first
        %  point on it does not neccessarily alternate like the rest.
        signDx = sign(diff(x(1:3)));
        if signDx(1) <= 0 % The first point is larger or equal to the second
            if signDx(1) == signDx(2) % Want alternating signs
                x(2) = [];
                ind(2) = [];
                len = len-1;
            end
        else % First point is smaller than the second
            if signDx(1) == signDx(2) % Want alternating signs
                x(1) = [];
                ind(1) = [];
                len = len-1;
            end
        end
    end
    % Skip the first point if it is smaller so we always start on a
    %   maxima
    if x(1) >= x(2)
        ii = 0;
    else
        ii = 1;
    end
    % Preallocate max number of maxima
    maxPeaks = ceil(len/2);
    peakLoc = zeros(maxPeaks,1);
    peakMag = zeros(maxPeaks,1);
    cInd = 1;
    % Loop through extrema which should be peaks and then valleys
    while ii < len
        ii = ii+1; % This is a peak
        % Reset peak finding if we had a peak and the next peak is bigger
        %   than the last or the left min was small enough to reset.
        if foundPeak
            tempMag = minMag;
            foundPeak = false;
        end
        % Found new peak that was lager than temp mag and selectivity larger
        %   than the minimum to its left.
        if x(ii) > tempMag && x(ii) > leftMin + sel
            tempLoc = ii;
            tempMag = x(ii);
        end
        % Make sure we don't iterate past the length of our vector
        if ii == len
            break; % We assign the last point differently out of the loop
        end
        ii = ii+1; % Move onto the valley
        % Come down at least sel from peak
        if ~foundPeak && tempMag > sel + x(ii)
            foundPeak = true; % We have found a peak
            leftMin = x(ii);
            peakLoc(cInd) = tempLoc; % Add peak to index
            peakMag(cInd) = tempMag;
            cInd = cInd+1;
        elseif x(ii) < leftMin % New left minima
            leftMin = x(ii);
        end
    end
    % Check end point
    if includeEndpoints
        if x(end) > tempMag && x(end) > leftMin + sel
            peakLoc(cInd) = len;
            peakMag(cInd) = x(end);
            cInd = cInd + 1;
        elseif ~foundPeak && tempMag > minMag % Check if we still need to add the last point
            peakLoc(cInd) = tempLoc;
            peakMag(cInd) = tempMag;
            cInd = cInd + 1;
        end
    elseif ~foundPeak
        if x(end) > tempMag && x(end) > leftMin + sel
            peakLoc(cInd) = len;
            peakMag(cInd) = x(end);
            cInd = cInd + 1;
        elseif tempMag > min(x0(end), x(end)) + sel
            peakLoc(cInd) = tempLoc;
            peakMag(cInd) = tempMag;
            cInd = cInd + 1;
        end
    end
    % Create output
    if cInd > 1
        peakInds = ind(peakLoc(1:cInd-1));
        peakMags = peakMag(1:cInd-1);
    else
        peakInds = [];
        peakMags = [];
    end
else % This is a monotone function where an endpoint is the only peak
    [peakMags,xInd] = max(x);
    if includeEndpoints && peakMags > minMag + sel
        peakInds = ind(xInd);
    else
        peakMags = [];
        peakInds = [];
    end
end
% Apply threshold value.  Since always finding maxima it will always be
%   larger than the thresh.
if ~isempty(thresh)
    m = peakMags>thresh;
    peakInds = peakInds(m);
    peakMags = peakMags(m);
end
if interpolate && ~isempty(peakMags)
    middleMask = (peakInds > 1) & (peakInds < len0);
    noEnds = peakInds(middleMask);
    magDiff = x0(noEnds + 1) - x0(noEnds - 1);
    magSum = x0(noEnds - 1) + x0(noEnds + 1)  - 2 * x0(noEnds);
    magRatio = magDiff ./ magSum;
    peakInds(middleMask) = peakInds(middleMask) - magRatio/2;
    peakMags(middleMask) = peakMags(middleMask) - magRatio .* magDiff/8;
end
% Rotate data if needed
if flipData
    peakMags = peakMags.';
    peakInds = peakInds.';
end
% Change sign of data if was finding minima
if extrema < 0
    peakMags = -peakMags;
    x0 = -x0;
end
% Plot if no output desired
if nargout == 0
    if isempty(peakInds)
        disp('No significant peaks found')
    else
        figure;
        plot(1:len0,x0,'.-',peakInds,peakMags,'ro','linewidth',2);
    end
else
    varargout = {peakInds,peakMags};
end

end
%% Function to find neibouring pixels 

function [Iadj , Radj, Nfound ] = neighbourND( index, sizeA, res )
% function  [Iadj , Radj, Nfound] = neighbour3D( index,  sizeA, res )
% Calculate the linear indices for neighboring points in a matrix 
% Second output is and array of distances based on an input resolution vector
% This resolution vector defaults to ones(1,ndims)
% The output Nfound reports the number of neighbours found in within the
% matrix. For 2D we expect up to 8, for 3D up to 26 etc...
% 
% Example 1:
% A is a 128x128x16 image data matrix with a spatial resolution of
% 0.1x 0.25x .375 mm^3 
% to get the neighbouring point linear indices for point 456 we do
% sizeA = [128 128 16]
% [ Iadj , Radj, Nfound] = neighbourND( 456, sizeA, [ .10 .25 .375] )
%
% NEW: now index can be a column array with linear indices
% Output Iadj will be Nx8 (2D) or Nx26 (3D) etc and Radj will be 
% a row array 1x8 or 1x26 etc... 
%
% Example 2:
% create points near the center of a 144x192x16 matrix
% spatial resolution .3 x .3x 5 mm^3
% idx = (-6:1:6)+((144*192*3)+144*96+76)
%[ Iadj , Radj, Nfound] = neighbourND( idx , [144,192, 32] , [.3, 0.3, 5])
% Results in 11x26 matrix Iadj, 
% 26 distances in Radj and Nfound is 26
%
% The neighbour indices outside the matrix will be zero!
% when a single index is entered the outside points are still removed so a
% point in a 3D matrix at the edge can sill return 17 neighbours or even less
% when it is a corner.
%==============================================

%==============================================
% Ronald Ouwerkerk 2010 NIH/NIDDK 
% New version: Now handles arrays of indices
% This script is made available on Matlab file exchange by the author 
% for use by other Matlab programmers.
% This script is not intended for commercial use.
% If used for published work a reference or acknowledgement is greatly 
% appreciated.
% The function was tested for several 1D(col and row), 2D, 3D and 4D cases
% I cannot be sure that it really works for all dimensionalities. 
% Let me know if you find a bug (and feel free to squash it for me)
%==============================================

% Set defaults and process input parameters
% first two are arbitary values for a demo
if nargin <1
    % default index [7,6,2]
    index = 128*128+ 128*5+7
end

if nargin < 2
    % default size 128x128xN with N big enough for the index
    i3 = floor( index /128/128);
    disp( 'Demo mode')
    sizeA =[128, 128, i3+2]
end

% Get dimensionality
ndimA = length( sizeA );

%Set default resolution to isotropic distances
if nargin < 3
    res =ones(1, length( sizeA) );
else
    if length(res) < ndimA;
        errstr = sprintf('\nError in %s.\n The length of the resolution array (%d) must equal the number of matrix dimensions (%d)\n', ...
                            mfilename,length(res)  , ndimA  );
        disp(errstr)
        help( mfilename)
        return
    else
        % reduce the resolution array, last digit is probably slice
        % thickness, irrelevant if we have one slice only
        res = res( 1:ndimA );
    end
end

% explicit version of ind2sub 
% ind2sub requires multiple output arguments, one for each dimension
ilin = index(:);
np = length( ilin );
imat = ones( np, ndimA);

for di = ndimA:-1:2    
    blocksize = prod( sizeA( 1:(di-1)  ) );
    ndi = 1+ floor( ( ilin-1) / blocksize );
    ilin = ilin- (ndi -1) *blocksize;
    imat(:,di) = ndi;
end
imat(:,1) = ilin;

% Find the indices of neighbours
% Get all the index permutations for neighbours ( -1, +1) over all
% dimensions. The total number of neighbours should be three  to the power Ndim
% minus one if we discard the original point itself

% initialize the shift index array
nneighb = 3^ndimA;
nbi = zeros( nneighb, ndimA);

di = ndimA;
while ( di ) 
    N = 3^(di-1);
    ni = 1:N;
    while( ni(end) < nneighb+1 )
        for val=[-1, 0, 1]
              nbi( ni ,di ) = val;
              ni = ni+ N;
        end
    end
    di = di-1;
end

% Create distance matrix
d = ones(nneighb, 1) * res;
d = d.*abs( nbi );
% create a row vector with distances
dvec = sqrt( sum( d.^2, 2))';
% Get index to exclude the original point: distance = 0
notorig = logical( dvec > 0 );

% Add the input index array to nbi to get all neighbours
% set up the array for neighbour indices
nd = length( index);
Iadj = zeros( nd, nneighb );
kdo = notorig(ones(nd,1), : ); 

for di = 1:ndimA
    indices = imat( :, di );
    shifts = nbi( :, di )';
    neighbindices = indices( :, ones( 1,nneighb)) +shifts( ones(nd, 1), : ) ;
    maxmat = sizeA( di );
    % set up mask matrix to keep indices within limits and excllude the original point
    s = logical( neighbindices <= maxmat );
    s =logical( neighbindices > 0 ) & s;
    kdo = kdo & s;
    % Calculate the linear index
    if di == 1       
        Iadj( kdo ) =  neighbindices( kdo );
    else
        blocksize = prod( sizeA( 1:(di-1)  ) );
        m = neighbindices-1;
        Iadj(kdo )  = Iadj(kdo )+ m(kdo)*blocksize;
    end
end

% Select only the sensible points for the neighbour index and distances matrices
% Remove columns that have no valid indices anywhere at all (e.g. origin)
% for shorter index lists with  all points near the edges more may be
% removed.
if nd == 1
    allkdo = any( kdo, 1);
    Iadj = Iadj( :, allkdo);
    Radj = dvec( allkdo );
    Nfound = length(  find( allkdo ) );
else
    Nfound = nneighb-1;
    Radj = dvec;
    iself = (Radj == 0);
    Iadj = Iadj(:,~iself);
    Radj = Radj(~iself);
end

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
    if ndims(data) >= 4, errcode = 2; assert(false); end;
    [height, width, depth] = size(data);
    tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
%     tagstruct.Photometric = Tiff.Photometric.MinIsWhite;
%     tagstruct.Photometric = Tiff.Photometric.Mask;
%     tagstruct.Photometric = Tiff.Photometric.Separated;
else
    if ndims(data) >= 5, errcode = 2; assert(false); end;
    [height, width, cc, depth] = size(data); % cc: color channels. 3: rgb, 4: rgb with alpha channel
    if cc ~= 3 && cc ~= 4, errcode = 3; assert(false); end;
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
                while ~tfile.lastDirectory(); % Append a new image to the last directory of an exiting file
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
    display(sprintf('The file was saved successfully. Elapsed time : %.3f s.', tElapsed));
end

catch exception
%% Exception management
    if exist('tfile', 'var'), tfile.close(); end
    switch errcode
        case 1
            if options.message, error '''data'' is empty.'; end;
        case 2
            if options.message, error 'Data dimension is too large.'; end;
        case 3
            if options.message, error 'Third dimesion (color depth) should be 3 or 4.'; end;
        case 4
            if options.message, error 'Color image cannot have int8, int16 or int32 format.'; end;
        case 5
            if options.message, error 'Unsupported Matlab data type. (char, logical, cell, struct, function_handle, class)'; end;
        case 6
            if options.message, error 'File already exists.'; end;
        case 7
            if options.message, error(['Failed to open the file ''' path '''.']); end;
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

