%% This is little script for the extraction of trace from the area that corresponds to the whiskers of animals being recorded while head restrained
% The Red channel values are extracted for a rectagular area
% USE

% (1) Choose the .avi file you want to process. 
% IMPORTANT
% only .avi RGB video files work for this

% (2) Choose the area you want to extract the trace for
% To do this draw the rectagular area by clicking and draging and then
% right-click and select crop-image

% The script will proceed from that point and may become a bit glitchy
% since it is using a parfor loop to go through all frames and extrace the
% mean value of pixels for that area for the red channel only

% The trace is then exported as a csv file in the same directory where the
% video file is.

% This trace can then be imported into the OxyDynamics_stats analysis script 

% From Antonis Asiminas, PhD
% Copenhagen, Denmark, February 2022,
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

% Choose what file do I want
[Videofile, Videopath]= uigetfile({'*.*','All Files (*.*)'},'File Selector');

%
%creating the object to read the avi file
mov=VideoReader([Videopath,'\',Videofile]);

%%
%using the first frame to get draw the area of interest
img = read(mov,1);
[~,rect] = imcrop(img);

rect=fix(rect);
RTrace=NaN(mov.NumFrames,1);
parfor i=1:mov.NumFrames

    img = read(mov,i); 
    RTrace(i)=mean(imcrop(img(:,:,1),rect),'all'); 

end


%% Smooth and detrend trace 

RTrace_smoothed = smoothdata(RTrace,'gaussian',20);
%% Quantifying variance with the sliding window
% Whisking = more variance
RTrace_var = movvar(RTrace_smoothed,100);

%%
TableOut=table(RTrace_smoothed,RTrace_var);
%% save the traces in a csv
writetable(TableOut,[Videopath,'\',Videofile(1:end-4),'_RTraces.csv']);
