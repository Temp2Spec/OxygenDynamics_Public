% A simple script and a few functions tomake alligned videos and plots for
% presentations
% Data is related to OXygen Dynamics project with Felix RM Beinlich

% Written by Antonis Asiminas, PhD
% June 2023, COpenhagen


%% Inputs

% Duration of recording in seconds
RecDur=1200;
% Acquisition frame rates in Hz
FR_Oxygen=1;
FR_beh=25;

%start and stop times in seconds
start=601;
stop=800;

% How many time to speed up
speedX=3;

%files to open
tracedata='CombDataExample_plotvideos.xlsx';
Oxygetif='ID_18_1_mcor.tif';
Pockettif='IM_OxySinks_BWID18_20201008.tif';
Behvideo='ID18_behavior_0DLC_resnet50_Posture4Sep13shuffle1_258000_filtered_labeled.mp4';


%% Getting the oxygen recording trimmed, sped up (3x) and saved as an mp4

trimTiffStack(Oxygetif, 'OxygenRec_video_601-800_3x.avi', start*FR_Oxygen, stop*FR_Oxygen, FR_Oxygen*speedX)

trimTiffStack(Pockettif, 'Pockets_video_601-800_3x.avi', start*FR_Oxygen, stop*FR_Oxygen, FR_Oxygen*speedX)


%% Getting the trace data from the xl file

dataTable = readtable(tracedata);

Mov=dataTable.Mov;
BLI=dataTable.BLI; BLI=BLI(~isnan(BLI));
HowMany=dataTable.HowMany; HowMany=HowMany(~isnan(HowMany));
MovStart=dataTable.MovStart; MovStart=MovStart(~isnan(MovStart));

%% Generating videos with line plots

generateLinePlotVideo(start:1/FR_beh:stop, Mov(start*FR_beh:stop*FR_beh), start*FR_beh, stop*FR_beh, FR_beh*speedX, 'Movement_601-800_3x.avi','Time(s)','Displacement(px)','Movement')
generateLinePlotVideo(start:1/FR_Oxygen:stop, BLI(start*FR_Oxygen:stop*FR_Oxygen), start*FR_Oxygen, stop*FR_Oxygen, FR_Oxygen*speedX, 'BLI_601-800_3x.avi','Time(s)','ΔBLI/BLI','Oxygen')
generateLinePlotVideo(start:1/FR_Oxygen:stop, HowMany(start*FR_Oxygen:stop*FR_Oxygen), start*FR_Oxygen, stop*FR_Oxygen, FR_Oxygen*speedX, 'HowMany_601-800_3x.avi','Time(s)','Number of Pockets','Number of Pockets')

%% Getting the behavioural video trimmed and sped up

% I am triming frames from the end of the video. This is after discussion
% with Felix
% I am saving the file with 3X the FR to speed it up. BE AWARE!! FR >240Hz
% is not permitted

trimMP4Video(Behvideo, 'Beh_video.avi', start*FR_beh, stop*FR_beh,FR_beh*speedX);

%% Testing
% x=start:1/FR_beh:stop;
% y=Mov(start*FR_beh:stop*FR_beh);
% startIdx=start*FR_beh;
% endIdx=stop*FR_beh;
% frameRate=FR_beh*3;
% outputFilename='Movement_601-800_3x.avi';
% xtitle='Time(s)';
% ytitle='Displacement(px)';
% plottitle='Movement';

%% FUNCTIONS %%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% A function that generates a video in MATLAB, showing a line plot progressively with all the points, 
% based on the given start and end points, and saves it as an MP4 video

% x: Vector of x-coordinates.
% y: Vector of y-coordinates.
% startIdx: Starting index of the points to include in the line plot.
% endIdx: Ending index of the points to include in the line plot.
% frameRate: Desired frame rate for the video.
% outputFilename: Output filename for the generated MP4 video
% xtitle: The x axis title
% ytitle: The y axis title
% plottitle: The plot title

function generateLinePlotVideo(x, y, startIdx, endIdx, frameRate, outputFilename,xtitle, ytitle, plottitle)
    

    % Compute the x-axis and y-axis limits using the calculated minimum and maximum values      
    xMin = min(x);
    xMax = max(x);
    yMin = min(y);
    yMax = max(y); 
    % Create figure and axes
    figure('visible','off');
    axes('NextPlot', 'replacechildren');

    % Initialize video object
    vidObj = VideoWriter(outputFilename);
    vidObj.FrameRate = frameRate;
    open(vidObj);

    % Iterate through the points
    for i = 1:endIdx-startIdx+1
        % Plot current set of points
        plot(x(1:i), y(1:i),'LineWidth',2);
        xlim([xMin, xMax]);
        ylim([yMin, yMax]);
        xlabel(xtitle);
        ylabel(ytitle);
        title(plottitle);

        % Capture frame and write to video
        frame = getframe(gcf);
        writeVideo(vidObj, frame);
    end

    % Finalize video
    close(vidObj);
end


%%  A function that opens a TIFF stack file, 
% trims frames based on specified start and end frame indices, 
% and saves the trimmed frames as an MP4 video

function trimTiffStack(inputFilename, outputFilename, startFrame, endFrame, outputFrameRate)
    
% Open the input TIFF stack
    [inputStack,~,~]=loadtiff(inputFilename);

    % Create a video writer object
    outputVideo = VideoWriter(outputFilename);
    outputVideo.FrameRate = outputFrameRate;
    open(outputVideo);

    % Read and write frames within the specified trim range
    for frameIndex = startFrame:endFrame
        % Read the current frame from the TIFF stack
        currentFrame = inputStack(:,:,frameIndex);

        % Write the current frame to the output video
        writeVideo(outputVideo, (currentFrame-min(currentFrame(:)))./max(currentFrame(:)));
    end

    % Finalize the video
    close(outputVideo);
end


%%  A function that opens an MP4 file, 
% trims frames based on specified start and end frame indices, 
% and saves the trimmed frames as a new MP4 video

function trimMP4Video(inputFilename, outputFilename, startFrame, endFrame, outputFrameRate)
    % Open the input MP4 file
    inputVideo = VideoReader(inputFilename);

    % Create a video writer object
    outputVideo = VideoWriter(outputFilename);
    outputVideo.FrameRate = outputFrameRate;
    open(outputVideo);

    % Read and write frames within the specified trim range
    for frameIndex = startFrame:endFrame
        % Read the current frame
        currentFrame = read(inputVideo, frameIndex);

        % Write the current frame to the output video
        writeVideo(outputVideo, currentFrame);
    end

    % Finalize the video
    close(inputVideo);
    close(outputVideo);
end






















