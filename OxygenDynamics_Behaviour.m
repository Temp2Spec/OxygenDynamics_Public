%% Script for the analysis of behavioural data collected during recording of oxygen dynamics in the mouse barrel cortex
% One custom function needed to open the .abf file. This does not work properly if I have it at the end of the script. 

% A few of the variable used here are coming from the OxygenDynamics_Master script. This is because this analysis script is meant to be called
% from the OxygenDynamics_Wrapper and eventually share workspace with PocketAnalysisWrapper and HypoxicPocket_Master
% If you want to run this individually, you need to modify.

% from Antonis Asiminas, PhD
% Center for Translational Neuromedicine
% University of Copenhagen, Dec 2021


%% Some parameters 
% Change these according to the set-up/experiment 
% The sample frequency for pupil and posture data is 25Hz
BehFs=25;
% the sampling freq for puffs is 1KHz
PuffsSFs=1000;
% The pixel ratio for pupil videos in mm/pxl
PupilPxlRatio = 0.014;
% The pixel ratio for body videos in mm/pxl
PostuPxlRatio = 0.043;
% RecDur in seconds is taken by the main analysis script 

%%
disp('Creating folder for behavioural analysis matrices');

 if not(isfolder('Behaviour_Output')) || strOW=='Y' %if there not a folder with previous oxygen sinks output or we want to overwrite previous data

     Behoutputfolder='Behaviour_Output';
    
 
 else %if there is a folder with previous output and we do not want to overwrite.

     Behoutputfolder=['Behaviour_Output_',datestr(now,30)]; %adds an identifier in the form of yyyymmddTHHMMSS
     
 end 
 
mkdir(Behoutputfolder)

%% Analysis of the posture data %%
% This data may not exist so the entire analyis is run on condition of the file name cell not being empty
% Data is stored in a .csv file 

if ~isempty(Posture) && all(~isnan(Posture))
    
    Posturedata=readtable([Posture,'.csv'],'Delimiter',',');
    PosturedataM=table2array(Posturedata); %converting to array to make life easier
    
    % Checking if the posture data array is a cell array filled with characters
% instead of double 2D array

    if iscellstr(PosturedataM) %if so
        
        PosturedataMtemp=cellfun(@str2double,PosturedataM); %convert the elements to double. Any text will converted to NaNs
        
        PosturedataM=PosturedataMtemp(all(~isnan(PosturedataMtemp),2),:); %remove the NaNs 
    
    end
    clear PosturedataMtemp Posturedata
    % Here I trim the tracking data to match the imaging data.
    % Of course if the beh data are shorter than the imaging data, we need to address that in the combination stage
    
    if size(PosturedataM, 1) > RecDur*BehFs %RecDur is from the master analysis script
        
        PosturedataM=PosturedataM(1:RecDur*BehFs+1,:); % +1 to end up with 30000 points after differentiation
 
    end
    
    % Remove and imputate datapoints with low tracking confidence
    % this is a logical vector with TRUE values in the positions that the
    % probability was less than 95%. 
    Rpawout=PosturedataM(:,4)<0.95;
    RpawinIdx=find(PosturedataM(:,4)>=0.95);
    %
    Lpawout=PosturedataM(:,7)<0.95;
    LpawinIdx=find(PosturedataM(:,7)>=0.95);
    %
    Noseout=PosturedataM(:,10)<0.95;
    NoseinIdx=find(PosturedataM(:,10)>=0.95);
    
    for c=1:length(Rpawout) % Lpawout and Noseout have the same length (the length of recording)
        
        if Rpawout(c) % if this is TRUE (tracking less than 95% confidence)
            
            % find the index of closest non-outlier point
            [~,I]=min(abs(RpawinIdx-c));
            % replace the coordinates
            PosturedataM(c,2:4) = PosturedataM(RpawinIdx(I),2:4);
        end
        
        % repeat for left paw
        if Lpawout(c)
            
            % find the index of closest non-outlier point
            [~,I]=min(abs(LpawinIdx-c));
            % replace the coordinates
            PosturedataM(c,5:7) = PosturedataM(LpawinIdx(I),5:7);
        end
        
        % repeat for nose
        if Noseout(c)
            
            % find the index of closest non-outlier point
            [~,I]=min(abs(NoseinIdx-c));
            % replace the coordinates
            PosturedataM(c,8:10) = PosturedataM(NoseinIdx(I),8:10);
        end
  
    end
   
    
    Right_Nose_Dist=sqrt(sum((PosturedataM(1:end-1,2:3)-PosturedataM(1:end-1,8:9)).^2,2)); %right paw distance from nose
    Left_Nose_Dist=sqrt(sum((PosturedataM(1:end-1,5:6)-PosturedataM(1:end-1,8:9)).^2,2)); %left paw distance from nose
    
    %creating a logical vector with true only when there is a putative
    %logical event (both paws are less than 10mm from nose)
    Groomlogical= Right_Nose_Dist<10/PostuPxlRatio & Left_Nose_Dist<10/PostuPxlRatio; 

    %Separating the different tracked body parts and using Pythagoras theorem to calculate the instateneous movement
    RpawMovpx=sqrt(sum((diff(PosturedataM(:,2:3))).^2, 2));
    LpawMovpx=sqrt(sum((diff(PosturedataM(:,5:6))).^2, 2));
    NoseMovpx=sqrt(sum((diff(PosturedataM(:,8:9))).^2, 2));
    
    % Convert to mm
    RpawMovmm=RpawMovpx*PostuPxlRatio;
    LpawMovmm=LpawMovpx*PostuPxlRatio;
    NoseMovmm=NoseMovpx*PostuPxlRatio;

    %finding global minimum and maximum movement to use in normalisation.
    %most likely from nose movement
    Minmove=min(horzcat(RpawMovmm,LpawMovmm,NoseMovmm),[],'all');
    Maxmove=max(horzcat(RpawMovmm,LpawMovmm,NoseMovmm),[],'all');

    % normalise the values to modulate between 0 and 1. The min and max values are from all data so the movement of different points can be comparable 
    RpawMovZ=(RpawMovmm-Minmove)./Maxmove;
    LpawMovZ=(LpawMovmm-Minmove)./Maxmove;
    NoseMovZ=(NoseMovmm-Minmove)./Maxmove; 
    
    
    %% OUTPUTS

% I save the individual metrics as .mat files to make my life easier when combining

%Tifffiles(1).folder is a variable from the master analysis script and it is the file path to the folder with data from this recording


save([Tifffiles(1).folder,'\',Behoutputfolder,'\','RpawMovmm.mat'], 'RpawMovmm')
save([Tifffiles(1).folder,'\',Behoutputfolder,'\','LpawMovmm.mat'], 'LpawMovmm')
save([Tifffiles(1).folder,'\',Behoutputfolder,'\','NoseMovmm.mat'], 'NoseMovmm')
save([Tifffiles(1).folder,'\',Behoutputfolder,'\','RpawMovZ.mat'], 'RpawMovZ')
save([Tifffiles(1).folder,'\',Behoutputfolder,'\','LpawMovZ.mat'], 'LpawMovZ')
save([Tifffiles(1).folder,'\',Behoutputfolder,'\','NoseMovZ.mat'], 'NoseMovZ')
save([Tifffiles(1).folder,'\',Behoutputfolder,'\','Groomlogical.mat'], 'Groomlogical') 
 
end

%% Analysis of the pupilometry data %%
% As previously, this data may not exist so the entire analyis is run on condition of the file name cell not being empty
% Data is stored in a .csv file 
% 

if ~isempty(Pupil) & all(~isnan(Pupil))
%%
    Pupildata=readtable([Pupil,'.csv'],'Delimiter',',');
    PupildataM=table2array(Pupildata); %converting to array to make life easier
    % As before 
    
    if iscellstr(PupildataM)
        
        PupildataMtemp=cellfun(@str2double,PupildataM); %convert the elements to double. Any text will converted to NaNs
        
        PupildataM=PupildataMtemp(all(~isnan(PupildataMtemp),2),:); %remove the NaNs 
    
    end
    
    if size(PupildataM, 1) > RecDur*BehFs
        
        PupildataM=PupildataM(1:RecDur*BehFs,:); 
 
    end

    % Remove and imputate datapoints with low tracking confidence
    LHorout=PupildataM(:,4)<0.5;
    LHorinIdx=find(~LHorout);
    %
    RHorout=PupildataM(:,7)<0.5;
    RHorinIdx=find(~RHorout);
    %
    Topout=PupildataM(:,10)<0.5;
    TopinIdx=find(~Topout);
    %
    Bottomout=PupildataM(:,13)<0.5;
    BottominIdx=find(~Bottomout);
    
   
    for c=1:length(LHorout) % RHorout,Topout and Bottomout have the same length
        
        if LHorout(c)       
            % find the index of closest non-outlier point
            [~,I]=min(abs(LHorinIdx-c));
            % replace the coordinates
            PupildataM(c,2:4) = PupildataM(LHorinIdx(I),2:4);
        end
        
        % repeat for right tracking point on the pupil
        if RHorout(c)   
            % find the index of closest non-outlier point
            [~,I]=min(abs(RHorinIdx-c));
            % replace the coordinates
            PupildataM(c,5:7) = PupildataM(RHorinIdx(I),5:7);
        end
        
        % repeat for top tracking point on pupil
        if Topout(c)
            
            % find the index of closest non-outlier point
            [~,I]=min(abs(TopinIdx-c));
            % replace the coordinates
            PupildataM(c,8:10) = PupildataM(TopinIdx(I),8:10);
        end
        
        % repeat for bottom tracking point on pupil
        if Bottomout(c)
            
            % find the index of closest non-outlier point
            [~,I]=min(abs(BottominIdx-c));
            % replace the coordinates
            PupildataM(c,11:13) = PupildataM(BottominIdx(I),11:13);
        end

    end

    % Calculate the diameter of the pupil along the horizontal and vertical axis
    HorzDiampx=sqrt(sum((PupildataM(:,2:3)-PupildataM(:,5:6)).^2,2));
    VertDiampx=sqrt(sum((PupildataM(:,8:9)-PupildataM(:,11:12)).^2,2));
    % Finding the outliers based on the speed of transition 
    Horzout=isoutlier(diff(HorzDiampx));
    Vertout=isoutlier(diff(VertDiampx));
    %this is to have the indeces of the non-outliers and help with the data imptutation
    HorzinIndx=find(~Horzout);
    VertinIndx=find(~Horzout);
    
    % Replacing the outliers with the closest non outlier value
    for c=1:length(Horzout) %both Horzout and Vertout have the same length
       
        
        if Horzout(c) % if it is an outlier
            % find the index of closest non-outlier 
            [~,I]=min(abs(HorzinIndx-c));
            % and replace the outlier value with the closest non-outlier
            
            HorzDiampx(c) = HorzDiampx(HorzinIndx(I));
        
        end
        % repeat the process for the vertical diameter
        if Vertout(c) % if it is an outlier
            % find the index of closest non-outlier 
            [~,I]=min(abs(VertinIndx-c));
            % and replace the outlier value with the closest non-outlier
            
            VertDiampx(c) = VertDiampx(VertinIndx(I));
        
        end

    end
    
    % Convert to mm
    HorzDiammm=HorzDiampx*PupilPxlRatio;
    VertDiammm=VertDiampx*PupilPxlRatio;
     
    % This is now the average diameter based on horizontal and vertical.
    PupilDiammm=(VertDiammm+HorzDiammm)./2;

    % Normalise the values
    PupilDiamZ=(PupilDiammm-min(PupilDiammm))./max(PupilDiammm);

    %% OUTPUTS

% I save the individual metrics as .mat files to make my life easier when combining

save([Tifffiles(1).folder,'\',Behoutputfolder,'\','PupilDiammm.mat'], 'PupilDiammm')
save([Tifffiles(1).folder,'\',Behoutputfolder,'\','PupilDiamZ.mat'], 'PupilDiamZ')    
    
end

%% Analysis of the whisker stimulation data %%
% As previously, this data may not exist so the entire analyis is run on condition of the file name cell not being empty
% Data is stored in a .abf file from Clampex

% Channel 1 is the imaging camera (1Hz)
% Channel 2 is noise
% Channel 3 is a behavioural camera (25Hz)
% Channel 4 is the puff 
% Channel 5 is a behavioural camera (25Hz)

if ~isempty(Puff) & all(~isnan(Puff))
%%

    if isa(Puff,'double')        
        Puff=num2str(Puff);
    end
    
    [Puffdata,SampF]=abfload([Puff,'.abf'],'start',0,'stop','e');
    
    %this is to label the pulses of the puffs as 1 
    Pufflogical=Puffdata(:,5)>1;
         
    if size(Pufflogical, 1) > RecDur*PuffsSFs
    
        Pufflogical=Pufflogical(1:RecDur*PuffsSFs,:);     
    end
    %this is to find indeces of the individual pulses within the puffs
    Pulses=find(Pufflogical);
    
    for p=1:length(Pulses)-1   
        %if the next pulse is less than 200 samples away then feel the samples between (if any) with true
        % this is a bit random but the puff pulses are 10sec apart and each stimulation is 150ms apart from the other so 200ms is ok
        if  Pulses(p+1)-Pulses(p)<200
            
             Pufflogical(Pulses(p):Pulses(p+1))=true;
             
        end
    end
    
    %at the end of this Pufflogical has information about the start and then end of the pulse train only
    %this is ok since the imaging sampling frequency is low so there no point looking at individual puffs  
    
    
        %% OUTPUTS

  save([Tifffiles(1).folder,'\',Behoutputfolder,'\','Pufflogical.mat'], 'Pufflogical')
  
end

