%% Script for the statistical analysis of oxygen dynamics as detected by the OxygenDynamics_Master script 
% The oxygen sinks need to be manually curated using the
% OxygenDynamics_Sinks_Curation app.

%What does the script do?
% It reads the same csv file as the OxygenDynamics Wrapper containing paths and metadata for the dataset
% IMPORTANT!!! **This analysis script should be run in the same directory where wrapper is called

% The script goes through the data folders:
% Finds the latest ManualCurOxySinksData OR OxySinksfolder folder (depending on user input) and
% Table_Out_Surges folder for each recording
% Collates the oxygen sink and surges data it finds into separate tables
% Also gets the traces of the bins in the recording area (Mean_ROI_TraceZ) and populates a N-row cell vector (N = number of data folders) 

% Finds the latest Behaviour_Output folder and collect the data for
% RpawMovmm LpawMovmm, NoseMovmm, Groomlogical, PupilDiammm, and Pufflogical
% and place them in a Nx6 cell array (N = number of data folders). 
% IMPORTANT!!! Given that some recordings do NOT have behavioural data or air puffs, some cells will be empty

% STEPS

% (1) Loads all the data and metadata 
% (2) Checks the metadata and asks for input to prints the levels for each grouping variable. 

% IMPORTANT!!
% This data structure will be assumed in the next steps
% The experimental design that Felix has used is the following 
% Between subject factor: genotype (usually there is only one genotype)
% Between subject factor: type of aneasthetic 
% Within subject factor (lvl1): Aneasthesia status (awake-aneasthetised)
% Within subject factor (lvl2): Recording type (baseline-stimulation). This
% can be derived from the existance or not of puff data

% The structure of the experiment dictates 8 separate groups that form the
% 8 columns in the filters cell arrays
%              ISO                           K/X                DrugID
%              /\                            /\
%             /  \                          /  \
%       AWAKE     ANESTHE             AWAKE     ANESTHE         Condition
%        /\          /\                /\          /\  
%       /  \        /  \              /  \        /  \ 
%     (1)  (2)    (3)  (4)          (5)  (6)    (7)  (8)
%  B/LINE PUFFS  B/LINE PUFFS    B/LINE PUFFS  B/LINE PUFFS     StimCond
%
%


% (3) Calculates a series of additional logical behavioural vectors that will connect the behavioural
% readout with Oxygen dynamics
% (4) Calculate some extra metrics for each oxygen sink and surge locus
% (5) Creates data filters according to the data structure found in the metadata 
% (6) Filter data filtering and prepare tables that can be readily used for statistical analysis

% (7) Exporting results. 

% (8) NOT READY YET! Performs generalised linear mixed-effects modelling analysis 
% (This step is not finished yet. I need to feed the data to R OR take sometime to implement a LME
% model for hypothesis testing in Matlab). For now analysis can be done in
% SPSS or Graphpad



% From Antonis Asiminas, PhD
% Copenhagen, Denmark, January 2022,

%Last updated 1 October 2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (1.1) Load the file with all data paths and metadata.
% The same file is used by the wrapper script and has provided all the
% metadata to the datatables produced by the master analysis script

Masterfolder=pwd;
InputD=readtable('ios_stack.csv','Delimiter',',');
if size(InputD,2)<2
    InputD=readtable('ios_stack.csv');
end
%%
Paths=table2cell(InputD(:,{'Paths'}));
Postures=table2cell(InputD(:,{'PostureFile'}));
Pupils=table2cell(InputD(:,{'PupilFile'}));
Puffs=table2cell(InputD(:,{'PuffsFile'}));
Whiskings=table2cell(InputD(:,{'WhiskingFile'}));
Mice=table2cell(InputD(:,{'Mouse'}));
Genotypes=table2cell(InputD(:,{'Genotype'}));
Conditions=table2cell(InputD(:,{'Condition'}));
DrugIDs=table2cell(InputD(:,{'DrugID'}));
Promoters=table2cell(InputD(:,{'Promoter'}));
SampleFs=table2cell(InputD(:,{'SampleF'}));
Pixelsizes=table2cell(InputD(:,{'Pixelsize'}));
Puffs_2Use=table2cell(InputD(:,{'Puff_2use'}));

%% Some inputs for behaviour
% The sample frequency for pupil and posture data is 25Hz
BehFs=25;
% the sampling freq for puffs is 1KHz
PuffsSFs=1000;

%% Inputs from the user to check if curated oxygensink data should be used and what folders

Currated = questdlg('Do you have curated Oxygen sinks data ready for analysis?',...
    'What data',...
    'Yes', 'No','No');

iOS = questdlg('Is this BLI or iOS imaging?',...
    'BLI or iOS',...
    'BLI', 'iOS','BLI');

%
Choosespecific = questdlg(['Do you you want to choose the specific output folders for each recording?' newline 'If not, you can still select whether the newest or oldest output will be used from all recordings'],...
    'What folders to use?',...
    'Yes', 'No','Yes');
if  strcmp(Choosespecific,'No')

    SinkFold_OldORNew = questdlg('Do you want to use the oldest or the most recent data for oxygen sinks?', ...
	'Which folder to use?', ...
	'Oldest','Recent','Recent');

    SurgeFold_OldORNew = questdlg('Do you want to use the oldest or the most recent data for oxygen surges?', ...
	'Which folder to use?', ...
	'Oldest','Recent','Recent');

    BehFold_OldORNew = questdlg('Do you want to use the oldest or the most recent data for behaviour?', ...
	'Which folder to use?', ...
	'Oldest','Recent','Recent');
end

%% Additional input

% This should be a vector with times (in sec from start of recording) that
% specific manual events happened. 
% For example, increase in CO2

% The CRITICAL thing here is that ALL recordings in the study analysed with
% this script are assumed to have the same timings!

ManualEvents=[];

%% (1.2) Collecting data from all data folders and adding the metadata from the source scv 

Table_OxygenSinks_OutCombo=[];
Table_OxygenSurges_OutCombo=[];
ROIs_Traces=cell(height(InputD),6);
Surges_Area=cell(height(InputD),6);
Sinks_Area=cell(height(InputD),6);
Sinks_Traces=cell(height(InputD),6);
Behaviour_data_combo=cell(height(InputD),10);
%%
for datai=1:length(Paths)
%
    cd(Paths{datai}); % go to the data folder of that session

    if isempty(strfind(Paths{1},'\'))
    
        DatafileID=Paths{datai};
    else
        
        DatafileID=Paths{datai}(strfind(Paths{datai},'\')+1:end);
    end 

    Mous = Mice{datai};
    %sometimes the mouse ids are numbers. In order for the rest of the
    %analysis to work I need to do this check and convert to strings
    if ~ischar(Mous)
        Mous=num2str(Mous);
    end
    Cond = Conditions{datai};
    %PiSz = Pixelsizes{datai};
    Gen = Genotypes{datai};
    Promo = Promoters{datai};
    Drug = DrugIDs{datai};
    Posture = Postures{datai};
    Pupil = Pupils{datai};
    Puff = Puffs{datai};
    Whisking = Whiskings{datai};
      
    % finding the folders that contain the sink outputs from master analysis script         
    list=dir;  
    list=list([list.isdir]);
    list=list(3:end);

    wrapper = @(x) contains(x,'OxygenSinks_Output'); 
    SinkFolders=list(cell2mat(cellfun(wrapper, {list(:).name}, 'UniformOutput', false))); 
%
    if  strcmp(Choosespecific,'No')
    
        if strcmp(SinkFold_OldORNew,'Oldest')
        
            % choose the older folder with that name       
            [~,I]=min([SinkFolders.datenum]); 
            SinksDataFolder= [SinkFolders(I).folder,'\',SinkFolders(I).name];
 
        else  
            % or choose the most recent
            [~,I]=max([SinkFolders.datenum]);
            SinksDataFolder= [SinkFolders(I).folder,'\',SinkFolders(I).name];

        end

    else

        % alternatively choose from a GUI prompt 
        [indx_xlsx,tf_xlsx] = listdlg('PromptString',{'Select the oxygen sinks output folder you want.',''},...
            'SelectionMode','single','ListSize',[250,100],'ListString',{SinkFolders.name,''});

        if ~(tf_xlsx) % if nothing was selected 
      
            error('No data folder was selected.')

        end

        SinksDataFolder= [SinkFolders(indx_xlsx).folder,'\',SinkFolders(indx_xlsx).name];


    end
    %%
    %repeat for surges folder
    wrapper = @(x) contains(x,'OxygenSurges');
    SurgeFolders=list(cell2mat(cellfun(wrapper, {list(:).name}, 'UniformOutput', false))); 
    

    if  strcmp(Choosespecific,'No')
    
        if strcmp(SurgeFold_OldORNew,'Oldest')
        
        [~,I]=min([SurgeFolders.datenum]);
        SurgesDataFolder= [SurgeFolders(I).folder,'\',SurgeFolders(I).name];
    
        else

        [~,I]=max([SurgeFolders.datenum]);
        SurgesDataFolder= [SurgeFolders(I).folder,'\',SurgeFolders(I).name];

        end

    else

        [indx_xlsx,tf_xlsx] = listdlg('PromptString',{'Select the oxygen surges output folder you want.',''},...
            'SelectionMode','single','ListSize',[250,100],'ListString',{SurgeFolders.name,''});

        if ~(tf_xlsx) % if nothing was selected 
      
            error('No data folder was selected.')

        end
 
        SurgesDataFolder= [SurgeFolders(indx_xlsx).folder,'\',SurgeFolders(indx_xlsx).name];

    end
%%
    %repeat for behaviour folder
    wrapper = @(x) contains(x,'Behaviour');
    BehFolders=list(cell2mat(cellfun(wrapper, {list(:).name}, 'UniformOutput', false))); 

  
    if  strcmp(Choosespecific,'No') 
    
        if strcmp(BehFold_OldORNew,'Oldest')
        
        [~,I]=min([BehFolders.datenum]);
        BehavDataFolder= [BehFolders(I).folder,'\',BehFolders(I).name];
    
        else

        [~,I]=max([BehFolders.datenum]);
        BehavDataFolder= [BehFolders(I).folder,'\',BehFolders(I).name];

        end

    else

        [indx_xlsx,tf_xlsx] = listdlg('PromptString',{'Select the behavioural output folder you want.',''},...
            'SelectionMode','single','ListSize',[250,100],'ListString',{BehFolders.name,''});

        if ~(tf_xlsx) % if nothing was selected 
      
            error('No data folder was selected.')
  
        end

        BehavDataFolder= [BehFolders(indx_xlsx).folder,'\',BehFolders(indx_xlsx).name];

    end
 %%
    cd(SinksDataFolder) %go in the selected manual curation folder
    d= dir('*mat'); %get the name of the mat file in this folder   
    %load the table with the info on the sinks
    wrapper = @(x) contains(x,'Urefined');
    load(d(cellfun(wrapper, {d(:).name})).name,"Table_OxygenSinks_Out");
    

    if strcmp(Currated, 'Yes')
        if ismember('Included_Pocs', who('-file', d.name)) %if I'm going to use currated data I need to get the logical vector to filter the oxygen sink data
           
            load(d.name,"Included_Pocs")
            Table_OxygenSinks_Out=Table_OxygenSinks_Out(Included_Pocs,:);

        else
            %a disclaimer if curated data were chosen to be analysed but the
            %filter from the curator is not present

            fprintf(['No Oxygen sink currated data were fround in data folder', Paths{datai},'\n'])

% I NEED TO REVISIT THE ANALYSIS AND THE ESTIMATION OF POTENTIAL NOISE BCZ IT THROWS OFF ALL THE SINKS IN SEVERAL CASES.            
            %then if the .mat file contains a vector indicating potential
            %noise sinks (if the master script is from after
            %February 2022)
%             if ismember('Trace_PotentialNoise', who('-file', d.name))
% 
%                 load(d.name,"Trace_PotentialNoise"); % load the filter
% 
%                 %filter the data based on this 
%                 Table_OxygenSinks_Out=Table_OxygenSinks_Out(Trace_PotentialNoise,:);
% 
%             end
            
        end

%    else %if the user indicated that they want non-curated data

%         if ismember('Trace_PotentialNoise', who('-file', d.name)) %if the .mat file contains a vector indicating potential
%             %noise sinks (if the master script is from after February 2022)
% 
%             load(d.name,"Trace_PotentialNoise"); % load the filter
% 
%             %filter the data based on that
%             Table_OxygenSinks_Out=Table_OxygenSinks_Out(Trace_PotentialNoise,:);
  
%        end

    end 

  
    % Checkpoint if there are NO sinks identified to skip to the next

    if size(Table_OxygenSinks_Out,1)>0
    
    
        %select only the columns with the data
    Table_OxygenSinks_Out=Table_OxygenSinks_Out(:,{'RecDuration','RecAreaSize',...
    'MeanCentroid_x','MeanCentroid_y','MeanOxySinkArea_um','MeanOxySinkFilledArea_um','MeanOxySinkDiameter_um','MeanOxySinkPerimeter_um','MeanCircularity',...
    'MeanOxySinkBoundingBox','NumOxySinkEvents','Start','Duration', 'NormOxySinkAmp','Size_modulation','OxySink_Pxls_all'});
    
    %getting the metadata from the input table
    Experiment=cell(height(Table_OxygenSinks_Out),1); Experiment(:)={DatafileID};
    Mouse=cell(height(Table_OxygenSinks_Out),1); Mouse(:)={Mous};    
    Condition=cell(height(Table_OxygenSinks_Out),1);Condition(:)={Cond};
    DrugID=cell(height(Table_OxygenSinks_Out),1);DrugID(:)={Drug};
    if ~isempty(Puff) && all(~isnan(Puff)) 
        PuffStim=true(height(Table_OxygenSinks_Out),1);
    else
        PuffStim=false(height(Table_OxygenSinks_Out),1);
    end
    Genotype=cell(height(Table_OxygenSinks_Out),1); Genotype(:)={Gen};
    Promoter=cell(height(Table_OxygenSinks_Out),1);Promoter(:)={Promo};
    
    %adding the metadata in the table
    Table_OxygenSinks_Out = addvars(Table_OxygenSinks_Out,Experiment,Mouse,Condition,DrugID,PuffStim,Genotype,Promoter,'Before','RecDuration');

    %add the current table to the combo table
    Table_OxygenSinks_OutCombo = [Table_OxygenSinks_OutCombo;Table_OxygenSinks_Out];

    
    %catching also the area of the pockets
    Sinks_Area{datai,1}=DatafileID;
    Sinks_Area{datai,2}=Mous;
    Sinks_Area{datai,3}=Cond;
    Sinks_Area{datai,4}=Drug;
    Sinks_Area{datai,5}=PuffStim(1); 
   
    load(d(cellfun(wrapper, {d(:).name})).name,'OxySinkArea_all');

    if exist('Included_Pocs', 'var') %this means the currated data have been loaded. no need to ask for the response to the prompt again
        OxySinkArea_all=OxySinkArea_all(Included_Pocs,:);
    %else % filter them with the logical vector for potential noise if not
        %OxySinkArea_all=OxySinkArea_all(Trace_PotentialNoise,:);
    end
    Sinks_Area{datai,6}=OxySinkArea_all;
    
    

    %catching the traces of the pockets 
    Sinks_Traces{datai,1}=DatafileID;
    Sinks_Traces{datai,2}=Mous;
    Sinks_Traces{datai,3}=Cond;
    Sinks_Traces{datai,4}=Drug;
    Sinks_Traces{datai,5}=PuffStim(1);
    %load the array with the sink traces
    load(d(cellfun(wrapper, {d(:).name})).name,"Mean_OxySink_Trace_Convo");
    
    
    if exist('Included_Pocs', 'var') %this means the currated data have been loaded. no need to ask for the response to the prompt again
        Mean_OxySink_Trace_Convo=Mean_OxySink_Trace_Convo(Included_Pocs,:);

    %else % filter them with the logical vector for potential noise if not
        %Mean_OxySink_Trace_Convo=Mean_OxySink_Trace_Convo(Trace_PotentialNoise,:);

    end
    Sinks_Traces{datai,6}=Mean_OxySink_Trace_Convo;

    end


    %% Now get the oxygensurge data
    cd(SurgesDataFolder) %go in the most recently modified oxygensurges data
    d= dir('*mat'); 
    wrapper = @(x) contains(x,'Surges');    
    load(d(cellfun(wrapper, {d(:).name})).name,"Table_OxygenSurges_Out") %load the data table 
    

    % Checkpoint if there are NO surges identified to skip to the next
    if size(Table_OxygenSurges_Out,1)>0

    Table_OxygenSurges_Out=Table_OxygenSurges_Out(:,{'RecDuration_Surge','RecAreaSize_Surge',...
    'MeanCentroid_Surge_x','MeanCentroid_Surge_y','MeanOxySurgeArea_um','MeanOxySurgeFilledArea_um','MeanOxySurgeDiameter_um','MeanOxySurgePerimeter_um','MeanCircularity_Surge',...
    'MeanOxySurgeBoundingBox','NumOxySurgeEvents','Start_Surge','Duration_Surge','NormOxySurgeAmp','Size_Surge_modulation','OxySurge_Pxls_all'});
    
    %getting the metadata from the input table
    Experiment=cell(height(Table_OxygenSurges_Out),1); Experiment(:)={DatafileID};
    Mouse=cell(height(Table_OxygenSurges_Out),1); Mouse(:)={Mous};    
    Condition=cell(height(Table_OxygenSurges_Out),1);Condition(:)={Cond};
    DrugID=cell(height(Table_OxygenSurges_Out),1);DrugID(:)={Drug};    
    if ~isempty(Puff) && all(~isnan(Puff)) 
        PuffStim=true(height(Table_OxygenSurges_Out),1);
    else
        PuffStim=false(height(Table_OxygenSurges_Out),1);
    end
    Genotype=cell(height(Table_OxygenSurges_Out),1); Genotype(:)={Gen};
    Promoter=cell(height(Table_OxygenSurges_Out),1);Promoter(:)={Promo};
    
    Table_OxygenSurges_Out = addvars(Table_OxygenSurges_Out,Experiment,Mouse,Condition,DrugID,PuffStim,Genotype,Promoter,'Before','RecDuration_Surge');

    Table_OxygenSurges_OutCombo=[Table_OxygenSurges_OutCombo;Table_OxygenSurges_Out];
    
    
    % Collecting the oxygen surge area
    Surges_Area{datai,1}=DatafileID;
    Surges_Area{datai,2}=Mous;
    Surges_Area{datai,3}=Cond;
    Surges_Area{datai,4}=Drug;
    Surges_Area{datai,5}=PuffStim(1);
    load(d(cellfun(wrapper, {d(:).name})).name,'OxySurgeArea_all');
    Surges_Area{datai,6}=OxySurgeArea_all;


    end
%%
    % in the same mat file there is also the ROI traces

    
    
        ROIs_Traces{datai,1}=DatafileID;  
        ROIs_Traces{datai,2}=Mous;    
        ROIs_Traces{datai,3}=Cond;   
        ROIs_Traces{datai,4}=Drug;    
        ROIs_Traces{datai,5}=PuffStim(1);
    if strcmp(iOS, 'BLI')
        load(d(cellfun(wrapper, {d(:).name})).name,'Mean_ROI_TraceZ');
        ROIs_Traces{datai,6}=Mean_ROI_TraceZ;

    end
    
    %% collecting behavioural data (if any)
    cd(BehavDataFolder)
  
    if ~isempty(Posture) && all(~isnan(Posture))
        
        load('LpawMovmm.mat');
        load('LpawMovZ.mat');
        load('RpawMovmm.mat');
        load('RpawMovZ.mat');
        load('NoseMovmm.mat');
        load('Groomlogical.mat');
        Behaviour_data_combo{datai,1}=LpawMovmm;
        Behaviour_data_combo{datai,2}=RpawMovmm;
        Behaviour_data_combo{datai,3}=NoseMovmm;
        Behaviour_data_combo{datai,4}=Groomlogical;
        Behaviour_data_combo{datai,5}=LpawMovZ;
        Behaviour_data_combo{datai,6}=RpawMovZ;
    end

    if ~isempty(Pupil) && all(~isnan(Pupil))
        load('PupilDiammm.mat');
        load('PupilDiamZ');
        Behaviour_data_combo{datai,7}=PupilDiammm;
        Behaviour_data_combo{datai,8}=PupilDiamZ;
    end
        
    if ~isempty(Puff) && all(~isnan(Puff)) 
        load('Pufflogical.mat')
        Behaviour_data_combo{datai,9}=Pufflogical;
    end

    
    if ~isempty(Whisking) && all(~isnan(Whisking)) 

        list=dir;    
        wrapper = @(x) contains(x,Whisking);

        % the whisking tracking file sometimes is in the behavioural_output
        % folder and sometime in the main data folder. Here I am checking
        if size(list(cell2mat(cellfun(wrapper, {list(:).name}, 'UniformOutput', false))),1)>0

            whisk=readtable(Whisking,'Delimiter',',');

        else
        
            cd ..
            whisk=readtable(Whisking,'Delimiter',',');

        end


        if size(InputD,2)<2    
            whisk=readtable(Whisking);
        end

        %store the whisking trace variance that will help me get the
        %whisking events 
        Behaviour_data_combo{datai,10}=whisk.RTrace_var;

        if size(Behaviour_data_combo{datai,10}, 1) > Table_OxygenSinks_OutCombo.RecDuration{end}*BehFs 
               
            % clipping the trace in case it is longer than the recording
            % *BehFs
            Behaviour_data_combo{datai,10}=Behaviour_data_combo{datai,10}(1:Table_OxygenSinks_OutCombo.RecDuration{end}*BehFs); 
     
        end

    end

    cd(Masterfolder) %return to the masterfolder before going to the next data folder
     
    clear Table_OxygenSinks_Out Included_Pocs Table_OxygenSurges_Out LpawMovmm RpawMovmm NoseMovmm Groomlogical PupilDiammm Pufflogical...
        LpawMovZ RpawMovZ PupilDiamZ Mean_ROI_TraceZ OxySurgeArea_all 
end

clear d I Experiment Mouse Condition DrugID PuffStim Genotype Promoter Table_OxygenSinks_Out Table_OxygenSurges_Out Included_Pocs...
    DatafileID Mous Cond PiSz Gen Promo Drug Posture Pupil Puff list wrapper SinkFolders SurgeFolders BehFolders BehFold_OldORNew...
    SinkFold_OldORNew SurgeFold_OldORNew OxySurgeArea_all OxySinkArea_all BehavDataFolder SinksDataFolder SurgesDataFolder ...
    indx_xlsx tf_xlsx Trace_PotentialNoise whisk Mean_OxySink_Trace_Convo                                 

%% (2) Check metadata and define levels for grouping variables

% I get the unique identifiers for each variable from the oxygen surges
% simply because they are not curated and we will get something out for
% every recording. I could easily get that from the metadata csv as well.
%Experiment_unique=uniqueStrCell(Table_OxygenSurges_OutCombo.Experiment);

%Do a check to see if the mouse IDs in this experiment are numbers or are
%there strings (for example contain letters and numbers)
if ~ischar(Table_OxygenSinks_OutCombo.Mouse{1})

    Table_OxygenSinks_OutCombo.Mouse=cellfun(@num2str,Table_OxygenSinks_OutCombo.Mouse,'UniformOutput',false); %converting to a string to

end

Mice_unique=uniqueStrCell(Table_OxygenSinks_OutCombo.Mouse);
Conditions_unique=uniqueStrCell(Table_OxygenSinks_OutCombo.Condition);
Drugs_unique=uniqueStrCell(Table_OxygenSinks_OutCombo.DrugID);
StimCond_unique=unique(Table_OxygenSinks_OutCombo.PuffStim); %with or without you

fprintf(['The dataset contains ', num2str(datai), ' recording sessions... \n']); 

fprintf(['From ', num2str(length(Mice_unique)), ' mice... \n']); 
for i=1:length(Mice_unique)
    fprintf([Mice_unique{i}, '\n']);
end

fprintf(['The between subjects grouping variable DrugID has ', num2str(length(Drugs_unique)), ' levels... \n']); 
for i=1:length(Drugs_unique)
    fprintf([Drugs_unique{i}, '\n']);
end

fprintf(['The within subjects grouping variable Condition has ', num2str(length(Conditions_unique)), ' levels... \n']); 
for i=1:length(Conditions_unique)
    fprintf([Conditions_unique{i}, '\n']);
end

fprintf(['The within subjects grouping variable Stimulation has ', num2str(length(StimCond_unique)), ' levels... \n']); 
for i=1:length(StimCond_unique)
    if StimCond_unique(i)
       fprintf('With stimulation \n');
    else
        fprintf('Without stimulation \n');
    end
end
clear i
%% (3) Creating some additional logical vectors for the movement of left paw and pupil dilation. These will be used for the probability calculations later
% I also collect previous calculated logical vectors in the same cell
% array. I should probably include this part in the behavioural analysis
% script to declutter this script

% At this point I do not use the logical vectors except from the puffs
% There is not a real hypothesis for either movement of pupil dialation
% events that can be testing with this analysis




% If there are behavioural data available for recordings or manual events, 
% indicate what time window should be used for calculating the influence of behaviour in hypoxic events
if any(~cellfun(@isempty, Behaviour_data_combo),'all') || ~isempty(ManualEvents)   
    prompt = {'Enter window for manual event:','Enter window for body movement/whisking:','Enter window for pupil:','Enter window for whisker stimulation:'};
    dlgtitle = 'Behaviour time window(s)';
    dims = [1 65];
    definput = {'3','3','3','3'};
    Twindows = inputdlg(prompt,dlgtitle,dims,definput);    
end
clear prompt dlgtitle dims definput %keeping the workspace clean
%%
% Cell array containing all logical vectors
Behaviouraldatalogical=cell(size(Behaviour_data_combo,1),6);
% Column 1: logical vector for left paw movement events (calculated here)
% Column 2: logical vector for right paw movement events (calculated here)
% Column 3: logical vector for grooming events (copied)
% Column 4: logical vector for pupil dialation events (calculated here)
% Column 5: logical vector for puffs (copied)
% Column 6: logical vector for whisking events based on the variance of red
% channel of the behavioural video

for dat_i=1:size(Behaviour_data_combo,1)

    if ~isempty(Behaviour_data_combo{dat_i, 1}) %start with LPaw movement because the recording happens at the right cortex
    
        %creating an "empty" vector where I will add the events of left paw movement
        Behaviouraldatalogical{dat_i, 1}=false(length(Behaviour_data_combo{dat_i, 1}),1);
        
        %finding the location of the movement peaks. The 0.5mm is a bit random but it is based on the movement profile when plotting data. 
        [~,locs] = findpeaks(Behaviour_data_combo{dat_i, 1},'MinPeakProminence',0.5); 

        %going through the identified peaks adding "true" in the false logical vector created earlier. For every peak I add points before and after the peak to be true 
        %This is bcz any movement realistically takes at least half second (SF for behaviour is 25Hz) 
        for peaks=1:length(locs)
        
            % a couple of conditionals to avoid crushing when the period
            % extends before or after the legth of recordings
            if locs(peaks)+fix(BehFs/2)<=length(Behaviouraldatalogical{dat_i, 1}) && locs(peaks)-fix(BehFs/2)>=1
     
                Behaviouraldatalogical{dat_i, 1}(locs(peaks)-fix(BehFs/2):locs(peaks)+fix(BehFs/2))=true;

            elseif locs(peaks)+fix(BehFs/2)>length(Behaviouraldatalogical{dat_i, 1})
                
                Behaviouraldatalogical{dat_i, 1}(locs(peaks)-fix(BehFs/2):end)=true;
            else
            
                Behaviouraldatalogical{dat_i, 1}(1:locs(peaks)+fix(BehFs/2))=true;
            end
 
        end
        
        % Now I need to refine these movement event windows
        LPawMovsprops=regionprops(Behaviouraldatalogical{dat_i, 1},'Area','PixelIdxList');
        
        for event_i=1:length(LPawMovsprops)
            
            if LPawMovsprops(event_i).Area<BehFs %if the movement event is smaller that a second
            
                %then filter out this event. 
                Behaviouraldatalogical{dat_i, 1}(LPawMovsprops(event_i).PixelIdxList)=false;
                
                %the reason is that temporal resolution for oxygen imaging is 1Hz so there is no point including events that are shorter than that 1sec.               
            end      
        end
            
        %copy the logical vectors for grooming to the Behaviouraldatalogical bcz it makes more sense.
        Behaviouraldatalogical{dat_i, 3}=Behaviour_data_combo{dat_i,4};
        
    end 

    %repeat for right paw
    if ~isempty(Behaviour_data_combo{dat_i, 2}) 
    
        %creating an "empty" vector where I will add the events of right paw movement
        Behaviouraldatalogical{dat_i, 2}=false(length(Behaviour_data_combo{dat_i, 2}),1);
        
        %finding the location of the movement peaks. The 0.5mm is a bit random but it is based on the movement profile when plotting data. 
        [~,locs] = findpeaks(Behaviour_data_combo{dat_i, 2},'MinPeakProminence',0.5); 

        %going through the identified peaks adding "true" in the false logical vector created earlier. For every peak I add points before and after the peak to be true 
        %This is bcz any movement realistically takes at least half second (SF for behaviour is 25Hz) 
        for peaks=1:length(locs)
        
            % a couple of conditionals to avoid crushing when the period
            % extends before or after the legth of recordings
            if locs(peaks)+fix(BehFs/2)<=length(Behaviouraldatalogical{dat_i, 2}) && locs(peaks)-fix(BehFs/2)>=1
     
                Behaviouraldatalogical{dat_i, 2}(locs(peaks)-fix(BehFs/2):locs(peaks)+fix(BehFs/2))=true;

            elseif locs(peaks)+fix(BehFs/2)>length(Behaviouraldatalogical{dat_i, 2})
                
                Behaviouraldatalogical{dat_i, 2}(locs(peaks)-fix(BehFs/2):end)=true;
            else
            
                Behaviouraldatalogical{dat_i, 2}(1:locs(peaks)+fix(BehFs/2))=true;
            end
 
        end
       
        % Now I need to refine these movement event windows
        LPawMovsprops=regionprops(Behaviouraldatalogical{dat_i, 2},'Area','PixelIdxList');
        
        for event_i=1:length(LPawMovsprops)
            
            if LPawMovsprops(event_i).Area<BehFs %if the movement event is smaller that a second
            
                %then filter out this event. 
                Behaviouraldatalogical{dat_i, 2}(LPawMovsprops(event_i).PixelIdxList)=false;
                
                %the reason is that temporal resolution for oxygen imaging is 1Hz so there is no point including events that are shorter than that 1sec.               
            end      
        end

    end 

    %repeat process for pupil dialation events
   
    if ~isempty(Behaviour_data_combo{dat_i, 7})

        %creating an "empty" vector where I will add the events of pupil dialation
        Behaviouraldatalogical{dat_i, 4}=false(length(Behaviour_data_combo{dat_i, 7}),1);
        
        %smoothing the pupil diameter data
        pupilsmooth_temp=smooth(Behaviour_data_combo{dat_i, 7}, 200);
        %calculating the change in pupil diameter
        pupilsmoothdiff_temp=diff(pupilsmooth_temp);
        
        [locs, ~] = peakfinder(pupilsmoothdiff_temp);
        
        for peaks=1:length(locs)

            if locs(peaks)+fix(BehFs/2)<=length(Behaviouraldatalogical{dat_i, 4}) && locs(peaks)-fix(BehFs/2)>=1
     
                Behaviouraldatalogical{dat_i, 4}(locs(peaks)-fix(BehFs/2):locs(peaks)+fix(BehFs/2))=true;

            elseif locs(peaks)+fix(BehFs/2)>length(Behaviouraldatalogical{dat_i, 4})
                
                Behaviouraldatalogical{dat_i, 4}(locs(peaks)-fix(BehFs/2):end)=true;
            else
            
                Behaviouraldatalogical{dat_i, 4}(1:locs(peaks)+fix(BehFs/2))=true;
            end
        end 
    end

    if ~isempty(Behaviour_data_combo{dat_i, 9})
        
        %copy the logical vectors for puffs to the Behaviouraldatalogical bcz it makes more sense.
        Behaviouraldatalogical{dat_i, 5}=Behaviour_data_combo{dat_i,9};

    end


    %and finally for the whisking events
    if ~isempty(Behaviour_data_combo{dat_i, 10})


        % The 5 for variance is a bit random but it is based on the movement profile when plotting data.
        Behaviouraldatalogical{dat_i, 6}=Behaviour_data_combo{dat_i, 10}>5;
        
        
        % Now I need to refine these movement event windows
        LPawMovsprops=regionprops(Behaviouraldatalogical{dat_i, 6},'Area','PixelIdxList');
        
        for event_i=1:length(LPawMovsprops)
            
            if LPawMovsprops(event_i).Area<25 %if the whisking event is smaller that a second
            
                %then filter out this event. 
                Behaviouraldatalogical{dat_i, 6}(LPawMovsprops(event_i).PixelIdxList)=false;
                
                %the reason is that temporal resolution for oxygen imaging is 1Hz so there is no point including events that are shorter than that 1sec.               
            end      
        end



    end


end

%keeping the workspace clean
clear locs dat_i event_i LPawMovsprops pupilsmooth_temp pupilsmoothdiff_temp PupilDilationprops peaks 
%% (4) Calculate additional metrics for each oxygen sink and surge locus

% These are metrics that can be calculated for each oxygen sink and surge.
% I can then append the Table and add more if I want to later
OxySinkArea_Norm=(Table_OxygenSinks_OutCombo.MeanOxySinkArea_um)./cell2mat(Table_OxygenSinks_OutCombo.RecAreaSize); %area of pocket normalised to the area of recording (sink area per um squared)
OxySinkArea_Filled_Norm=(Table_OxygenSinks_OutCombo.MeanOxySinkFilledArea_um)./cell2mat(Table_OxygenSinks_OutCombo.RecAreaSize); %area of pocket normalised to the area of recording
NumOxySinkEvents_Norm=((Table_OxygenSinks_OutCombo.NumOxySinkEvents)*60)./cell2mat(Table_OxygenSinks_OutCombo.RecDuration); %number of events normalised by the time of recording in seconds. Multiplying by 60 to have events per minute
MeanOxySinkEvent_Duration=cellfun(@mean,Table_OxygenSinks_OutCombo.Duration);
%POOL ALL EVENTS DURATIONS FOR DISTRIBUTION PER GROUP id
MeanOxySinkEvent_NormAmp=cellfun(@mean,Table_OxygenSinks_OutCombo.NormOxySinkAmp); %the amplitude of signal change relative to the signal for 20sec before the event (mean across all events of a oxygen sink event)
MeanOxySinkEvent_SizeModulation=cellfun(@mean,Table_OxygenSinks_OutCombo.Size_modulation); % mean size modulation of oxygen sinks during events. This size end-size start/size start. Negative means the sink get smaller

AdditionalOxySinkMetrics=table(OxySinkArea_Norm,OxySinkArea_Filled_Norm,NumOxySinkEvents_Norm,MeanOxySinkEvent_Duration,MeanOxySinkEvent_NormAmp,MeanOxySinkEvent_SizeModulation);

OxySurgeArea_Norm=(Table_OxygenSurges_OutCombo.MeanOxySurgeArea_um)./cell2mat(Table_OxygenSurges_OutCombo.RecAreaSize_Surge); %area of surge locus normalised to the area of recording (sink area per um squared)
OxySurgeArea_Filled_Norm=(Table_OxygenSurges_OutCombo.MeanOxySurgeFilledArea_um)./cell2mat(Table_OxygenSurges_OutCombo.RecAreaSize_Surge); %filled area of surge locus normalised to the area of recording
NumOxySurgeEvents_Norm=(Table_OxygenSurges_OutCombo.NumOxySurgeEvents)./cell2mat(Table_OxygenSurges_OutCombo.RecDuration_Surge); %number of events normalised by the time of recording in second (events per second)
MeanOxySurgeEvent_Duration=cellfun(@mean,Table_OxygenSurges_OutCombo.Duration_Surge);
MeanOxySurgeEvent_NormAmp=cellfun(@mean,Table_OxygenSurges_OutCombo.NormOxySurgeAmp); %the amplitude of signal change relative to the signal for 20sec before the event (mean across all events of a oxygen sink event)
MeanOxySurgeEvent_SizeModulation=cellfun(@mean,Table_OxygenSurges_OutCombo.Size_Surge_modulation); % mean size modulation of oxygen sinks during events. This size end-size start/size start. Negative means the sink get smaller

AdditionalOxySurgeMetrics=table(OxySurgeArea_Norm,OxySurgeArea_Filled_Norm,NumOxySurgeEvents_Norm,MeanOxySurgeEvent_Duration,MeanOxySurgeEvent_NormAmp,MeanOxySurgeEvent_SizeModulation);

clear OxySinkArea_Norm OxySinkArea_Filled_Norm NumOxySinkEvents_Norm MeanOxySinkEvent_Duration MeanOxySinkEvent_NormAmp MeanOxySinkEvent_SizeModulation ...
    OxySurgeArea_Norm OxySurgeArea_Filled_Norm NumOxySurgeEvents_Norm MeanOxySurgeEvent_Duration MeanOxySurgeEvent_NormAmp MeanOxySurgeEvent_SizeModulation

%%
if strcmp(iOS, 'BLI')
% Here I calculate a number of recording specific time series that can be
% plotted with behaviour and puffs (if present).
% These are, of course, recording specific.

%ROI mean signal for every sampling point
%ROIs entropy for every sampling point
%ROIs covariance coefficient

%ROIs mean signal differential
%ROIs signal differential entropy
%ROIs signal differential covariance coeffient
%These last three will tell me if the increase or decrease in oxygen (dynamics) exhibit spatial organisation (high entropy/cov coef) or not 

%A variable to capture correlations between traces
TraceCorrs=cell(size(ROIs_Traces,1),11);
TraceCorrs(:,1:5)=ROIs_Traces(:,1:5);

for i=1:size(ROIs_Traces,1) %go through all the recordings 
    %create temporary vectors to capture ROI entropy and mean ROI intensity
    %for each sampling point
    
    tempROImean=zeros(1,size(ROIs_Traces{i,6},2));
    tempCovCoef=zeros(1,size(ROIs_Traces{i,6},2));
    tempEntropyvec=zeros(1,size(ROIs_Traces{i,6},2));



    %the differential for each ROI
    tempROIdifs=diff(ROIs_Traces{i,6},1,2);

    tempROIdiffmean=zeros(1,size(ROIs_Traces{i,6},2));
    tempdiffCovCoef=zeros(1,size(ROIs_Traces{i,6},2));
    tempdiffEntropyvec=zeros(1,size(ROIs_Traces{i,6},2));


    for k=1:size(ROIs_Traces{i,6},2) % go through each sampling point 

        tempROImean(1,k)=mean(ROIs_Traces{i,6}(:,k)); %calculate the mean ROI intensity across the ROIs
        tp = mat2gray(ROIs_Traces{i,6}(:,k)); %converting to greyscale to avoid the negative values
        tempCovCoef(1,k)=std(tp)/mean(tp);
        tp=tp./sum(tp); %calculate the probability for each value 
        tempEntropyvec(1,k)=info_entropy(tp, 'bit'); % calculate Shannon's entropy for the ROI signal distribution
        

    end


    for k=1:size(tempROIdifs,2) % go through each sampling point for the differential 

        %notice that I store the values at k+1 position since the first
        %sampling does not have a differential
        tempROIdiffmean(1,k+1)=mean(tempROIdifs(:,k)); %calculate the mean ROI differential intensity across the ROIs
        tp = mat2gray(tempROIdifs(:,k)); %converting to greyscale to avoid the negative values
        tempdiffCovCoef(1,k+1)=std(tp)/mean(tp);
        tp=tp./sum(tp); %calculate the probability for each value 
        tempdiffEntropyvec(1,k+1)=info_entropy(tp, 'bit'); % calculate Shannon's entropy for the ROI signal distribution
        

    end

    TraceCorrs{i,6}=corr(zscore(tempROImean)',tempCovCoef');
    TraceCorrs{i,7}=corr(zscore(tempROImean)',tempEntropyvec');
    TraceCorrs{i,8}=corr(tempCovCoef',tempEntropyvec');

    TraceCorrs{i,9}=corr(zscore(tempROIdiffmean)',tempdiffCovCoef');
    TraceCorrs{i,10}=corr(zscore(tempROIdiffmean)',tempdiffEntropyvec');
    TraceCorrs{i,11}=corr(tempdiffCovCoef',tempdiffEntropyvec');


    %capture the calculate traces
    ROIs_Traces{i,7}=zscore(tempROImean);
    ROIs_Traces{i,8}=tempCovCoef;
    ROIs_Traces{i,9}=tempEntropyvec;

    ROIs_Traces{i,10}=zscore(tempROIdiffmean);
    ROIs_Traces{i,11}=tempdiffCovCoef;
    ROIs_Traces{i,12}=tempdiffEntropyvec;



end

clear i k tempEntropyvec tempROImean tempCovCoef tempROIdiffmean tempdiffCovCoef tempdiffEntropyvec tempROIdifs tp
end

%Another thing I can calculate is the deferential of each ROI and then the mean differential
%Then I can calculate the cov coef and entropy of that differential
%%

% Number of sink events at any given time and overall area
% Number of surge events at any given time and overall area
% It will be so much easier to get that from Surges and Sinks area
% variables. 
% I will update in the next weeks.

NumOngoingOxysinks=cell(datai,6);
NumOngoingOxysinks(:,1:5)=Sinks_Traces(:,1:5);
TotalSinkArea_Norm=cell(datai,6);
TotalSinkArea_Norm(:,1:5)=Sinks_Traces(:,1:5);
TotalSinkArea_um=cell(datai,6);
TotalSinkArea_um(:,1:5)=Sinks_Traces(:,1:5);
NumOngoingOxysurges=cell(datai,6);
NumOngoingOxysurges(:,1:5)=Surges_Area(:,1:5);
TotalSurgeArea=cell(datai,6);
TotalSurgeArea(:,1:5)=Surges_Area(:,1:5);
SinksRaster=cell(datai,6);
SinksRaster(:,1:5)=Sinks_Traces(:,1:5);
SurgesRaster=cell(datai,6);
SurgesRaster(:,1:5)= Surges_Area(:,1:5);


for i=1:datai %go through the different recordings

    %creating a filter to get only the sinks for that recording
    wrapper = @(x) strcmp(x, Sinks_Traces{i,1});
    TempIndex=cell2mat(cellfun(wrapper, Table_OxygenSinks_OutCombo.Experiment, 'UniformOutput', false)); 

    Firstsinkofrec=find(TempIndex,1,'first');

    tempvec_totalsinks=zeros(1,Table_OxygenSinks_OutCombo.RecDuration{Firstsinkofrec}); %the temporary vector that will have the number of sinks at every frame will be as long as the recording  
    tempvec_totalsinkarea_norm=zeros(1,Table_OxygenSinks_OutCombo.RecDuration{Firstsinkofrec}); %the temporary vector that will have the total area of sinks at every frame will be as long as the recording  
    tempvec_totalsinkarea_um=zeros(1,Table_OxygenSinks_OutCombo.RecDuration{Firstsinkofrec});   
    %get the starts,duration,amp, sizemod,pixels for the sinks in that recording
    thisrecording=Table_OxygenSinks_OutCombo(TempIndex,9:end);
 

    % create a raster for that recording to have all the hypoxic events
    % from the sink areas
    SinksRasterRec=false(height(thisrecording),Table_OxygenSinks_OutCombo.RecDuration{Firstsinkofrec});

    for k=1:height(thisrecording) %go through the sinks

        for g=1:length(thisrecording.Start{k}) %and for every sink go through the event starts

            %add a value in the vector from the start to start+duration
            tempvec_totalsinks(thisrecording.Start{k}(g):thisrecording.Start{k}(g)+thisrecording.Duration{k}(g)-1)=...
            tempvec_totalsinks(thisrecording.Start{k}(g):thisrecording.Start{k}(g)+thisrecording.Duration{k}(g)-1)+1;

            %add the size of the sink for that period normalised by the
            %area of recording
            tempvec_totalsinkarea_norm(thisrecording.Start{k}(g):thisrecording.Start{k}(g)+thisrecording.Duration{k}(g)-1)=...
                tempvec_totalsinkarea_norm(thisrecording.Start{k}(g):thisrecording.Start{k}(g)+thisrecording.Duration{k}(g)-1)+...
                (length(thisrecording.OxySink_Pxls_all{k})/thisrecording.RecAreaSize{k});


            tempvec_totalsinkarea_um(thisrecording.Start{k}(g):thisrecording.Start{k}(g)+thisrecording.Duration{k}(g)-1)=...
                tempvec_totalsinkarea_um(thisrecording.Start{k}(g):thisrecording.Start{k}(g)+thisrecording.Duration{k}(g)-1)+...
                (length(thisrecording.OxySink_Pxls_all{k})*Pixelsizes{i}^2);




             SinksRasterRec(k,thisrecording.Start{k}(g):thisrecording.Start{k}(g)+thisrecording.Duration{k}(g)-1)=true;

        end

    end

    NumOngoingOxysinks{i,6}=tempvec_totalsinks;
    TotalSinkArea_Norm{i,6}=tempvec_totalsinkarea_norm;
    TotalSinkArea_um{i,6}=tempvec_totalsinkarea_um;
    SinksRaster{i,6}=SinksRasterRec;

    if strcmp(iOS, 'BLI')
    %adding a few more correlation values between sink area and the covcoef/entropy measures
    TraceCorrs{i,12}=corr(tempvec_totalsinkarea_norm',ROIs_Traces{i,8}');
    TraceCorrs{i,13}=corr(tempvec_totalsinkarea_norm',ROIs_Traces{i,9}');

    end
    
    %repeat the process for surges
    tempvec_totalsurges=zeros(1,Table_OxygenSinks_OutCombo.RecDuration{Firstsinkofrec});
    tempvec_totalsurgearea=zeros(1,Table_OxygenSinks_OutCombo.RecDuration{Firstsinkofrec});
    %creating a filter to get only the surges for that recording
    wrapper = @(x) strcmp(x, ROIs_Traces{i,1});
    TempIndex=cell2mat(cellfun(wrapper, Table_OxygenSurges_OutCombo.Experiment, 'UniformOutput', false)); 

    %get the starts and durations for the sinks in that recording
    thisrecording=Table_OxygenSurges_OutCombo(TempIndex,9:end);


    % create a raster for that recording to have all the oxygen surge events
    % from the surge areas
    SurgesRasterRec=false(height(thisrecording),Table_OxygenSinks_OutCombo.RecDuration{Firstsinkofrec});

    for k=1:height(thisrecording) %go through the surges

        for g=1:length(thisrecording.Start_Surge{k}) %and for every surge go through the events

            %add a value in the vector from the start to start+duration
            tempvec_totalsurges(thisrecording.Start_Surge{k}(g):thisrecording.Start_Surge{k}(g)+thisrecording.Duration_Surge{k}(g)-1)=...
                tempvec_totalsurges(thisrecording.Start_Surge{k}(g):thisrecording.Start_Surge{k}(g)+thisrecording.Duration_Surge{k}(g)-1)+1;


            tempvec_totalsurgearea(thisrecording.Start_Surge{k}(g):thisrecording.Start_Surge{k}(g)+thisrecording.Duration_Surge{k}(g)-1)=...
                tempvec_totalsurgearea(thisrecording.Start_Surge{k}(g):thisrecording.Start_Surge{k}(g)+thisrecording.Duration_Surge{k}(g)-1)+...
                (length(thisrecording.OxySurge_Pxls_all{k})/thisrecording.RecAreaSize_Surge{k});

            SurgesRasterRec(k,thisrecording.Start_Surge{k}(g):thisrecording.Start_Surge{k}(g)+thisrecording.Duration_Surge{k}(g)-1)=true;

        end

    end

    NumOngoingOxysurges{i,6}=tempvec_totalsurges;
    TotalSurgeArea{i,6}=tempvec_totalsurgearea;
    SurgesRaster{i,6}=SurgesRasterRec;
    
    if strcmp(iOS, 'BLI')
    %adding a few more correlation values between surge area and the covcoef/entropy measures
    TraceCorrs{i,14}=corr(tempvec_totalsurgearea',ROIs_Traces{i,8}');
    TraceCorrs{i,15}=corr(tempvec_totalsurgearea',ROIs_Traces{i,9}');

    end

end

clear i k g wrapper TempIndex tempvec_totalsinks tempvec_totalsinkarea_norm tempvec_totalsurges tempvec_totalsurgearea thisrecording SinksRasterRec SurgesRasterRec tempvec_totalsinkarea_um


%% Connecting physiology with behaviour 
% Bin the movement and the pupil size
% in 10bins and then see what the different traces values are during those
% periods. The problem here is that I need to average on seconds to make
% sense for the slow bioiluminescence data
% Behaviour_data_combo;
% Reminder
% Column 1: LpawMovmm;
% Column 2: RpawMovmm;
% Column 3: NoseMovmm;
% Column 4: Groomlogical;
% Column 5: LpawMovZ;
% Column 6: RpawMovZ;
% Column 7: PupilDiammm;
% Column 8: PupilDiamZ;
% Column 9: Puffs;

%Create the cell arrays that will capture the binned trace data for each
%behavioural metric
BinsLPaw_10percentiles=cell(size(InputD,1),11);
BinsRPaw_10percentiles=cell(size(InputD,1),11);
BinsPupil_10percentiles=cell(size(InputD,1),11);
% Structure for every one
% each column is a recording
% Column 1:  ROIs mean signal 
% Column 2:  ROIs cov coef
% Column 3:  ROIs entropy
% Column 4:  NumOngoingOxysinks;
% Column 5:  TotalSinkArea;
% Column 6:  NumOngoingOxysurges;
% Column 7: TotalSurgeArea;
% Column 8: ROIs diff mean signal 
% Column 9: ROIs diff cov coef
% Column 10: ROIs diff entropy

for i=1:size(Behaviour_data_combo,1) %go through each recording
    if ~isempty(Behaviour_data_combo{i,1}) %if there are behavioural tracking data for this recording

        %resize the behaviour to have the same length with imaging data
        Lpaw_resize=reshape(Behaviour_data_combo{i,1},size(Behaviour_data_combo{i,1},1)/size(ROIs_Traces{i,6},2),[]);
        Lpaw_resize_mean=mean(Lpaw_resize,1);

        Rpaw_resize=reshape(Behaviour_data_combo{i,2},size(Behaviour_data_combo{i,2},1)/size(ROIs_Traces{i,6},2),[]);
        Rpaw_resize_mean=mean(Rpaw_resize,1);

        %find the boundaries of 10 percentiles for each behavioural measure  
        Lpaw_binedges = prctile(Lpaw_resize_mean,[0 10 20 30 40 50 60 70 80 90 100]);
        Rpaw_binedges = prctile(Rpaw_resize_mean,[0 10 20 30 40 50 60 70 80 90 100]);

        tempvec_ROImean=nan(1,length(Lpaw_binedges)-1);
        tempvec_ROIscovcoef=nan(1,length(Lpaw_binedges)-1);
        tempvec_ROIsentropy=nan(1,length(Lpaw_binedges)-1);
        tempvec_NumOngoingOxysinks=nan(1,length(Lpaw_binedges)-1);
        tempvec_TotalSinkArea_norm=nan(1,length(Lpaw_binedges)-1);
        tempvec_TotalSinkArea_um=nan(1,length(Lpaw_binedges)-1);
        tempvec_NumOngoingOxysurges=nan(1,length(Lpaw_binedges)-1);
        tempvec_TotalSurgeArea=nan(1,length(Lpaw_binedges)-1);
        tempvec_ROIsdiffmean=nan(1,length(Lpaw_binedges)-1);
        tempvec_ROIsdiffcovcoef=nan(1,length(Lpaw_binedges)-1);
        tempvec_ROIsdiffentropy=nan(1,length(Lpaw_binedges)-1);

        for k=2:length(Lpaw_binedges)

            % find the elements in Lpaw_resize_mean  and Rpaw_resize_mean
            % that belong in the current percentile bin
            % and extract the values from each trace
            if strcmp(iOS, 'BLI')
            tempvec_ROImean(k-1)=mean(ROIs_Traces{i,7}(Lpaw_resize_mean>Lpaw_binedges(k-1) & Lpaw_resize_mean<=Lpaw_binedges(k)));
            tempvec_ROIscovcoef(k-1)=mean(ROIs_Traces{i,8}(Lpaw_resize_mean>Lpaw_binedges(k-1) & Lpaw_resize_mean<=Lpaw_binedges(k)));
            tempvec_ROIsentropy(k-1)=mean(ROIs_Traces{i,9}(Lpaw_resize_mean>Lpaw_binedges(k-1) & Lpaw_resize_mean<=Lpaw_binedges(k)));
            
            tempvec_ROIsdiffmean(k-1)=mean(ROIs_Traces{i,10}(Lpaw_resize_mean>Lpaw_binedges(k-1) & Lpaw_resize_mean<=Lpaw_binedges(k)));
            tempvec_ROIsdiffcovcoef(k-1)=mean(ROIs_Traces{i,11}(Lpaw_resize_mean>Lpaw_binedges(k-1) & Lpaw_resize_mean<=Lpaw_binedges(k)));
            tempvec_ROIsdiffentropy(k-1)=mean(ROIs_Traces{i,12}(Lpaw_resize_mean>Lpaw_binedges(k-1) & Lpaw_resize_mean<=Lpaw_binedges(k)));

            end

            tempvec_NumOngoingOxysinks(k-1)=mean(NumOngoingOxysinks{i,6}(Lpaw_resize_mean>Lpaw_binedges(k-1) & Lpaw_resize_mean<=Lpaw_binedges(k)));
            tempvec_TotalSinkArea_norm(k-1)=mean(TotalSinkArea_Norm{i,6}(Lpaw_resize_mean>Lpaw_binedges(k-1) & Lpaw_resize_mean<=Lpaw_binedges(k)));

            tempvec_TotalSinkArea_um(k-1)=mean(TotalSinkArea_um{i,6}(Lpaw_resize_mean>Lpaw_binedges(k-1) & Lpaw_resize_mean<=Lpaw_binedges(k)));
           
            tempvec_NumOngoingOxysurges(k-1)=mean(NumOngoingOxysurges{i,6}(Lpaw_resize_mean>Lpaw_binedges(k-1) & Lpaw_resize_mean<=Lpaw_binedges(k)));
            tempvec_TotalSurgeArea(k-1)=mean(TotalSurgeArea{i,6}(Lpaw_resize_mean>Lpaw_binedges(k-1) & Lpaw_resize_mean<=Lpaw_binedges(k)));
            


        end  
        BinsLPaw_10percentiles{i,1}=tempvec_ROImean;
        BinsLPaw_10percentiles{i,2}=tempvec_ROIscovcoef;
        BinsLPaw_10percentiles{i,3}=tempvec_ROIsentropy;
        BinsLPaw_10percentiles{i,4}=tempvec_NumOngoingOxysinks;
        BinsLPaw_10percentiles{i,5}=tempvec_TotalSinkArea_norm;
        BinsLPaw_10percentiles{i,6}=tempvec_NumOngoingOxysurges;
        BinsLPaw_10percentiles{i,7}=tempvec_TotalSurgeArea;
        BinsLPaw_10percentiles{i,8}=tempvec_ROIsdiffmean;
        BinsLPaw_10percentiles{i,9}=tempvec_ROIsdiffcovcoef;
        BinsLPaw_10percentiles{i,10}=tempvec_ROIsdiffentropy;
        BinsLPaw_10percentiles{i,11}=tempvec_TotalSinkArea_um;

        %repeat the process for the Rpaw tracking

        tempvec_ROImean=nan(1,length(Rpaw_binedges)-1);
        tempvec_ROIscovcoef=nan(1,length(Rpaw_binedges)-1);
        tempvec_ROIsentropy=nan(1,length(Rpaw_binedges)-1);
        tempvec_NumOngoingOxysinks=nan(1,length(Rpaw_binedges)-1);
        tempvec_TotalSinkArea_norm=nan(1,length(Rpaw_binedges)-1);
        tempvec_TotalSinkArea_um=nan(1,length(Lpaw_binedges)-1);
        tempvec_NumOngoingOxysurges=nan(1,length(Rpaw_binedges)-1);
        tempvec_TotalSurgeArea=nan(1,length(Rpaw_binedges)-1);
        tempvec_ROIsdiffmean=nan(1,length(Rpaw_binedges)-1);
        tempvec_ROIsdiffcovcoef=nan(1,length(Rpaw_binedges)-1);
        tempvec_ROIsdiffentropy=nan(1,length(Rpaw_binedges)-1);

        for k=2:length(Rpaw_binedges)

            if strcmp(iOS, 'BLI')
            tempvec_ROImean(k-1)=mean(ROIs_Traces{i,7}(Rpaw_resize_mean>Rpaw_binedges(k-1) & Rpaw_resize_mean<=Rpaw_binedges(k)));
            tempvec_ROIscovcoef(k-1)=mean(ROIs_Traces{i,8}(Rpaw_resize_mean>Rpaw_binedges(k-1) & Rpaw_resize_mean<=Rpaw_binedges(k)));
            tempvec_ROIsentropy(k-1)=mean(ROIs_Traces{i,9}(Rpaw_resize_mean>Rpaw_binedges(k-1) & Rpaw_resize_mean<=Rpaw_binedges(k)));
            tempvec_ROIsdiffmean(k-1)=mean(ROIs_Traces{i,10}(Rpaw_resize_mean>Rpaw_binedges(k-1) & Rpaw_resize_mean<=Rpaw_binedges(k)));
            tempvec_ROIsdiffcovcoef(k-1)=mean(ROIs_Traces{i,11}(Rpaw_resize_mean>Rpaw_binedges(k-1) & Rpaw_resize_mean<=Rpaw_binedges(k)));
            tempvec_ROIsdiffentropy(k-1)=mean(ROIs_Traces{i,12}(Rpaw_resize_mean>Rpaw_binedges(k-1) & Rpaw_resize_mean<=Rpaw_binedges(k)));
            end

            tempvec_NumOngoingOxysinks(k-1)=mean(NumOngoingOxysinks{i,6}(Rpaw_resize_mean>Rpaw_binedges(k-1) & Rpaw_resize_mean<=Rpaw_binedges(k)));
            tempvec_TotalSinkArea_norm(k-1)=mean(TotalSinkArea_Norm{i,6}(Rpaw_resize_mean>Rpaw_binedges(k-1) & Rpaw_resize_mean<=Rpaw_binedges(k)));

            tempvec_TotalSinkArea_um(k-1)=mean(TotalSinkArea_um{i,6}(Rpaw_resize_mean>Rpaw_binedges(k-1) & Rpaw_resize_mean<=Rpaw_binedges(k)));

            tempvec_NumOngoingOxysurges(k-1)=mean(NumOngoingOxysurges{i,6}(Rpaw_resize_mean>Rpaw_binedges(k-1) & Rpaw_resize_mean<=Rpaw_binedges(k)));
            tempvec_TotalSurgeArea(k-1)=mean(TotalSurgeArea{i,6}(Rpaw_resize_mean>Rpaw_binedges(k-1) & Rpaw_resize_mean<=Rpaw_binedges(k)));
            

        end
        BinsRPaw_10percentiles{i,1}=tempvec_ROImean;
        BinsRPaw_10percentiles{i,2}=tempvec_ROIscovcoef;
        BinsRPaw_10percentiles{i,3}=tempvec_ROIsentropy;
        BinsRPaw_10percentiles{i,4}=tempvec_NumOngoingOxysinks;
        BinsRPaw_10percentiles{i,5}=tempvec_TotalSinkArea_norm;
        BinsRPaw_10percentiles{i,6}=tempvec_NumOngoingOxysurges;
        BinsRPaw_10percentiles{i,7}=tempvec_TotalSurgeArea;
        BinsRPaw_10percentiles{i,8}=tempvec_ROIsdiffmean;
        BinsRPaw_10percentiles{i,9}=tempvec_ROIsdiffcovcoef;
        BinsRPaw_10percentiles{i,10}=tempvec_ROIsdiffentropy;
        BinsRPaw_10percentiles{i,11}=tempvec_TotalSinkArea_um;
      

    else


        %if there is no behaviour for that recording, then add a NaN
        %vector. I have hard coded a ten element one based on 10
        %percentiles but I should revisit this to make it more robust
        BinsLPaw_10percentiles{i,1}=nan(1,10);
        BinsLPaw_10percentiles{i,2}=nan(1,10);
        BinsLPaw_10percentiles{i,3}=nan(1,10);
        BinsLPaw_10percentiles{i,4}=nan(1,10);
        BinsLPaw_10percentiles{i,5}=nan(1,10);
        BinsLPaw_10percentiles{i,6}=nan(1,10);
        BinsLPaw_10percentiles{i,7}=nan(1,10);
        BinsLPaw_10percentiles{i,8}=nan(1,10);
        BinsLPaw_10percentiles{i,9}=nan(1,10);
        BinsLPaw_10percentiles{i,10}=nan(1,10);
        BinsLPaw_10percentiles{i,11}=nan(1,10);



        BinsRPaw_10percentiles{i,1}=nan(1,10);
        BinsRPaw_10percentiles{i,2}=nan(1,10);
        BinsRPaw_10percentiles{i,3}=nan(1,10);
        BinsRPaw_10percentiles{i,4}=nan(1,10);
        BinsRPaw_10percentiles{i,5}=nan(1,10);
        BinsRPaw_10percentiles{i,6}=nan(1,10);
        BinsRPaw_10percentiles{i,7}=nan(1,10);
        BinsRPaw_10percentiles{i,8}=nan(1,10);
        BinsRPaw_10percentiles{i,9}=nan(1,10);
        BinsRPaw_10percentiles{i,10}=nan(1,10);
        BinsRPaw_10percentiles{i,11}=nan(1,10);


    end

    %repeat for the pupil tracking
    if ~isempty(Behaviour_data_combo{i,7}) %if there are pupil tracking data for this recording


        %resize the behaviour to have the same length with imaging data
        Pupil_resize=reshape(Behaviour_data_combo{i,7},size(Behaviour_data_combo{i,7},1)/size(ROIs_Traces{i,6},2),[]);
        Pupil_resize_mean=mean(Pupil_resize,1);

        %find the boundaries of 10 percentiles for each behavioural measure  
        Pupil_binedges = prctile(Pupil_resize_mean,[0 10 20 30 40 50 60 70 80 90 100]);

        tempvec_ROImean=nan(1,length(Pupil_binedges)-1);
        tempvec_ROIscovcoef=nan(1,length(Pupil_binedges)-1);
        tempvec_ROIsentropy=nan(1,length(Pupil_binedges)-1);
        tempvec_NumOngoingOxysinks=nan(1,length(Pupil_binedges)-1);
        tempvec_TotalSinkArea_norm=nan(1,length(Pupil_binedges)-1);
        tempvec_TotalSinkArea_um=nan(1,length(Pupil_binedges)-1);
        tempvec_NumOngoingOxysurges=nan(1,length(Pupil_binedges)-1);
        tempvec_TotalSurgeArea=nan(1,length(Pupil_binedges)-1);
        tempvec_ROIsdiffmean=nan(1,length(Pupil_binedges)-1);
        tempvec_ROIsdiffcovcoef=nan(1,length(Pupil_binedges)-1);
        tempvec_ROIsdiffentropy=nan(1,length(Pupil_binedges)-1);

        for k=2:length(Pupil_binedges)

            if strcmp(iOS, 'BLI')
            tempvec_ROImean(k-1)=mean(ROIs_Traces{i,7}(Pupil_resize_mean>Pupil_binedges(k-1) & Pupil_resize_mean<=Pupil_binedges(k)));
            tempvec_ROIscovcoef(k-1)=mean(ROIs_Traces{i,8}(Pupil_resize_mean>Pupil_binedges(k-1) & Pupil_resize_mean<=Pupil_binedges(k)));
            tempvec_ROIsentropy(k-1)=mean(ROIs_Traces{i,9}(Pupil_resize_mean>Pupil_binedges(k-1) & Pupil_resize_mean<=Pupil_binedges(k)));

            tempvec_ROIsdiffmean(k-1)=mean(ROIs_Traces{i,10}(Pupil_resize_mean>Pupil_binedges(k-1) & Pupil_resize_mean<=Pupil_binedges(k)));
            tempvec_ROIsdiffcovcoef(k-1)=mean(ROIs_Traces{i,11}(Pupil_resize_mean>Pupil_binedges(k-1) & Pupil_resize_mean<=Pupil_binedges(k)));
            tempvec_ROIsdiffentropy(k-1)=mean(ROIs_Traces{i,12}(Pupil_resize_mean>Pupil_binedges(k-1) & Pupil_resize_mean<=Pupil_binedges(k)));

            end


            tempvec_NumOngoingOxysinks(k-1)=mean(NumOngoingOxysinks{i,6}(Pupil_resize_mean>Pupil_binedges(k-1) & Pupil_resize_mean<=Pupil_binedges(k)));
            tempvec_TotalSinkArea_norm(k-1)=mean(TotalSinkArea_Norm{i,6}(Pupil_resize_mean>Pupil_binedges(k-1) & Pupil_resize_mean<=Pupil_binedges(k)));

            tempvec_TotalSinkArea_um(k-1)=mean(TotalSinkArea_um{i,6}(Pupil_resize_mean>Pupil_binedges(k-1) & Pupil_resize_mean<=Pupil_binedges(k)));

            tempvec_NumOngoingOxysurges(k-1)=mean(NumOngoingOxysurges{i,6}(Pupil_resize_mean>Pupil_binedges(k-1) & Pupil_resize_mean<=Pupil_binedges(k)));
            tempvec_TotalSurgeArea(k-1)=mean(TotalSurgeArea{i,6}(Pupil_resize_mean>Pupil_binedges(k-1) & Pupil_resize_mean<=Pupil_binedges(k)));


        end


        BinsPupil_10percentiles{i,1}=tempvec_ROImean;          
        BinsPupil_10percentiles{i,2}=tempvec_ROIscovcoef;            
        BinsPupil_10percentiles{i,3}=tempvec_ROIsentropy;            
        BinsPupil_10percentiles{i,4}=tempvec_NumOngoingOxysinks;            
        BinsPupil_10percentiles{i,5}=tempvec_TotalSinkArea_norm;            
        BinsPupil_10percentiles{i,6}=tempvec_NumOngoingOxysurges;            
        BinsPupil_10percentiles{i,7}=tempvec_TotalSurgeArea;            
        BinsPupil_10percentiles{i,8}=tempvec_ROIsdiffmean;            
        BinsPupil_10percentiles{i,9}=tempvec_ROIsdiffcovcoef;        
        BinsPupil_10percentiles{i,10}=tempvec_ROIsdiffentropy;
        BinsPupil_10percentiles{i,11}=tempvec_TotalSinkArea_um;


    end


end

clear i k Lpaw_resize Lpaw_resize_mean Rpaw_resize Rpaw_resize_mean Lpaw_binedges Rpaw_binedges  tempvec_ROImean tempvec_ROIscovcoef tempvec_ROIsentropy tempvec_NumOngoingOxysinks tempvec_TotalSinkArea_norm...
    tempvec_NumOngoingOxysurges tempvec_TotalSurgeArea tempvec_ROIsdiffmean tempvec_ROIsdiffcovcoef tempvec_ROIsdiffentropy Pupil_resize Pupil_resize_mean Pupil_binedges tempvec_TotalSinkArea_um


%% Event-based analysis
% For now I do this only for mean ROIsignal but I can repeat for all other
% trace types I have calculated (e.g cov coef, entropy, oxygen sink area etc)

% start with the manual events 
if strcmp(iOS, 'BLI')
if ~isempty(ManualEvents) %check if there are any in this set of recordings

    %initiate the variable that will capture the mean snipped trace for
    %each recording
    ROI_MeanTraceSnip_ManEvents=cell(1,length(ManualEvents));
    

    for k=1:length(ManualEvents) % and go through the different events 

        
        %create a cell array that will be populated with the ROI mean signal snipets from
        %the manual event k
        SeparateEvent=cell(datai,6);
        SeparateEvent(:,1:5)=ROIs_Traces(:,1:5);
        
        for i=1:size(ROIs_Traces,1) %now go through the recordings one by one  


            SeparateEvent{i,6}=ROIs_Traces{i,7}(ManualEvents(k)*SampleFs{i}-str2double(Twindows{1})*SampleFs{i}:ManualEvents(k)*SampleFs{i}+str2double(Twindows{1})*SampleFs{i});
         
        end

        ROI_MeanTraceSnip_ManEvents{k}=SeparateEvent;

    end
end

clear i k ROItraceSnip

% Now Lpaw movement events

if any(~cellfun(@isempty,Behaviouraldatalogical(:,1))) %check if there are any recordings with logical vectors for left paw movement events


    %initiate the variable that will capture the mean snipped trace for
    %each recording
    ROI_MeanTraceSnip_LPawEvents=cell(size(ROIs_Traces,1),6);
    ROI_MeanTraceSnip_LPawEvents(:,1:5)=ROIs_Traces(:,1:5);
    
    for i=1:size(ROIs_Traces,1) %now go through the recordings one by one 

        if ~isempty(Behaviouraldatalogical(i,1)) % second check if there is a logical vector for that specific recording

            Events=bwconncomp(Behaviouraldatalogical{i,1}); %find the events

            
            ROItraceSnip=nan(length(Events.PixelIdxList), str2double(Twindows{2})*SampleFs{i}*2+1); %preallocate the space for all the signal snippets

            for k=1:length(Events.PixelIdxList)


                %index the ROItraces based on the start of each event
                %(corrected for sampling frequency differences) and the
                %timewindows selected
                %I need to include some conditionals in case the event
                %start is too close to the start or the end of the
                %recording

                if fix(Events.PixelIdxList{k}(1)*SampleFs{i}/BehFs)-str2double(Twindows{2})*SampleFs{i}<1 %if the event -window extends beyond the start.

                    diff=1-(fix(Events.PixelIdxList{k}(1)*SampleFs{i}/BehFs)-str2double(Twindows{2})*SampleFs{i}); % find the difference

                    %and then add the smaller sniper in the right way
                    ROItraceSnip(k,1+diff:end)=ROIs_Traces{i,7}(1:fix(Events.PixelIdxList{k}(1)*SampleFs{i}/BehFs)+str2double(Twindows{2})*SampleFs{i});

                elseif fix(Events.PixelIdxList{k}(1)*SampleFs{i}/BehFs)+str2double(Twindows{2})*SampleFs{i}>length(ROIs_Traces{i,7}) %if the event +window extends past the end
                    
                    diff= (fix(Events.PixelIdxList{k}(1)*SampleFs{i}/BehFs)+str2double(Twindows{2})*SampleFs{i})-length(ROIs_Traces{i,7}); % find the difference

                    ROItraceSnip(k,1:end-diff)=ROIs_Traces{i,7}(fix(Events.PixelIdxList{k}(1)*SampleFs{i}/BehFs)-str2double(Twindows{2})*SampleFs{i}:end);

                else

                    ROItraceSnip(k,:)=ROIs_Traces{i,7}(fix(Events.PixelIdxList{k}(1)*SampleFs{i}/BehFs)-str2double(Twindows{2})*SampleFs{i}:fix(Events.PixelIdxList{k}(1)*SampleFs{i}/BehFs)+str2double(Twindows{2})*SampleFs{i});

                end

            end

            ROI_MeanTraceSnip_LPawEvents{i,6}=mean(ROItraceSnip,1,'omitnan');

        end

    end

end


clear i k ROItraceSnip Events diff

% Now repeat for Rpaw movement events

if any(~cellfun(@isempty,Behaviouraldatalogical(:,2))) %check if there are any recordings with logical vectors for left paw movement events


    %initiate the variable that will capture the mean snipped trace for
    %each recording
    ROI_MeanTraceSnip_RPawEvents=cell(size(ROIs_Traces,1),6);
    ROI_MeanTraceSnip_RPawEvents(:,1:5)=ROIs_Traces(:,1:5);
    
    for i=1:size(ROIs_Traces,1) %now go through the recordings one by one 

        if ~isempty(Behaviouraldatalogical(i,2)) % second check if there is a logical vector for that specific recording

            Events=bwconncomp(Behaviouraldatalogical{i,2}); %find the events

            
            ROItraceSnip=nan(length(Events.PixelIdxList), str2double(Twindows{2})*SampleFs{i}*2+1); %preallocate the space for all the signal snippets

            for k=1:length(Events.PixelIdxList)

                %index the ROItraces based on the start of each event
                %(corrected for sampling frequency differences) and the
                %timewindows selected

                if fix(Events.PixelIdxList{k}(1)*SampleFs{i}/BehFs)-str2double(Twindows{2})*SampleFs{i}<1

                    diff=1-(fix(Events.PixelIdxList{k}(1)*SampleFs{i}/BehFs)-str2double(Twindows{2})*SampleFs{i});

                    ROItraceSnip(k,1+diff:end)=ROIs_Traces{i,7}(1:fix(Events.PixelIdxList{k}(1)*SampleFs{i}/BehFs)+str2double(Twindows{2})*SampleFs{i});

                elseif fix(Events.PixelIdxList{k}(1)*SampleFs{i}/BehFs)+str2double(Twindows{2})*SampleFs{i}>length(ROIs_Traces{i,7})

                    diff= (fix(Events.PixelIdxList{k}(1)*SampleFs{i}/BehFs)+str2double(Twindows{2})*SampleFs{i})-length(ROIs_Traces{i,7});

                    ROItraceSnip(k,1:end-diff)=ROIs_Traces{i,7}(fix(Events.PixelIdxList{k}(1)*SampleFs{i}/BehFs)-str2double(Twindows{2})*SampleFs{i}:end);

                else

                    ROItraceSnip(k,:)=ROIs_Traces{i,7}(fix(Events.PixelIdxList{k}(1)*SampleFs{i}/BehFs)-str2double(Twindows{2})*SampleFs{i}:fix(Events.PixelIdxList{k}(1)*SampleFs{i}/BehFs)+str2double(Twindows{2})*SampleFs{i});
                end
                

            end

            ROI_MeanTraceSnip_RPawEvents{i,6}=mean(ROItraceSnip,1,'omitnan');

        end

    end

end


clear i k ROItraceSnip Events diff


% Repeat for grooming events

if any(~cellfun(@isempty,Behaviouraldatalogical(:,3))) %check if there are any recordings with logical vectors for left paw movement events


    %initiate the variable that will capture the mean snipped trace for
    %each recording
    ROI_MeanTraceSnip_GroomEvents=cell(size(ROIs_Traces,1),6);
    ROI_MeanTraceSnip_GroomEvents(:,1:5)=ROIs_Traces(:,1:5);
    
    for i=1:size(ROIs_Traces,1) %now go through the recordings one by one 

        if ~isempty(Behaviouraldatalogical(i,3)) % second check if there is a logical vector for that specific recording

            Events=bwconncomp(Behaviouraldatalogical{i,3}); %find the events

            
            ROItraceSnip=nan(length(Events.PixelIdxList), str2double(Twindows{2})*SampleFs{i}*2+1); %preallocate the space for all the signal snippets

            for k=1:length(Events.PixelIdxList)


                %index the ROItraces based on the start of each event
                %(corrected for sampling frequency differences) and the
                %timewindows selected
                if fix(Events.PixelIdxList{k}(1)*SampleFs{i}/BehFs)-str2double(Twindows{2})*SampleFs{i}<1

                    diff=1-(fix(Events.PixelIdxList{k}(1)*SampleFs{i}/BehFs)-str2double(Twindows{2})*SampleFs{i});

                    ROItraceSnip(k,1+diff:end)=ROIs_Traces{i,7}(1:fix(Events.PixelIdxList{k}(1)*SampleFs{i}/BehFs)+str2double(Twindows{2})*SampleFs{i});

            
                elseif fix(Events.PixelIdxList{k}(1)*SampleFs{i}/BehFs)+str2double(Twindows{2})*SampleFs{i}>length(ROIs_Traces{i,7})

                    diff=(fix(Events.PixelIdxList{k}(1)*SampleFs{i}/BehFs)+str2double(Twindows{2})*SampleFs{i})-length(ROIs_Traces{i,7});


                    ROItraceSnip(k,1:end-diff)=ROIs_Traces{i,7}(fix(Events.PixelIdxList{k}(1)*SampleFs{i}/BehFs)-str2double(Twindows{2})*SampleFs{i}:end);

            
                else
                
                    ROItraceSnip(k,:)=ROIs_Traces{i,7}(fix(Events.PixelIdxList{k}(1)*SampleFs{i}/BehFs)-str2double(Twindows{2})*SampleFs{i}:fix(Events.PixelIdxList{k}(1)*SampleFs{i}/BehFs)+str2double(Twindows{2})*SampleFs{i});
            
                end
            end

            ROI_MeanTraceSnip_GroomEvents{i,6}=mean(ROItraceSnip,1,'omitnan');

        end

    end

end


clear i k ROItraceSnip Events diff


% Repeat for pupil dialation events

if any(~cellfun(@isempty,Behaviouraldatalogical(:,4))) %check if there are any recordings with logical vectors for left paw movement events


    %initiate the variable that will capture the mean snipped trace for
    %each recording
    ROI_MeanTraceSnip_PupilEvents=cell(size(ROIs_Traces,1),6);
    ROI_MeanTraceSnip_PupilEvents(:,1:5)=ROIs_Traces(:,1:5);
    
    for i=1:size(ROIs_Traces,1) %now go through the recordings one by one 

        if ~isempty(Behaviouraldatalogical(i,4)) % second check if there is a logical vector for that specific recording

            Events=bwconncomp(Behaviouraldatalogical{i,4}); %find the events

            
            ROItraceSnip=nan(length(Events.PixelIdxList), str2double(Twindows{3})*SampleFs{i}*2+1); %preallocate the space for all the signal snippets

            for k=1:length(Events.PixelIdxList)


                %index the ROItraces based on the start of each event
                %(corrected for sampling frequency differences) and the
                %timewindows selected

                if fix(Events.PixelIdxList{k}(1)*SampleFs{i}/BehFs)-str2double(Twindows{3})<1

                    diff=1-(fix(Events.PixelIdxList{k}(1)*SampleFs{i}/BehFs)-str2double(Twindows{3}));
                
                    ROItraceSnip(k,1+diff:end)=ROIs_Traces{i,7}(1:fix(Events.PixelIdxList{k}(1)*SampleFs{i}/BehFs)+str2double(Twindows{3})*SampleFs{i});
                
                elseif fix(Events.PixelIdxList{k}(1)*SampleFs{i}/BehFs)+str2double(Twindows{3})>length(ROIs_Traces{i,7})

                
                    diff=(fix(Events.PixelIdxList{k}(1)*SampleFs{i}/BehFs)+str2double(Twindows{3}))-length(ROIs_Traces{i,7});

                    ROItraceSnip(k,1:end-diff)=ROIs_Traces{i,7}(fix(Events.PixelIdxList{k}(1)*SampleFs{i}/BehFs)-str2double(Twindows{3})*SampleFs{i}:end);
                
                
                else

                     ROItraceSnip(k,:)=ROIs_Traces{i,7}(fix(Events.PixelIdxList{k}(1)*SampleFs{i}/BehFs)-str2double(Twindows{3})*SampleFs{i}:fix(Events.PixelIdxList{k}(1)*SampleFs{i}/BehFs)+str2double(Twindows{3})*SampleFs{i});

                end

            end

            ROI_MeanTraceSnip_PupilEvents{i,6}=mean(ROItraceSnip,1,'omitnan');

        end

    end

end


clear i k ROItraceSnip Events diff


% For puffs events

if any(~cellfun(@isempty,Behaviouraldatalogical(:,5))) %check if there are any recordings with logical vectors for puff events


    %initiate the variable that will capture the mean snipped trace for
    %each recording
    ROI_MeanTraceSnip_PuffEvents=cell(size(ROIs_Traces,1),6);
    ROI_MeanTraceSnip_PuffEvents(:,1:5)=ROIs_Traces(:,1:5);
    
    for i=1:size(ROIs_Traces,1) %now go through the recordings one by one 

        if ~isempty(Behaviouraldatalogical(i,5)) % second check if there is a logical vector for that specific recording

            Events=bwconncomp(Behaviouraldatalogical{i,5}); %find the events

            
            ROItraceSnip=nan(length(Events.PixelIdxList), str2double(Twindows{4})*SampleFs{i}*2+1); %preallocate the space for all the signal snippets

            for k=1:length(Events.PixelIdxList)

                %index the ROItraces based on the start of each event
                %(corrected for sampling frequency differences) and the
                %timewindows selected

                if fix(Events.PixelIdxList{k}(1)*SampleFs{i}/PuffsSFs)-str2double(Twindows{4})*SampleFs{i}<1

                    diff=1-(fix(Events.PixelIdxList{k}(1)*SampleFs{i}/PuffsSFs)-str2double(Twindows{4})*SampleFs{i});

                    ROItraceSnip(k,1+diff:end)=ROIs_Traces{i,7}(1:fix(Events.PixelIdxList{k}(1)*SampleFs{i}/PuffsSFs)+str2double(Twindows{4})*SampleFs{i});

                elseif fix(Events.PixelIdxList{k}(1)*SampleFs{i}/PuffsSFs)+str2double(Twindows{4})*SampleFs{i}>length(ROIs_Traces{i,7})

                    diff=(fix(Events.PixelIdxList{k}(1)*SampleFs{i}/PuffsSFs)+str2double(Twindows{4})*SampleFs{i})-length(ROIs_Traces{i,7});

                    ROItraceSnip(k,1:end-diff)=ROIs_Traces{i,7}(fix(Events.PixelIdxList{k}(1)*SampleFs{i}/PuffsSFs)-str2double(Twindows{4})*SampleFs{i}:end);

                else

                    ROItraceSnip(k,:)=ROIs_Traces{i,7}(fix(Events.PixelIdxList{k}(1)*SampleFs{i}/PuffsSFs)-str2double(Twindows{4})*SampleFs{i}:fix(Events.PixelIdxList{k}(1)*SampleFs{i}/PuffsSFs)+str2double(Twindows{4})*SampleFs{i});

                end
                

            end

            ROI_MeanTraceSnip_PuffEvents{i,6}=mean(ROItraceSnip,1,'omitnan');

        end

    end

end


% And finally for whisking events

if any(~cellfun(@isempty,Behaviouraldatalogical(:,6))) %check if there are any recordings with logical vectors for whisking events


    %initiate the variable that will capture the mean snipped trace for
    %each recording
    ROI_MeanTraceSnip_WhiskEvents=cell(size(ROIs_Traces,1),6);
    ROI_MeanTraceSnip_WhiskEvents(:,1:5)=ROIs_Traces(:,1:5);
    
    for i=1:size(ROIs_Traces,1) %now go through the recordings one by one 

        if ~isempty(Behaviouraldatalogical(i,6)) % second check if there is a logical vector for that specific recording

            Events=bwconncomp(Behaviouraldatalogical{i,6}); %find the events
            
            ROItraceSnip=nan(length(Events.PixelIdxList), str2double(Twindows{2})*SampleFs{i}*2+1); %preallocate the space for all the signal snippets

            for k=1:length(Events.PixelIdxList)


                %index the ROItraces based on the start of each event
                %(corrected for sampling frequency differences) and the
                %timewindows selected
                %I need to include some conditionals in case the event
                %start is too close to the start or the end of the
                %recording

                if fix(Events.PixelIdxList{k}(1)*SampleFs{i}/BehFs)-str2double(Twindows{2})*SampleFs{i}<1 %if the event -window extends beyond the start.

                    diff=1-(fix(Events.PixelIdxList{k}(1)*SampleFs{i}/BehFs)-str2double(Twindows{2})*SampleFs{i}); % find the difference

                    %and then add the smaller sniper in the right way
                    ROItraceSnip(k,1+diff:end)=ROIs_Traces{i,7}(1:fix(Events.PixelIdxList{k}(1)*SampleFs{i}/BehFs)+str2double(Twindows{2})*SampleFs{i});

                elseif fix(Events.PixelIdxList{k}(1)*SampleFs{i}/BehFs)+str2double(Twindows{2})*SampleFs{i}>length(ROIs_Traces{i,7}) %if the event +window extends past the end
                    
                    diff= (fix(Events.PixelIdxList{k}(1)*SampleFs{i}/BehFs)+str2double(Twindows{2})*SampleFs{i})-length(ROIs_Traces{i,7}); % find the difference

                    ROItraceSnip(k,1:end-diff)=ROIs_Traces{i,7}(fix(Events.PixelIdxList{k}(1)*SampleFs{i}/BehFs)-str2double(Twindows{2})*SampleFs{i}:end);

                else

                    ROItraceSnip(k,:)=ROIs_Traces{i,7}(fix(Events.PixelIdxList{k}(1)*SampleFs{i}/BehFs)-str2double(Twindows{2})*SampleFs{i}:fix(Events.PixelIdxList{k}(1)*SampleFs{i}/BehFs)+str2double(Twindows{2})*SampleFs{i});

                end

            end

            ROI_MeanTraceSnip_WhiskEvents{i,6}=mean(ROItraceSnip,1,'omitnan');

        end

    end

end

clear i k ROItraceSnip Events diff

end
%% Alligning the traces of all oxygen sink events for each recording

for i = 1:size(Sinks_Traces,1) % go through the recordings

    %get the info for the sinks in a given recording
    SinksTemp=Table_OxygenSinks_OutCombo(logical(strcmp(Table_OxygenSinks_OutCombo.Experiment(:),Sinks_Traces{i,1})...
        .*strcmp(Table_OxygenSinks_OutCombo.Mouse(:),Sinks_Traces{i,2})),:);

    %% If the sinks are curated this is redundant. But if there are not,
    % crap can go through and this analysis crashes. So I need to filter
    % some crap out.
    for j = 1:height(SinksTemp)

        if any(SinksTemp.Duration{j}>150) || any(SinksTemp.Duration{j}<3)

            
            SinksTemp.NumOxySinkEvents(j)=SinksTemp.NumOxySinkEvents(j)-sum(SinksTemp.Duration{j}>150|SinksTemp.Duration{j}<3);
            SinksTemp.Start{j}=SinksTemp.Start{j}(SinksTemp.Duration{j}<=150&SinksTemp.Duration{j}>=3);
            SinksTemp.NormOxySinkAmp{j}=SinksTemp.NormOxySinkAmp{j}(SinksTemp.Duration{j}<=150&SinksTemp.Duration{j}>=3);
            SinksTemp.Size_modulation{j}=SinksTemp.Size_modulation{j}(SinksTemp.Duration{j}<=150&SinksTemp.Duration{j}>=3);
            SinksTemp.Duration{j}=SinksTemp.Duration{j}(SinksTemp.Duration{j}<=150&SinksTemp.Duration{j}>=3);
            
        end
 
    end

    SinksTemp=SinksTemp(SinksTemp.NumOxySinkEvents~=0,:);

    %%

    %finding what is the longest duration of a sink across all loci for this recording
    maxdursink=max(cellfun(@max, SinksTemp.Duration));
    totalsinkevents=sum(cellfun(@length, SinksTemp.Duration));
    %initiate the array that will capture all the sinks events separate
    %traces for that recording
    SinkAlligned_temp = nan(totalsinkevents,maxdursink);

    counter=1; %initiate a counter
    for j = 1:height(SinksTemp) %each loci

        for k = 1:length(SinksTemp.Start{j}) %each oxygen sink event

            %get the trace for the sink event from the loci
            tempeventtrace=Sinks_Traces{i,6}(j,SinksTemp.Start{j}(k):SinksTemp.Start{j}(k)+SinksTemp.Duration{j}(k)-1);

            %adjust it's values so start is 0 and the rest of the values
            %are relative to that
            tempeventtrace=tempeventtrace-tempeventtrace(1);

            SinkAlligned_temp(counter,1:length(tempeventtrace))=tempeventtrace;

            counter=counter+1;

        end

    end

    Sinks_Traces{i,7}=SinkAlligned_temp;

end

clear i j k SinksTemp tempeventtrace maxdursink SinkAlligned_temp totalsinkevents
%% (5)  Creating filters for the different conditions

% The structure of the experiment dictates 8 separate groups that form the
% 8 columns in the filters cell arrays
%              ISO                           K/X                DrugID
%              /\                            /\
%             /  \                          /  \
%       AWAKE     ANESTHE             AWAKE     ANESTHE         Condition
%        /\          /\                /\          /\  
%       /  \        /  \              /  \        /  \ 
%     (1)  (2)    (3)  (4)          (5)  (6)    (7)  (8)
%  B/LINE PUFFS  B/LINE PUFFS    B/LINE PUFFS  B/LINE PUFFS     StimCond
%
%

% Cell arrays where I will store the filters
FiltersOxySinksMetrics=cell(4,length(Drugs_unique)*length(Conditions_unique)*length(StimCond_unique));
FiltersOxySurgesMetrics=cell(4,length(Drugs_unique)*length(Conditions_unique)*length(StimCond_unique));
% Column 1 has awake/baseline recording for ISO
% Column 2 has awake/puffs recording for ISO
% Column 3 has anesthetised/baseline recording for ISO
% Column 4 has anesthetised/puffs recording for ISO
% Column 5-8 have the repeat for KX
counter=1;
for i=1:length(Drugs_unique) %go through DrugIDs first (Between subjects variable)
    % a temp logical vector for each drug for the sinktable/and metrics
    % a temp logical vector for each drug for the surgetable/and metrics
    wrapper = @(x) strcmp(x, Drugs_unique{i});
    TempOxySinkIndexDrugID=cell2mat(cellfun(wrapper,Table_OxygenSinks_OutCombo.DrugID(:), 'UniformOutput', false));
    TempOxySurgeIndexDrugID=cell2mat(cellfun(wrapper,Table_OxygenSurges_OutCombo.DrugID(:), 'UniformOutput', false));
    %these indeces are meaningful for the additional metrics tables

    for k=1:length(Conditions_unique) %for each drug ID go through the conciousness condition of mice (this is a within-subject variable and therefor nested in the DrugID loop)
 
        % a temp logical vector for each condtion for the sinktable/and metrics
        % a temp logical vector for each condtion for the surgetable/and metrics
        wrapper = @(x) strcmp(x, Conditions_unique{k});
        TempOxySinkIndexCondition=cell2mat(cellfun(wrapper,Table_OxygenSinks_OutCombo.Condition(:), 'UniformOutput', false));
        TempOxySurgeIndexCondition=cell2mat(cellfun(wrapper,Table_OxygenSurges_OutCombo.Condition(:), 'UniformOutput', false));
                
        for g=1:length(StimCond_unique) %for type of stimulation (with/without). This is a second level nested within-subject variable and therefore nested in the Condition loop

            % a temp logical vector for each stim condition for the sinktable/and metrics
            % a temp logical vector for each stim condition for the surgetable/and metrics
            TempOxySinkIndexPuffs=Table_OxygenSinks_OutCombo.PuffStim(:)==StimCond_unique(g);
            TempOxySurgeIndexPuffs=Table_OxygenSurges_OutCombo.PuffStim(:)==StimCond_unique(g);

            %multiply all temp logical vectors for sinks elementwise to get one and store
            %it in the FiltersOxysinksMetrics
            FiltersOxySinksMetrics{4,counter}=logical(TempOxySinkIndexDrugID.*TempOxySinkIndexCondition.*TempOxySinkIndexPuffs);

            %multiply all logical vectors for surgess elementwise to get one and store
            %it in the FiltersOxysurgesMetrics
            FiltersOxySurgesMetrics{4,counter}=logical(TempOxySurgeIndexDrugID.*TempOxySurgeIndexCondition.*TempOxySurgeIndexPuffs);

            FiltersOxySinksMetrics{1,counter}= Drugs_unique{i}; 
            FiltersOxySinksMetrics{2,counter}= Conditions_unique{k};    
               
            FiltersOxySurgesMetrics{1,counter}= Drugs_unique{i};
            FiltersOxySurgesMetrics{2,counter}= Conditions_unique{k};

            if StimCond_unique(g)                   
                FiltersOxySinksMetrics{3,counter}='With stimulation';
                FiltersOxySurgesMetrics{3,counter}='With stimulation';
            else
                FiltersOxySinksMetrics{3,counter}='Without stimulation';
                FiltersOxySurgesMetrics{3,counter}='Without stimulation';
               
            end

            
            FiltersOxySinksMetrics{5,counter}=unique(Table_OxygenSinks_OutCombo.Mouse(FiltersOxySinksMetrics{4,counter}));    
            FiltersOxySurgesMetrics{5,counter}=unique(Table_OxygenSurges_OutCombo.Mouse(FiltersOxySurgesMetrics{4,counter}));


            counter=counter+1;
        end
    end
end

clear i k g counter wrapper TempOxySinkIndexDrugID TempOxySurgeIndexDrugID TempOxySinkIndexCondition TempOxySurgeIndexCondition TempOxySinkIndexPuffs TempOxySurgeIndexPuffs


% Create filters for the ROIs data and the Num of Sinks and surges time-series 
% These filter give me the row in ROIs_Traces and NumOngoing (Sinks and surges) that belong to one of the 8 groups listed above
Filters_ROIsandEvents=cell(4,length(Drugs_unique)*length(Conditions_unique)*length(StimCond_unique));
counter=1;
for i=1:length(Drugs_unique)

    wrapper = @(x) strcmp(x, Drugs_unique{i});
    TempROIsandEventsIndexDrugID=cell2mat(cellfun(wrapper,ROIs_Traces(:,4), 'UniformOutput', false));
    
    for k=1:length(Conditions_unique)

        wrapper = @(x) strcmp(x, Conditions_unique{k});
        TempROIsandEventsIndexCondition=cell2mat(cellfun(wrapper,ROIs_Traces(:,3), 'UniformOutput', false));
       
        for g=1:length(StimCond_unique)

            TempROIsandEventsIndexPuffs=cell2mat(ROIs_Traces(:,5))==StimCond_unique(g);

            Filters_ROIsandEvents{4,counter}=logical(TempROIsandEventsIndexDrugID.*TempROIsandEventsIndexCondition.*TempROIsandEventsIndexPuffs);

            Filters_ROIsandEvents{1,counter}= Drugs_unique{i};
            Filters_ROIsandEvents{2,counter}= Conditions_unique{k};

            if StimCond_unique(g)                   
                Filters_ROIsandEvents{3,counter}='With_stimulation';
            else
                Filters_ROIsandEvents{3,counter}='Without_stimulation';        
            end

            counter=counter+1;
        end

    end

end

clear i k g counter wrapper TempROIsandEventsIndexDrugID TempROIsandEventsIndexCondition TempROIsandEventsIndexPuffs

%% (6) Now using the filters generated in (5) divide the dataset and prepare the tables for export

FilteredSinkData_mice=cell(1,size(FiltersOxySinksMetrics,2));
FilteredSurgeData_mice=cell(1,size(FiltersOxySinksMetrics,2));

FilteredSinkData=cell(10,size(FiltersOxySinksMetrics,2));
FilteredSurgeData=cell(10,size(FiltersOxySinksMetrics,2));

for i=1:size(FiltersOxySinksMetrics,2)

    FilteredSinkData_mice{1,i}=Table_OxygenSinks_OutCombo.Mouse(FiltersOxySinksMetrics{4,i});
    
    FilteredSinkData{1,i}=Table_OxygenSinks_OutCombo.MeanOxySinkArea_um(FiltersOxySinksMetrics{4,i});
    FilteredSinkData{2,i}=Table_OxygenSinks_OutCombo.MeanOxySinkFilledArea_um(FiltersOxySinksMetrics{4,i});
    FilteredSinkData{3,i}=Table_OxygenSinks_OutCombo.MeanOxySinkDiameter_um(FiltersOxySinksMetrics{4,i});
    FilteredSinkData{4,i}=Table_OxygenSinks_OutCombo.MeanOxySinkPerimeter_um(FiltersOxySinksMetrics{4,i});
    FilteredSinkData{5,i}=AdditionalOxySinkMetrics.OxySinkArea_Norm(FiltersOxySinksMetrics{4,i});
    FilteredSinkData{6,i}=AdditionalOxySinkMetrics.OxySinkArea_Filled_Norm(FiltersOxySinksMetrics{4,i});
    FilteredSinkData{7,i}=AdditionalOxySinkMetrics.NumOxySinkEvents_Norm(FiltersOxySinksMetrics{4,i});
    FilteredSinkData{8,i}=AdditionalOxySinkMetrics.MeanOxySinkEvent_Duration(FiltersOxySinksMetrics{4,i});
    FilteredSinkData{9,i}=AdditionalOxySinkMetrics.MeanOxySinkEvent_NormAmp(FiltersOxySinksMetrics{4,i});
    FilteredSinkData{10,i}=AdditionalOxySinkMetrics.MeanOxySinkEvent_SizeModulation(FiltersOxySinksMetrics{4,i});

    FilteredSurgeData_mice{1,i}=Table_OxygenSurges_OutCombo.Mouse(FiltersOxySurgesMetrics{4,i});

    FilteredSurgeData{1,i}=Table_OxygenSurges_OutCombo.MeanOxySurgeArea_um(FiltersOxySurgesMetrics{4,i});
    FilteredSurgeData{2,i}=Table_OxygenSurges_OutCombo.MeanOxySurgeFilledArea_um(FiltersOxySurgesMetrics{4,i});
    FilteredSurgeData{3,i}=Table_OxygenSurges_OutCombo.MeanOxySurgeDiameter_um(FiltersOxySurgesMetrics{4,i});
    FilteredSurgeData{4,i}=Table_OxygenSurges_OutCombo.MeanOxySurgePerimeter_um(FiltersOxySurgesMetrics{4,i});
    FilteredSurgeData{5,i}=AdditionalOxySurgeMetrics.OxySurgeArea_Norm(FiltersOxySurgesMetrics{4,i});
    FilteredSurgeData{6,i}=AdditionalOxySurgeMetrics.OxySurgeArea_Filled_Norm(FiltersOxySurgesMetrics{4,i});
    FilteredSurgeData{7,i}=AdditionalOxySurgeMetrics.NumOxySurgeEvents_Norm(FiltersOxySurgesMetrics{4,i});
    FilteredSurgeData{8,i}=AdditionalOxySurgeMetrics.MeanOxySurgeEvent_Duration(FiltersOxySurgesMetrics{4,i});
    FilteredSurgeData{9,i}=AdditionalOxySurgeMetrics.MeanOxySurgeEvent_NormAmp(FiltersOxySurgesMetrics{4,i});
    FilteredSurgeData{10,i}=AdditionalOxySurgeMetrics.MeanOxySurgeEvent_SizeModulation(FiltersOxySurgesMetrics{4,i});

end

clear i
%%
GroupHeaders=FiltersOxySinksMetrics(1:3,:);

ExportReady=cell(10,2);
for i=1:size(FilteredSinkData,1) %for every metric

    %create temp cell arrays that have the length of the recording with the most sinks/surges and coloumn equal with the number of groups*the maximum number of mice that exist in a group    
    %this prealocation ensures that there will definitelly be enough space
    %in all dimensions to wrangle and rearrange data
    tempSinkmetric=cell(max(cellfun(@length,FilteredSinkData),[],'all'),length(Drugs_unique)*length(Conditions_unique)*length(StimCond_unique)*max(cellfun(@length,FiltersOxySinksMetrics(5,:)),[],'all')); 
    
    tempSurgemetric=cell(max(cellfun(@length,FilteredSurgeData),[],'all'),length(Drugs_unique)*length(Conditions_unique)*length(StimCond_unique)*max(cellfun(@length,FiltersOxySurgesMetrics(5,:)),[],'all'));
    
    counter=1;
    for k=1:size(FilteredSinkData,2) %go through the different groups 
    
        if ~isempty(FilteredSinkData{i,k}) %if sink data for that group exists
                    
            for jj=1:length(FiltersOxySinksMetrics{5,k})
            
                wrapper = @(x) strcmp(x, FiltersOxySinksMetrics{5,k}{jj});
                tempfiltermouse=cell2mat(cellfun(wrapper,FilteredSinkData_mice{1,k}, 'UniformOutput', false));
                
                tempSinkmetric(1,counter)= FiltersOxySinksMetrics(1,k);
                tempSinkmetric(2,counter)= FiltersOxySinksMetrics(2,k);
                tempSinkmetric(3,counter)= FiltersOxySinksMetrics(3,k);
                tempSinkmetric(4,counter)= FiltersOxySinksMetrics{5,k}(jj);                                          
                tempSinkmetric(5:4+length(FilteredSinkData{i,k}(tempfiltermouse)),counter)=num2cell(FilteredSinkData{i,k}(tempfiltermouse));

                counter=counter+1;
            end

        end

    end

    %clean up empty cells 
    tempSinkmetric=tempSinkmetric(any(~(cellfun(@isempty,tempSinkmetric)),2),:);
    tempSinkmetric=tempSinkmetric(:,any(~(cellfun(@isempty,tempSinkmetric)),1));

    %repreat for surges
    counter=1;
    for k=1:size(FilteredSurgeData,2) %go through the different groups

        if ~isempty(FilteredSurgeData{i,k}) %if surge data for that group exists   
            
            for jj=1:length(FiltersOxySurgesMetrics{5,k})

                wrapper = @(x) strcmp(x, FiltersOxySurgesMetrics{5,k}{jj});
                tempfiltermouse=cell2mat(cellfun(wrapper,FilteredSurgeData_mice{1,k}, 'UniformOutput', false));
                            
                tempSurgemetric(1,counter)= FiltersOxySurgesMetrics(1,k);
                tempSurgemetric(2,counter)= FiltersOxySurgesMetrics(2,k);
                tempSurgemetric(3,counter)= FiltersOxySurgesMetrics(3,k);
                tempSurgemetric(4,counter)= FiltersOxySurgesMetrics{5,k}(jj);                                          
                tempSurgemetric(5:4+length(FilteredSurgeData{i,k}(tempfiltermouse)),counter)=num2cell(FilteredSurgeData{i,k}(tempfiltermouse));
               
                counter=counter+1;
            end
        
        end
 
    end

    %clean up empty cells 
    tempSurgemetric=tempSurgemetric(any(~(cellfun(@isempty,tempSurgemetric)),2),:);
    tempSurgemetric=tempSurgemetric(:,any(~(cellfun(@isempty,tempSurgemetric)),1));

    ExportReady{i,1}=cell2table(tempSinkmetric);    
    ExportReady{i,2}=cell2table(tempSurgemetric);

end
clear i k tempSinkmetric tempSurgemetric jj tempfiltermouse

%%
if strcmp(iOS, 'BLI')
Pooled_Traces=cell(14,size(Filters_ROIsandEvents,2));
Pooled_Traces(1:3,:)=Filters_ROIsandEvents(1:3,:);
Pooled_Traces_mice=cell(14,size(Filters_ROIsandEvents,2));
Pooled_Traces_mice(1:3,:)=Filters_ROIsandEvents(1:3,:);

for g=1:size(Filters_ROIsandEvents,2) %go through the different experimental groups
    %index the behaviour with the logical vector indacating which recordings belong to that group
    tempbehcombo=Behaviour_data_combo(Filters_ROIsandEvents{4,g},:); %get the behavioural tracking for the recordings of that group
    tempPuffs2use=Puffs_2Use(Filters_ROIsandEvents{4,g});
    tempSampleFs=SampleFs(Filters_ROIsandEvents{4,g},:); %get the imaging sampling frequencies for the recordings of that group
    tempMice=Mice(Filters_ROIsandEvents{4,g},:);    
    %creating temporary data based on the filterring
    tempROIs_Traces=ROIs_Traces(Filters_ROIsandEvents{4,g},:);
    tempNumOngoingOxysinks=NumOngoingOxysinks(Filters_ROIsandEvents{4,g},:);          
    tempTotalSinkArea_Norm=TotalSinkArea_Norm(Filters_ROIsandEvents{4,g},:);                   
    tempNumOngoingOxysurges=NumOngoingOxysurges(Filters_ROIsandEvents{4,g},:);
    tempTotalSurgeArea=TotalSurgeArea(Filters_ROIsandEvents{4,g},:);     
    tempTotalSinkArea_um=TotalSinkArea_um(Filters_ROIsandEvents{4,g},:);
    
for i=1:size(tempbehcombo,1)

    if ~isempty(tempbehcombo{i,9})

        %get the logical vector in a temp variable
        pfflogical=tempbehcombo{i,9};
        %label all the individual puffs    
        puff_num=bwlabel(pfflogical);
        
        %creating the temp arrays that have rows equal to the number of
        %puffs and columns equal to the number of sampling points to
        %capture 30s before the start of puffs 60s of puffs and 30s after
        %there is one more column

        %this need to change, I need to take that from the Twindows rather
        %than making this here. I also think that this has been coded again
        %in the event-based analysis

        pooled_ROImean=nan(max(puff_num), length(-30*tempSampleFs{i}:60*tempSampleFs{i}));
        pooled_ROIcovcoef=nan(max(puff_num), length(-30*tempSampleFs{i}:60*tempSampleFs{i}));
        pooled_ROIentro=nan(max(puff_num), length(-30*tempSampleFs{i}:60*tempSampleFs{i}));
        pooled_NumOxySi=nan(max(puff_num), length(-30*tempSampleFs{i}:60*tempSampleFs{i}));
        pooled_AreaOxySi=nan(max(puff_num), length(-30*tempSampleFs{i}:60*tempSampleFs{i}));
        pooled_NumOxySu=nan(max(puff_num), length(-30*tempSampleFs{i}:60*tempSampleFs{i}));
        pooled_AreaOxySu=nan(max(puff_num), length(-30*tempSampleFs{i}:60*tempSampleFs{i}));
        pooled_ROIdiffmean=nan(max(puff_num), length(-30*tempSampleFs{i}:60*tempSampleFs{i}));
        pooled_ROIdiffcovcoef=nan(max(puff_num), length(-30*tempSampleFs{i}:60*tempSampleFs{i}));
        pooled_ROIdiffentro=nan(max(puff_num), length(-30*tempSampleFs{i}:60*tempSampleFs{i}));
        pooled_AreaOxySi_um=nan(max(puff_num), length(-30*tempSampleFs{i}:60*tempSampleFs{i}));


        %%%%%

        %this can be only one columns table that will be used 100 lines
        %below when converting to table...in this case I will covert the
        %traces to table and ADD this column table with them

        %in this column I add the mouse number to make things easier for
        %plotting and analysis in xl

        Title="Mouse";
        sz=[max(puff_num),1];
        varTypes="string";
        temp_mice=table('Size',sz,'VariableTypes',varTypes,'VariableNames',Title);

        temp_mice.Mouse(:)=tempMice{i};
        %%%%%%

        %go through each puff and add it to the plot using the patch function 
        for k=1:max(puff_num)
            if ~isempty(intersect(k,str2num(tempPuffs2use{i}))) %if that puff is in the list of included puffs
                
            %find the start of the puff and convert to index that makes                
            %sense for the sampling rate of imaging data/traces               
            PuffStart=nearest(find(puff_num==k,1,'first')/PuffsSFs)*tempSampleFs{i}; 
            %go through the different traces and collect the trace segments that correspond to the window around the start of the stimulation       
            
            if PuffStart-30*tempSampleFs{i} >0 && PuffStart+60*tempSampleFs{i}<length(tempROIs_Traces{i,7})+1
            
            pooled_ROImean(k,1:end)=tempROIs_Traces{i,7}(PuffStart-30*tempSampleFs{i}:PuffStart+60*tempSampleFs{i});
            pooled_ROIcovcoef(k,1:end)=tempROIs_Traces{i,8}(PuffStart-30*tempSampleFs{i}:PuffStart+60*tempSampleFs{i});
            pooled_ROIentro(k,1:end)=tempROIs_Traces{i,9}(PuffStart-30*tempSampleFs{i}:PuffStart+60*tempSampleFs{i});
            pooled_NumOxySi(k,1:end)=tempNumOngoingOxysinks{i,6}(PuffStart-30*tempSampleFs{i}:PuffStart+60*tempSampleFs{i});
            pooled_AreaOxySi(k,1:end)=tempTotalSinkArea_Norm{i,6}(PuffStart-30*tempSampleFs{i}:PuffStart+60*tempSampleFs{i});            
            pooled_NumOxySu(k,1:end)=tempNumOngoingOxysurges{i,6}(PuffStart-30*tempSampleFs{i}:PuffStart+60*tempSampleFs{i});
            pooled_AreaOxySu(k,1:end)=tempTotalSurgeArea{i,6}(PuffStart-30*tempSampleFs{i}:PuffStart+60*tempSampleFs{i});
            pooled_ROIdiffmean(k,1:end)=tempROIs_Traces{i,10}(PuffStart-30*tempSampleFs{i}:PuffStart+60*tempSampleFs{i});
            pooled_ROIdiffcovcoef(k,1:end)=tempROIs_Traces{i,11}(PuffStart-30*tempSampleFs{i}:PuffStart+60*tempSampleFs{i});
            pooled_ROIdiffentro(k,1:end)=tempROIs_Traces{i,12}(PuffStart-30*tempSampleFs{i}:PuffStart+60*tempSampleFs{i});
            pooled_AreaOxySi_um(k,1:end)=tempTotalSinkArea_um{i,6}(PuffStart-30*tempSampleFs{i}:PuffStart+60*tempSampleFs{i});

            elseif PuffStart-30*tempSampleFs{i} <1

            pooled_ROImean(k,end-(length(tempROIs_Traces{i,7}(1:PuffStart+60*tempSampleFs{i}))-1):end)=tempROIs_Traces{i,7}(1:PuffStart+60*tempSampleFs{i});
            pooled_ROIcovcoef(k,end-(length(tempROIs_Traces{i,8}(1:PuffStart+60*tempSampleFs{i}))-1):end)=tempROIs_Traces{i,8}(1:PuffStart+60*tempSampleFs{i});
            pooled_ROIentro(k,end-(length(tempROIs_Traces{i,9}(1:PuffStart+60*tempSampleFs{i}))-1):end)=tempROIs_Traces{i,9}(1:PuffStart+60*tempSampleFs{i});
            pooled_NumOxySi(k,end-(length(tempNumOngoingOxysinks{i,6}(1:PuffStart+60*tempSampleFs{i}))-1):end)=tempNumOngoingOxysinks{i,6}(1:PuffStart+60*tempSampleFs{i});
            pooled_AreaOxySi(k,end-(length(tempTotalSinkArea_Norm{i,6}(1:PuffStart+60*tempSampleFs{i}))-1):end)=tempTotalSinkArea_Norm{i,6}(1:PuffStart+60*tempSampleFs{i});            
            pooled_NumOxySu(k,end-(length(tempNumOngoingOxysurges{i,6}(1:PuffStart+60*tempSampleFs{i}))-1):end)=tempNumOngoingOxysurges{i,6}(1:PuffStart+60*tempSampleFs{i});
            pooled_AreaOxySu(k,end-(length(tempTotalSurgeArea{i,6}(1:PuffStart+60*tempSampleFs{i}))-1):end)=tempTotalSurgeArea{i,6}(1:PuffStart+60*tempSampleFs{i});
            pooled_ROIdiffmean(k,end-(length(tempROIs_Traces{i,10}(1:PuffStart+60*tempSampleFs{i}))-1):end)=tempROIs_Traces{i,10}(1:PuffStart+60*tempSampleFs{i});
            pooled_ROIdiffcovcoef(k,end-(length(tempROIs_Traces{i,11}(1:PuffStart+60*tempSampleFs{i}))-1):end)=tempROIs_Traces{i,11}(1:PuffStart+60*tempSampleFs{i});
            pooled_ROIdiffentro(k,end-(length(tempROIs_Traces{i,12}(1:PuffStart+60*tempSampleFs{i}))-1):end)=tempROIs_Traces{i,12}(1:PuffStart+60*tempSampleFs{i});
            pooled_AreaOxySi_um(k,end-(length(tempTotalSinkArea_um{i,6}(1:PuffStart+60*tempSampleFs{i}))-1):end)=tempTotalSinkArea_um{i,6}(1:PuffStart+60*tempSampleFs{i});    

            elseif PuffStart+60*tempSampleFs{i}>length(tempROIs_Traces{i,7})

            pooled_ROImean(k,1:length(tempROIs_Traces{i,7}(PuffStart-30*tempSampleFs{i}:end)))=tempROIs_Traces{i,7}(PuffStart-30*tempSampleFs{i}:end);
            pooled_ROIcovcoef(k,1:length(tempROIs_Traces{i,8}(PuffStart-30*tempSampleFs{i}:end)))=tempROIs_Traces{i,8}(PuffStart-30*tempSampleFs{i}:end);
            pooled_ROIentro(k,1:length(tempROIs_Traces{i,9}(PuffStart-30*tempSampleFs{i}:end)))=tempROIs_Traces{i,9}(PuffStart-30*tempSampleFs{i}:end);
            pooled_NumOxySi(k,1:length(tempNumOngoingOxysinks{i,6}(PuffStart-30*tempSampleFs{i}:end)))=tempNumOngoingOxysinks{i,6}(PuffStart-30*tempSampleFs{i}:end);
            pooled_AreaOxySi(k,1:length(tempTotalSinkArea_Norm{i,6}(PuffStart-30*tempSampleFs{i}:end)))=tempTotalSinkArea_Norm{i,6}(PuffStart-30*tempSampleFs{i}:end);            
            pooled_NumOxySu(k,1:length(tempNumOngoingOxysurges{i,6}(PuffStart-30*tempSampleFs{i}:end)))=tempNumOngoingOxysurges{i,6}(PuffStart-30*tempSampleFs{i}:end);
            pooled_AreaOxySu(k,1:length(tempTotalSurgeArea{i,6}(PuffStart-30*tempSampleFs{i}:end)))=tempTotalSurgeArea{i,6}(PuffStart-30*tempSampleFs{i}:end);
            pooled_ROIdiffmean(k,1:length(tempROIs_Traces{i,10}(PuffStart-30*tempSampleFs{i}:end)))=tempROIs_Traces{i,10}(PuffStart-30*tempSampleFs{i}:end);
            pooled_ROIdiffcovcoef(k,1:length(tempROIs_Traces{i,11}(PuffStart-30*tempSampleFs{i}:end)))=tempROIs_Traces{i,11}(PuffStart-30*tempSampleFs{i}:end);
            pooled_ROIdiffentro(k,1:length(tempROIs_Traces{i,12}(PuffStart-30*tempSampleFs{i}:end)))=tempROIs_Traces{i,12}(PuffStart-30*tempSampleFs{i}:end);
            pooled_AreaOxySi_um(k,1:length(tempTotalSinkArea_um{i,6}(PuffStart-30*tempSampleFs{i}:end)))=tempTotalSinkArea_um{i,6}(PuffStart-30*tempSampleFs{i}:end);
            end
            
            end

        end


        %removing the rows for the Puffs that were not included in the
        %pooling
        pooled_ROImean= pooled_ROImean(ismember(1:max(puff_num),str2num(tempPuffs2use{i})),:);       
        pooled_ROIcovcoef= pooled_ROIcovcoef(ismember(1:max(puff_num),str2num(tempPuffs2use{i})),:);         
        pooled_ROIentro= pooled_ROIentro(ismember(1:max(puff_num),str2num(tempPuffs2use{i})),:);       
        pooled_NumOxySi= pooled_NumOxySi(ismember(1:max(puff_num),str2num(tempPuffs2use{i})),:);         
        pooled_AreaOxySi= pooled_AreaOxySi(ismember(1:max(puff_num),str2num(tempPuffs2use{i})),:);          
        pooled_NumOxySu= pooled_NumOxySu(ismember(1:max(puff_num),str2num(tempPuffs2use{i})),:);            
        pooled_AreaOxySu= pooled_AreaOxySu(ismember(1:max(puff_num),str2num(tempPuffs2use{i})),:);          
        pooled_ROIdiffmean= pooled_ROIdiffmean(ismember(1:max(puff_num),str2num(tempPuffs2use{i})),:);             
        pooled_ROIdiffcovcoef= pooled_ROIdiffcovcoef(ismember(1:max(puff_num),str2num(tempPuffs2use{i})),:);            
        pooled_ROIdiffentro= pooled_ROIdiffentro(ismember(1:max(puff_num),str2num(tempPuffs2use{i})),:);
        pooled_AreaOxySi_um= pooled_AreaOxySi_um(ismember(1:max(puff_num),str2num(tempPuffs2use{i})),:);

        temp_mice=temp_mice.Mouse((ismember(1:max(puff_num),str2num(tempPuffs2use{i})))');

        Pooled_Traces{4,g}=vertcat(Pooled_Traces{4,g},pooled_ROImean);
        Pooled_Traces{5,g}=vertcat(Pooled_Traces{5,g},pooled_ROIcovcoef);
        Pooled_Traces{6,g}=vertcat(Pooled_Traces{6,g},pooled_ROIentro);
        Pooled_Traces{7,g}=vertcat(Pooled_Traces{7,g},pooled_NumOxySi);
        Pooled_Traces{8,g}=vertcat(Pooled_Traces{8,g},pooled_AreaOxySi);
        Pooled_Traces{9,g}=vertcat(Pooled_Traces{9,g},pooled_NumOxySu);
        Pooled_Traces{10,g}=vertcat(Pooled_Traces{10,g},pooled_AreaOxySu);
        Pooled_Traces{11,g}=vertcat(Pooled_Traces{11,g},pooled_ROIdiffmean);
        Pooled_Traces{12,g}=vertcat(Pooled_Traces{12,g},pooled_ROIdiffcovcoef);
        Pooled_Traces{13,g}=vertcat(Pooled_Traces{13,g},pooled_ROIdiffentro);
        Pooled_Traces{14,g}=vertcat(Pooled_Traces{14,g},pooled_AreaOxySi_um);

        Pooled_Traces_mice{4,g}=vertcat(Pooled_Traces_mice{4,g},temp_mice);
        Pooled_Traces_mice{5,g}=vertcat(Pooled_Traces_mice{5,g},temp_mice);
        Pooled_Traces_mice{6,g}=vertcat(Pooled_Traces_mice{6,g},temp_mice);
        Pooled_Traces_mice{7,g}=vertcat(Pooled_Traces_mice{7,g},temp_mice);
        Pooled_Traces_mice{8,g}=vertcat(Pooled_Traces_mice{8,g},temp_mice);
        Pooled_Traces_mice{9,g}=vertcat(Pooled_Traces_mice{9,g},temp_mice);
        Pooled_Traces_mice{10,g}=vertcat(Pooled_Traces_mice{10,g},temp_mice);
        Pooled_Traces_mice{11,g}=vertcat(Pooled_Traces_mice{11,g},temp_mice);
        Pooled_Traces_mice{12,g}=vertcat(Pooled_Traces_mice{12,g},temp_mice);
        Pooled_Traces_mice{13,g}=vertcat(Pooled_Traces_mice{13,g},temp_mice);
        Pooled_Traces_mice{14,g}=vertcat(Pooled_Traces_mice{14,g},temp_mice);

    end

end
end

end
clear g i k tempbehcombo tempSampleFs tempMice pfflogical puff_num PuffStart pooled_ROImean pooled_ROIcovcoef pooled_ROIentro pooled_NumOxySi pooled_AreaOxySi pooled_NumOxySu pooled_AreaOxySu...
  pooled_ROIdiffmean pooled_ROIdiffcovcoef pooled_ROIdiffentro tempPuffs2use tempROIs_Traces tempNumOngoingOxysinks tempTotalSinkArea_Norm tempNumOngoingOxysurges tempTotalSurgeArea tempotalSinkArea_um...
  temp_mice Title sz varTypes

%%
if strcmp(iOS, 'BLI')
for i=1:size(Pooled_Traces,2) %for each experimental group
    for j=4:size(Pooled_Traces,1) %for each pooled metric

        if ~isempty(Pooled_Traces{j,i})

            % divide the values with the average of baseline -30s to puff
            % start
            Pooled_Traces{j,i}=Pooled_Traces{j,i}(:,1:end)./mean(Pooled_Traces{j,i}(:,1:31),2,'omitnan');

            %adjusting the value at the start of the puff to 0
            Pooled_Traces{j,i}=Pooled_Traces{j,i}(:,1:end)-Pooled_Traces{j,i}(:,31);

            %converting to table to make it easier to export
            Pooled_Traces{j,i}=array2table(Pooled_Traces{j,i});
  
            Title=["Mouse"];
            varTypes="string";  
            sz=[size(Pooled_Traces_mice{j,i},1),1];
            temp_mice=table('Size',sz,'VariableTypes',varTypes,'VariableNames',Title);
            temp_mice.Mouse(:)=Pooled_Traces_mice{j,i};

            Pooled_Traces{j,i}=[temp_mice Pooled_Traces{j,i}];
            % this is where I need to bring the mice numbers in from the
            % table

            newNames = append(string(-30:60),"sec");
            allVars = 2:width(Pooled_Traces{j,i});
            Pooled_Traces{j,i}=renamevars(Pooled_Traces{j,i},allVars,newNames);

        end

    end
end
end
clear i j newNames allVars temp_mice

%%
Titles1={'DrugID';'Condition';'Stimulation';'ROI Mean';'ROI CovCoef';'ROI Entropy';'OxySinks HowMany';'OxySinksArea_Norm';'OxySurgesHowMany';'OxySurgesArea';...
    'ROIdiff Mean';'ROIdiff CovCoef';'ROIdiff Entropy';'OxySinksArea_um'};
Titles2={'DrugID';'Condition';'Stimulation';'ROI Mean vs CovCoef';'ROI Mean vs Entropy';'CovCoef vs Entropy';'ROIdiff Mean vs diffCovCoef';'ROIdiff Mean vs diffEntropy';...
    'diffCovCoef vs diffEntropy';'OxySinksArea vs CovCoef';'OxySinksArea vs Entropy';'OxySurgesArea vs CovCoef';'OxySurgesArea vs Entropy'};
%%
if strcmp(iOS, 'BLI')
%I add the names of the trace types to make my life easier when exporting
Pooled_Traces=[Titles1,Pooled_Traces];
end
%%
ExportTraces=GroupHeaders;
ExportTraceCorrs=GroupHeaders;
for i=1:size(ExportTraces,2) %go through the different expeimental groups

    if any(Filters_ROIsandEvents{4,i}) %if there are recordings that fall under that experimental group (one of the eight)

        %the number of traces for each experimental group should be the same for all metrics because they are dependent on how many recordings belong to each group

        tempcell=cell(4,sum(Filters_ROIsandEvents{4,i})); %making a headers cell that will use for the rest
        tempcell(1,:)=ExportTraces(1,i);
        tempcell(2,:)=ExportTraces(2,i);
        tempcell(3,:)=ExportTraces(3,i);

        tempcell(4,:)=ROIs_Traces(Filters_ROIsandEvents{4,i},2)';

        %as you can see need to do some acrobatics to rearrange the data
        %into tables that will alow me easy data exporting    
        
        % I first find what is the longer recording
        wrapper = @(x) size(x, 2);     

        if strcmp(iOS, 'BLI')
        maxcols = max(cellfun(wrapper,ROIs_Traces(Filters_ROIsandEvents{4,i},7)));  %get the number of columns of the widest array
        padded = cellfun(@(m) [m, nan(size(m, 1), maxcols - size(m, 2))], ROIs_Traces(Filters_ROIsandEvents{4,i},7), 'UniformOutput', false);  %pad each array
        ExportTraces{4,i}=cell2table([tempcell;num2cell((cell2mat(padded))')]); 
        
        maxcols = max(cellfun(wrapper,ROIs_Traces(Filters_ROIsandEvents{4,i},8)));  %get the number of columns of the widest array
        padded = cellfun(@(m) [m, nan(size(m, 1), maxcols - size(m, 2))], ROIs_Traces(Filters_ROIsandEvents{4,i},8), 'UniformOutput', false);  %pad each array    
        ExportTraces{5,i}=cell2table([tempcell;num2cell((cell2mat(padded))')]);     
        
        maxcols = max(cellfun(wrapper,ROIs_Traces(Filters_ROIsandEvents{4,i},9)));  %get the number of columns of the widest array
        padded = cellfun(@(m) [m, nan(size(m, 1), maxcols - size(m, 2))], ROIs_Traces(Filters_ROIsandEvents{4,i},9), 'UniformOutput', false);  %pad each array        
        ExportTraces{6,i}=cell2table([tempcell;num2cell((cell2mat(padded))')]);
        

        maxcols = max(cellfun(wrapper,ROIs_Traces(Filters_ROIsandEvents{4,i},10)));  %get the number of columns of the widest array
        padded = cellfun(@(m) [m, nan(size(m, 1), maxcols - size(m, 2))], ROIs_Traces(Filters_ROIsandEvents{4,i},10), 'UniformOutput', false);  %pad each array        
        ExportTraces{11,i}=cell2table([tempcell;num2cell((cell2mat(padded))')]);
        
        maxcols = max(cellfun(wrapper,ROIs_Traces(Filters_ROIsandEvents{4,i},11)));  %get the number of columns of the widest array
        padded = cellfun(@(m) [m, nan(size(m, 1), maxcols - size(m, 2))], ROIs_Traces(Filters_ROIsandEvents{4,i},11), 'UniformOutput', false);  %pad each array        
        ExportTraces{12,i}=cell2table([tempcell;num2cell((cell2mat(padded))')]);
        
        maxcols = max(cellfun(wrapper,ROIs_Traces(Filters_ROIsandEvents{4,i},12)));  %get the number of columns of the widest array
        padded = cellfun(@(m) [m, nan(size(m, 1), maxcols - size(m, 2))], ROIs_Traces(Filters_ROIsandEvents{4,i},12), 'UniformOutput', false);  %pad each array        
        ExportTraces{13,i}=cell2table([tempcell;num2cell((cell2mat(padded))')]);

        end


        maxcols = max(cellfun(wrapper,NumOngoingOxysinks(Filters_ROIsandEvents{4,i},6)));  %get the number of columns of the widest array
        padded = cellfun(@(m) [m, nan(size(m, 1), maxcols - size(m, 2))], NumOngoingOxysinks(Filters_ROIsandEvents{4,i},6), 'UniformOutput', false);  %pad each array        
        ExportTraces{7,i}=cell2table([tempcell;num2cell((cell2mat(padded))')]);
        
        maxcols = max(cellfun(wrapper,TotalSinkArea_Norm(Filters_ROIsandEvents{4,i},6)));  %get the number of columns of the widest array
        padded = cellfun(@(m) [m, nan(size(m, 1), maxcols - size(m, 2))], TotalSinkArea_Norm(Filters_ROIsandEvents{4,i},6), 'UniformOutput', false);  %pad each array        
        ExportTraces{8,i}=cell2table([tempcell;num2cell((cell2mat(padded))')]);
        
        maxcols = max(cellfun(wrapper,NumOngoingOxysurges(Filters_ROIsandEvents{4,i},6)));  %get the number of columns of the widest array
        padded = cellfun(@(m) [m, nan(size(m, 1), maxcols - size(m, 2))], NumOngoingOxysurges(Filters_ROIsandEvents{4,i},6), 'UniformOutput', false);  %pad each array        
        ExportTraces{9,i}=cell2table([tempcell;num2cell((cell2mat(padded))')]);
        
        maxcols = max(cellfun(wrapper,TotalSurgeArea(Filters_ROIsandEvents{4,i},6)));  %get the number of columns of the widest array
        padded = cellfun(@(m) [m, nan(size(m, 1), maxcols - size(m, 2))], TotalSurgeArea(Filters_ROIsandEvents{4,i},6), 'UniformOutput', false);  %pad each array        
        ExportTraces{10,i}=cell2table([tempcell;num2cell((cell2mat(padded))')]);
        
        

        maxcols = max(cellfun(wrapper,TotalSinkArea_um(Filters_ROIsandEvents{4,i},6)));  %get the number of columns of the widest array
        padded = cellfun(@(m) [m, nan(size(m, 1), maxcols - size(m, 2))], TotalSinkArea_um(Filters_ROIsandEvents{4,i},6), 'UniformOutput', false);  %pad each array        
        ExportTraces{14,i}=cell2table([tempcell;num2cell((cell2mat(padded))')]);

     
        if strcmp(iOS, 'BLI')
        tempcell=cell(3,1); %making a headers cell that will use for the rest of the correlations 
        tempcell(1,1)=ExportTraceCorrs(1,i);
        tempcell(2,1)=ExportTraceCorrs(2,i);
        tempcell(3,1)=ExportTraceCorrs(3,i);
   
        ExportTraceCorrs{4,i}=cell2table([tempcell;num2cell(cell2mat(TraceCorrs(Filters_ROIsandEvents{4,i},6)))]);
        ExportTraceCorrs{5,i}=cell2table([tempcell;num2cell(cell2mat(TraceCorrs(Filters_ROIsandEvents{4,i},7)))]); 
        ExportTraceCorrs{6,i}=cell2table([tempcell;num2cell(cell2mat(TraceCorrs(Filters_ROIsandEvents{4,i},8)))]); 
        ExportTraceCorrs{7,i}=cell2table([tempcell;num2cell(cell2mat(TraceCorrs(Filters_ROIsandEvents{4,i},9)))]); 
        ExportTraceCorrs{8,i}=cell2table([tempcell;num2cell(cell2mat(TraceCorrs(Filters_ROIsandEvents{4,i},10)))]); 
        ExportTraceCorrs{9,i}=cell2table([tempcell;num2cell(cell2mat(TraceCorrs(Filters_ROIsandEvents{4,i},11)))]); 
        ExportTraceCorrs{10,i}=cell2table([tempcell;num2cell(cell2mat(TraceCorrs(Filters_ROIsandEvents{4,i},12)))]); 
        ExportTraceCorrs{11,i}=cell2table([tempcell;num2cell(cell2mat(TraceCorrs(Filters_ROIsandEvents{4,i},13)))]); 
        ExportTraceCorrs{12,i}=cell2table([tempcell;num2cell(cell2mat(TraceCorrs(Filters_ROIsandEvents{4,i},14)))]); 
        ExportTraceCorrs{13,i}=cell2table([tempcell;num2cell(cell2mat(TraceCorrs(Filters_ROIsandEvents{4,i},15)))]); 

        end
    else %if no recordings for that category, populate a trace with nans based on the length of the first recroding trace


        tempcell=cell(3,1);
        tempcell(1,1)=ExportTraces(1,i);
        tempcell(2,1)=ExportTraces(2,i);
        tempcell(3,1)=ExportTraces(3,i);

        if strcmp(iOS, 'BLI')
        ExportTraces{4,i}=cell2table([tempcell;num2cell(nan(size(ROIs_Traces{1,7},2),1))]);
        ExportTraces{5,i}=cell2table([tempcell;num2cell(nan(size(ROIs_Traces{1,8},2),1))]);
        ExportTraces{6,i}=cell2table([tempcell;num2cell(nan(size(ROIs_Traces{1,9},2),1))]);
        ExportTraces{11,i}=cell2table([tempcell;num2cell(nan(size(ROIs_Traces{1,10},2),1))]);
        ExportTraces{12,i}=cell2table([tempcell;num2cell(nan(size(ROIs_Traces{1,11},2),1))]);
        ExportTraces{13,i}=cell2table([tempcell;num2cell(nan(size(ROIs_Traces{1,12},2),1))]);
        end

        ExportTraces{7,i}=cell2table([tempcell;num2cell(nan(size(NumOngoingOxysinks{1,6},2),1))]);
        ExportTraces{8,i}=cell2table([tempcell;num2cell(nan(size(TotalSinkArea_Norm{1,6},2),1))]);
        ExportTraces{9,i}=cell2table([tempcell;num2cell(nan(size(NumOngoingOxysurges{1,6},2),1))]);
        ExportTraces{10,i}=cell2table([tempcell;num2cell(nan(size(TotalSurgeArea{1,6},2),1))]);

        ExportTraces{14,i}=cell2table([tempcell;num2cell(nan(size(TotalSinkArea_um{1,6},2),1))]);

        if strcmp(iOS, 'BLI')
        tempcell=cell(3,1); %making a headers cell that will use for the rest
        tempcell(1,1)=ExportTraceCorrs(1,i);
        tempcell(2,1)=ExportTraceCorrs(2,i);
        tempcell(3,1)=ExportTraceCorrs(3,i);


        ExportTraceCorrs{4,i}=cell2table([tempcell;num2cell(nan(size(TraceCorrs(1,6),1),1))]);
        ExportTraceCorrs{5,i}=cell2table([tempcell;num2cell(nan(size(TraceCorrs(1,7),1),1))]); 
        ExportTraceCorrs{6,i}=cell2table([tempcell;num2cell(nan(size(TraceCorrs(1,8),1),1))]); 
        ExportTraceCorrs{7,i}=cell2table([tempcell;num2cell(nan(size(TraceCorrs(1,9),1),1))]); 
        ExportTraceCorrs{8,i}=cell2table([tempcell;num2cell(nan(size(TraceCorrs(1,10),1),1))]); 
        ExportTraceCorrs{9,i}=cell2table([tempcell;num2cell(nan(size(TraceCorrs(1,11),1),1))]); 
        ExportTraceCorrs{10,i}=cell2table([tempcell;num2cell(nan(size(TraceCorrs(1,12),1),1))]); 
        ExportTraceCorrs{11,i}=cell2table([tempcell;num2cell(nan(size(TraceCorrs(1,13),1),1))]); 
        ExportTraceCorrs{12,i}=cell2table([tempcell;num2cell(nan(size(TraceCorrs(1,14),1),1))]); 
        ExportTraceCorrs{13,i}=cell2table([tempcell;num2cell(nan(size(TraceCorrs(1,15),1),1))]); 

        end

    end

end
ExportTraces=[Titles1,ExportTraces];
if strcmp(iOS, 'BLI')
ExportTraceCorrs=[Titles2,ExportTraceCorrs];
end

clear i tempcell wrapper maxcols padded


%%
%Here I filter the binned traces and prepare them for export
ExportLPawBinnedTraces=GroupHeaders;
ExportRPawBinnedTraces=GroupHeaders;
ExportPupilBinnedTraces=GroupHeaders;
if strcmp(iOS, 'BLI')
for i=1:size(ExportLPawBinnedTraces,2)
    %there have to be recordings which beloging to this group and all of them
    %should have LPaw tracking (this may cause problems when I have
    %recordings from older experiments that did not have behaviour
    if any(Filters_ROIsandEvents{4,i}) && all(~(cellfun(@isempty,BinsLPaw_10percentiles(Filters_ROIsandEvents{4,i},1)))) 

        tempcell=cell(4,10); %making a headers cell that will use for the rest
        tempcell(1,:)=ExportLPawBinnedTraces(1,i);
        tempcell(2,:)=ExportLPawBinnedTraces(2,i);
        tempcell(3,:)=ExportLPawBinnedTraces(3,i);
        tempcell(4,:)=num2cell(10:10:100); %names for the bins
   
        ExportLPawBinnedTraces{4,i}=cell2table([tempcell;num2cell(cell2mat(BinsLPaw_10percentiles(Filters_ROIsandEvents{4,i},1)))]); 
        ExportLPawBinnedTraces{5,i}=cell2table([tempcell;num2cell(cell2mat(BinsLPaw_10percentiles(Filters_ROIsandEvents{4,i},2)))]);     
        ExportLPawBinnedTraces{6,i}=cell2table([tempcell;num2cell(cell2mat(BinsLPaw_10percentiles(Filters_ROIsandEvents{4,i},3)))]);
        ExportLPawBinnedTraces{7,i}=cell2table([tempcell;num2cell(cell2mat(BinsLPaw_10percentiles(Filters_ROIsandEvents{4,i},4)))]);
        ExportLPawBinnedTraces{8,i}=cell2table([tempcell;num2cell(cell2mat(BinsLPaw_10percentiles(Filters_ROIsandEvents{4,i},5)))]);
        ExportLPawBinnedTraces{9,i}=cell2table([tempcell;num2cell(cell2mat(BinsLPaw_10percentiles(Filters_ROIsandEvents{4,i},6)))]);
        ExportLPawBinnedTraces{10,i}=cell2table([tempcell;num2cell(cell2mat(BinsLPaw_10percentiles(Filters_ROIsandEvents{4,i},7)))]);
        ExportLPawBinnedTraces{11,i}=cell2table([tempcell;num2cell(cell2mat(BinsLPaw_10percentiles(Filters_ROIsandEvents{4,i},8)))]);
        ExportLPawBinnedTraces{12,i}=cell2table([tempcell;num2cell(cell2mat(BinsLPaw_10percentiles(Filters_ROIsandEvents{4,i},9)))]);
        ExportLPawBinnedTraces{13,i}=cell2table([tempcell;num2cell(cell2mat(BinsLPaw_10percentiles(Filters_ROIsandEvents{4,i},10)))]);
        ExportLPawBinnedTraces{14,i}=cell2table([tempcell;num2cell(cell2mat(BinsLPaw_10percentiles(Filters_ROIsandEvents{4,i},11)))]);

    else

        tempcell=cell(4,10); %making a headers cell that will use for the rest
        tempcell(1,:)=ExportLPawBinnedTraces(1,i);
        tempcell(2,:)=ExportLPawBinnedTraces(2,i);
        tempcell(3,:)=ExportLPawBinnedTraces(3,i);
        tempcell(4,:)=num2cell(10:10:100); %names for the bins


        ExportLPawBinnedTraces{4,i}=cell2table([tempcell;num2cell(nan(1,10))]); 
        ExportLPawBinnedTraces{5,i}=cell2table([tempcell;num2cell(nan(1,10))]);    
        ExportLPawBinnedTraces{6,i}=cell2table([tempcell;num2cell(nan(1,10))]);
        ExportLPawBinnedTraces{7,i}=cell2table([tempcell;num2cell(nan(1,10))]);
        ExportLPawBinnedTraces{8,i}=cell2table([tempcell;num2cell(nan(1,10))]);
        ExportLPawBinnedTraces{9,i}=cell2table([tempcell;num2cell(nan(1,10))]);
        ExportLPawBinnedTraces{10,i}=cell2table([tempcell;num2cell(nan(1,10))]);
        ExportLPawBinnedTraces{11,i}=cell2table([tempcell;num2cell(nan(1,10))]);
        ExportLPawBinnedTraces{12,i}=cell2table([tempcell;num2cell(nan(1,10))]);
        ExportLPawBinnedTraces{13,i}=cell2table([tempcell;num2cell(nan(1,10))]);
        ExportLPawBinnedTraces{14,i}=cell2table([tempcell;num2cell(nan(1,10))]);

    end

    %repeat for Rpaw binned trace data
     if any(Filters_ROIsandEvents{4,i}) && all(~(cellfun(@isempty,BinsRPaw_10percentiles(Filters_ROIsandEvents{4,i},1)))) 

        tempcell=cell(4,10); %making a headers cell that will use for the rest
        tempcell(1,:)=ExportRPawBinnedTraces(1,i);
        tempcell(2,:)=ExportRPawBinnedTraces(2,i);
        tempcell(3,:)=ExportRPawBinnedTraces(3,i);
        tempcell(4,:)=num2cell(10:10:100); %names for the bins
   
        ExportRPawBinnedTraces{4,i}=cell2table([tempcell;num2cell(cell2mat(BinsRPaw_10percentiles(Filters_ROIsandEvents{4,i},1)))]); 
        ExportRPawBinnedTraces{5,i}=cell2table([tempcell;num2cell(cell2mat(BinsRPaw_10percentiles(Filters_ROIsandEvents{4,i},2)))]);     
        ExportRPawBinnedTraces{6,i}=cell2table([tempcell;num2cell(cell2mat(BinsRPaw_10percentiles(Filters_ROIsandEvents{4,i},3)))]);
        ExportRPawBinnedTraces{7,i}=cell2table([tempcell;num2cell(cell2mat(BinsRPaw_10percentiles(Filters_ROIsandEvents{4,i},4)))]);
        ExportRPawBinnedTraces{8,i}=cell2table([tempcell;num2cell(cell2mat(BinsRPaw_10percentiles(Filters_ROIsandEvents{4,i},5)))]);
        ExportRPawBinnedTraces{9,i}=cell2table([tempcell;num2cell(cell2mat(BinsRPaw_10percentiles(Filters_ROIsandEvents{4,i},6)))]);
        ExportRPawBinnedTraces{10,i}=cell2table([tempcell;num2cell(cell2mat(BinsRPaw_10percentiles(Filters_ROIsandEvents{4,i},7)))]);
        ExportRPawBinnedTraces{11,i}=cell2table([tempcell;num2cell(cell2mat(BinsRPaw_10percentiles(Filters_ROIsandEvents{4,i},8)))]);
        ExportRPawBinnedTraces{12,i}=cell2table([tempcell;num2cell(cell2mat(BinsRPaw_10percentiles(Filters_ROIsandEvents{4,i},9)))]);
        ExportRPawBinnedTraces{13,i}=cell2table([tempcell;num2cell(cell2mat(BinsRPaw_10percentiles(Filters_ROIsandEvents{4,i},10)))]);
        ExportRPawBinnedTraces{14,i}=cell2table([tempcell;num2cell(cell2mat(BinsRPaw_10percentiles(Filters_ROIsandEvents{4,i},11)))]);

    else

        tempcell=cell(4,10); %making a headers cell that will use for the rest
        tempcell(1,:)=ExportRPawBinnedTraces(1,i);
        tempcell(2,:)=ExportRPawBinnedTraces(2,i);
        tempcell(3,:)=ExportRPawBinnedTraces(3,i);
        tempcell(4,:)=num2cell(10:10:100); %names for the bins


        ExportRPawBinnedTraces{4,i}=cell2table([tempcell;num2cell(nan(1,10))]); 
        ExportRPawBinnedTraces{5,i}=cell2table([tempcell;num2cell(nan(1,10))]);    
        ExportRPawBinnedTraces{6,i}=cell2table([tempcell;num2cell(nan(1,10))]);
        ExportRPawBinnedTraces{7,i}=cell2table([tempcell;num2cell(nan(1,10))]);
        ExportRPawBinnedTraces{8,i}=cell2table([tempcell;num2cell(nan(1,10))]);
        ExportRPawBinnedTraces{9,i}=cell2table([tempcell;num2cell(nan(1,10))]);
        ExportRPawBinnedTraces{10,i}=cell2table([tempcell;num2cell(nan(1,10))]);
        ExportRPawBinnedTraces{11,i}=cell2table([tempcell;num2cell(nan(1,10))]);
        ExportRPawBinnedTraces{12,i}=cell2table([tempcell;num2cell(nan(1,10))]);
        ExportRPawBinnedTraces{13,i}=cell2table([tempcell;num2cell(nan(1,10))]);
        ExportRPawBinnedTraces{14,i}=cell2table([tempcell;num2cell(nan(1,10))]);

    end

    %repeat for pupil binned trace data
     if any(Filters_ROIsandEvents{4,i}) && all(~(cellfun(@isempty,BinsPupil_10percentiles(Filters_ROIsandEvents{4,i},1))))

        tempcell=cell(4,10); %making a headers cell that will use for the rest
        tempcell(1,:)=ExportPupilBinnedTraces(1,i);
        tempcell(2,:)=ExportPupilBinnedTraces(2,i);
        tempcell(3,:)=ExportPupilBinnedTraces(3,i);
        tempcell(4,:)=num2cell(10:10:100); %names for the bins
   
        ExportPupilBinnedTraces{4,i}=cell2table([tempcell;num2cell(cell2mat(BinsPupil_10percentiles(Filters_ROIsandEvents{4,i},1)))]); 
        ExportPupilBinnedTraces{5,i}=cell2table([tempcell;num2cell(cell2mat(BinsPupil_10percentiles(Filters_ROIsandEvents{4,i},2)))]);     
        ExportPupilBinnedTraces{6,i}=cell2table([tempcell;num2cell(cell2mat(BinsPupil_10percentiles(Filters_ROIsandEvents{4,i},3)))]);
        ExportPupilBinnedTraces{7,i}=cell2table([tempcell;num2cell(cell2mat(BinsPupil_10percentiles(Filters_ROIsandEvents{4,i},4)))]);
        ExportPupilBinnedTraces{8,i}=cell2table([tempcell;num2cell(cell2mat(BinsPupil_10percentiles(Filters_ROIsandEvents{4,i},5)))]);
        ExportPupilBinnedTraces{9,i}=cell2table([tempcell;num2cell(cell2mat(BinsPupil_10percentiles(Filters_ROIsandEvents{4,i},6)))]);
        ExportPupilBinnedTraces{10,i}=cell2table([tempcell;num2cell(cell2mat(BinsPupil_10percentiles(Filters_ROIsandEvents{4,i},7)))]);
        ExportPupilBinnedTraces{11,i}=cell2table([tempcell;num2cell(cell2mat(BinsPupil_10percentiles(Filters_ROIsandEvents{4,i},8)))]);
        ExportPupilBinnedTraces{12,i}=cell2table([tempcell;num2cell(cell2mat(BinsPupil_10percentiles(Filters_ROIsandEvents{4,i},9)))]);
        ExportPupilBinnedTraces{13,i}=cell2table([tempcell;num2cell(cell2mat(BinsPupil_10percentiles(Filters_ROIsandEvents{4,i},10)))]);
        ExportPupilBinnedTraces{14,i}=cell2table([tempcell;num2cell(cell2mat(BinsPupil_10percentiles(Filters_ROIsandEvents{4,i},11)))]);

    else

        tempcell=cell(4,10); %making a headers cell that will use for the rest
        tempcell(1,:)=ExportPupilBinnedTraces(1,i);
        tempcell(2,:)=ExportPupilBinnedTraces(2,i);
        tempcell(3,:)=ExportPupilBinnedTraces(3,i);
        tempcell(4,:)=num2cell(10:10:100); %names for the bins


        ExportPupilBinnedTraces{4,i}=cell2table([tempcell;num2cell(nan(1,10))]); 
        ExportPupilBinnedTraces{5,i}=cell2table([tempcell;num2cell(nan(1,10))]);    
        ExportPupilBinnedTraces{6,i}=cell2table([tempcell;num2cell(nan(1,10))]);
        ExportPupilBinnedTraces{7,i}=cell2table([tempcell;num2cell(nan(1,10))]);
        ExportPupilBinnedTraces{8,i}=cell2table([tempcell;num2cell(nan(1,10))]);
        ExportPupilBinnedTraces{9,i}=cell2table([tempcell;num2cell(nan(1,10))]);
        ExportPupilBinnedTraces{10,i}=cell2table([tempcell;num2cell(nan(1,10))]);
        ExportPupilBinnedTraces{11,i}=cell2table([tempcell;num2cell(nan(1,10))]);
        ExportPupilBinnedTraces{12,i}=cell2table([tempcell;num2cell(nan(1,10))]);
        ExportPupilBinnedTraces{13,i}=cell2table([tempcell;num2cell(nan(1,10))]);
        ExportPupilBinnedTraces{14,i}=cell2table([tempcell;num2cell(nan(1,10))]);

    end
 
end


ExportLPawBinnedTraces=[Titles1,ExportLPawBinnedTraces];
ExportRPawBinnedTraces=[Titles1,ExportRPawBinnedTraces];
ExportPupilBinnedTraces=[Titles1,ExportPupilBinnedTraces];
end

% Given that I am using percentiles for each recording the
% different bins will not correspond to the same movement speed/pupil size
% for the different recordings. Therefore when I average data for each of
% the 10 bins this maybe wrong.
% I can do a control analysis with the mm movement to see if mice move
% wildly differently. If not, binning the zscored data makes sense (..ish)
% Another way may be to bin the data for each recording at a mm and then
% pile up all sampling points from all animals and calculate
% average....This will be a nightmare to analyse statistically


%% Prepare the event-based analysis for extraction
if strcmp(iOS, 'BLI')
if exist('ROI_MeanTraceSnip_ManEvents', 'var')

   
    Export_ROI_MeanTraceSnip_ManEvents=[];
    
    for i=1:length(ROI_MeanTraceSnip_ManEvents)

        tempheads=(ROI_MeanTraceSnip_ManEvents{i}(:,1:5))';
    
        tempcell=(num2cell(cell2mat(ROI_MeanTraceSnip_ManEvents{i}(:,6))))';
        
        temp=[tempheads;tempcell];

        for g=1:size(ROI_MeanTraceSnip_ManEvents{i},1)
    
            stringVector{g} = strcat('Event',int2str(i),'Recording',int2str(g));

        end
        
        temp=cell2table(temp,"VariableNames",stringVector);

        Export_ROI_MeanTraceSnip_ManEvents=[Export_ROI_MeanTraceSnip_ManEvents,temp];

    end

end


clear i tempheads temp stringVector


if exist('ROI_MeanTraceSnip_LPawEvents', 'var')
        
    Export_ROI_MeanTraceSnip_LPawEvents=ROI_MeanTraceSnip_LPawEvents(:,1:5);

    tempcell=num2cell(cell2mat(ROI_MeanTraceSnip_LPawEvents(:,6)));

    Export_ROI_MeanTraceSnip_LPawEvents=[Export_ROI_MeanTraceSnip_LPawEvents,tempcell];

    Export_ROI_MeanTraceSnip_LPawEvents=cell2table(Export_ROI_MeanTraceSnip_LPawEvents);

end

if exist('ROI_MeanTraceSnip_RPawEvents', 'var')
      
    Export_ROI_MeanTraceSnip_RPawEvents=ROI_MeanTraceSnip_RPawEvents(:,1:5);

    tempcell=num2cell(cell2mat(ROI_MeanTraceSnip_RPawEvents(:,6)));

    Export_ROI_MeanTraceSnip_RPawEvents=[Export_ROI_MeanTraceSnip_RPawEvents,tempcell];

    Export_ROI_MeanTraceSnip_RPawEvents=cell2table(Export_ROI_MeanTraceSnip_RPawEvents);

end

if exist('ROI_MeanTraceSnip_GroomEvents', 'var')
      
    Export_ROI_MeanTraceSnip_GroomEvents=ROI_MeanTraceSnip_GroomEvents(:,1:5);

    tempcell=num2cell(cell2mat(ROI_MeanTraceSnip_GroomEvents(:,6)));

    Export_ROI_MeanTraceSnip_GroomEvents=[Export_ROI_MeanTraceSnip_GroomEvents,tempcell];

    Export_ROI_MeanTraceSnip_GroomEvents=cell2table(Export_ROI_MeanTraceSnip_GroomEvents);

end

if exist('ROI_MeanTraceSnip_PupilEvents', 'var')
      
    Export_ROI_MeanTraceSnip_PupilEvents=ROI_MeanTraceSnip_PupilEvents(:,1:5);

    tempcell=num2cell(cell2mat(ROI_MeanTraceSnip_PupilEvents(:,6)));

    Export_ROI_MeanTraceSnip_PupilEvents=[Export_ROI_MeanTraceSnip_PupilEvents,tempcell];

    Export_ROI_MeanTraceSnip_PupilEvents=cell2table(Export_ROI_MeanTraceSnip_PupilEvents);

end

if exist('ROI_MeanTraceSnip_PuffEvents', 'var')
     
    Export_ROI_MeanTraceSnip_PuffEvents=ROI_MeanTraceSnip_PuffEvents(:,1:5);

    tempcell=num2cell(cell2mat(ROI_MeanTraceSnip_PuffEvents(:,6)));

    Export_ROI_MeanTraceSnip_PuffEvents=[Export_ROI_MeanTraceSnip_PuffEvents,tempcell];

    Export_ROI_MeanTraceSnip_PuffEvents=cell2table(Export_ROI_MeanTraceSnip_PuffEvents);

end

if exist('ROI_MeanTraceSnip_WhiskEvents', 'var')
     
    Export_ROI_MeanTraceSnip_WhiskEvents=ROI_MeanTraceSnip_WhiskEvents(:,1:5);

    tempcell=num2cell(cell2mat(ROI_MeanTraceSnip_WhiskEvents(:,6)));

    Export_ROI_MeanTraceSnip_WhiskEvents=[Export_ROI_MeanTraceSnip_WhiskEvents,tempcell];

    Export_ROI_MeanTraceSnip_WhiskEvents=cell2table(Export_ROI_MeanTraceSnip_WhiskEvents);

end


clear tempcell


%% Prepare the alligned sink traces for export

ExportreadySinksalligned=cell(2,size(Filters_ROIsandEvents,2));
for i=1:size(Filters_ROIsandEvents,2) 

    ExportreadySinksalligned{1,i}=[Filters_ROIsandEvents{1,i},Filters_ROIsandEvents{2,i},Filters_ROIsandEvents{3,i}];

end

for i=1:size(Filters_ROIsandEvents,2) %go thought the different groups

    %getting only the cells that belong to recordings in that experimental
    %group
    tempcell=Sinks_Traces(Filters_ROIsandEvents{4,i},7);

    tempMetadata=Sinks_Traces(Filters_ROIsandEvents{4,i},1:5);

    wrapper = @(x) size(x, 2); 
    maxcols = max(cellfun(wrapper,tempcell));  %get the max number of columns of the widest array (longer duration sink across all recordings for that group)   
    temparray=nan(sum(cellfun(@height, tempcell)),maxcols);
    tempMetadata_pool=cell(sum(cellfun(@height, tempcell)),5);

    counter=1;
    for j=1:size(tempcell,1) %go through the recordings
      
        temparray(counter:counter+size(tempcell{j},1)-1,1:size(tempcell{j},2))=tempcell{j}; %combine the traces

        tempMetadata_pool(counter:counter+size(tempcell{j},1)-1,1)={tempMetadata{j,1}};
        tempMetadata_pool(counter:counter+size(tempcell{j},1)-1,2)={tempMetadata{j,2}};
        tempMetadata_pool(counter:counter+size(tempcell{j},1)-1,3)={tempMetadata{j,3}};
        tempMetadata_pool(counter:counter+size(tempcell{j},1)-1,4)={tempMetadata{j,4}};
        tempMetadata_pool(counter:counter+size(tempcell{j},1)-1,5)={Filters_ROIsandEvents{3,i}};

        counter=counter+size(tempcell{j},1);
    end


    ExportreadySinksalligned{2,i}=cell2table([tempMetadata_pool, num2cell(temparray)]);
 
end

clear i j tempcell wrapper maxcols temparray tempMetadata_pool tempMetadata counter
%%
end
Table_OxygenSinks_OutCombo=[Table_OxygenSinks_OutCombo,AdditionalOxySinkMetrics];
Table_OxygenSurges_OutCombo=[Table_OxygenSurges_OutCombo,AdditionalOxySurgeMetrics];

%% (7) Exporting results

disp('Creating folder for sorted data and statistical analysis matrices');
Statsoutputfolder=['Stats_Output_',datestr(now,30)]; %adds an identifier in the form of yyyymmddTHHMMSS
Figuresoutputfolder=['Figures_Output_',datestr(now,30)];
mkdir(Statsoutputfolder)
mkdir(Figuresoutputfolder)
%saving almost everything in a mat file
if strcmp(iOS, 'BLI')
save([Masterfolder,'\',Statsoutputfolder,'\','DataOutput.mat'],...
    'Table_OxygenSurges_OutCombo','Table_OxygenSinks_OutCombo','FiltersOxySinksMetrics','FiltersOxySurgesMetrics','Filters_ROIsandEvents','ExportTraces', 'Behaviour_data_combo','Behaviouraldatalogical','ROIs_Traces','NumOngoingOxysinks',...
    'TotalSinkArea_Norm','NumOngoingOxysurges','TotalSurgeArea','ExportTraceCorrs','TraceCorrs','SinksRaster','SurgesRaster');
else
save([Masterfolder,'\',Statsoutputfolder,'\','DataOutput.mat'],...
    'Table_OxygenSurges_OutCombo','Table_OxygenSinks_OutCombo','FiltersOxySinksMetrics','FiltersOxySurgesMetrics','Filters_ROIsandEvents','ExportTraces','NumOngoingOxysinks',...
    'TotalSinkArea_Norm','NumOngoingOxysurges','TotalSurgeArea','SinksRaster','SurgesRaster');
end
%saving oxygen surges table for LME fitting
save([Masterfolder,'\',Statsoutputfolder,'\','SurgeTable4LME.mat'],'Table_OxygenSurges_OutCombo');
%saving oxygen sinks table for LME fitting
save([Masterfolder,'\',Statsoutputfolder,'\','SinkTable4LME.mat'],'Table_OxygenSinks_OutCombo');


cd(Statsoutputfolder)
Output_xlsx=sprintf('FilteredData.xlsx');
writetable(ExportReady{1,1},Output_xlsx,'Sheet','OxySinkArea_um','WriteVariableNames',0);
writetable(ExportReady{2,1},Output_xlsx,'Sheet','OxySinkFilledArea_um','WriteVariableNames',0);
writetable(ExportReady{3,1},Output_xlsx,'Sheet','OxySinkDiameter_um','WriteVariableNames',0);
writetable(ExportReady{4,1},Output_xlsx,'Sheet','OxySinkPerimeter_um','WriteVariableNames',0);
writetable(ExportReady{5,1},Output_xlsx,'Sheet','OxySinkArea_Norm','WriteVariableNames',0);
writetable(ExportReady{6,1},Output_xlsx,'Sheet','OxySinkArea_Filled_Norm','WriteVariableNames',0);
writetable(ExportReady{7,1},Output_xlsx,'Sheet','NumOxySinkEvents_Norm','WriteVariableNames',0);
writetable(ExportReady{8,1},Output_xlsx,'Sheet','MeanOxySinkEvent_Duration','WriteVariableNames',0);
writetable(ExportReady{9,1},Output_xlsx,'Sheet','MeanOxySinkEvent_NormAmp','WriteVariableNames',0);
writetable(ExportReady{10,1},Output_xlsx,'Sheet','MeanOxySinkEvent_SizeMod','WriteVariableNames',0);

writetable(ExportReady{1,2},Output_xlsx,'Sheet','OxySurgeArea_um','WriteVariableNames',0);
writetable(ExportReady{2,2},Output_xlsx,'Sheet','OxySurgeFilledArea_um','WriteVariableNames',0);
writetable(ExportReady{3,2},Output_xlsx,'Sheet','OxySurgeDiameter_um','WriteVariableNames',0);
writetable(ExportReady{4,2},Output_xlsx,'Sheet','OxySurgePerimeter_um','WriteVariableNames',0);
writetable(ExportReady{5,2},Output_xlsx,'Sheet','OxySurgeArea_Norm','WriteVariableNames',0);
writetable(ExportReady{6,2},Output_xlsx,'Sheet','OxySurgeArea_Filled_Norm','WriteVariableNames',0);
writetable(ExportReady{7,2},Output_xlsx,'Sheet','NumOxySurgeEvents_Norm','WriteVariableNames',0);
writetable(ExportReady{8,2},Output_xlsx,'Sheet','MeanOxySurgeEvent_Duration','WriteVariableNames',0);
writetable(ExportReady{9,2},Output_xlsx,'Sheet','MeanOxySurgeEvent_NormAmp','WriteVariableNames',0);
writetable(ExportReady{10,2},Output_xlsx,'Sheet','MeanOxySurgeEvent_SizeMod','WriteVariableNames',0);
%%
if strcmp(iOS, 'BLI')
%exporting alligned sinks
for i=1:size(ExportreadySinksalligned,2)

    findW=strfind(ExportreadySinksalligned{1,i},'W');
    writetable(ExportreadySinksalligned{2,i},Output_xlsx,'Sheet',...
        ['SinksLevel',ExportreadySinksalligned{1,i}(1:7),ExportreadySinksalligned{1,i}(findW(end):strfind(ExportreadySinksalligned{1,i},'_')-1),'stim'],'WriteVariableNames',0);

end

clear i findW
end

counter=1;
for i=4:size(ExportTraces,1) %for every trace type

    for k=2:size(ExportTraces,2) %for every experimental group
    
        %export the table on a new sheet.
        %the sheet is the same for all k experimental groups of the same
        %trace type
        %using xlscol I am placing the traces one next to the other
        if ~isempty(ExportTraces{i,k})
        writetable(ExportTraces{i,k},Output_xlsx,'Sheet',ExportTraces{i,1},'Range', [xlscol(counter),'1'],'WriteVariableNames',0);

        counter=counter+ size(ExportTraces{i,k},2);
        end
    end

    counter=1; %reset the counter before going to the next trace type
end

if strcmp(iOS, 'BLI')
counter=1;
for i=4:size(ExportTraceCorrs,1) %for every trace type

    for k=2:size(ExportTraceCorrs,2) %for every experimental group
    
        %export the table on a new sheet.
        %the sheet is the same for all k experimental groups of the same
        %trace type
        %using xlscol I am placing the traces one next to the other
        writetable(ExportTraceCorrs{i,k},Output_xlsx,'Sheet',ExportTraceCorrs{i,1},'Range', [xlscol(counter),'1'],'WriteVariableNames',0);

        counter=counter+ size(ExportTraceCorrs{i,k},2);
    end

    counter=1; %reset the counter before going to the next trace type
end

clear i k
end

%this is likely to cause problems because not all recordings that have
%pupil data have Lpaw and Rpaw tracking data!
%I need to separate
if strcmp(iOS, 'BLI')
counter=1;
for i=4:size(ExportLPawBinnedTraces,1) %for every trace type

    for k=2:size(ExportLPawBinnedTraces,2) %for every experimental group
    
        %export the tables on  new sheets.
        %the sheet is the same for all k experimental groups of the same
        %trace type
        %using xlscol I am placing the traces one next to the other
        
        writetable(ExportLPawBinnedTraces{i,k},Output_xlsx,'Sheet',[ExportLPawBinnedTraces{i,1},'LPawbin'],'Range', [xlscol(counter),'1']);
        writetable(ExportRPawBinnedTraces{i,k},Output_xlsx,'Sheet',[ExportLPawBinnedTraces{i,1},'RPawbin'],'Range', [xlscol(counter),'1']);
        writetable(ExportPupilBinnedTraces{i,k},Output_xlsx,'Sheet',[ExportPupilBinnedTraces{i,1},'Pupilbin'],'Range', [xlscol(counter),'1']);

        counter=counter+ size(ExportPupilBinnedTraces{i,k},2);
    end

    counter=1; %reset the counter before going to the next trace type
end

clear i k
end
if strcmp(iOS, 'BLI')
counter=1;
for i=4:size(Pooled_Traces,1) %for every trace type

    for k=2:size(Pooled_Traces,2) %for every experimental group

        if ~isempty(Pooled_Traces{i,k})

            tempcell=cell(3,width(Pooled_Traces{i,k})); 
            tempcell(1,:)={Pooled_Traces{1,k}};
            tempcell(2,:)={Pooled_Traces{2,k}};
            tempcell(3,:)={Pooled_Traces{3,k}};

            writetable(cell2table(tempcell),Output_xlsx,'Sheet',[Pooled_Traces{i,1},'Pooled_Traces'],'Range', [xlscol(counter),'1']);
            
            writetable(Pooled_Traces{i,k},Output_xlsx,'Sheet',[Pooled_Traces{i,1},'Pooled_Traces'],'Range', [xlscol(counter),'5']);


            counter=counter+ width(Pooled_Traces{i,k});
        end

    end

    counter=1; %reset the counter before going to the next trace type
end

clear i k tempcell
end

%exporting the Event-based analysis
if exist('Export_ROI_MeanTraceSnip_ManEvents','var')

    writetable(Export_ROI_MeanTraceSnip_ManEvents,Output_xlsx,'Sheet','ROI_MeanTraceSnip_ManEvents');

end

if exist('Export_ROI_MeanTraceSnip_LPawEvents','var')

    writetable(Export_ROI_MeanTraceSnip_LPawEvents,Output_xlsx,'Sheet','ROI_MeanTraceSnip_LPawEvents');

end

if exist('Export_ROI_MeanTraceSnip_RPawEvents','var')

    writetable(Export_ROI_MeanTraceSnip_RPawEvents,Output_xlsx,'Sheet','ROI_MeanTraceSnip_RPawEvents');

end

if exist('Export_ROI_MeanTraceSnip_GroomEvents','var')

    writetable(Export_ROI_MeanTraceSnip_GroomEvents,Output_xlsx,'Sheet','ROI_MeanTraceSnip_GroomEvents');

end

if exist('Export_ROI_MeanTraceSnip_PupilEvents','var')

    writetable(Export_ROI_MeanTraceSnip_PupilEvents,Output_xlsx,'Sheet','ROI_MeanTraceSnip_PupilEvents');

end

if exist('Export_ROI_MeanTraceSnip_PuffEvents','var')

    writetable(Export_ROI_MeanTraceSnip_PuffEvents,Output_xlsx,'Sheet','ROI_MeanTraceSnip_PuffEvents');

end

if exist('Export_ROI_MeanTraceSnip_WhiskEvents','var')

    writetable(Export_ROI_MeanTraceSnip_WhiskEvents,Output_xlsx,'Sheet','ROI_MeanTraceSnip_WhiskEvents');

end





cd(Masterfolder)

%% Plotting traces

% Variables I will be using in the trace plotting
% Behaviour_data_combo;
% Reminder
% Column 1: LpawMovmm;
% Column 2: RpawMovmm;
% Column 3: NoseMovmm;
% Column 4: Groomlogical;
% Column 5: LpawMovZ;
% Column 6: RpawMovZ;
% Column 7: PupilDiammm;
% Column 8: PupilDiamZ;
% Column 9: Puffs;

% Behaviouraldatalogical;
% Reminder
% Column 1: logical vector for left paw movement events (calculated here)
% Column 2: logical vector for right paw movement events (calculated here)
% Column 3: logical vector for grooming events (copied)
% Column 4: logical vector for pupil dialation events (calculated here)
% Column 5: logical vector for puffs (copied)

% Oxygen data
% ExportTraces;
% Row 4:  ROIs mean signal 
% Row 5:  ROIs cov coef
% Row 6:  ROIs entropy
% Row 7:  NumOngoingOxysinks;
% Row 8:  TotalSinkArea;
% Row 9:  NumOngoingOxysurges;
% Row 10: TotalSurgeArea;
% Row 11: ROIs diff mean signal 
% Row 12: ROIs diff cov coef
% Row 13: ROIs diff entropy

%% PLot the traces with the puffs as semitransparent ribons (awake or aneasthetised)
%the plotting needs to be on the average of each group
%this makes sense in the puffs sessions only
%instead of looping through the recordings I need to loop through the filtered group data

%ask if there are any sessions with puffs
if any(~cellfun(@isempty,Behaviouraldatalogical(:,5)))

for i=2:size(ExportTraces,2) %loop through the groups
    
    %
    
    % the puff protocol is the same for all recordings so I can get the puff logical vector from the fist recording that has it        
    pfflogical=Behaviouraldatalogical{find(~cellfun(@isempty,Behaviouraldatalogical(:,5)),1,'first'),5};

    %check if there are puffs for that group based on the tag on row 3 "With stimulation" 
    %also check that the data in the table are not NaNs (I just check the
    %first column(mouse)
    if strcmp(ExportTraces{3,i},'With stimulation') && ~isnan(ExportTraces{4,i}.Var1{5})

        for data_trace=4:size(ExportTraces,1) %go through the different traces
 
            % Convert the 5:end rows of the table to matrix and calculate
            % the animal mean
            tracemean=mean(cell2mat(ExportTraces{data_trace,i}{5:end,:}),2);
            % plot the animalmeanExportTraces with animalstemExportTraces as shadow 
            Figure1=figure('visible','off');
            options.handle=Figure1;
            options.color_area=hex2rgb('#BFBFBF'); %grey
            options.color_line=hex2rgb('#000000'); %black
            %I take the hexcodes from here https://davidmathlogic.com/colorblind/
            options.alpha=0.7; %transparency of the shaded area
            options.line_width=1;
            options.error='std';
            plot_areaerrorbar(cell2mat(ExportTraces{data_trace,i}{5:end,:})',options);
            % use the ExportTraces(1:3,i) and ExportTraces(data_trace,1) info to create a title for the graph
            title([ExportTraces{1,i},ExportTraces{2,i},ExportTraces{3,i},' ',ExportTraces{data_trace,1}]);
            xlabel('Time(sec)') 
            ylabel('Z')
            % adapt the y axis limits to extend 10% of the max value at either end
            ylim([min(tracemean)-max(tracemean)*0.1 max(tracemean)*1.1])
 
            hold on
            %label all the individual puffs
            puff_num=bwlabel(pfflogical);

            %go through each puff and add it to the plot using the patch function
            for k=1:max(puff_num)
                %find the start of the puff and convert to index that makes
                %sense for the sampling rate of imaging data/traces
                boxX(1:2)=(find(puff_num==k,1,'first')/PuffsSFs)*SampleFs{find(~cellfun(@isempty,Behaviouraldatalogical(:,5)),1,'first')}; 
                %find the end of the puff and do the same conversion
                boxX(3:4)=(find(puff_num==k,1,'last')/PuffsSFs)*SampleFs{find(~cellfun(@isempty,Behaviouraldatalogical(:,5)),1,'first')};
                %this was to define the start and the end of each stimulation period   
                %(top left, bottom left, top right, bottom right)
                %for the boxY I will get the max and min of the y axis    
                boxY=[min(tracemean)-max(tracemean)*0.1 max(tracemean)*1.1 max(tracemean)*1.1 min(tracemean)-max(tracemean)*0.1];

                %plottinga each stimulation as a rectangular with transparency
                patch(boxX,boxY,hex2rgb('#FFC20A'),'EdgeColor','none','FaceAlpha',0.3)

            end

            hold off

            exportgraphics(gcf,[Masterfolder,'\',Figuresoutputfolder,'\',ExportTraces{1,i},...
                ExportTraces{2,i},ExportTraces{3,i},' ',ExportTraces{data_trace,1},'.pdf'],'ContentType','vector')

            close
        end
    end
end

clear i data_trace pfflogical k tracemean  options Figure1

end




%DEFINITELLY add

%get for each puff a window around and plot
% window -30s 30s puffs +30s after
% use -30s average as baseline to normlise the rest
% timepoint zero is start of the puffing 
% y-shift for the signal to be 0 at timepoint 0

%
% I also need to decide what pair of metrics makes sense to plot together
% Maybe all the ones that I have calculated the cor coef for?

%example double trace plotting that may have different sampling frequency
%create the linear spaced x axes for the behaviour and oxygen data 
% x1=linspace(0,600,length(behtest));
% x2=linspace(0,600,length(ROIs_Traces{4,7}));
% x3=linspace(0,600,length(ROIs_Traces{4,7}));
% 
% fig=figure;
% left_color = hex2rgb('#1A85FF'); % Y labels are blue.
% right_color = hex2rgb('#D41159'); % X labels are red.
% set(fig,'defaultAxesColorOrder',[left_color; right_color]);
% yyaxis left
% if min(behtest)<0 
%    
%     ylim([min(behtest)*1.1 max(behtest)*1.1])
% else 
%     ylim([min(behtest)-max(behtest)*0.1 max(behtest)*1.1])
% 
% end
% plot(x1,(behtest)','Color','#1A85FF');
% 
% yyaxis right
% if min(ROIs_Traces{4,7})<0 
%    
%     ylim([min(ROIs_Traces{4,7})*1.1 max(ROIs_Traces{4,7})*1.1])
% else 
%     ylim([min(ROIs_Traces{4,7})-max(ROIs_Traces{4,7})*0.1 max(ROIs_Traces{4,7})*1.1])
% 
% end
% plot(x2,(ROIs_Traces{4,7})','Color','#D41159');

%% (8) Linear Mix effects modeling analysis 

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Function to find peaks in noisy data
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

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

%%

function latestfile = getlatestfile(directory,type, filetype, partname)
%% This function returns the latest file or folder from the directory passsed as input argument

% Inputs:
% DIRECTORY is the *FULL* path of a folder (required!!)
% TYPE is a string that is either 'file' or 'folder' (required if additional arguments are used!!) 
% FILETYPE is a string that indicates the file extention e.g. '*.csv'. Of course this is not needed when the type input is 'folder'
% PARTNAME is a string containing part of the name of folder or file we are after.

% Examples:
% This will return the most recent file of any type in the specified directory
% latestfile = getlatestfile(directory)
% This will return the most recent folder in the specified directory
% latestfile = getlatestfile(directory,'folder')
% This will return the most recent .mat file in the specified directory
% latestfile = getlatestfile(directory,'file','*.mat')
% This will return the most recent .mat file that contains 'Combo' in its name in the specified directory
% latestfile = getlatestfile(directory,'file','*.mat', '*Combo')
% This will return the most recent file of any type that contains 'Combo' in its name in the specified directory
% latestfile = getlatestfile(directory,'file',[], '*Combo') 
% This will return the most recent folder that contains 'Combo' in its name in the specified directory
% latestfile = getlatestfile(directory,'folder',[], '*Combo') 

% Output:
% The function returns the full path of the most recent folder or file of the type indicated

% Written by Antonis Asiminas, PhD
% Copenhagen, Denmark, August 2021,
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
if nargin < 2
    % default mode is to look for the later file 
    disp('The most recent file of any type will be returned')
    type ='file';
end    

if nargin < 3
    
    filetype=[]; 
    
end

if nargin < 4
    
    partname=[];
    
end


switch type              
    case 'file'
        if ~isempty(filetype)                  
            if ~isempty(partname)
                %Go into the folder of interest
                cd(directory)
                %Get the files of the filetype of interest that contain the specified string in their name
                dirc = dir([partname,filetype]);
                %I contains the index to the biggest number which is the latest file 
                [~,I] = max([dirc(:).datenum]); 
            else
                %Go into the folder of interest
                cd(directory)
                %Get the files of the filetype of interest
                dirc = dir(filetype);
                %I contains the index to the biggest number which is the latest file 
                [~,I] = max([dirc(:).datenum]); 
            end
    
        else
            if ~isempty(partname)
                %Go into the folder of interest
                cd(directory)
                %Get the directory contents that contain the specified string in their name
                dirc = dir([partname,'*']);
                %Get only the files            
                dirc = dirc(find(~cellfun(@isdir,{dirc(:).name})));
                %I contains the index to the biggest number which is the latest file
                [~,I] = max([dirc(:).datenum]);  
            else
                %Go into the folder of interest            
                cd(directory)            
                %Get the directory contents            
                dirc = dir;              
                %Get only the files            
                dirc = dirc(find(~cellfun(@isdir,{dirc(:).name})));            
                %I contains the index to the biggest number which is the latest file
                [~,I] = max([dirc(:).datenum]);            
            end
        end
             
    case 'folder'        
        if ~isempty(partname)
            %Go into the folder of interest       
            cd(directory)
            %Get the directory contents that contain the specified string in their name    
            dirc = dir([partname,'*']);
            %Get only the folders            
            dirc = dirc(find(cellfun(@isdir,{dirc(:).name})));  
            %I contains the index to the biggest number which is the latest file  
            [~,I] = max([dirc(:).datenum]);         
        else
            %Go into the folder of interest       
            cd(directory)        
            %Get the directory contents          
            dirc = dir;         
            %Get only the folders            
            dirc = dirc(find(cellfun(@isdir,{dirc(:).name})));            
            % Remove the first two entries in dirc are '.' and '..' which are classified as folders         
            dirc = dirc(3:end);         
            %I contains the index to the biggest number which is the latest file  
            [~,I] = max([dirc(:).datenum]); 
        end     
end
    

if ~isempty(I)
    latestfile = [directory,'\',dirc(I).name];
end

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

%%
function Entropy = info_entropy (input, b)
% In information theory, entropy is a measure of the uncertainty associated 
% with a random variable. In this context, the term usually refers to the 
% Shannon entropy, which quantifies the expected value of the information 
% contained in a message, usually in units such as bits. In this context, a 
% 'message' means a specific realization of the random variable.
% Shannon denoted the entropy H of a discrete random variable X with possible 
% values {x1, ..., xn} as,
%
%     H(X) = E(I(X)).
% 
% E is the expected value,
% I is the information content of X
% 
% I(X) is itself a random variable. If p denotes the probability mass function 
% of X then the entropy can explicitly be written as,
%             n                n                         n
%     H(X) = SUM p(xi)I(xi) = SUM p(xi)logb(1/p(xi)) = -SUM p(xi)logb(p(xi))
%            i=1              i=1                       i=1
%
% Example usage
% -------------
% in = [.25 .25 .25 .25];
% b = 'bit';
% Entropy = info_entropy (in, b)
Entropy = 1;
if b == 'bit'
    info_sz = size (input);
    n = info_sz (1, 2);
    info_sum = sum (input);
    tol = 5 * eps(info_sum); %this is the allowed tolerance based on floating-point relative accuracy (eps)
    if single(info_sum) > 1 + tol
        error ('Sum of probabilities cannot be greater than 1');
    end
    Entropy = 0;
    for i=1:n
        if isequal(input(1,i),0)
            tmp = 0;
        else
            tmp = (-input(1,i)*log2(input(1,i)));
        end
        Entropy = Entropy + tmp;
    end
        
elseif b == 'nat'
    info_sz = size (input);
    n = info_sz (1, 2);
    info_sum = sum (input);
    tol = 5 * eps(info_sum); %this is the allowed tolerance based on floating-point relative accuracy (eps)
    if single(info_sum) > 1 + tol 
        error ('Sum of probabilities cannot be greater than 1');
    end
    Entropy = 0;
    for i=1:n
        if isequal(input(1,i),0)
            tmp = 0;
        else
            tmp = (-input(1,i)*log(input(1,i)));
        end
        Entropy = Entropy + tmp;
    end
        
elseif b == 'dit'
    info_sz = size (input);
    n = info_sz (1, 2);
    info_sum = sum (input);
    if single(info_sum) > 1
        error ('Sum of probabilities cannot be greater than 1');
    end
    Entropy = 0;
    for i=1:n
        if isequal(input(1,i),0)
            tmp = 0;
        else
            tmp = (-input(1,i)*log10(input(1,i)));
        end
        Entropy = Entropy + tmp;
    end
        
end
end

%% This function takes an array and returns the Shannon Entropy of an array.

function [H] = Entropy_Array(X, varargin)

% A probability density function is derived from the array.
%
% By default the program selects a number of bins proportional to the square
% root of the number of points in the array, however the Freedman?Diaconis rule ('fd') 
% and Sturge's rule ('sturges') can be selected.
%
%
% Ex: X = randn(1000,1);
%
% H_default = Entropy_Array(X);
% H_fd = Entropy_Array(X,'fd');
%
%
% 3/30/16 - Ben Lucas
% Create an empirical probability distribution p for the values of X:
% Select number of bins - This uses a square root method by default,
if isempty(varargin)
    num_bins = ceil(sqrt(length(X)));
else 
    method = varargin{1};
    if strcmpi(method,'fd')
        num_bins = ceil(range(X)./(2*iqr(X)*length(X)^(-1/3)));
    elseif strcmpi(method, 'sturges')
        num_bins = ceil(1 + log2(length(X)));
    elseif strcmpi(method, 'sqrt')
        num_bins = ceil(sqrt(length(X)));
    end
end
bins = linspace(min(X), max(X), num_bins);
bins(end) = bins(end) + 1; % eliminates the edge case by setting max(X) < bins(end)
X = sort(X);
p = zeros(1,length(bins)-1);
for i = 1:(length(bins)-1)
    p(i) = sum((X>=bins(i)) & (X<bins(i+1)));
end
p = p/sum(p);   % normalize to PDF
p(p==0) = [];   % eliminate empty bins
H = sum(-p.*log2(p));



end

%% XLSCOL Convert Excel column letters to numbers or vice versa.
%   B = XLSCOL(A) takes input A, and converts to corresponding output B.
%   The input may be a number, a string, an array or matrix, an Excel
%   range, a cell, or a combination of each within a cell, including nested
%   cells and arrays. The output maintains the shape of the input and
%   attempts to "flatten" the cell to remove nesting.  Numbers and symbols
%   within strings or Excel ranges are ignored.
%
%   Examples
%   --------
%       xlscol(256)   % returns 'IV'
%
%       xlscol('IV')  % returns 256
%
%       xlscol([405 892])  % returns {'OO' 'AHH'}
%
%       xlscol('A1:IV65536')  % returns [1 256]
%
%       xlscol({8838 2430; 253 'XFD'}) % returns {'MAX' 'COL'; 'IS' 16384}
%
%       xlscol(xlscol({8838 2430; 253 'XFD'})) % returns same as input
%
%       b = xlscol({'A10' {'IV' 'ALL34:XFC66'} {'!@#$%^&*()'} '@#$' ...
%         {[2 3]} [5 7] 11})
%       % returns {1 [1x3 double] 'B' 'C' 'E' 'G' 'K'}
%       %   with b{2} = [256 1000 16383]
%
%   Notes
%   -----
%       CELLFUN and ARRAYFUN allow the program to recursively handle
%       multiple inputs.  An interesting side effect is that mixed input,
%       nested cells, and matrix shapes can be processed.
%
%   See also XLSREAD, XLSWRITE.
%
%   Version 1.1 - Kevin Crosby
% DATE      VER  NAME          DESCRIPTION
% 07-30-10  1.0  K. Crosby     First Release
% 08-02-10  1.1  K. Crosby     Vectorized loop for numerics.
% Contact: Kevin.L.Crosby@gmail.com

function [b] = xlscol(a)

base = 26;
if iscell(a)
  b = cellfun(@xlscol, a, 'UniformOutput', false); % handles mixed case too
elseif ischar(a)
  if ~isempty(strfind(a, ':')) % i.e. if is a range
    b = cellfun(@xlscol, regexp(a, ':', 'split'));
  else % if isempty(strfind(a, ':')) % i.e. if not a range
    b = a(isletter(a));        % get rid of numbers and symbols
    if isempty(b)
      b = {[]};
    else % if ~isempty(a);
      b = double(upper(b)) - 64; % convert ASCII to number from 1 to 26
      n = length(b);             % number of characters
      b = b * base.^((n-1):-1:0)';
    end % if isempty(a)
  end % if ~isempty(strfind(a, ':')) % i.e. if is a range
elseif isnumeric(a) && numel(a) ~= 1
  b = arrayfun(@xlscol, a, 'UniformOutput', false);
else % if isnumeric(a) && numel(a) == 1
  n = ceil(log(a)/log(base));  % estimate number of digits
  d = cumsum(base.^(0:n+1));   % offset
  n = find(a >= d, 1, 'last'); % actual number of digits
  d = d(n:-1:1);               % reverse and shorten
  r = mod(floor((a-d)./base.^(n-1:-1:0)), base) + 1;  % modulus
  b = char(r+64);  % convert number to ASCII
end % if iscell(a)
% attempt to "flatten" cell, by removing nesting
if iscell(b) && (iscell([b{:}]) || isnumeric([b{:}]))
  b = [b{:}];
end % if iscell(b) && (iscell([b{:}]) || isnumeric([ba{:}]))

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

%% Function to convert hex colour code to rgb vector
function [ rgb ] = hex2rgb(hex,range)
% hex2rgb converts hex color values to rgb arrays on the range 0 to 1. 
% 
% 
% * * * * * * * * * * * * * * * * * * * * 
% SYNTAX:
% rgb = hex2rgb(hex) returns rgb color values in an n x 3 array. Values are
%                    scaled from 0 to 1 by default. 
%                    
% rgb = hex2rgb(hex,256) returns RGB values scaled from 0 to 255. 
% 
% 
% * * * * * * * * * * * * * * * * * * * * 
% EXAMPLES: 
% 
% myrgbvalue = hex2rgb('#334D66')
%    = 0.2000    0.3020    0.4000
% 
% 
% myrgbvalue = hex2rgb('334D66')  % <-the # sign is optional 
%    = 0.2000    0.3020    0.4000
% 
%
% myRGBvalue = hex2rgb('#334D66',256)
%    = 51    77   102
% 
% 
% myhexvalues = ['#334D66';'#8099B3';'#CC9933';'#3333E6'];
% myrgbvalues = hex2rgb(myhexvalues)
%    =   0.2000    0.3020    0.4000
%        0.5020    0.6000    0.7020
%        0.8000    0.6000    0.2000
%        0.2000    0.2000    0.9020
% 
% 
% myhexvalues = ['#334D66';'#8099B3';'#CC9933';'#3333E6'];
% myRGBvalues = hex2rgb(myhexvalues,256)
%    =   51    77   102
%       128   153   179
%       204   153    51
%        51    51   230
% 
% HexValsAsACharacterArray = {'#334D66';'#8099B3';'#CC9933';'#3333E6'}; 
% rgbvals = hex2rgb(HexValsAsACharacterArray)
% 
% * * * * * * * * * * * * * * * * * * * * 
% Chad A. Greene, April 2014
%
% Updated August 2014: Functionality remains exactly the same, but it's a
% little more efficient and more robust. Thanks to Stephen Cobeldick for
% the improvement tips. In this update, the documentation now shows that
% the range may be set to 256. This is more intuitive than the previous
% style, which scaled values from 0 to 255 with range set to 255.  Now you
% can enter 256 or 255 for the range, and the answer will be the same--rgb
% values scaled from 0 to 255. Function now also accepts character arrays
% as input. 
% 
% * * * * * * * * * * * * * * * * * * * * 
% See also rgb2hex, dec2hex, hex2num, and ColorSpec. 
% 
%% Input checks:
assert(nargin>0&nargin<3,'hex2rgb function must have one or two inputs.') 
if nargin==2
    assert(isscalar(range)==1,'Range must be a scalar, either "1" to scale from 0 to 1 or "256" to scale from 0 to 255.')
end
%% Tweak inputs if necessary: 
if iscell(hex)
    assert(isvector(hex)==1,'Unexpected dimensions of input hex values.')
    
    % In case cell array elements are separated by a comma instead of a
    % semicolon, reshape hex:
    if isrow(hex)
        hex = hex'; 
    end
    
    % If input is cell, convert to matrix: 
    hex = cell2mat(hex);
end
if strcmpi(hex(1,1),'#')
    hex(:,1) = [];
end
if nargin == 1
    range = 1; 
end
%% Convert from hex to rgb: 
switch range
    case 1
        rgb = reshape(sscanf(hex.','%2x'),3,[]).'/255;
    case {255,256}
        rgb = reshape(sscanf(hex.','%2x'),3,[]).';
    
    otherwise
        error('Range must be either "1" to scale from 0 to 1 or "256" to scale from 0 to 255.')
end
end


%%
% ----------------------------------------------------------------------- %
% Function plot_areaerrorbar plots the mean and standard deviation of a   %
% set of data filling the space between the positive and negative mean    %
% error using a semi-transparent background, completely customizable.     %
%                                                                         %
%   Input parameters:                                                     %
%       - data:     Data matrix, with rows corresponding to observations  %
%                   and columns to samples.                               %
%       - options:  (Optional) Struct that contains the customized params.%
%           * options.handle:       Figure handle to plot the result.     %
%           * options.color_area:   RGB color of the filled area.         %
%           * options.color_line:   RGB color of the mean line.           %
%           * options.alpha:        Alpha value for transparency.         %
%           * options.line_width:   Mean line width.                      %
%           * options.x_axis:       X time vector.                        %
%           * options.error:        Type of error to plot (+/-).          %
%                   if 'std',       one standard deviation;               %
%                   if 'sem',       standard error mean;                  %
%                   if 'var',       one variance;                         %
%                   if 'c95',       95% confidence interval.              %
% ----------------------------------------------------------------------- %
%   Example of use:                                                       %
%       data = repmat(sin(1:0.01:2*pi),100,1);                            %
%       data = data + randn(size(data));                                  %
%       plot_areaerrorbar(data);                                          %
% ----------------------------------------------------------------------- %
%   Author:  Victor Martinez-Cagigal                                      %
%   Date:    30/04/2018                                                   %
%   E-mail:  vicmarcag (at) gmail (dot) com                               %
% ----------------------------------------------------------------------- %
function plot_areaerrorbar(data, options)
    % Default options
    if(nargin<2)
        options.handle     = figure(1);
        options.color_area = [128 193 219]./255;    % Blue theme
        options.color_line = [ 52 148 186]./255;
        %options.color_area = [243 169 114]./255;    % Orange theme
        %options.color_line = [236 112  22]./255;
        options.alpha      = 0.5;
        options.line_width = 2;
        options.error      = 'std';
    end
    if(isfield(options,'x_axis')==0), options.x_axis = 1:size(data,2); end
    options.x_axis = options.x_axis(:);
    
    % Computing the mean and standard deviation of the data matrix
    data_mean = mean(data,1);
    data_std  = std(data,0,1);
    
    % Type of error plot
    switch(options.error)
        case 'std', error = data_std;
        case 'sem', error = (data_std./sqrt(size(data,1)));
        case 'var', error = (data_std.^2);
        case 'c95', error = (data_std./sqrt(size(data,1))).*1.96;
    end
    
    % Plotting the result
    figure(options.handle);
    x_vector = [options.x_axis', fliplr(options.x_axis')];
    patch = fill(x_vector, [data_mean+error,fliplr(data_mean-error)], options.color_area);
    set(patch, 'edgecolor', 'none');
    set(patch, 'FaceAlpha', options.alpha);
    hold on;
    plot(options.x_axis, data_mean, 'color', options.color_line, ...
        'LineWidth', options.line_width);
    hold off;
    
end

