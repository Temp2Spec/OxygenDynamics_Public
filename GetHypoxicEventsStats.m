%% Script that is a variation of the stats script that gets the information for individual hypoxic events

% Antonis Asiminas, PhD
% Center for Translational Neuromedicine
% University of Copenhagen, 28 Sep 2023

% Last updated 2 October 2023


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (1) Load the file with all data paths and metadata.
% The same file is used by the wrapper script and has provided all the
% metadata to the datatables produced by the master analysis script
clear;


Masterfolder=pwd;
InputD=readtable('ios_stack.csv','Delimiter',',');
if size(InputD,2)<2
    InputD=readtable('ios_stack.csv' );
end

%%

Paths=table2cell(InputD(:,{'Paths'}));
Mice=table2cell(InputD(:,{'Mouse'}));
Puffs=table2cell(InputD(:,{'PuffsFile'}));
Genotypes=table2cell(InputD(:,{'Genotype'}));
Conditions=table2cell(InputD(:,{'Condition'}));
DrugIDs=table2cell(InputD(:,{'DrugID'}));
Promoters=table2cell(InputD(:,{'Promoter'}));
SampleFs=table2cell(InputD(:,{'SampleF'}));
Pixelsizes=table2cell(InputD(:,{'Pixelsize'}));

%% Do you wnat the oldest or newest output from Master

SinkFold_OldORNew = questdlg('Do you want to use the oldest or the most recent data for oxygen sinks?', ...
	'Which folder to use?', ...
	'Oldest','Recent','Recent');

%% Ask if the recordings had epochs to be able to separate events based on the epoch

Distinct = questdlg('Do the recordings in this dataset contain distinct behavioural epochs?',...
    'Ys or No',...
    'Yes', 'No','No');

switch Distinct
    case 'Yes'

        prompt = 'How many epochs?';
        dlgTitle = 'Number of epochs';
        numLines = 1;
        defaultInput = {'1'};  % Default value in the input field

        userInput = inputdlg(prompt, dlgTitle, numLines, defaultInput);

        % Check if the user clicked Cancel or entered an empty value
        if isempty(userInput)
            disp('Operation canceled or no value entered.');
        else
            % Convert the user input to a numerical value
            Epochs = str2double(userInput{1});
    
            % Check if the input is a valid numerical value
            if ~isnan(Epochs)
                disp(['You entered: ', num2str(Epochs)]);
            else
                disp('Invalid input. Please enter a numerical value.');
            end
        end


        prompt = ['Enter the start and end of all ',num2str(Epochs),' epochs in seconds. This will applied to all recordings (e.g. 1,600,601,1200,1201,1800).'];
        dlgTitle = 'Epoch starts';
        numLines = 2;
        defaultInput = {'0,0,0,0'};  % Default values in the input fields

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

                Epochs_Starts=Epochs_Limits(1:2:end);
                Epochs_Ends=Epochs_Limits(2:2:end);

            else
               
                disp('Invalid input. Please enter numerical values.');
            end
  
        end

end

%% (2) Collecting data from all data folders and adding the metadata from the source scv 

Table_OxygenSinks_OutCombo=[];
Sinks_Area=cell(height(InputD),10);
%
for datai=1:length(Paths)

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
    PxlSz = Pixelsizes{datai};
    Puff = Puffs{datai};

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
    wrapper = @(x) contains(x,'Urefined');
    load(d(cellfun(wrapper, {d(:).name})).name,"Table_OxygenSinks_Out");
    DIR = fullfile(Tifffile.folder,'\', Tifffile.name);

    [IM_BW,~,~]=loadtiff(DIR);

    %converting to logical to make my life easier with size
    IM_BW=logical(IM_BW);

     %% 
    if ismember('Included_Pocs', who('-file', d(cellfun(wrapper, {d(:).name})).name)) %if I'm going to use currated data I need to get the logical vector to filter the oxygen sink data
                     
        load(d(cellfun(wrapper, {d(:).name})).name,"Included_Pocs")           
        Table_OxygenSinks_Out=Table_OxygenSinks_Out(Included_Pocs,:);
        
       
    else
            %a disclaimer if curated data were chosen to be analysed but the
            %filter from the curator is not present
          
            fprintf(['No Oxygen sink currated data were fround in data folder ', DatafileID ,'\n'])   
    end
   
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
        Pixelsize=cell(height(Table_OxygenSinks_Out),1);Pixelsize(:)={PxlSz};

        %adding the metadata in the table
        Table_OxygenSinks_Out = addvars(Table_OxygenSinks_Out,Experiment,Mouse,Condition,DrugID,PuffStim,Genotype,Promoter,Pixelsize,'Before','RecDuration');


        %add the current table to the combo table
        Table_OxygenSinks_OutCombo = [Table_OxygenSinks_OutCombo;Table_OxygenSinks_Out];

        wrapper = @(x) contains(x,'Urefined');
        %catching also the area of the pockets
        Sinks_Area{datai,1}=DatafileID;
        Sinks_Area{datai,2}=Mous;
        Sinks_Area{datai,3}=Cond;
        Sinks_Area{datai,4}=Drug;
        Sinks_Area{datai,5}=PuffStim(1);
        Sinks_Area{datai,6}=Promo;
        Sinks_Area{datai,7}=PxlSz;               
        load(d(cellfun(wrapper, {d(:).name})).name,'OxySink_Map');
        Sinks_Area{datai,8}=size(OxySink_Map);
        load(d(cellfun(wrapper, {d(:).name})).name,'OxySinkArea_all');
        if exist('Included_Pocs', 'var') %this means the currated data have been loaded. no need to ask for the response to the prompt again
            OxySinkArea_all=OxySinkArea_all(Included_Pocs,:);
            %else % filter them with the logical vector for potential noise if not
            %OxySinkArea_all=OxySinkArea_all(Trace_PotentialNoise,:);
        end
        Sinks_Area{datai,9}=OxySinkArea_all;
        
        Sinks_Area{datai,10}=IM_BW;

    end

    cd(Masterfolder)

end

    clear Included_Pocs d I datai Experiment Mouse Condition DrugID Genotype Promoter Pixelsize...
        DatafileID Mous Cond PxlSz Gen Promo Drug list wrapper SinkFolders...
        SinkFold_OldORNew OxySinkArea_all SinksDataFolder OxySink_Map Table_OxygenSinks_Out IM_BW
         
cd(Masterfolder)

%getting rid of any empty rows from directories that did not contain data
Sinks_Area=Sinks_Area(~cellfun(@isempty,Sinks_Area(:,1)),:);
%% (3) Get the hypoxic even specific metrics
OxygenSinksInfo_all=cell(size(Sinks_Area,1),1);
for datai=1:size(Sinks_Area,1)
    OxygenSinksInfo_rec=cell(size(Sinks_Area{datai,10},3),1);

    for i=1:size(Sinks_Area{datai,10},3)

        OxygenSinksInfo_rec{i}=regionprops(Sinks_Area{datai,10}(:,:,i),'Area','PixelIdxList');

    end

    OxygenSinksInfo_all{datai}=OxygenSinksInfo_rec;
end

clear OxygenSinksInfo_rec
%%
%this will contain putative unique pockets(rows) and the pixels belonging to events for every frame of the recording.
Overall_OxygenSinks_Pxllist_all=cell(1,length(OxygenSinksInfo_all));
sink_overlap_thres=0.6; %percentage of overlap that is needed between sinks in two different frames

for datai=1:length(OxygenSinksInfo_all)

    %creating a filter to get only the sinks for that recording
    wrapper = @(x) strcmp(x,  Sinks_Area{datai,1});
    TempIndex=cell2mat(cellfun(wrapper, Table_OxygenSinks_OutCombo.Experiment, 'UniformOutput', false)); 
    thisrecording=Table_OxygenSinks_OutCombo(TempIndex,9:end);

    Overall_OxygenSinks_Pxllist_rec=cell(size(Sinks_Area{datai,9}));
    
    for i=1:size(Sinks_Area{datai,9},1) %for every pocket

        for ij=1:size(Sinks_Area{datai,9},2) %%for every frame

            
           if ~isempty(Sinks_Area{datai,9}{i,ij}) %check if the pocket had an ongoing event in that frame


               % comparing with the size of the identified pockets in the
               % same frame.
               %if there is only pocket with the same size
               if sum([OxygenSinksInfo_all{datai}{ij}(:).Area]==Sinks_Area{datai,9}{i,ij})==1
               
               
                   %just capture the pixel indices for this frame
                   Overall_OxygenSinks_Pxllist_rec{i,ij}=OxygenSinksInfo_all{datai}{ij}([OxygenSinksInfo_all{datai}{ij}(:).Area]==Sinks_Area{datai,9}{i,ij}).PixelIdxList;


                   % replace the area with nan so it cannot be added to
                   % another pocket
                   OxygenSinksInfo_all{datai}{ij}([OxygenSinksInfo_all{datai}{ij}(:).Area]==Sinks_Area{datai,9}{i,ij}).Area=nan;

    
               %if there are more pixel lists than pockets in that frame
               elseif length(OxygenSinksInfo_all{datai}{ij})>sum(~cellfun(@isempty,Sinks_Area{datai,9}(:,ij)))
                   
                   %create all possible combos to add
                   allcombos=nchoosek(1:length(OxygenSinksInfo_all{datai}{ij}),2);

                   allcomboAreas=[OxygenSinksInfo_all{datai}{ij}(allcombos(:,1)).Area]+[OxygenSinksInfo_all{datai}{ij}(allcombos(:,2)).Area];

                   if sum(allcomboAreas==Sinks_Area{datai,9}{i,ij})==1 % check if there is a combo that gives the area of Sink I'm asfter


                       constituentpockets=allcombos(allcomboAreas==Sinks_Area{datai,9}{i,ij},:);

                       for t=constituentpockets
                       
                           %capture the pixel indices from all constituents of
                           %this frame into one
                           Overall_OxygenSinks_Pxllist_rec{i,ij}=[Overall_OxygenSinks_Pxllist_rec{i,ij};[OxygenSinksInfo_all{datai}{ij}(t).PixelIdxList]];

                       end

                       % replace the area of all constituent pockets with nans so it cannot be added to
                       % another pocket
                       for t=constituentpockets
                       
                           OxygenSinksInfo_all{datai}{ij}(t).Area=nan;

                       end


                   end
               
               
               else %if neither a single area or a combination matched the area from Sinks_Area, find the one with the most overlap


                   tempintersect=nan(length(OxygenSinksInfo_all{datai}{ij}),1);
                   
                   for k=1:length(OxygenSinksInfo_all{datai}{ij}) %go through the pixles lists of the pockets of this frame


                       %calculate the intersection between the pixels of
                       %the pockets for that frame and the unique pixels
                       %for that pocket across recording.
                       tempintersect(k)=length(intersect([OxygenSinksInfo_all{datai}{ij}(k).PixelIdxList],thisrecording.OxySink_Pxls_all{i}));


                   end



                   if sum(logical(tempintersect))==1


                        %and get the pixels form the pocket with the highest
                        %intersection             
                        Overall_OxygenSinks_Pxllist_rec{i,ij}=OxygenSinksInfo_all{datai}{ij}(tempintersect==max(tempintersect)).PixelIdxList;

                        % replace the area with nan so it cannot be added to
                        % another pocket
                        OxygenSinksInfo_all{datai}{ij}(tempintersect==max(tempintersect)).Area=nan;

                   end


               end
           end

        end

        
    end

    Overall_OxygenSinks_Pxllist_all{datai}=Overall_OxygenSinks_Pxllist_rec;

end

clear OxygenSinksInfo_all counter pixl_temp tempintersect Overall_OxygenSinks_Pxllist_rec
%%
Eventspecificmetrics=[];

for datai=1:size(Sinks_Area,1) %for each experiment

    %creating a filter to get only the sinks for that recording
    wrapper = @(x) strcmp(x,  Sinks_Area{datai,1});
    TempIndex=cell2mat(cellfun(wrapper, Table_OxygenSinks_OutCombo.Experiment, 'UniformOutput', false)); 
    thisrecording=Table_OxygenSinks_OutCombo(TempIndex,:);
  

    %getting the metadata from the input table       
    Experiment=cell(sum(thisrecording.NumOxySinkEvents(:)),1); Experiment(:)={thisrecording.Experiment{1}};       
    Mouse=cell(sum(thisrecording.NumOxySinkEvents(:)),1); Mouse(:)={thisrecording.Mouse{1}};           
    Condition=cell(sum(thisrecording.NumOxySinkEvents(:)),1);Condition(:)={thisrecording.Condition{1}};        
    DrugID=cell(sum(thisrecording.NumOxySinkEvents(:)),1);DrugID(:)={thisrecording.DrugID{1}};       
    Genotype=cell(sum(thisrecording.NumOxySinkEvents(:)),1); Genotype(:)={thisrecording.Genotype{1}};
    PuffStim=cell(sum(thisrecording.NumOxySinkEvents(:)),1); PuffStim(:)={thisrecording.PuffStim(1)};
    Promoter=cell(sum(thisrecording.NumOxySinkEvents(:)),1);Promoter(:)={thisrecording.Promoter{1}};     
    Pixelsize=cell(sum(thisrecording.NumOxySinkEvents(:)),1);Pixelsize(:)={thisrecording.Pixelsize{1}};

    PocketID=cell(sum(thisrecording.NumOxySinkEvents(:)),1);
    counter=1;
    for g=1:height(thisrecording)

        for gg=1:thisrecording.NumOxySinkEvents(g)

            PocketID(counter)={[thisrecording.Mouse{g},'_',num2str(g)]};

            counter=counter+1;
        end

    end
    EventID=cell(sum(thisrecording.NumOxySinkEvents(:)),1);
    counter=1;
    for g=1:height(thisrecording)

        for gg=1:thisrecording.NumOxySinkEvents(g)

            EventID(counter)={[PocketID{counter},'_',num2str(gg)]};

            counter=counter+1;
        end

    end


    %Add the information about amplitude and duration that is already in the Table_OxygenSinks_OutCombo
    
    NormOxySinkAmp=nan(sum(thisrecording.NumOxySinkEvents(:)),1); 
    counter=1;
    for g=1:height(thisrecording)

        for gg=1:thisrecording.NumOxySinkEvents(g)

             NormOxySinkAmp(counter)=thisrecording.NormOxySinkAmp{g}(gg);

            counter=counter+1;
        end

    end
    
    Start=nan(sum(thisrecording.NumOxySinkEvents(:)),1); 
    counter=1;
    for g=1:height(thisrecording)

        for gg=1:thisrecording.NumOxySinkEvents(g)

            Start(counter)=thisrecording.Start{g}(gg);

            counter=counter+1;
        end

    end

    Duration=nan(sum(thisrecording.NumOxySinkEvents(:)),1); 
    counter=1;
    for g=1:height(thisrecording)

        for gg=1:thisrecording.NumOxySinkEvents(g)

            Duration(counter)=thisrecording.Duration{g}(gg);

            counter=counter+1;
        end

    end
    
    Size_modulation=nan(sum(thisrecording.NumOxySinkEvents(:)),1); 
    counter=1;
    for g=1:height(thisrecording)

        for gg=1:thisrecording.NumOxySinkEvents(g)

            Size_modulation(counter)=thisrecording.Size_modulation{g}(gg);

            counter=counter+1;
        end

    end


    Area_um=nan(sum(thisrecording.NumOxySinkEvents(:)),1);
    Area_norm=nan(sum(thisrecording.NumOxySinkEvents(:)),1);
    FilledArea_um=nan(sum(thisrecording.NumOxySinkEvents(:)),1);
    FilledArea_norm=nan(sum(thisrecording.NumOxySinkEvents(:)),1);
    Centroid_x=nan(sum(thisrecording.NumOxySinkEvents(:)),1);
    Centroid_y=nan(sum(thisrecording.NumOxySinkEvents(:)),1);
    Circularity=nan(sum(thisrecording.NumOxySinkEvents(:)),1);
    Perimeter_um=nan(sum(thisrecording.NumOxySinkEvents(:)),1);
    Diameter_um=nan(sum(thisrecording.NumOxySinkEvents(:)),1);

    counter=1;
    for k=1:height(thisrecording) %for every pocket of this recording

        %create a cell vector with all the pixels of hypoxic events of this pocket
        temp_pocket=Overall_OxygenSinks_Pxllist_all{datai}(k,:);
        %create a logical vector that has tru values in only the frame that
        %this pocket had an event
        temp_pocket_bool= ~cellfun(@isempty,temp_pocket);

        % finding the individual events
        hpevents_temp=regionprops(temp_pocket_bool,'Area','PixelIdxList');

        if length(hpevents_temp)>thisrecording.NumOxySinkEvents(k)

            hpevents_temp=[hpevents_temp(1:thisrecording.NumOxySinkEvents(k))];

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
                Single_BW=false(Sinks_Area{datai,8});


                Single_BW(Overall_OxygenSinks_Pxllist_all{datai}{k,hpevents_temp(j).PixelIdxList(f)})=true;

                props_eventframe_single=regionprops(Single_BW,'Area','PixelIdxList','Centroid','Circularity',"FilledArea","EquivDiameter","Perimeter");


                if isempty([props_event_temp.Area])
                    
                    props_event_temp=[props_eventframe_single(:)];

                else

                    props_event_temp=cat(1, [props_event_temp(:)], [props_eventframe_single(:)]);
                end
            
            end


            centroidxy=[props_event_temp.Centroid];
            

            % capturing the results
            Area_um(counter)=mean([props_event_temp.Area])*thisrecording.Pixelsize{k}^2;
            Area_norm(counter)=(mean([props_event_temp.Area])*thisrecording.Pixelsize{k}^2)/thisrecording.RecAreaSize{k};
            FilledArea_um(counter)=mean([props_event_temp.FilledArea])*thisrecording.Pixelsize{k}^2;
            FilledArea_norm(counter)=(mean([props_event_temp.FilledArea])*thisrecording.Pixelsize{k}^2)/thisrecording.RecAreaSize{k};
            Centroid_x(counter)=mean(centroidxy(1:2:end));
            Centroid_y(counter)=mean(centroidxy(2:2:end));
            Circularity(counter)= mean([props_event_temp.Circularity]);
            Perimeter_um(counter)=mean([props_event_temp.Perimeter])*thisrecording.Pixelsize{k};
            Diameter_um(counter)=mean([props_event_temp.EquivDiameter])*thisrecording.Pixelsize{k};


            counter=counter+1;

        end

    end

    Eventspecificmetrics_rec=table(Experiment,Mouse,Condition,DrugID,Genotype,PuffStim,Promoter,Pixelsize,PocketID,EventID,...
        NormOxySinkAmp, Start, Duration, Size_modulation, Area_um,Area_norm,FilledArea_um,FilledArea_norm,Centroid_x,Centroid_y,...
        Circularity,Perimeter_um,Diameter_um);


    %add the current table to the combo table
    Eventspecificmetrics = [Eventspecificmetrics;Eventspecificmetrics_rec];
    

end

clear  temp_pocket datai wrapper TempIndex thisrecording temp_pocket_bool Single_BW props_eventframe_single props_event_temp centroidxy counter


%% Check metadata and define levels for grouping variables

if ~ischar(Eventspecificmetrics.Mouse{1})

    Eventspecificmetrics.Mouse=cellfun(@num2str,Eventspecificmetrics.Mouse,'UniformOutput',false); %converting to a string to

end

Mice_unique=uniqueStrCell(Eventspecificmetrics.Mouse);
Conditions_unique=uniqueStrCell(Eventspecificmetrics.Condition);
Drugs_unique=uniqueStrCell(Eventspecificmetrics.DrugID);
StimCond_unique=unique(Table_OxygenSinks_OutCombo.PuffStim); %with or without you

fprintf(['The dataset contains ', num2str(size(Sinks_Area,1)), ' recording sessions... \n']); 

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
% Column 1 has awake/baseline recording for ISO
% Column 2 has awake/puffs recording for ISO
% Column 3 has anesthetised/baseline recording for ISO
% Column 4 has anesthetised/puffs recording for ISO
% Column 5-8 have the repeat for KX

counter=1;
for i=1:length(Drugs_unique) %go through DrugIDs first (Between subjects variable)
    % a temp logical vector for each drug for the sinktable/and metrics

    wrapper = @(x) strcmp(x, Drugs_unique{i});
    TempOxySinkIndexDrugID=cell2mat(cellfun(wrapper,Eventspecificmetrics.DrugID(:), 'UniformOutput', false));
   
    %these indeces are meaningful for the additional metrics tables

    for k=1:length(Conditions_unique) %for each drug ID go through the conciousness condition of mice (this is a within-subject variable and therefor nested in the DrugID loop)
 
        % a temp logical vector for each condtion for the sinktable/and metrics

        wrapper = @(x) strcmp(x, Conditions_unique{k});
        TempOxySinkIndexCondition=cell2mat(cellfun(wrapper,Eventspecificmetrics.Condition(:), 'UniformOutput', false));
        
                
        for g=1:length(StimCond_unique) %for type of stimulation (with/without). This is a second level nested within-subject variable and therefore nested in the Condition loop

            % a temp logical vector for each stim condition for the sinktable/and metrics

            TempOxySinkIndexPuffs=[Eventspecificmetrics.PuffStim{:}]'==StimCond_unique(g);
      
            %multiply all temp logical vectors for sinks elementwise to get one and store
            %it in the FiltersOxysinksMetrics
            FiltersOxySinksMetrics{4,counter}=logical(TempOxySinkIndexDrugID.*TempOxySinkIndexCondition.*TempOxySinkIndexPuffs);

            %multiply all logical vectors for surgess elementwise to get one and store
            %it in the FiltersOxysurgesMetrics

            FiltersOxySinksMetrics{1,counter}= Drugs_unique{i}; 
            FiltersOxySinksMetrics{2,counter}= Conditions_unique{k};    


            if StimCond_unique(g)                   
                FiltersOxySinksMetrics{3,counter}='With stimulation';
                
            else
                FiltersOxySinksMetrics{3,counter}='Without stimulation';
                
            end

            
            FiltersOxySinksMetrics{5,counter}=unique(Eventspecificmetrics.Mouse(FiltersOxySinksMetrics{4,counter}));    

            counter=counter+1;
        end
    end
end

clear i k g counter wrapper TempOxySinkIndexDrugID TempOxySurgeIndexDrugID TempOxySinkIndexCondition TempOxySurgeIndexCondition TempOxySinkIndexPuffs TempOxySurgeIndexPuffs



%% Now using the filters generated to divide the dataset and prepare the tables for export

FilteredSinkData_mice=cell(1,size(FiltersOxySinksMetrics,2));

FilteredSinkData=cell(11,size(FiltersOxySinksMetrics,2));

for i=1:size(FiltersOxySinksMetrics,2)

    FilteredSinkData_mice{1,i}=Eventspecificmetrics.Mouse(FiltersOxySinksMetrics{4,i});
    
    FilteredSinkData{1,i}=Eventspecificmetrics.NormOxySinkAmp(FiltersOxySinksMetrics{4,i});
    FilteredSinkData{2,i}=Eventspecificmetrics.Start(FiltersOxySinksMetrics{4,i});
    FilteredSinkData{3,i}=Eventspecificmetrics.Duration(FiltersOxySinksMetrics{4,i});
    FilteredSinkData{4,i}=Eventspecificmetrics.Size_modulation(FiltersOxySinksMetrics{4,i});
    FilteredSinkData{5,i}=Eventspecificmetrics.Area_um(FiltersOxySinksMetrics{4,i});
    FilteredSinkData{6,i}=Eventspecificmetrics.Area_norm(FiltersOxySinksMetrics{4,i});
    FilteredSinkData{7,i}=Eventspecificmetrics.FilledArea_um(FiltersOxySinksMetrics{4,i});
    FilteredSinkData{8,i}=Eventspecificmetrics.FilledArea_norm(FiltersOxySinksMetrics{4,i});
    FilteredSinkData{9,i}=Eventspecificmetrics.Circularity(FiltersOxySinksMetrics{4,i});
    FilteredSinkData{10,i}=Eventspecificmetrics.Perimeter_um(FiltersOxySinksMetrics{4,i});
    FilteredSinkData{11,i}=Eventspecificmetrics.Diameter_um(FiltersOxySinksMetrics{4,i});

    

end

clear i


%% Getting data ready for exporting

if exist('Epochs') % if the dataset contains recordings with defined epochs
    
    ExportReady=cell(size(FilteredSinkData,1),1);
    for i=1:size(FilteredSinkData,1) %for every metric

        %create temp cell arrays that have the length of the recording with the most sinks/surges and coloumn equal with the number of groups*the maximum number of mice that exist in a group*the number of epochs    
        %this prealocation ensures that there will definitelly be enough space
        %in all dimensions to wrangle and rearrange data
        tempSinkmetric=cell(max(cellfun(@length,FilteredSinkData),[],'all'),length(Drugs_unique)*length(Conditions_unique)*length(StimCond_unique)*max(cellfun(@length,FiltersOxySinksMetrics(5,:)),[],'all')*Epochs); 
    
        counter=1;
        for k=1:size(FilteredSinkData,2) %go through the different groups 
    
            if ~isempty(FilteredSinkData{i,k}) %if sink data for that group exists
                    
                for jj=1:length(FiltersOxySinksMetrics{5,k}) %go the different mice of that group
            
                    for epch=1:length(Epochs_Starts) %for each distinct epoch
                        wrapper = @(x) strcmp(x, FiltersOxySinksMetrics{5,k}{jj});
                        %create a filter for the data of that mouse for the epoch of interest in the
                        %group of this iteration
                        tempfiltermouse_epoch=logical(cell2mat(cellfun(wrapper,FilteredSinkData_mice{1,k}, 'UniformOutput', false)).*(FilteredSinkData{2,k}>=Epochs_Starts(epch)&FilteredSinkData{2,k}<Epochs_Ends(epch))); 
                
                        tempSinkmetric(1,counter)= FiltersOxySinksMetrics(1,k);
                        tempSinkmetric(2,counter)= FiltersOxySinksMetrics(2,k);
                        tempSinkmetric(3,counter)= FiltersOxySinksMetrics(3,k);
                        tempSinkmetric(4,counter)= {[FiltersOxySinksMetrics{5,k}{jj},'_','epoch_', num2str(epch)]};                                          
                        tempSinkmetric(5:4+length(FilteredSinkData{i,k}(tempfiltermouse_epoch)),counter)=num2cell(FilteredSinkData{i,k}(tempfiltermouse_epoch));

                        counter=counter+1;
                    end
                end

            end

        end

        %clean up empty cells 
        tempSinkmetric=tempSinkmetric(any(~(cellfun(@isempty,tempSinkmetric)),2),:);
        tempSinkmetric=tempSinkmetric(:,any(~(cellfun(@isempty,tempSinkmetric)),1));


        ExportReady{i,1}=cell2table(tempSinkmetric);    


    end

else
    ExportReady=cell(size(FilteredSinkData,1),1);
    for i=1:size(FilteredSinkData,1) %for every metric

        %create temp cell arrays that have the length of the recording with the most sinks/surges and coloumn equal with the number of groups*the maximum number of mice that exist in a group    
        %this prealocation ensures that there will definitelly be enough space
        %in all dimensions to wrangle and rearrange data
        tempSinkmetric=cell(max(cellfun(@length,FilteredSinkData),[],'all'),length(Drugs_unique)*length(Conditions_unique)*length(StimCond_unique)*max(cellfun(@length,FiltersOxySinksMetrics(5,:)),[],'all')); 
    
        counter=1;
        for k=1:size(FilteredSinkData,2) %go through the different groups 
    
            if ~isempty(FilteredSinkData{i,k}) %if sink data for that group exists
                    
                for jj=1:length(FiltersOxySinksMetrics{5,k}) %go the different mice of that group
            
                    wrapper = @(x) strcmp(x, FiltersOxySinksMetrics{5,k}{jj});
                    tempfiltermouse=cell2mat(cellfun(wrapper,FilteredSinkData_mice{1,k}, 'UniformOutput', false)); %create a filet for the data of that mouse in the group of this iteration
                
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


        ExportReady{i,1}=cell2table(tempSinkmetric);    


    end

end


clear i k tempSinkmetric tempSurgemetric jj tempfiltermouse

%% (7) Exporting results

disp('Creating folder for sorted data and statistical analysis matrices');
Statsoutputfolder=['HypoxicEventSpecificMetrics_Output_',datestr(now,30)]; %adds an identifier in the form of yyyymmddTHHMMSS
mkdir(Statsoutputfolder)

%saving almost everything in a mat file

save([Masterfolder,'\',Statsoutputfolder,'\','HypoxicEventSpecificMetrics4LME.mat'],...
    'Table_OxygenSinks_OutCombo','FiltersOxySinksMetrics','Eventspecificmetrics');

cd(Statsoutputfolder)
Output_xlsx=sprintf('SortedSinkEventMetrics.xlsx');
writetable(ExportReady{1,1},Output_xlsx,'Sheet','NormOxySinkAmp','WriteVariableNames',0);
writetable(ExportReady{3,1},Output_xlsx,'Sheet','Duration','WriteVariableNames',0);
writetable(ExportReady{4,1},Output_xlsx,'Sheet','Size_modulation','WriteVariableNames',0);
writetable(ExportReady{5,1},Output_xlsx,'Sheet','Area_um','WriteVariableNames',0);
writetable(ExportReady{6,1},Output_xlsx,'Sheet','Area_norm','WriteVariableNames',0);
writetable(ExportReady{7,1},Output_xlsx,'Sheet','FilledArea_um','WriteVariableNames',0);
writetable(ExportReady{8,1},Output_xlsx,'Sheet','FilledArea_norm','WriteVariableNames',0);
writetable(ExportReady{9,1},Output_xlsx,'Sheet','Circularity','WriteVariableNames',0);
writetable(ExportReady{10,1},Output_xlsx,'Sheet','Perimeter_um','WriteVariableNames',0);
writetable(ExportReady{11,1},Output_xlsx,'Sheet','Diameter_um','WriteVariableNames',0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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

%% not used bits of code

% %this will contain putative unique pockets(rows) and the pixels belonging to events for every frame of the recording.
% Overall_OxygenSinks_Pxllist_all=cell(1,length(OxygenSinksInfo_all));
% sink_overlap_thres=0.6; %percentage of overlap that is needed between sinks in two different frames
% 
% for datai=1:length(OxygenSinksInfo_all)
% 
%     Overall_OxygenSinks_Pxllist_rec=cell(1,length(OxygenSinksInfo_all{datai}));
%     counter=1;
% 
%     for i=1:length(OxygenSinksInfo_all{datai})-3 %for every video frame 
% 
%         for ij=1:length(OxygenSinksInfo_all{datai}{i}) %%for every pocket in that frame
% 
%             if ~isnan(OxygenSinksInfo_all{datai}{i}(ij).PixelIdxList)
% 
%                 %create a temporary variable containing the pixels (linear
%                 %indexing) of this oxygen sink frame
%             
%                 pixl_temp=OxygenSinksInfo_all{datai}{i}(ij).PixelIdxList;
%                 
%                 Overall_OxygenSinks_Pxllist_rec{counter,i}=OxygenSinksInfo_all{datai}{i}(ij).PixelIdxList; %start an entry on that identified pocket-frame
% 
%                 %and then...
%                 for f=i+1:length(OxygenSinksInfo_all{datai})%...sweep forward to all the subsequent frames
%                         
%                     for fj=1:length(OxygenSinksInfo_all{datai}{f}) % and check the different oxygen sinks in those frames. Again going from smaller to larger
%             
%                         if ~isnan(OxygenSinksInfo_all{datai}{f}(fj).PixelIdxList) %if a subsequent oxygen sink frame has not been removed from the pool 
%                         
%                             %get the number of pixels belonging to the start pocket that correspond to the overlap thershold
%                             F_overlap=length(pixl_temp)*sink_overlap_thres; 
%                             %get the number of pixels belonging to the current pocket that correspond to the overlap thershold
%                             NextF_overlap=length(OxygenSinksInfo_all{datai}{f}(fj).PixelIdxList)*sink_overlap_thres; 
%                         
%                             if length(intersect(pixl_temp,OxygenSinksInfo_all{datai}{f}(fj).PixelIdxList))>NextF_overlap || length(intersect(pixl_temp,OxygenSinksInfo_all{datai}{f}(fj).PixelIdxList))>F_overlap    
%                                 %and if there is an significant overlap between oxygen sink frames of interest and subsequent oxygen sink frame... 
% 
%                                 %add the pixel list of that pocket-frame in the list
%                                 % I use vercat in the eventuality the preceding
%                                 % sink frame(s) overlaps with more than one sink
%                                 % area (fj) in the current frame (f)
%                                 Overall_OxygenSinks_Pxllist_rec{counter,f}=vertcat(Overall_OxygenSinks_Pxllist_rec{counter,f},OxygenSinksInfo_all{datai}{f}(fj).PixelIdxList);                                   
%                                                                             
%                                 %Now I get the unique pixels of frames that this oxygen sink exists and check the overlap with oxygen sinks in subsequent frames                                                     
%                                 %The reason I do that is that an oxygen sink in                                                      
%                                 %frames n and n+3 might not have overlap above                                                       
%                                 %the minimum but n+2 and n+3 may have. There is a                              
%                                 %drift in the centroid of the pocket.                    
%                                 pixl_temp=unique(vertcat(pixl_temp,OxygenSinksInfo_all{datai}{f}(fj).PixelIdxList));
%   
%                                 OxygenSinksInfo_all{datai}{f}(fj).PixelIdxList=NaN; %removing that pocket-frame from the pool
% 
%                             end
%                         end
%                     
%                     end
%                 end
%   
%                 OxygenSinksInfo_all{datai}{i}(ij).PixelIdxList=NaN; %remove that pocket-frame from the pool before going to the next
%             
%                 counter=counter+1; %count up for my indexing of the identified pockets
% 
%             end
%         end
% 
%         Overall_OxygenSinks_Pxllist_all{datai}=Overall_OxygenSinks_Pxllist_rec;
%     end
% 
% 
% end
% 
% clear OxygenSinksInfo_all counter pixl_temp
