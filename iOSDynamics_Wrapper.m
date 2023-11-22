%% Wrapper script to run the iOS Dynamics analysis in all datafolders with iOS data. 
% This is a modified version of the Oxygen Dynamics analysis 

% The input is a CSV file with the relative paths of the data
% folders. For example for current folder 'C:\User\Alldata\' that contains folders 'Mouse1' and 'Mouse2' etc
% the paths should be 'Mouse1\Recording1' 'Mouse1\Recording2' 'Mouse2\Recording1' etc
 
% The same CSV file can contain other metadata related to each recording,
% including, the animal genotype the type of recording (e.g awake vs aneasthetised), the sampling frequency or the type of lense (that will affect pixel size) etc.
% Check the example file 

% The wrapper script goes through the data and calls the main OxygenDynamics analysis script to analyse sinks and surges of oxygen
% If a recording comes with behavioural data, the wrapper also calls a behavioural analysis script 


% The wrapper and OxygenDynamics share a workspace, in that way the metadata can be included in the output table of
% each data session and make statistical analysis easier. 

% IMPORTANT!!
% Please use underscores instead of spaces in folder and file names of all data

% The names of the three behavioural data files should be stated without their extention. 
% It is expected that they will be csv files for pupil and movement and a abf file for the air puffs
% These files MUST be in the same folder where the imaging data file is. 

% As always, remember to add the folder with the scripts in the path.

% From Antonis Asiminas, PhD
% Update 26 Sep 2023

%%
InputD=readtable('iOS_DataPathsExample.csv','Delimiter',',');

%For some reason sometimes the csv is not being read correctly and it InputD ends up being
%a single columns of all info This looks like it is an Excel thing..
if size(InputD,2)<2
    InputD=readtable('iOS_DataPathsExample.csv');
end

%%
answer1 = questdlg('Do you want to analyse data you have analysed before?',...
    'What data',...
    'Yes', 'No','No');

switch answer1

    case 'Yes'
    
        strAgain = 'Y';
        answer2 = questdlg('Do you want to overwrite previous analysis (if any)?','Yes', 'No');

        switch answer2
            case 'Yes'
                strOW = 'Y';
            case 'No'
                strOW = 'N';
        end

    case 'No'
        
        strAgain = 'N';
        strOW = 'N';
                      
end

clear answer1 answer2

%% This is another prompt to ask if you just want the itf conversioning to deltaF/F process only

answer1 = questdlg('Do you ONLY want to convert the raw tif flies to dF/F tifs?',...
    'What data',...
    'Only df/f tifs', 'All analysis','All analysis');

%% Email alerts info
% INFO!!
% the initial idea was to send an email though a gmail account I generated
% but the KU firewall settings did not allow it. I found a work around that
% involves the use of outlook. Thus 

% User input
%prompt = {'Enter email alerts recipient:','Enter password for matlabalertctn email:'};
%dlgtitle = 'Input';
%dims = [1 35];
%definput = {'felix.beinlich@sund.ku.dk','PW'};
%source = 'matlabalertctn@gmail.com';  %from address (gmail)
% destination = answer{1};              %to address (any mail service)
% myEmailPassword = answer{2};          %the password to the 'from' account

prompt = {'Enter email alerts recipient:'};
dlgtitle = 'Input';
dims = [1 35];
definput = {'felix.beinlich@sund.ku.dk'};
answer = inputdlg(prompt,dlgtitle,dims,definput);
destination = answer{1};
%the the function at the bottom of the script if you want to specify an
%account you want to send the email from. This needs to be already
%connected to your outlook!

% %set up SMTP service for Gmail
% setpref('Internet','E_mail',source);
% setpref('Internet','SMTP_Server','smtp.gmail.com');
% setpref('Internet','SMTP_Username',source);
% setpref('Internet','SMTP_Password',myEmailPassword);
% 
% % Gmail server.
% props = java.lang.System.getProperties;
% props.setProperty('mail.smtp.auth','true');
% props.setProperty('mail.smtp.socketFactory.class', 'javax.net.ssl.SSLSocketFactory');
% props.setProperty('mail.smtp.socketFactory.port','465');

% For more info see 
% https://se.mathworks.com/matlabcentral/answers/346491-use-of-sendmail-function-to-send-an-email-from-a-gmail-account
% https://se.mathworks.com/matlabcentral/answers/433031-sending-e-mail-error#answer_349841
%% Separate the information contained in the table

Paths=table2cell(InputD(:,{'Paths'}));
Postures=table2cell(InputD(:,{'PostureFile'}));
Pupils=table2cell(InputD(:,{'PupilFile'}));
Puffs=table2cell(InputD(:,{'PuffsFile'}));
Mice=table2cell(InputD(:,{'Mouse'}));
Genotypes=table2cell(InputD(:,{'Genotype'}));
Conditions=table2cell(InputD(:,{'Condition'}));
DrugIDs=table2cell(InputD(:,{'DrugID'}));
Promoters=table2cell(InputD(:,{'Promoter'}));
SampleFs=table2cell(InputD(:,{'SampleF'}));
Pixelsizes=table2cell(InputD(:,{'Pixelsize'}));

%% Creating an empty table to get all the tables from each recording created in the master script.
% this is only for the oxygen surges since the oxygen sink data need to be refined 

Masterfolder=pwd;

% Create two structures to catch errors from the oxygen dymanics and
% behavioural metrics' analyses
iOSDynamicsfailures = struct('FileName', {}, 'ERROR', {}, 'MESSAGE', {});    
Behaviourfailures = struct('FileName', {}, 'ERROR', {}, 'MESSAGE', {});

for datai=1:length(Paths)
     
    cd(Paths{datai}); % go through the different paths and get the info for the session of the iteration
    SFs = SampleFs{datai};
    Mous = Mice{datai};
    Cond = Conditions{datai};
    PiSz = Pixelsizes{datai};
    Gen = Genotypes{datai};
    Promo = Promoters{datai};
    Drug = DrugIDs{datai};
    Posture = Postures{datai};
    Pupil = Pupils{datai};
    Puff = Puffs{datai};


   
    switch answer1

                
        case 'Only df/f tifs'

            if not(isfolder('OxygenSinks_Output')) || strAgain == 'Y' || strOW == 'Y' %the last one is redundant. There is a check in the master analysis script
            
                run iOS_Tiffout % call the analysis script

            end

            
        case 'All analysis'
    
        
            if not(isfolder('OxygenSinks_Output')) || strAgain == 'Y' || strOW == 'Y' %the last one is redundant. There is a check in the master analysis script

            
                try

                    run iOSDynamics_Master % call the analysis script        

                catch ME %catch the error message
              
            
                    errorMessage = sprintf('Error in OxygenDynamics_Master.m.\nThe error reported by MATLAB is:\n\n%s', ME.message);      

                    subj = ['Error alert in OxygenDynamics_Master when analysing ',Paths{datai}];  % subject line    
          
                    %FIREWALL SETTINGs DO NOT PERMIT THIS IN DESKTOP COMPUTERS PLUGGED IN THE NETWORK (comment when if this produces an ERROR)                         
                    sendolmail(destination,subj,errorMessage);                                               
                    iOSDynamicsfailures(end + 1).FileName = [Paths{datai}];                                         
                    iOSDynamicsfailures(end).ERROR  = getReport(ME);         
                    iOSDynamicsfailures(end).MESSAGE  = [ME.message];              
                    uiwait(warndlg(errorMessage));
            
        
                end
        
    
            end

          
    
            if not(isfolder('Behaviour_Output')) || strAgain == 'Y' 
                     
        
                try
                                     
                    run OxygenDynamics_Behaviour % call the behavioural data analysis script
                   
        
                catch ME %catch the error message
                             
                    errorMessage = sprintf('Error in OxygenDynamics_Behaviour.m.\nThe error reported by MATLAB is:\n\n%s', ME.message);                                                 
                           
                    subj = 'Error alert for OxygenDynamics_Behaviour';  % subject line                       
                            
                    %FIREWALL SETTINGs DO NOT PERMIT THIS IN DESKTOP COMPUTERS PLUGGED IN THE NETWORK (comment when if this produces an ERROR)                                
                    sendolmail(destination,subj,errorMessage);
                
                    Behaviourfailures(end + 1).FileName = [Paths{datai}];                             
                    Behaviourfailures(end).ERROR  = getReport(ME);                            
                    Behaviourfailures(end).MESSAGE  = [ME.message];
                                 
                    uiwait(warndlg(errorMessage));
                       
        
                end
     
            end

    end
    cd(Masterfolder) %return to the masterfolder before going to the next data folder
             
 
    clearvars -except datai Masterfolder Paths Mice DrugIDs Conditions Genotypes Promoters SampleFs Pixelsizes Postures Pupils Puffs strOW strAgain InputD source...
        destination myEmailPassword setpref props answer1 iOSDynamicsfailures Behaviourfailures
end


% If there were any failures, save the files that contain the info
if ~isempty(iOSDynamicsfailures)

    save(iOSDynamicsfailures.mat','iOSDynamicsfailures');
end

if ~isempty(Behaviourfailures)
    save('Behaviourfailures.mat','Behaviourfailures');
end

subj = 'Message from Matlab';  % subject line
msg = 'iOSDynamics_Wrapper finished running successfully';     % main body of email.   
sendolmail(destination,subj,msg);
%FIREWALL SETTINGs DO NOT PERMIT THIS IN DESKTOP COMPUTERS PLUGGED IN THE NETWORK (uncomment when that has been addressed)
%sendmail(destination,subj,msg); 
% % Remove the preferences (for privacy reasons)
% setpref('Internet','E_mail','');
% setpref('Internet','SMTP_Server','''');
% setpref('Internet','SMTP_Username','');
% setpref('Internet','SMTP_Password','');

%% Function that opens outlook and sends email
% for more info see
% https://se.mathworks.com/matlabcentral/answers/94446-can-i-send-e-mail-through-matlab-using-microsoft-outlook
function sendolmail(to,subject,body,attachments)
%Sends email using MS Outlook. The format of the function is 
%Similar to the SENDMAIL command.
% Create object and set parameters.
h = actxserver('outlook.Application');
mail = h.CreateItem('olMail');
mail.Subject = subject;
mail.To = to;
mail.BodyFormat = 'olFormatHTML';
mail.HTMLBody = body;
%mail.SentOnBehalfOfName = "example@email.com"; %change to the email you
%want send from. For example if you have a matlab alerts account
% Add attachments, if specified.
if nargin == 4
    for i = 1:length(attachments)
        mail.attachments.Add(attachments{i});
    end
end
% Send message and release object.
mail.Send;
h.release;
end
