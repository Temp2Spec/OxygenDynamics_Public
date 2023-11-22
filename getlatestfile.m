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