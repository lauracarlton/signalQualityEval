

% go to the NN22 Data Transfer/date folder
% Specify the path to the directory
directory_path = '/Users/lauracarlton/Library/CloudStorage/GoogleDrive-lcarlton@bu.edu/.shortcut-targets-by-id/1Qja7XXGiJbseXCER2cQ-c7FviFH3swGF/NN22 Data Transfer/240329';

cd /Users/lauracarlton/Documents/Homer3/
setpaths
cd '/Users/lauracarlton/Library/CloudStorage/GoogleDrive-lcarlton@bu.edu/.shortcut-targets-by-id/1Qja7XXGiJbseXCER2cQ-c7FviFH3swGF/NN22 Data Transfer/convertBin2Snirf'
path(path,cd)
cd(directory_path)



%%
desired_extension = '.bin'; 

% Get a list of files in the directory
file_list = dir(fullfile(directory_path, ['*' desired_extension]));

% loop through each file and convert bin to snirf 
for f = 1:length(file_list)
    
    [~, name, ~] = fileparts(file_list(f).name);
    disp(name)

    snirf = convertBintoSnirfv2(name, 1)

end







