function [info, config, data] = Read_data(data)

%% array configuration and measurement info
config.dir = fullfile(data.dir, 'configuration/config.txt');
if ~exist(config.dir, 'file');
    config.dir = fullfile(data.dir, 'array.txt');
end
info.dir = fullfile(data.dir, 'info.txt');

%% Read in info file
% fprintf('Read in info file\n')
fid = fopen(info.dir); % open file
info_file = textscan(fid, '%s %s %s');  % It basically reads the info.txt document and stores 3 strings per line in a cell array in this case.
fclose(fid);

% organize data.
index = find(strcmp([info_file{1,1}], 'acoustic_sample_frequency'));
acoustic_sample_frequency = str2num(info_file{1,2}{index,1});
info.sf = acoustic_sample_frequency;

index = find(strcmp([info_file{1,1}], 'number_of_microphones'));
info.N = str2double(info_file{1,2}{index,1});

index = find(strcmp([info_file{1,1}], 'amplification'));
info.amplification = [info_file{1,2}{index,1} ' ' info_file{1,3}{index,1}];

index = find(strcmp([info_file{1,1}], 'acoustic_camera'));
info.acoustic_camera = [info_file{1,2}{index,1}];

index = find(strcmp([info_file{1,1}], 'optical_camera'));
info.optical_camera = info_file{1,2}{index,1};

index = find(strcmp([info_file{1,1}], 'start_time'));
info.start_time = info_file{1,2}{index,1}; % Reference time of the start of the recording.

index = find(strcmp([info_file{1,1}], 'start_timestamp'));
info.start_timestamp = [info_file{1,3}{index,1} ' ' info_file{1,2}{index,1}];

if strcmp(info.optical_camera, 'TRUE')
    index = find(strcmp([info_file{1,1}], 'optical_frame_rate'));
    info.optical_frame_rate = str2num(info_file{1,2}{index,1});
    
    fid = fopen(char(fullfile(data.dir, 'optical_times')),'r');
    data.optical_times = fread(fid,'uint32');
    fclose(fid);
    
    index = find(strcmp([info_file{1,1}], 'picture_width'));
    info.picture_width = str2double(info_file{1,2}{index,1});
    
    index = find(strcmp([info_file{1,1}], 'picture_height'));
    info.picture_height = str2double(info_file{1,2}{index,1});
    
    index = find(strcmp([info_file{1,1}], 'optical_exposure_time'));
    info.optical_exposure_time  = str2num(info_file{1,2}{index,1}); %adapted for comma as decimal separation.
    
end

fid = fopen(info.dir); % open file
info.comments = fgetl(fid);
while ~strncmp(char(info.comments),'comments',1)
    info.comments = fgetl(fid); %
end
fclose(fid);

info.comments = char(info.comments);
info.comments = info.comments(1,10:end); %Starts storing the comments after the word 'comments'

%% Read in array configuration
% fprintf('Read in array configuration\n')

fid = fopen(config.dir);
if strcmp(info.acoustic_camera, 'camera1') || strcmp(info.acoustic_camera, 'CRIO') % For the 115 fly-over analysis the CRIO camera was used.
    C = textscan(fid, '%d%f%f');
    config.channel(1:info.N) = C{1,1}; 
    config.x(1:info.N) = C{1,2};                    % microphones x-position
    config.y(1:info.N) = C{1,3};                    % microphones y-position
else
    C = textscan(fid, '%s%d%f%f');
    config.index = C{1,1};
    config.channel(1:info.N) = C{1,2}; 
    config.x(1:info.N) = C{1,3};                    % microphones x-position
    config.y(1:info.N) = C{1,4};                    % microphones y-position    
%     config.x(1:info.N) = C{1,4};                    % microphones x-position
%     config.y(1:info.N) = -C{1,3};                    % microphones y-position    

end
fclose(fid);

config.z(1:info.N) = 0;                         % z-position (assumed 0)

%% Read in measurement data
% fprintf('Read in measurement data\n')

% check whether file is calibrated
file = char(fullfile(data.dir, 'acoustic_data.cal'));
file_info = dir(file);
if isempty(file_info) || file_info.bytes == 0 && exist(data.calibration_folder,'dir')   %In case there is not calibration from LabVIEW.
    file = char(fullfile(data.dir, 'acoustic_data'));
    fid = fopen(file,'r');
    data.file_full = fread(fid,[info.N inf],'int16',0,'l').'; % ACOUSTIC DATA.
    fclose(fid);
    
    % 32 mic camera. For the 115 fly-over analysis the CRIO camera was used.
    if strcmp(info.acoustic_camera, 'camera1') || strcmp(info.acoustic_camera, 'CRIO')
        import = importdata(fullfile(data.calibration_folder,'ADCdata.txt'),'\t',0);
        scalar = import(:,1);
        bias = import(:,2);
        
        % apply ADC calibration 
        for i = 1:8
            j = (i-1)*4 + 1;
            data.file_full(:,j:j+3) = 2*(data.file_full(:,j:j+3)*scalar(i) + ones(size(data.file_full,1),4)*bias(i))/1e9; 
        end
        
        import = importdata(fullfile(data.calibration_folder,'Microphone calibration','miccalib.txt'),'\t',0);
        Responses = import(:,2);

        % scale microphone responce and amplification
        for i = 1 : info.N
            data.file_full(:,i) = data.file_full(:,i)/Responses(config.channel(i))/13.03;
        end
        
    else    % camera2-camera3
        data.file_full = data.file_full*2.5/32768; % B to V (bits to volts)

        % apply amplification calibration
        if strcmp(info.amplification,'Low Amplification')
            disp('Low Amplification')
            data.file_full = data.file_full/3.985;
        else
%             disp('High Amplification')
            data.file_full = data.file_full/28.03;
        end

        % scale microphone response
        import = importdata(fullfile(data.calibration_folder,'Microphone calibration','miccalibBundle.txt'),'\t',0);
        Responses = import(:,2);
        k=1;
        for i = 0:7 % letter            
            for j = 0:7 % channel                
                channel = config.channel(j*8+i+1)-16; % CHANGED! IN ORDER TO READ SOME CALIBRATION DATA
                data.file_full(:,i*8+j+1) = data.file_full(:,i*8+j+1)/Responses((channel-1)*8+i+1);
                data.file_full(:,i*8+j+1) = data.file_full(:,i*8+j+1) - mean(data.file_full(:,i*8+j+1));
                channels(k,:) = [(channel-1)*8+i+1, i*8+j+1];
                k = k+1;
            end
        end
    end
    
    info.calibrated = 0;
%     disp('Warning: data calibrated using local Matlab linear correction, not LabVIEW frequency dependant version')
else %In case there is calibration from LabVIEW.
    fid = fopen(file,'r');
    if strcmp(info.acoustic_camera, 'CRIO') %For the 115 fly-over analysis the CRIO camera was used.
        data.file_full = fread(fid, [info.N inf], 'double', 0, 'l').';
    else
        data.file_full = fread(fid, [info.N inf], 'single', 0, 'l').';
    end
    fclose(fid);    
    info.calibrated = 1;
end
% 
% % Adjust dataset to chosen percentage of data
% data.length = length(data.file_full(:,1));
% data.start_sample = data.time_offset * acoustic_sample_frequency + 1; %First sample to consider (t0), selected by user.
% data.stop_sample = data.start_sample + data.time_range * acoustic_sample_frequency; %Last sample to consider (tfinal), selected by user.
% 
% if (data.stop_sample > data.length) 
%     error('Start time and duration exceeds data length.'); 
% end
% 
% data.file_sample = data.file_full(data.start_sample:data.stop_sample,:); %consider only the data from time interval of interest.
% data.ns = length(data.file_sample(:,1)); % number of samples
% 
% % fprintf('Chosen complete data subset has %.0f samples\n', data.ns)

end
