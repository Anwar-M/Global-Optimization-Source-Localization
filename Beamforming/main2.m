% Script for tests and single images
clear all; clc;
c = 343;
pref = 2e-5;

file_path{1} = 'O:\Anechoic Chamber Experiment\DATA\2016-06-21_16-50-56';
data.dir = file_path{1};
data.calibration_folder = 'O:\Anechoic Chamber Experiment\Anechoic Chamber\calibration';
[info, config, data] = Read_data(data);
data.sf = info.sf;
data.x = config.x; data.y = config.y; data.z = config.z;
clear config info;

data_check(:,1) = 1:64;
data_check(:,2) = mean(data.file_full);
data_check(:,3) = abs(data_check(:,2)) > 1e-14;
data_check(:,4) = rms(data.file_full);
data_check(:,5) = (data_check(:,4) > 0.1)+(data_check(:,4) < 0.02);
data_check(:,6) = data_check(:,3)+data_check(:,5);
disp(['Defect mics: ' num2str(sum(data_check(:,6)>0))]);
data_check(:,7) = 20*log10(data_check(:,4)/2e-5);

remove_mics = (data_check(:,1).*(data_check(:,6)>0));
select_mics = data_check(:,1).*(remove_mics == 0);
remove_mics(remove_mics == 0) = [];
select_mics(select_mics == 0) = [];

N_mics = length(select_mics);

minval = min(data.file_full(:));
maxval = max(data.file_full(:));

data.x(remove_mics) = [];
data.y(remove_mics) = [];
data.z(remove_mics) = [];
data.z = -data.y*sind(4);
data.file_full(:,remove_mics) = [];


t_start = 2;
t_end = 30;

fcenter = 6100;
f_low = fcenter-100;%/2^(1/6);
f_high = fcenter+100;%*2^(1/6);
[CSM, freqs] = developCSM(data.file_full, f_low, f_high, data.sf, 0.01, .5, t_start, t_end);

%% Beamforming
method = 3;
xmin = -.8; xmax = .8;
ymin = -.8; ymax = .8;
resolution = 0.01; % 1 cm
z = 1.87;

switch method
    case 1
        [x_t_x, x_t_y, A] = FastBeamforming3(CSM, z, freqs, ...
           [xmin xmax ymin ymax], resolution, [data.x; data.y; data.z], c);
    case 2
        [x_t_x, x_t_y, A] = CleanPSF(CSM, z, freqs, ...
            [xmin xmax ymin ymax], resolution, [data.x; data.y; data.z], c);
    case 3
        [x_t_x, x_t_y, A] = CleanSC(CSM, z, freqs, ...
            [xmin xmax ymin ymax], resolution, [data.x; data.y; data.z], c);
    case 4
        
        x_t_x = xmin:resolution:xmax;
        x_t_y = ymin:resolution:ymax;
        settings.g_exponent = 4;
        [S,~,~,~] = BF_functional(data.file_full, data.sf, c, ...
                                  x_t_x, x_t_y, z, f, 100, ...
                                  [data.x; data.y; data.z].', settings);
    otherwise
        error('wtf man..?');
end

SPL = 20*log10( sqrt(2*0.5*real(A)) / 2e-5 );
    
maxval = max(SPL(:));
%% Finish Images and overlays

dynamic_range = 2*12/2; % dB

maxval = ceil(max(SPL(:)));
minval = maxval - dynamic_range + 1;
figure;
SPL(SPL < minval - .5) = 0;
imagesc(x_t_x, x_t_y, SPL, [minval-.5 round(maxval)]);

set(0,'defaulttextinterpreter','latex');
hXLabel = xlabel('$x$ [m]');
hYLabel = ylabel('$y$ [m]');
set([hXLabel, hYLabel], ...
    'FontName', 'AvantGarde', 'FontSize', 14);

set(gca, ...
'YDir','normal', ...
'XTick', linspace(xmin, xmax, 5), ...
'YTick', linspace(ymin, ymax, 5), ...
'XColor'      , [.3 .3 .3], ...
'YColor'      , [.3 .3 .3], ...
'Fontsize'    , 16, ...
'FontName'   , 'Helvetica' );
set(gca,'xgrid', 'on', 'ygrid', 'on', 'gridlinestyle', ':', 'xcolor', 'k', 'ycolor', 'k');

axis equal; axis([xmin xmax ymin ymax]);
cb = colorbar;
jetmod = jet;
jetmod(1,:) = 1;
colormap(jetmod);
cb.Limits = [minval maxval];
