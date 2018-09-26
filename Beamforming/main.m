% Script for tests and single images
clear all; clc;

addpath('..\');
CsmReaderb7a;
c = csound;

save_data = 0;
pref = 2e-5;

data.x = mic_info(1,:); data.y = mic_info(2,:); 
N_mics = size(mic_info,2);
data.z = mic_info(3,:); 

fcenter = 10000;
f_low = fcenter-50;%-3e2;%/2^(1/6);
f_high = fcenter+50;%+3e2;%*2^(1/6);
fsel = (freqs>=f_low).*(freqs<=f_high);
fselind = find(fsel);
freqs = freqs(boolean(fsel));
CSM = zeros(N_mics,N_mics,numel(freqs));
for I = 1:numel(freqs)
    CSM(:,:,I) = squeeze(cpreal(fselind(I),:,:) + 1i*cpimag(fselind(I),:,:));
end

%% Beamforming
method = 3;
xmin = -.25; xmax = .25;
ymin = -.25; ymax = .25;
resolution = 0.01; % 1 cm
z = 0.75;

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
minval = maxval - dynamic_range;
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
