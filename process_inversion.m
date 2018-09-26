function process_inversion
close all;

data_path = '.\DATA';
flist = dir([data_path '\*.mat']);
% flist = dir('inv_*-01.mat');

save_img = 0;

n_src = 4;
xmin = -.25; xmax = .25;
ymin = -.25; ymax = .25;
n_grid = 50;
x = linspace(xmin, xmax, n_grid);
y = linspace(ymin, ymax, n_grid);
z = 0.75;

for j = 99:99%1:length(flist)
    SPL = zeros(n_grid, n_grid);
    tic;
    load([data_path '\' flist(j).name]);
    toc
    
%     load(['inversion_' num2str(f(j)) 'Hz_128.mat']);
    
    A = zeros(n_src*Nruns, 1);
    for i = 1:Nruns
        
        for k = 1:n_src
            x_val = F1((k-1)*4+1, i);
            y_val = F1((k-1)*4+2, i);

            if (x_val < xmin)||(x_val > xmax)||(y_val < ymin)||(y_val > ymax)
                continue
            end
            
            x_i = round((x_val - xmin)*(n_grid - 1)/(xmax - xmin)) + 1;
            y_i = round((y_val - ymin)*(n_grid - 1)/(ymax - ymin)) + 1;
            
            SPL(y_i, x_i) = 20*log10( sqrt(0.5*F1(k*4, i)) / (4*pi*1*2e-5) );
        end
            
    end
%     maxval = max(SPL(:));
    maxval = 60;
    minval = 54;
    figure;
    h2 = imagesc(x, y, SPL, [minval round(maxval)]);
%     selection = (SPL ~= 0);
%     set(h2, 'AlphaData', selection*0.625);
%     imagesc(x, y, SPL); 
    plot_settings('$x$ [m]', '$y$ [m]', ['$f =$ ' num2str(f) ' Hz'], ...
        [-.2 .2], [-.2 .2], 1, [1 minval maxval], save_img, ...
        ['INV-BF' num2str(f) ' Hz']);
    
    figure;
    rand_ind = randi(Nruns,10,1);
    plot(10*G2(rand_ind,:).', 'LineWidth', 1);
    plot_settings('Generation', '$E \times 10$', ['$f =$ ' num2str(f) ' Hz'], ...
        [0 Ng], 10*[0 .15], 0, 0, save_img, ['INV-E' num2str(f) ' Hz']);    
end




end

function plot_settings(xlab, ylab, titlab, xlim, ylim, ...
    ax_equal, clr_bar, save_imgs, filename)

set(0,'defaulttextinterpreter','latex');
    
hXLabel = xlabel(xlab);
hYLabel = ylabel(ylab);
hTitle = title(titlab);

if ax_equal axis equal; end;
set( gca                       , ...
    'FontName'   , 'Helvetica' );
set([hXLabel, hYLabel, hTitle], ...
    'FontName'   , 'AvantGarde');
set([hXLabel, hYLabel, hTitle]  , ...
    'FontSize'   , 14          );
set(gca, ...
  'YDir'        , 'normal' , ...
  'Box'         , 'on'     , ...
  'TickDir'     , 'in'     , ...
  'TickLength'  , [.01 .01] , ...
  'YGrid'       , 'off'      , ...
  'XLim'        , [xlim(1)-.05 xlim(2)+.05], ...
  'YLim'        , [ylim(1)-.05 ylim(2)+.05], ...
  'XTick'       , linspace(xlim(1), xlim(2), 5), ...
  'YTick'       , linspace(ylim(1), ylim(2), 5), ...
  'XColor'      , [.3 .3 .3], ...
  'YColor'      , [.3 .3 .3], ...
  'Fontsize'    , 16, ...
  'LineWidth'   , 1         );
set(gca,'xgrid', 'on', 'ygrid', 'on', 'gridlinestyle', ':', 'xcolor', 'k', 'ycolor', 'k');

if clr_bar(1)
    cb = colorbar;
    colormap(flipud(hot));
    title(cb, 'SPL (dB)', 'Fontsize', 12, 'FontName', 'Helvetica');
    set(cb, 'Fontsize', 16, 'FontName', 'Helvetica', 'YTick', linspace(clr_bar(2),clr_bar(3),7));
end

scalefactor = .85;
g = get(gca,'Position');
g(1:2) = g(1:2) + (1-scalefactor)/2*g(3:4);
g(3:4) = scalefactor*g(3:4);
set(gca,'Position',g);
    
if save_imgs
    set(gcf, 'PaperPosition', [0 0 12 10]);
    set(gcf, 'PaperSize', [12 10]);
    print('-dpng', '-r600', ['.\imgs\' filename]);
    close;
end

end