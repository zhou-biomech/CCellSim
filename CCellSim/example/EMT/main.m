close all; clear;

% Add PATH reference in order to run solver
addpath('../../src/');

% folder and path
path = './data';
if ~exist(path); mkdir(path); end

% cells and moleculars
totgap = 0.4;
[cells, molnames, bdrys] = init(totgap);

clear cells
load([path,'/','600.mat']);
delete([cells.h]); delete([cells.hc]);

figure('position',[50,50,1000,1000])
bdrys.plot; 

for i = 1:length(molnames); eval([molnames{i}, '=', num2str(i),';']); end
cells.plot('pip3'); colorbar;
axis image;axis([-100,100,-100,100]); drawnow

% max time and time step
tmax = 1200;
dt = min([cells(1).mols.dt])/5; dt = 1/ceil(1/dt);

% Initialize progress bar:
upd = progbar(tmax);

for t = 0:dt:tmax
    
    if t>100
        for c = cells
            c.mols.updligdsrc(dt);
            c.mols.updnoligsrc(dt,t);
            c.mols.diffusion(dt);
            c.csk2mol(dt)
        end
    end
    
    for c = cells
        c.force();
    end

    if t<=200
        for c = cells
            e = mean(c.p) - [0,0];
            e = e/vecnorm(e);
            c.f = c.f - e*20;
        end
    end

    for c = cells
        c.cellbdry(bdrys);
        c.update(dt);
    end
    
    cells.cellcell();
    
    for c = cells
        c.cellcsk();
    end
    
    
    if mod(t,10)==0
        cells.plot(pip3); colorbar; drawnow
        frame = getframe(gcf);
        save([path, '/', num2str(t), '.mat'], 'cells');
        saveas(gcf, [path, '/', num2str(t)], 'jpg')
        upd(t);
    end
end