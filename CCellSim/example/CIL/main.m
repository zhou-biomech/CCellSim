close all; clear;

% Add PATH reference in order to run solver
addpath('../../src/');

% folder and path
path = './data';
if ~exist(path); mkdir(path); end

% cells and moleculars
totgap = 0.1;
[cells, molnames, bdrys] = init(totgap);

figure('position',[50,50,1000,1000])
bdrys.plot;

for i = 1:length(molnames); eval([molnames{i}, '=', num2str(i),';']); end
cells.plot('pip3'); colorbar; 
axis image;axis([-100,100,-100,100]); drawnow


% max time and time step
tmax = 2000;
t0 = 100;
dt = min([cells(1).mols.dt])/5;
dt = 1/ceil(1/dt);

% Initialize progress bar:
upd = progbar(tmax);


for t = 0:dt:tmax
    
    for c = cells
        if t>200
            r = sqrt(sum(cells(j).p.^2,2));
            c.mols(fMLP).conc = 0.04 + psrcdiff(r);
        else
            c.mols(fMLP).conc = 0.04;
        end
        c.mols.updligdsrc(dt);
        c.mols.updnoligsrc(dt,t);
        c.mols.diffusion(dt);
        c.csk2mol(dt)
    end

    for c = cells
        c.force();
    end
    
    for c = cells
        c.cellbdry(bdrys);
    end
    
    for c = cells
        if t<=100
            c.update(dt, true);
        else
            c.update(dt);
        end
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
