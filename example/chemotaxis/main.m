close all; clear;


% Add PATH reference in order to run solver
addpath('../../src/');

% folder and path
path = ['./data'];
if ~exist(path); mkdir(path); end 

% cells and moleculars
totgap = 0.4;
[cells, molnames, bdrys] = init(totgap);


figure('position',[50,50,1000,1000])
bdrys.plot;

for i = 1:length(molnames); eval([molnames{i}, '=', num2str(i),';']); end
axis image; axis([-100,100,-100,100]); drawnow

% max time and time step
tmax = 2000;
t_fMLP = 200;
dt = min([cells(1).mols.dt])/10;
dt = 1/ceil(1/dt);


% Initialize progress bar:
upd = progbar(tmax);

for t = 0:dt:tmax
    for c = cells
        c.mols(fMLP).conc = 0.04;   % uniform
        
        if t>=500
            r = sqrt(sum(cells(j).p.^2,2));
            c.mols(fMLP).conc = 0.04 + psrcdiff(r);
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
