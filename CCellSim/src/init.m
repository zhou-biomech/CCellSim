function [cells, molnames, bdrys] = init(totgap)
mol = config('molbasic');

molnames = mol(:,2);   domain = mol(:,3);     ligd = mol(:,4);       
    rcpt = mol(:,5);      kon = mol(:,6);     koff = mol(:,7);   
     a2f = mol(:,8); initconc = mol(:,9); dfucoeff = mol(:,10); 

nmols = length(molnames);
for i = 1:length(molnames); eval([molnames{i}, '=', num2str(i),';']); end

% -------------------------------------------------------------------------

cel = config('celbasic');

for j = 1:size(cel,1)
    cells(j) = Cell(cel(j,1), cel(j,2), cel(j,3), cel(j,4), totgap);
    
    for i = 1:nmols
        cells(j).mols(i) = Molecule(molnames{i}, domain{i}, initconc{i},...
                                    dfucoeff{i},  cells(j), a2f{i});
    end
    
    for i = 1:nmols
        cells(j).mols(i).ligd = cells(j).mols(eval(ligd{i}));
        cells(j).mols(i).rcpt = cells(j).mols(eval(rcpt{i}));
        cells(j).mols(i).kon  = kon{i};
        cells(j).mols(i).koff = koff{i};
    end
end

% -------------------------------------------------------------------------

bdrys = config('bdryset');