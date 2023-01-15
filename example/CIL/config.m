% /home/zhou/Dropbox/software2Dcells/multicells20220629/data

function out = config(type, varargin)
switch type
    case 'celbasic'; out = celbasic;
    case 'molbasic'; out = molbasic;
    case 'bdryset' ; out = bdryset;
    case 'rcptsrc' ; [mols, dt] = varargin{:}; out = rcptsrc(mols, dt);
    case 'noligsrc'; [mols, dt, t] = varargin{:}; out = noligsrc(mols,dt,t);
    case 'csk2mol' ; [ cel, dt] = varargin{:}; out = csk2mol(cel, dt);
end

%% ------------------------------------------------------------------------
function [cel, nndes] = celbasic()
nndes = 45;
% cell position, radius and number of nodes
% cel = [-5, -5, 5, nndes
%         5,  5, 5, nndes];

n = 15;
r1 = 30;
r2 = 60;

t1 = [0:360/6:360-360/6]' +  360/18;
t2 = [0:360/9:360-360/9]';

cel = [[r1.*cosd(t1) r1.*sind(t1); r2.*cosd(t2) r2.*sind(t2)], ones(n,2).*[5, nndes]];

%% ------------------------------------------------------------------------

function bdrys = bdryset()
bdrys = Bdry.empty;

%% ------------------------------------------------------------------------

function mol = molbasic()
[~, nndes] = celbasic();
GPCRtot = 4000;     % number of GPCRs
GPCRn= round(GPCRtot/nndes)/100;
mol = {...
%id molnames     domain   ligd    rcpt           kon  koff  a2f  initconc dfucoeff
 1 'pip2'         'mebre'  'pten'  '[]'            10    1  -10     10    1
 2 'pip3'         'mebre'  'pi3k'  '[]'            10    1   10   0.20    1
 3 'pi3k'         'plasm'  '[]'    'pip3'         NaN  NaN    0   0.06   10
 4 'pten'         'plasm'  '[]'    'pip2'         NaN  NaN    0   0.07   10
 5 'Rac'          'mebre'  '[]'    '[]'           NaN  NaN    0      3  0.1
 6 'Rho'          'mebre'  '[]'    '[]'           NaN  NaN    0   1.25  0.1
 7 'FilGAP'       'plasm'  '[]'    '[]'           NaN  NaN    0      0   10
 8 'Rac_in'       'mebre'  '[]'    '[]'           NaN  NaN    0    7.5   10
 9 'Rho_in'       'mebre'  '[]'    '[]'           NaN  NaN    0      3   10
10 'activator'    'mebre'  '[]'    '[]'            20    1    0      0    1
11 'inhibitor_m'  'mebre'  '[]'    '[]'            20    1    0      0    1
12 'inhibitor'    'mebre'  '[]'    '[]'           NaN  NaN    0      0   10
13 'fMLP'         'mebre'  '[]'    'GPCR'         NaN  NaN    0      0    0
14 'GPCR'         'mebre'  'fMLP'  '[]'        1/3e-2    1    0  GPCRn    0
15 'GEF_Rac'      'plasm'  '[]'    'activator'    NaN  NaN    0   0.07   10
16 'GEF_Rho'      'plasm'  '[]'    'inhibitor_m'  NaN  NaN    0   0.07   10
};

%% ------------------------------------------------------------------------

function mols = noligsrc(mols, dt, t)
isadh = mols(1).cell.isadh;
for i = 1:length(mols); eval([mols(i).name,'= mols(',num2str(i),');']); end

kpe = 15; kpt = 15;
ke = kpe*pip2.conc./(6+pip2.conc);
kp = kpt*pip3.conc./(6+pip3.conc);
ka = Rac.conc/6;
kr = Rho.conc/1.5;

% pip2
pip2.src = -ke.*ka.*pi3k.eff + kp.*kr.*pten.eff;
pip2.src = pip2.src./pip2.a2n*dt;

% pip3
pip3.src = -kp.*kr.*pten.eff + ke.*ka.*pi3k.eff;
pip3.src = pip3.src./pip3.a2n*dt;

% ----------------------

% inhibitor
kI = 3;% associate rate 3 uM-1.s-1
k_source = 1;%s-1
inhibitor.src = k_source*fMLP.eff-kI*inhibitor.conc;
inhibitor.src = inhibitor.src./inhibitor.a2n*dt;

% activator
delt_a=0.2;     % s-1 0.2
k_AI=100;       %uM-1.S-1              100 is ref
activator.src = k_source*fMLP.eff - delt_a*activator.amt*activator.a2n -...
                k_AI*inhibitor_m.conc.*activator.amt*activator.a2n;%s-1 mols
activator.src = activator.src./activator.a2n*dt;
        
% inhibitor
delt_I = 0.2;%s-1                             %100 is ref
inhibitor_m.src = -k_AI*activator.conc.*inhibitor_m.amt*inhibitor_m.a2n+kI*inhibitor.conc- ...
                   delt_I.*inhibitor_m.amt*inhibitor_m.a2n; %s-1
inhibitor_m.src = inhibitor_m.src./inhibitor_m.a2n*dt;
        
% Rac
k_GAP_Rac=1;                        % s-1 ref
k_Rac=3*k_GAP_Rac/7.5;              %S-1 const            
Rac1=k_Rac.*Rac_in.amt*Rac.a2n;     %S-1      self production
A0=100; b=40;% 30 is good

Ang=1./(1+A0*exp(-b*pip3.conc.*FilGAP.conc));
Rac2=-Ang.*Rac.amt*Rac.a2n*0.3;     % 0.3
Rac3=-Rac.amt*Rac.a2n*k_GAP_Rac;    %s-1   self delay % k_GAP_Rac 负责 Rac 下调
%Rac4=0;
% Rac4=0.02*GEF_Rac.eff.*Rac_in.amt*Rac.a2n;
Rac4=0.02*GEF_Rac.eff.*Rac_in.amt*Rac.a2n; % oct 14 0.05
Rac5 = -0.5.*isadh.*Rac.amt*Rac.a2n; % for CIL


Rac.src=Rac1+Rac2+Rac3+Rac4+Rac5;
Rac.src = Rac.src./Rac.a2n*dt;

% Rac_in
Rac_in.src=-(Rac1+Rac2+Rac3+Rac4+Rac5);
Rac_in.src = Rac_in.src./Rac_in.a2n*dt;

% Rho
k_GAP_Rho=1;                                        % s-1
Rho_ini=1.25;                                       %uM
LL_Rho_c=3;                                         %uM
k_Rho=(Rho_ini)*k_GAP_Rho/LL_Rho_c;                 %S-1
Rho1=k_Rho*Rho_in.amt*Rho.a2n;                      %s-1
Rho2=-k_GAP_Rho*Rho.amt*Rho.a2n;
Rho3=Ang.*Rho_in.amt*Rho.a2n*0.3; % 0.3
%Rho4=0;
% Rho4=0.02*GEF_Rho.eff.*Rho_in.amt*Rho.a2n;
Rho4=0.02*GEF_Rho.eff.*Rho_in.amt*Rho.a2n; % oct 14; 0.05
Rho5 = 0.5.*isadh.*Rho_in.amt*Rho.a2n; % for CIL

Rho.src=Rho1+Rho2+Rho3+Rho4+Rho5;
Rho.src = Rho.src./Rho.a2n*dt;
        
% Rho_in
Rho_in.src=-(Rho1+Rho2+Rho3+Rho4+Rho5);
Rho_in.src = Rho_in.src./Rho_in.a2n*dt;

% -------------------------------------------------------------------------

function cel = csk2mol(cel, dt)
mols = cel.mols;
for i = 1:length(mols); eval([mols(i).name,'= mols(',num2str(i),');']); end

k_on_GAP = 50;
k_slow = 0; %4.0669
k_fast = 4; %0.1876

ind = find(cel.rcsk<0.9);
is = ismember(cel.tcsk, ind);
po_is=is(:);
deg = abs(cel.dAcsk(:));
cel.isopen = (cel.isopen | (deg>15)) & (deg>5);%10
ind = cel.isopen.*po_is;

k_gap = k_fast.* ind+k_slow.*ones(size(deg,1),1);
cel.res = cel.res-cel.res.* k_gap*dt; % GAP remain

c_gap = cel.tot_gap - sum(cel.res);
res_n = cel.res0-cel.res;  %avilable site
cel.res = cel.res + k_on_GAP*dt.*res_n*c_gap;
FilGAP.conc=c_gap; % GAP in cyto
