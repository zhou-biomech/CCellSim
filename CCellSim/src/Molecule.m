classdef Molecule < handle
   properties
       name                        % name
       domain                      % diffusion domain
       s                           % Length or area
       D                           % diffusion coefficient
       dt
       amt                         % amount = conc .* s (Length or area)
       num                         % molecular number
       conc                        % concentration
       cell = Cell.empty           % cell
       
       src                         % source term
       ligd                        % ligand molecule
       rcpt                        % receptor molecule
       
       kon 
       koff
       eff  = 0                    % molecule on membrane caculated by mc
       
       a2f                         % amout -> force
   end
   
   properties (Constant)
       a2n = 100;
   end
   
   methods
       function obj = Molecule(name, domain, conc, D, cell, a2f)
           if nargin<=5; a2f = 0; end
           obj.name = name;
           obj.domain = domain;
           obj.D = D;
           obj.cell = cell;
           obj.amt = conc.*obj.s;
           
           obj.src = zeros(size(cell.x));
           obj.a2f = a2f;
           obj.dt = 0.25*round(min(obj.cell.l).^2/obj.D, 1, 'significant');
       end
       
       % ------------------------------------------------------------------
       
       function obj = set.conc(obj, conc)
          obj.amt =  conc .* obj.s;
       end
       
       function conc = get.conc(obj)
           conc = obj.amt./obj.s .* ones(size(obj.cell.l));
       end
       
       % ------------------------------------------------------------------
       
       function num = get.num(obj)
           num = floor(obj.a2n*obj.amt);
       end
       % ------------------------------------------------------------------
       
       function s = get.s(obj)
           switch obj.domain               
               case {'plasm','cytoplasm',1}
                   s = obj.cell.A;
               case {'mebre','membrane',2}
                   s = obj.cell.lc;
           end
       end
       
       % ------------------------------------------------------------------
       
       function diffusion(objs,dt)
           for k = 1:length(objs)
               obj = objs(k);
               switch obj.domain
                   case {'plasm','cytoplasm',1}
                       obj.amt = obj.amt + sum(obj.src); %%??
                   case {'mebre','membrane',2}
                       if obj.D<=0
                           obj.amt = obj.amt + obj.src; continue; %%??
                       end
                       
                       C = obj.conc;
                       N = length(C);

                       i = [N, 1:N-1];
                       j = [1, 2:N  ];
                       
                       J = obj.D * (C(j) - C(i))./obj.cell.l;
                       obj.amt(j) = obj.amt(j) - J*dt;
                       obj.amt(i) = obj.amt(i) + J*dt;
                       
                       obj.amt = obj.amt + obj.src;
               end
           end
       end
       
       % ------------------------------------------------------------------
       
       function updligdsrc(objs, dt)
           for k = 1:length(objs)
               obj = objs(k);
               rcpt = obj.rcpt;
               if isempty(rcpt); continue; end
               obj = mcasdissn(obj, rcpt, rcpt.kon, rcpt.koff, dt);
               obj.src = obj.src./obj.a2n;
           end
       end
       
      
       % ------------------------------------------------------------------  
       
       function updnoligsrc(objs,dt,t)
           config('noligsrc', objs, dt,t);
       end

       % ------------------------------------------------------------------
   end
    
end

%% ------------------------------------------------------------------------


function ligd = mcasdissn(ligd, rcpt, kon, koff, dt)
% MC for Association and dissociation kinetics
% rcpt receptor   pip3    on mebre
% ligd ligand     pi3k    on plasm

n = rcpt.num - ligd.eff;
p = 1 - exp(-kon*ligd.conc*dt);

on = binornd(n,p);
ligd.eff = ligd.eff + on;

n = ligd.eff;
p = ones(size(n)) - exp(-koff*dt);
off = binornd(n,p);
ligd.eff = ligd.eff-off;

ligd.src = off - on;
end

% -------------------------------------------------------------------------

function r = binornd(n,p)
r = zeros(size(n));
for i = 1:length(n)
    r(i) = sum(rand(n(i),1)<p(i));
end
end
