function p = get_grid(trial, forward_in, params,p);
% Stub function!

if ~isfield(p, 'avals')
    disp('resetting avals')
    keyboard
    p.da_grid =1;
    p.da = 1;
    p.avals =-100:p.da:100; 
end
if ~isfield(forward_in, 'pdf')
    forward.ma = forward_in.ma(end);
    forward.va = forward_in.va(end);
    forward.T = forward_in.T(end);
    forward = compute_pdf(forward, p.avals, p);
else
    forward = forward_in;
end

if trial.pokedR == 0
    xint    = cumsum(forward.pdf(end,:).*p.da);
    dex     = find(xint > p.error_tolerance, 1)-1;
    if dex <= 0; dex = 1; end;
    p.gridb = p.avals(dex);
    p.a_grid= p.gridb:p.da_grid:0 + params(7);
%    p.a_grid= [p.gridb:2:-10 -9:0.5:0] + params(7);

else
    xint    = 1-cumsum(forward.pdf(end,:).*p.da);
    dex     = find(xint < p.error_tolerance, 1);
    if isempty(dex); dex = length(p.avals); end;
    p.gridb = p.avals(dex);
    p.a_grid= 0:p.da_grid:p.gridb + params(7);

end

% simple version
%p.a_grid    = -p.b:p.da_grid:p.b + params(7);

% set scale factors for non-uniform grid
p.a_grid_s  = [p.da_grid (diff(p.a_grid(1:end-1)) + diff(p.a_grid(2:end)))/2 p.da_grid];

if p.da_grid > p.da
    error('delta mode spacing larger than evaluation spacing, will create funky results. Proceed with caution!');
end



% I could do a couple things here:
% I could expand avals - by recursively calling with larger and larger p.avals
% I could trim a_grid to stay inside avals
% will need to handle the pokedR = 0 case as well
if max(abs(p.a_grid)) > max(abs(p.avals))
%    disp('You need to expand the size of avals, because the probability mass is hitting the edge')
    p.avals = (p.avals(1)-1):p.da:(p.avals(end)+1);
    forward = rmfield(forward, 'pdf');
    p = get_grid(trial, forward, params,p);   
end
