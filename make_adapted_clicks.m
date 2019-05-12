function [L R Lnclicks Rnclicks extras] = make_adapted_clicks(leftbups, rightbups, phi, tau_phi, psi, tau_psi, varargin)
% function [L R] = make_adapted_clicks(leftbups, rightbups, phi, tau_phi, psi, tau_psi, varargin)
%
% compute the effectives sizes of each click in leftbups and rightbups
% considering the sensory adaptation parameters
% 
% returns L and R, the same sizes as leftbups and rightbups.
% if there's no appreciable adaptation, returns all ones

pairs = {
    'cross_side_suppression'  0.002; ... % clicks on different sides closer than this duration (in sec) don't count for nclicks
}; parseargs(varargin, pairs);


Lsame = ones(size(leftbups)); 
Rsame = ones(size(rightbups)); 

% magnitude of stereo clicks set to zero
if ~isempty(leftbups) && ~isempty(rightbups) && abs(leftbups(1)-rightbups(1)) < eps,
    Lsame(1) = 0;
    Rsame(1) = 0;
end;

% if there's appreciable same-side adaptation
if abs(phi - 1) > eps, 
    % inter-click-intervals
    ici_L = diff(leftbups);
    ici_R = diff(rightbups);

    for i = 2:numel(leftbups),
        last_L = tau_phi*log(1-Lsame(i-1)*phi);
        Lsame(i) = 1 - exp((-ici_L(i-1) + last_L)/tau_phi);
    end;
    
    for i = 2:numel(rightbups),
        last_R = tau_phi*log(1-Rsame(i-1)*phi);
        Rsame(i) = 1 - exp((-ici_R(i-1) + last_R)/tau_phi);
    end;
    
    Lsame = real(Lsame);
    Rsame = real(Rsame);
end;

Lother = ones(size(leftbups));
Rother = ones(size(rightbups));

% if there's appreciable across-side adaptation
if abs(psi - 1) > eps, 
    lefts  = [leftbups(:)  -ones(numel(leftbups),1)];
    rights = [rightbups(:) +ones(numel(rightbups),1)];
    allbups = sortrows([lefts; rights])'; % one bup in each col, second row has side bup was on
    
    adapted = ones(1,size(allbups,2)); 
    nclicks = ones(size(adapted)); 
    
    % now let's go through and figure all the across-side adaptive effects
    for c = 1:size(allbups,2)-1,
        if allbups(2,c)~=allbups(2,c+1), % if this bup and the next are on opposite sides
            dt = allbups(1,c+1) - allbups(1,c);
            if dt <= cross_side_suppression,
                adapted(c) = 0;
                adapted(c+1) = 0;
                nclicks(c) = 0.5;
                nclicks(c+1) = 0.5;
            else
                % strength of the cross-side adaptation is weighed by the
                % magnitude of the preceeding click
                adapted(c+1) = 1 - exp(-dt/tau_psi + log(1 - (adapted(c)*(psi-1)+1)));
            end;
        end;
    end;
    
    
    Lother = real(adapted(allbups(2,:)==-1));
    Rother = real(adapted(allbups(2,:)==+1));
    Lnclicks = nclicks(allbups(2,:)==-1);
    Rnclicks = nclicks(allbups(2,:)==+1);
else
    Lnclicks = ones(size(leftbups));
    Rnclicks = ones(size(rightbups));
end;

% now take the product of the two effects:
L = Lsame .* Lother;
R = Rsame .* Rother;

if nargout > 4,
    extras.Lsame = Lsame;
    extras.Rsame = Rsame;
    extras.Lother = Lother;
    extras.Rother = Rother;
end;