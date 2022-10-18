function [D] = compute_pdf(D,avals,p,varargin)
% Computes the PDF at each timepoint. 
% D is a distribution structure with the following fields
%   D.ma, a vector of mean values for the distribution at each time point
%   D.va, a vector of variance values for the distribution at each time point
%   D.T,  a vector of time points to evalute the distribution at
%
% avals is a vector of accumulation values that determines the grid the distribution is evaluated on
%
% Returns D, with the following fields added
%   D.pdf, the pdf at each timepoint
%   D.avals, the same vector avals passed in, is now stored in the structure

if length(D.ma) ~= length(D.va)
    error('mean and variance vectors must be the same length')
end

if length(avals) < 1
    error('Vector of A values must be greater than 1')
end

if length(unique(diff(avals))) > 1
    if sum(abs(diff(unique(diff(avals))))) > 1e-8
    error('p.avals is not a uniform grid! I cannot verify solution for a non uniform evaluation grid. Proceed with EXTREME caution!');
    end
end

if p.da > 1
    error('p.da > 1 may lead to unstable behavior at some point. Proceed with caution!');
end

% work through logic on inputs
single_gauss = true;
modedex = 1;
if length(varargin) > 0
    option = varargin{1};
    if ischar(option)
        if strcmp(option, 'mode') & size(D.ma,1) > 1
            modedex = varargin{2};
        elseif strcmp(option, 'mixture') & size(D.ma,1) > 1
            single_gauss = false;
        end
    end
end
if isfield(D, 'pdf')
    D = rmfield(D,'pdf');
end

if single_gauss
    % Compute PDF for each time point
    for i=1:length(D.ma)
        D.pdf(i,:) = normpdf(avals, D.ma(modedex,i),sqrt(D.va(modedex,i)));
        % normpdf doesn't normalize things correctly if all the mass is in one bin
        % Here I check for that case, and normalize.
        % Note that the distribution wont integrate to 1 if avals is too short, hence the 1-sided equality check. 
        da = min(diff(avals));
        if sum(da.*D.pdf(i,:)) - 1 > 0.01
            if length(unique(D.pdf(i,:))) > 2 
            % this was happening sometimes in the low-variance case....
            % If you remove this error, things should run to completetion,...
            % but the variance of the posterior distribution is sometimes ...
            % wrong. Proceed carefully!
                warning(['hitting weird case with integral summation '...
                    'and not just a 1-bin thing. '...
                    'Are you in the low-variance case?'])
            end
            D.pdf(i,:) = D.pdf(i,:)./sum(da.*D.pdf(i,:));
        end        
        
        if sum(da.*D.pdf(i,:)) < 0.1 & abs(D.ma(modedex,i)) < max(abs(avals))
%            keyboard
            x = diff(normcdf(avals, D.ma(modedex,i), sqrt(D.va(modedex,i))));
            z = [0 (x(1:end-1) + x(2:end))./2 0]; 
            D.pdf(i,:) = z./sum(da.*z);
        end

    end
else
    % need to compute mixture model of n modes 
    % We now assume that D.ma and D.va are (n x timesteps)    
    for i=1:size(D.ma,2)
        sigma(1,1,:) = D.va(:,i);
        if isfield(D, 's')
            M{i} = gmdistribution(D.ma(:,i), sigma,D.s(:,i));
        else
            M{i} = gmdistribution(D.ma(:,i), sigma);
        end
        PDF{i} = pdf(M{i}, avals');
        % Again the weird case with delta-functions is coming into play
        da = min(diff(avals));
        if sum(PDF{i}.*da) -1 > 0.01;
        % I should cache this in p (fine_grid, that is);
        % This is also a poor-mans solution. I really should just do adaptive checks.
            if max(abs(p.avals)) < max(abs(p.a_grid));
                disp('evaluation bins are too small for grid, type dbcont to continue with a hack solution');
                keyboard
                fine_grid = -100:da:100;
                temp = pdf(M{i}, fine_grid');
                temp_sum = sum(temp.*da);
                PDF{i} = PDF{i}./temp_sum;
            else
                PDF{i} = PDF{i}./sum(PDF{i}.*da);
            end
        end
        if sum(D.va(:,i)) < 0.1 & i > size(D.ma,2)/2
            y= diff(cdf(M{i},p.avals'));
            z = [0; (y(1:end-1) + y(2:end))./2; 0];
            PDF{i} = z./sum(da.*z);
        end
        pdfc(i,:) = PDF{i};
    end
%    D.M = M; % I Dont know why I would ever need this, so removing for space
%    D.PDF = PDF; % This is redundant with pdf, so I won't return it
    D.pdf =pdfc;
end

% Store a value grid
D.avals = avals;

