function adapted = adapt_vectorized_clicks(buptimes, nantimes, phi, tau_phi, varargin)
    % buptimes is a matrix of click times of size m X n where m is the number
    % of clicks and n is the number of trials. The array is NaN-padded for
    % all trials except the ones with the longest number of clicks.
    % adaptation occurrs between clicks independent of which side they
    % occurred on.
    % rewritten by Adrian, 2017 to be vectorized
    adapted=buptimes;
    adapted(~nantimes)=1;
    phi=max(0,phi);
    tau_phi = max(0,tau_phi);
    ici = diff(buptimes);
    for i = 2:size(buptimes,1)
        last = tau_phi * log(1 - adapted(i-1,:)*phi);
        %check = 1 - exp((-ici(i-1,:) + last)/tau_phi);
        adapted(i,:) = 1 + exp(-ici(i-1,:)/tau_phi).*(adapted(i-1,:)*phi-1);
        
        adapted(i,ici(i-1,:)<=0)=0;
        adapted(i-1,ici(i-1,:)<=0)=0;    
    end
    adapted = real(adapted);
end