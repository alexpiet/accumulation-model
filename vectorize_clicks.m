
function [buptimes,nantimes,streamIdx] = vectorize_clicks(leftbups,rightbups)
% function [buptimes,nantimes,streamIdx] = vectorize_clicks(leftbups,rightbups)
% code written by Adrian Bondy, originally as processVectorizedData
    nTrials=length(leftbups);
    nleft = cellfun('length',leftbups);
    nright = cellfun('length',rightbups);
    maxBups = max([nleft(:)+nright(:)]);
    buptimes=NaN(maxBups,nTrials);
    streamIdx = buptimes;
    for i=1:nTrials
        % concatenate left and right bups and keep track of which was which
        buptimes(1:nleft(i),i) = leftbups{i};
        buptimes(nleft(i)+1:(nright(i)+nleft(i)),i) = rightbups{i};
        streamIdx(1:nleft(i),i) = -1;
        streamIdx(nleft(i)+1:(nright(i)+nleft(i)),i) = 1;
    end
    maxNeeded = sum(all(isnan(buptimes),2)==0);
    buptimes = buptimes(1:maxNeeded,:);
    nantimes = isnan(buptimes);
    streamIdx = streamIdx(1:maxNeeded,:);
    % sort buptimes and stream idx
    [buptimes,sortIdx] = sort(buptimes);
    tmp = meshgrid(1:nTrials,1:maxNeeded);
    streamIdx = streamIdx(sub2ind(size(streamIdx),sortIdx(:),tmp(:)));
    streamIdx = reshape(streamIdx,size(buptimes));
end
