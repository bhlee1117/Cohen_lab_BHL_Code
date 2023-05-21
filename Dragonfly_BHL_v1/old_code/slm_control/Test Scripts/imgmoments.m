function mom = imgmoments(imgstack)
% imgmoments Calculate image moments
%
% 2016-2018 Vicente Parot
% Cohen Lab - Harvard University
%
    pk_an_r = (size(imgstack,1)-1)/2; % peak analysis radius
    ar = -pk_an_r:pk_an_r;
    assert(isequal(numel(ar),size(imgstack,1)),'imgs should have odd dimensions')
    assert(isequal(size(imgstack,1),size(imgstack,2)),'imgs should be square')
    [ay, ax] = ndgrid(ar,ar); % grid
%         for ib = 1:nc
%             rspots(:,:,ib) = rim(crow(ib)+ar,ccol(ib)+ar,ir); % crop in corrected rows and cols
%             nspots(:,:,ib) = rspots(:,:,ib) - min(min(rspots(:,:,ib))); % remove offset
%             nspots(:,:,ib) = nspots(:,:,ib)./max(max(nspots(:,:,ib))); % remove offset
%         end
        raw_ms = [ % raw area, raw 1st moments
            sum(bsxfun(@times,tovec(imgstack),ax(:).^0.*ay(:).^0)) % output 1
            sum(bsxfun(@times,tovec(imgstack),ax(:).^1.*ay(:).^0)) % output 2
            sum(bsxfun(@times,tovec(imgstack),ax(:).^0.*ay(:).^1)) % output 3
            ];
        newcoids = [ % centroid
            raw_ms(2,:)./raw_ms(1,:) % output 4
            raw_ms(3,:)./raw_ms(1,:) % output 5
            ];
        stdxy = sqrt([ % central std dev x, y
            sum(bsxfun(@times,tovec(imgstack),bsxfun(@minus,ax(:),newcoids(1,:)).^2)) % output 6
            sum(bsxfun(@times,tovec(imgstack),bsxfun(@minus,ay(:),newcoids(2,:)).^2)) % output 7
            ]./([1;1]*raw_ms(1,:)));
        stdxy = abs(stdxy);
        skxy = [ % central skewness, x, y, normalized to std dev
            sum(bsxfun(@times,tovec(imgstack),bsxfun(@minus,ax(:),newcoids(1,:)).^3)) % output 8
            sum(bsxfun(@times,tovec(imgstack),bsxfun(@minus,ay(:),newcoids(2,:)).^3)) % output 9
            ]./stdxy.^3./([1;1]*raw_ms(1,:));
        kurxy = [ % central kurtosis, x, y, normalized to std dev
            sum(bsxfun(@times,tovec(imgstack),bsxfun(@minus,ax(:),newcoids(1,:)).^4)) % output 10
            sum(bsxfun(@times,tovec(imgstack),bsxfun(@minus,ay(:),newcoids(2,:)).^4)) % output 11
            ]./stdxy.^4./([1;1]*raw_ms(1,:));
mom = [
    raw_ms
    newcoids
    stdxy
    skxy
    kurxy
    ];
