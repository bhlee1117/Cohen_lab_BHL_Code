function mom = imgmoments(imgstack)
    pk_an_r = (size(imgstack,1)-1)/2; % peak analysis radius
    ar = -pk_an_r:pk_an_r;
    assert(isequal(numel(ar),size(imgstack,1)),'check imgs')
    assert(isequal(size(imgstack,1),size(imgstack,2)),'imgs should be square')
    [ay, ax] = ndgrid(ar,ar); % grid
%         for ib = 1:nc
%             rspots(:,:,ib) = rim(crow(ib)+ar,ccol(ib)+ar,ir); % crop in corrected rows and cols
%             nspots(:,:,ib) = rspots(:,:,ib) - min(min(rspots(:,:,ib))); % remove offset
%             nspots(:,:,ib) = nspots(:,:,ib)./max(max(nspots(:,:,ib))); % remove offset
%         end
        raw_ms = [ % raw area, raw 1st moment
            sum(bsxfun(@times,tovec(imgstack),ax(:).^0.*ay(:).^0))
            sum(bsxfun(@times,tovec(imgstack),ax(:).^1.*ay(:).^0))
            sum(bsxfun(@times,tovec(imgstack),ax(:).^0.*ay(:).^1))
            ];
%         newcoids = [ % centroid
%             raw_ms(2,:)./raw_ms(1,:)
%             raw_ms(3,:)./raw_ms(1,:)
%             ];
        [~, idx] = max(tovec(imgstack));
        [i1, i2] = ind2sub([size(imgstack,1) size(imgstack,2)],idx);
        newcoids = [
            ar(i1)
            ar(i2)
            ];
        stdxy = sqrt([ % central std dev x, y
            sum(bsxfun(@times,tovec(imgstack),bsxfun(@minus,ax(:),newcoids(1,:)).^2))
            sum(bsxfun(@times,tovec(imgstack),bsxfun(@minus,ay(:),newcoids(2,:)).^2))
            ]./([1;1]*raw_ms(1,:)));
        stdxy = real(stdxy);
        skxy = [ % central skewness, x, y, normalized to std dev
            sum(bsxfun(@times,tovec(imgstack),bsxfun(@minus,ax(:),newcoids(1,:)).^3))
            sum(bsxfun(@times,tovec(imgstack),bsxfun(@minus,ay(:),newcoids(2,:)).^3))
            ]./stdxy.^3./([1;1]*raw_ms(1,:));
        kurxy = [ % central kurtosis, x, y, normalized to std dev
            sum(bsxfun(@times,tovec(imgstack),bsxfun(@minus,ax(:),newcoids(1,:)).^4))
            sum(bsxfun(@times,tovec(imgstack),bsxfun(@minus,ay(:),newcoids(2,:)).^4))
            ]./stdxy.^4./([1;1]*raw_ms(1,:));
mom = [
    raw_ms
    newcoids
    stdxy
    skxy
    kurxy
    ];
