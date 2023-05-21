% Hamamatsu camera noise removal
% Eli N. Weinstein, Cohen Lab, 5/5/2015

% Input:
% thsto -> Trace to be corrected
% Output:
% varargout{1} (thstc) -> Trace corrected based on average saturation of noise
% frames
% varargout{2} (noiset) -> Location of noise frames (in case you want them to be
% corrected by a set value, not the value measured here). 
% varargout{3} (satcoeff) -> Saturation of noise frames (exposure of noise
% frames/exposure of non-noise frames)
function varargout = HamNoiseRemoval500(thsto)
% Settings
% If you're debugging or want to understand what's going on, set to 1 and figures will pop up.
dsp = 0; 
% Window to average over. # of frames in between when noise frames are 3
% apart rather than 4.
wind = 4*30/2;
% Error adjustment steps. Determines the balance of importance between
% local maxima (which can be thrown off by spikes) and global averages (which
% can be thrown off by swamping the local) in refining the location of
% 3-spikes. Higher indicates greater local importance. 0 indicates no
% correction.
adj3iter = 2;
% Secondary smoothing scale (after averaging over phase in window). With
% current error correction methods, I've found it works best to turn this
% off (smwind = 1).
smwind = 1;
%% Local smoothing
bound1 = @(x) max(1, x);
bound2 = @(x) min(length(thsto), x);
% Smooth in a small local window (that can include max. 1 noise frame) by
% scaling min in that frame to 0 and max in that frame to 1.
thst = zeros(size(thsto));
for j = 1:length(thsto)
    thsseg = thsto(bound1(j-2):bound2(j+2));
    thst(j) = (thsto(j)-min(thsseg))/(max(thsseg)-min(thsseg));
end
%% Wide smoothing according to offset
% Identify each frame with its offset
indx4 = zeros(1, length(thst));
indx4(1:4:length(thst)) = 1*ones(1, length(1:4:length(thst)));
indx4(2:4:length(thst)) = 2*ones(1, length(2:4:length(thst)));
indx4(3:4:length(thst)) = 3*ones(1, length(3:4:length(thst)));
indx4(4:4:length(thst)) = 4*ones(1, length(4:4:length(thst)));

% Look in a window of size 30*4 around each frame and compute an average
% intensity for each offset within that window
wsc = zeros(4, length(thst));
for j = 1:length(thst)
    for jj = 1:4
        if bound1(j-wind+1) ~= j-wind+1
            r1 = bound1(j-wind+1);
            r2 = bound2(j + max([j - r1, 4]));
        elseif bound2(j+wind) ~= j+wind
            r2 = bound2(j+wind);
            r1 = bound1(j - max([r2 - j], 4));
        else
            r1 = j - wind+1;
            r2 = j + wind;
        end
        jjind = find(indx4(r1:r2) == jj);
        locw = thst(r1:r2);
        wsc(jj, j) = mean(locw(jjind));
    end
end

% Smooth within each offset; should function as a form of error correction.
wscs = zeros(size(wsc));
for jj = 1:4
    wscs(jj, :) = smooth(wsc(jj, :), smwind);
end

if dsp
    figure;
    subplot(3, 1, 1)
    plot(wsc')
    hold on
    plot(wscs', ':')
    hold off
    subplot(3, 1, 2)
    plot(thsto)
    subplot(3, 1, 3)
    plot(thst)
end
%% Assign noise frames
noiset = zeros(size(thst));
for j = 1:length(thst)
    [~, mi] = max(wscs(:, j));
    if mi == indx4(j)
        noiset(j) = 1;
    end
end
%% Error correction
[~, pkloc] = findpeaks(noiset);
dfpkloc = diff(pkloc);
for cc = 1:length(dfpkloc)
    % Sometimes a 3-spike is missed. We can only decide between the two 
    % middle frames based on which has a higher intensity
    if dfpkloc(cc) == 7 
        if thst(pkloc(cc)+3) > thst(pkloc(cc)+4) %thsto(pkloc(cc)+3) > thsto(pkloc(cc)+4)
            noiset(pkloc(cc)+3) = 1;
        else
            noiset(pkloc(cc)+4) = 1;
        end
    end
end
for ccc = 1:adj3iter
    [~, pkloc] = findpeaks(noiset);
    dfpkloc = diff(pkloc);
    for cc = 1:length(dfpkloc)
        % At every 3 spike, we may have to adjust either side
        if dfpkloc(cc) == 3
            % Choose side to be adjusted based on local peak size
            gp1 =  thst(pkloc(cc)-1) - thst(pkloc(cc)); 
            gp2 =  thst(pkloc(cc+1)+1) - thst(pkloc(cc+1));
            if gp2 > gp1
                if thst(pkloc(cc+1)+1) > thst(pkloc(cc+1))
                    noiset(pkloc(cc+1)) = 0;
                    noiset(pkloc(cc+1)+1) = 1;
                end
            else
                if thst(pkloc(cc)-1) > thst(pkloc(cc))
                    noiset(pkloc(cc)) = 0;
                    noiset(pkloc(cc)-1) = 1;
                end
            end
            
        end
    end
end
%% Estimate saturation of noise frames and correct original trace
% Re-smooth using our knowledge of which peaks are noise
[~, pkloc] = findpeaks(noiset);
noisev = noiset(pkloc);
pkloc = [0, pkloc, length(noiset)+1];
for k = 2:(length(pkloc)-1)
    lbase = mean(thsto([(pkloc(k)-1):-1:(pkloc(k-1)+1), (pkloc(k)+1):1:(pkloc(k+1)-1)]));
    noisev(k) = thsto(pkloc(k))/lbase;
end
% Average noise frame saturation
satcoeff = mean(noisev);
% Correct original trace
thstc = thsto.*(~noiset + (noiset./satcoeff));

if dsp
    figure
    subplot(2, 1, 1)
    plot(thsto)
    subplot(2, 1, 2)
    plot(thstc)
end
%% Output
varargout{1} = thstc;
varargout{2} = noiset;
varargout{3} = satcoeff;