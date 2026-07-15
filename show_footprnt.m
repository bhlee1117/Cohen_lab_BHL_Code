function mergedImg = show_footprnt(c_ftprnt, mov_mc, colr)
    [H, W, N] = size(c_ftprnt);

    if nargin < 3 || isempty(colr)
        colr = max(turbo(N), 0);
        colr(colr>1) = 1;
    end

    % Normalize footprints
    c_ftprnt = toimg(rescale2(tovec(c_ftprnt),1), H, W);
    c_ftprnt = mat2gray(c_ftprnt);
    mov_mc   = mat2gray(mov_mc);

    % ---- Fast colored composite: matrix multiply instead of broadcast+sum ----
    % (H*W x N) * (N x 3) -> (H*W x 3) -> H x W x 3
    Cvec      = reshape(c_ftprnt, H*W, N);
    mergedImg = reshape(Cvec * colr(1:N,:), H, W, 3);

    ax1 = nexttile([1 1]);
    imshow2(mergedImg, []);

    ax2 = nexttile([1 1]);
    imshow2(imfuse(mean(mov_mc,3), sum(c_ftprnt,3)), []);
    linkaxes([ax1 ax2],'xy')
end
