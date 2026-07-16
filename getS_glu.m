function S = getS_glu(G)
% GETS_GLU  Return the H x W x Nsyn glutamate footprint stack from a Glu_Result
% struct, whether S_glu is stored dense (legacy) or sparse [H*W x Nsyn] (compacted).
% For the sparse form G must also carry G.sz_glu = [H W].
S = G.S_glu;
if issparse(S)
    S = reshape(full(S), G.sz_glu(1), G.sz_glu(2), []);
end
end
