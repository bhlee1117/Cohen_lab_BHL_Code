function C = getGluCoord(G)
% GETGLUCOORD  Synapse centroids [x y] (Nsyn x 2) from a Glu_Result struct.
% Uses the precomputed G.GluCoord if present (compacted files), otherwise
% computes it from the footprints (legacy dense files).
if isfield(G,'GluCoord') && ~isempty(G.GluCoord)
    C = G.GluCoord;
else
    C = get_coord(getS_glu(G));
end
end
