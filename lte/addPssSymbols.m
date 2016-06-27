function [ grid ] = addPssSymbols(grid, lte)
% Insert PSS symbols in a LTE grid

    cellid = lte.cellid;
    pss_indices = lte.nUsedSubcarriers/2-31:lte.nUsedSubcarriers/2+30;

    %% Generate Zadoff-Chu Sequence
    n = [-31:-1 1:31];
    zc = exp((-j*cellid.*n.*(n+1))/63);
    
    %% Map to subcarriers
    grid(pss_indices, 7)  = zc;
    grid(pss_indices, 77) = zc;
end