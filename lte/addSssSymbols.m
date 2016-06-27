function [ grid ] = addSssSymbols(grid, lte)
% Insert SSS symbols in a LTE grid
% Note: this function is using Matlab's lteSSS function because
% SSS symbols generation is based in a series of gold sequence
% generation and rotations, which are note of interest to the 
% course.

    sss_from_matlab = lteSSS(struct('NCellID',lte.cellid));
    sss_indices = lte.nUsedSubcarriers/2-31:lte.nUsedSubcarriers/2+30;
    
    grid(sss_indices, 6)  = sss_from_matlab;
    grid(sss_indices, 76) = sss_from_matlab;
end