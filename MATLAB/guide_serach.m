% Guided Motif Discovery for Multidimensional Matrix Profile
% Chin-Chia Michael Yeh
%
% [motif_idx, motif_dim] = guide_serach(...
%     data, sub_len, pro_mul, pro_idx, n_dim)
%
% Output:
%     motif_idx: the index for the founded motifs (matrix)
%     motif_dim: the dimensions spanned by the found motifs (cell)
% Input:
%     data: input time series (matrix)
%     sub_len: interested subsequence length (scalar)
%     pro_mul: multidimensional matrix profile (matrix)
%     pro_idx: matrix profile index (matrix)
%     n_dim: the dimensionality of the motif that you wish to find (scalar)
%
% C.-C. M. Yeh, N. Kavantzas, and E. Keogh, "Matrix Profile VI: Meaningful
% Multidimensional Motif Discovery," IEEE ICDM 2017.
% https://sites.google.com/view/mstamp/
% http://www.cs.ucr.edu/~eamonn/MatrixProfile.html
%

function [motif_idx, motif_dim] = guide_serach(...
    data, sub_len, pro_mul, pro_idx, n_dim)
pro_mul = pro_mul(:, n_dim);
pro_idx = pro_idx(:, n_dim);
[~, motif_idx] = min(pro_mul);
motif_idx = sort([motif_idx, pro_idx(motif_idx)]);

motif_1 = data(motif_idx(1):motif_idx(1)+sub_len - 1, :);
motif_2 = data(motif_idx(2):motif_idx(2)+sub_len - 1, :);

[~, motif_dim] = sort(sum(abs(motif_1 - motif_2), 1), 'ascend');
motif_dim = sort(motif_dim(1:n_dim));
motif_dim = {motif_dim; motif_dim;};