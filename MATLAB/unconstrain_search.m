% MDL Based Motif Discovery for Multidimensional Matrix Profile
% Chin-Chia Michael Yeh
%
% [motif_idx, motif_dim] = unconstrain_search(...
%     data, sub_len, pro_mul, pro_idx, n_bit, k)
%
% Output:
%     motif_idx: the index for the founded motifs (matrix)
%     motif_dim: the dimensions spanned by the found motifs (cell)
% Input:
%     data: input time series (matrix)
%     sub_len: interested subsequence length (scalar)
%     pro_mul: multidimensional matrix profile (matrix)
%     pro_idx: matrix profile index (matrix)
%     n_bit: number of bit for discretization (scalar)
%     k: number of motif wish to retrieve, set to inf for retrieving
%        all possible k-motifs (scalar)
%
% C.-C. M. Yeh, N. Kavantzas, and E. Keogh, "Matrix Profile VI: Meaningful
% Multidimensional Motif Discovery," IEEE ICDM 2017.
% https://sites.google.com/view/mstamp/
% http://www.cs.ucr.edu/~eamonn/MatrixProfile.html
%

function [motif_idx, motif_dim] = unconstrain_search(...
    data, sub_len, pro_mul, pro_idx, n_bit, k)
exc_zone = round(0.5 * sub_len);
tot_dim = size(data, 2);
if isinf(k)
    k = size(pro_mul, 1);
end
motif_idx = zeros(k, 1);
motif_dim = cell(k, 1);
base_bit = n_bit * tot_dim * sub_len * 2;
for i = 1:k
    fprintf('finding motif %d ... \n', i);
    [val, idx_1] = min(pro_mul, [], 1);
    if any(isinf(val))
        motif_idx = motif_idx(1:k-1);
        motif_dim = motif_dim(1:k-1);
        break;
    end

    bit_sz = zeros(tot_dim, 1);
    idx_2 = zeros(tot_dim, 1);
    dim = cell(tot_dim, 1);
    for j = 1:tot_dim
        idx_2(j) = pro_idx(idx_1(j), j);
        motif_1 = data(idx_1(j):idx_1(j) + sub_len - 1, :);
        motif_2 = data(idx_2(j):idx_2(j) + sub_len - 1, :);
        [bit_sz(j), dim{j}] = get_bit_save(motif_1, motif_2, j, n_bit);
    end
    [best_bit, min_idx] = min(bit_sz);
    if best_bit > base_bit
        motif_idx = motif_idx(1:k-1);
        motif_dim = motif_dim(1:k-1);
        break;
    end
    motif_idx(i, 1) = idx_1(min_idx);
    motif_dim{i} = dim{min_idx};

    st_idx = max(1, motif_idx(i, 1) - exc_zone);
    ed_idx = min(size(pro_mul, 1), motif_idx(i, 1) + exc_zone);
    pro_mul(st_idx:ed_idx, :) = inf;
end
motif_dim = motif_dim(motif_idx ~= 0);
motif_idx = motif_idx(motif_idx ~= 0);


function [bit_sz, dim_id] = get_bit_save(motif_1, motif_2, n_dim, n_bit)
tot_dim = size(motif_1, 2);
sub_len = size(motif_1, 1);
split_pt = get_desc_split_pt(n_bit);
disc_1 = discretization(motif_1, split_pt);
disc_2 = discretization(motif_2, split_pt);

[~, dim_id] = sort(sum(abs(disc_1 - disc_2), 1), 'ascend');
dim_id = dim_id(1:n_dim);
motif_diff = disc_1(:, dim_id) - disc_2(:, dim_id);
n_val = length(unique(motif_diff));

bit_sz = n_bit * (tot_dim * sub_len * 2 - n_dim * sub_len);
bit_sz = bit_sz + n_dim * sub_len * log2(n_val) + n_val * n_bit;


function disc = discretization(motif, split_pt)
for i = 1:size(motif, 2)
    motif(:, i) = (motif(:, i) - mean(motif(:, i))) / ...
        std(motif(:, i), 1);
end
disc = zeros(size(motif));
for i = 1:length(split_pt)
    disc(motif < split_pt(i) & disc == 0) = i;
end
disc(disc == 0) = length(split_pt) + 1;


function split_pt = get_desc_split_pt(n_bit)
split_pt = norminv((1:(2^n_bit)-1)/(2^n_bit), 0, 1);