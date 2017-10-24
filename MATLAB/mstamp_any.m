% STAMP Based mSTAMP Implemented as an Anytime Algorithm
% Chin-Chia Michael Yeh
%
% [pro_mul, pro_idx] = mstamp_any(data, sub_len, pct_stop)
%
% Output:
%     pro_mul: multidimensional matrix profile (matrix)
%     pro_idx: matrix profile index (matrix)
% Input:
%     data: input time series (matrix)
%     sub_len: interested subsequence length (scalar)
%     pct_stop: stop percentage, a number from 0 to 1 (scalar)
%
% C.-C. M. Yeh, N. Kavantzas, and E. Keogh, "Matrix Profile VI: Meaningful
% Multidimensional Motif Discovery," IEEE ICDM 2017.
% https://sites.google.com/view/mstamp/
% http://www.cs.ucr.edu/~eamonn/MatrixProfile.html
%

function [pro_mul, pro_idx] = mstamp_any(data, sub_len, pct_stop)
%% get various length
exc_zone = round(sub_len / 2);
data_len = size(data, 1);
n_dim = size(data, 2);
pro_len = data_len - sub_len + 1;
itr_stop = round(pro_len * pct_stop);
if itr_stop < 1
    itr_stop = 1;
end
if itr_stop > pro_len
    itr_stop = pro_len;
end

%% check input
if sub_len > data_len / 2
    error(['Error: Time series is too short relative to desired ' ...
        'subsequence length']);
end
if sub_len < 4
    error('Error: Subsequence length must be at least 4');
end
if pct_stop < 0
    error('Error: Stop percentage must be greater than 0');
end
if pct_stop > 1
    error('Error: Stop percentage must be less than or equal to 1');
end

%% check skip position
skip_loc = false(pro_len, 1);
for i = 1:pro_len
    if any(isnan(reshape(data(i:i+sub_len-1, :), 1, []))) ...
            || any(isinf(reshape(data(i:i+sub_len-1, :), 1, [])))
        skip_loc(i) = true;
    end
end
data(isnan(data)) = 0;
data(isinf(data)) = 0;

%% initialization
data_freq = zeros((sub_len + data_len), n_dim);
data_mu = zeros(pro_len, n_dim);
data_sig = zeros(pro_len, n_dim);
first_prod = zeros(pro_len, n_dim);
for i = 1:n_dim
    [data_freq(:, i), data_mu(:, i), data_sig(:, i)] = ...
        mass_pre(data(:, i), data_len, sub_len);
    [~, first_prod(:, i)] = mass(...
        data_freq(:, i), data(1:sub_len, i), data_len, ...
        sub_len, data_mu(:, i), data_sig(:, i), ...
        data_mu(1, i), data_sig(1, i));
end

%% compute the matrix profile
dist_pro = zeros(pro_len, n_dim);
last_prod = zeros(pro_len, n_dim);
pro_mul = inf(pro_len, n_dim);
pro_idx = zeros(pro_len, n_dim);
idxs = randperm(pro_len);
idxs = idxs(1:itr_stop);
for j = 1:length(idxs)
    idx = idxs(j);
    fprintf('%d %d\n', j, pro_len);
    query = data(idx:idx+sub_len-1, :);
    for k = 1:n_dim
        [dist_pro(:, k), last_prod(:, k)] = ...
            mass(data_freq(:, k), query(:, k), ...
            data_len, sub_len, data_mu(:, k), ...
            data_sig(:, k), data_mu(idx, k), ...
            data_sig(idx, k));
    end
    dist_pro = real(dist_pro);
    dist_pro = max(dist_pro, 0);

    % apply exclusion zone
    exc_zone_st = max(1, idx - exc_zone);
    exc_zone_ed = min(pro_len, idx + exc_zone);
    dist_pro(exc_zone_st:exc_zone_ed, :) = inf;
    dist_pro(data_sig < eps) = inf;
    if skip_loc(idx)
        dist_pro = inf(size(dist_pro));
    end
    dist_pro(skip_loc, :) = inf;

    % figure out and store the nearest neighbor
    dist_pro_sort = sort(dist_pro, 2);
    dist_pro_cum = zeros(pro_len, 1);
    dist_pro_merg = zeros(pro_len, 1);
    for k = 1:n_dim
        dist_pro_cum = dist_pro_cum + dist_pro_sort(:, k);
        dist_pro_merg(:) = dist_pro_cum / k;
        update_idx = dist_pro_merg < pro_mul(:, k);
        pro_mul(update_idx, k) = dist_pro_merg(update_idx);
        pro_idx(update_idx, k) = idx;
    end
end
pro_mul = sqrt(pro_mul);


%% The following two functions are modified from the code provided in the following URL
%  http://www.cs.unm.edu/~mueen/FastestSimilaritySearch.html
function [data_freq, data_mu, data_sig] = mass_pre(data, data_len, sub_len)
data(data_len+1:(sub_len+data_len)) = 0;
data_freq = fft(data);
data_cum = cumsum(data);
data2_cum =  cumsum(data.^2);
data2_sum = data2_cum(sub_len:data_len) - ...
    [0; data2_cum(1:data_len-sub_len)];
data_sum = data_cum(sub_len:data_len) - ...
    [0; data_cum(1:data_len-sub_len)];
data_mu = data_sum./sub_len;
data_sig2 = (data2_sum./sub_len)-(data_mu.^2);
data_sig2 = real(data_sig2);
data_sig2 = max(data_sig2, 0);
data_sig = sqrt(data_sig2);

function [dist_pro, last_prod] = mass(data_freq, query, ...
    data_len, sub_len, data_mu, data_sig, query_mu, query_sig)
% pre-process query for fft
query = query(end:-1:1);
query(sub_len+1:(sub_len+data_len)) = 0;

% compute the product
query_freq = fft(query);
product_freq = data_freq.*query_freq;
product = ifft(product_freq);

% compute the distance profile
dist_pro = 2 * (sub_len - ...
    (product(sub_len:data_len) - sub_len*data_mu*query_mu)./...
    (data_sig * query_sig));
last_prod = real(product(sub_len:data_len));