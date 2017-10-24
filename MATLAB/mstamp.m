% STOMP Based mSTAMP with Constrained Search Implemented
% Chin-Chia Michael Yeh
%
% [pro_mul, pro_idx] = mstamp(data, sub_len, must_dim, exc_dim)
%
% Output:
%     pro_mul: multidimensional matrix profile (matrix)
%     pro_idx: matrix profile index (matrix)
% Input:
%     data: input time series (matrix)
%     sub_len: interested subsequence length (scalar)
%     must_dim: the dimension which must be included (vector)
%     exc_dim: the dimension which must be excluded (vector)
%
% C.-C. M. Yeh, N. Kavantzas, and E. Keogh, "Matrix Profile VI: Meaningful
% Multidimensional Motif Discovery," IEEE ICDM 2017.
% https://sites.google.com/view/mstamp/
% http://www.cs.ucr.edu/~eamonn/MatrixProfile.html
%

function [pro_mul, pro_idx] = ...
    mstamp(data, sub_len, must_dim, exc_dim)
%% get various length
exc_zone = round(sub_len / 2);
data_len = size(data, 1);
pro_len = data_len - sub_len + 1;
n_dim = size(data, 2);

%% check input
if sub_len > data_len / 2
    error(['Error: Time series is too short relative ', ...
        'to desired subsequence length']);
end
if sub_len < 4
    error('Error: Subsequence length must be at least 4');
end
if any(must_dim > n_dim)
    error(['Error: The must have dimension must be less ', ...
        'then the total dimension']);
end
if any(exc_dim > n_dim)
    error(['Error: The exclusion dimension must be less ', ...
        'then the total dimension']);
end
if ~isempty(intersect(must_dim, exc_dim))
    error(['Error: The same dimension is presented in both ', ...
        'the exclusion dimension and must have dimension']);
end

%% check skip position
n_exc = length(exc_dim);
n_must = length(must_dim);
mask_exc = false(n_dim, 1);
mask_exc(exc_dim) = true;
skip_loc = false(pro_len, 1);
for i = 1:pro_len
    if any(isnan(reshape(data(i:i+sub_len-1, ~mask_exc), 1, []))) ...
            || any(isinf(reshape(data(i:i+sub_len-1, ~mask_exc), 1, [])))
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
pro_mul = zeros(pro_len, n_dim);
pro_idx = zeros(pro_len, n_dim);
dist_pro = zeros(pro_len, n_dim);
last_prod = zeros(pro_len, n_dim);
drop_val = zeros(1, n_dim);
for i = 1:pro_len
    % compute the distance profile
    fprintf('%d %d\n', i, pro_len);
    query = data(i:i+sub_len-1, :);
    if i==1
        for j = 1:n_dim
            [dist_pro(:, j), last_prod(:, j)] = ...
                mass(data_freq(:, j), query(:, j), ...
                data_len, sub_len, data_mu(:, j), ...
                data_sig(:, j), data_mu(i, j), ...
                data_sig(i, j));
        end
    else
        last_prod(2:data_len - sub_len + 1, :) = ...
            last_prod(1:data_len - sub_len, :) ...
            - data(1:data_len - sub_len, :) ...
            .* repmat(drop_val, pro_len - 1, 1) ...
            + data(sub_len + 1:data_len, :) ...
            .* repmat(query(sub_len, :), pro_len - 1, 1);
        last_prod(1, :) = first_prod(i, :);
        dist_pro = 2 * (sub_len - (last_prod ...
            - sub_len * data_mu .* repmat(data_mu(i, :), pro_len, 1)) ...
            ./ (data_sig .* repmat(data_sig(i, :), pro_len, 1)));
    end
    dist_pro = real(dist_pro);
    drop_val(:) = query(1, :);

    % apply exclusion zone
    exc_st = max(1, i - exc_zone);
    exc_ed = min(pro_len, i+exc_zone);
    dist_pro(exc_st:exc_ed, :) = inf;
    dist_pro(data_sig < eps) = inf;
    if skip_loc(i) || any(data_sig(i, ~mask_exc) < eps)
        dist_pro = inf(size(dist_pro));
    end
    dist_pro(skip_loc, :) = inf;

    % apply dimension "must have" and "exclusion"
    dist_pro(:, exc_dim) = inf;
    mask_must = false(n_must, 1);
    mask_must(must_dim) = true;
    dist_pro_must = dist_pro(:, mask_must);
    dist_pro(:, mask_must) = -inf;
    dist_pro_sort = sort(dist_pro, 2);
    dist_pro_sort(:, 1:n_must) = dist_pro_must;

    % figure out and store the nearest neighbor
    dist_pro_cum = zeros(pro_len, 1);
    dist_pro_merg = zeros(pro_len, 1);
    for j = max(1, n_must):(n_dim - n_exc)
        dist_pro_cum = dist_pro_cum + dist_pro_sort(:, j);
        dist_pro_merg(:) = dist_pro_cum / j;
        [min_val, min_idx] = min(dist_pro_merg);
        pro_mul(i, j) = min_val;
        pro_idx(i, j) = min_idx;
    end
end

%% remove bad k setting in the returned matrix
pro_mul = sqrt(pro_mul);
pro_mul(:, 1:(n_must - 1)) = nan;
pro_mul(:, (n_dim - n_exc + 1):end) = nan;
pro_idx(:, 1:(n_must - 1)) = nan;
pro_idx(:, (n_dim - n_exc + 1):end) = nan;


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