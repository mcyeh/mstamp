% Plot the Motifs on the data
% Chin-Chia Michael Yeh
%
% plot_motif_on_data(data, sub_len, motif_idx, motif_dim)
%
% Input:
%     data: input time series (matrix)
%     sub_len: interested subsequence length (scalar)
%     motif_idx: the index for the founded motifs (matrix)
%     motif_dim: the dimensions spanned by the found motifs (cell)
%
% C.-C. M. Yeh, N. Kavantzas, and E. Keogh, "Matrix Profile VI: Meaningful
% Multidimensional Motif Discovery," IEEE ICDM 2017.
% https://sites.google.com/view/mstamp/
% http://www.cs.ucr.edu/~eamonn/MatrixProfile.html
%

function plot_motif_on_data(data, sub_len, motif_idx, motif_dim)
figure();
ax = axes();
hold(ax, 'on');

%% plot the data
for i = 1:size(data, 2)
    data(:, i) = data(:, i) - min(data(:, i));
    data(:, i) = data(:, i) / max(data(:, i));
    data(:, i) = data(:, i) + (i - 1) * 1.1;
    plot(data(:, i), 'color', 'k');
end


for i = 1:length(motif_idx)
    for k = 1:length(motif_dim{i})
        motif_location = motif_idx(i):motif_idx(i) + sub_len - 1;
        motif = data(motif_location, motif_dim{i}(k));
        plot(motif_location, motif, 'color', 'r');
    end
end

hold(ax, 'off');