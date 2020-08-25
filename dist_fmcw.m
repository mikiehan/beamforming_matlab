function dist = dist_fmcw(Y, params)

Fs = params.Fs;
c3 = params.c3; % multiplier for converting frequency to dist
dist_offset = params.dist_offset; % distance offset introduced by bela recording delay
dist_offset2 = params.dist_offset2;
enable_plot = params.enable_plot; % if plotting the profile

[N, L] = size(Y);
w = 0.5 - 0.5 * cos(2 * pi / L * (0 : L - 1));
start_idx = 150; % distance range
end_idx = 260; % distance range

dist = zeros(1, N);
padding = zeros(1, Fs - L);
for i = 1 : N
    Yext = [Y(i, :) .* w, padding];
    S = fft(Yext);
    S = abs(S(start_idx : end_idx));
    [~, idx] = max(S);
    %if (idx < numel(S)) && (idx > 1)
        %f = localRegressionComp(S, idx) + idx - 1 + start_idx;
    %else
        f = idx - 1 + start_idx;
    %end
        dist(i) = (f / c3 - 2 * dist_offset - i * dist_offset2) / 2;
    if (enable_plot)
        plot(S);
    end
end