function print_progress(done, total, e, n, t, r, num_runs)
    percent  = 100 * done / total;
    bar_len  = 30;
    filled   = round(bar_len * done / total);
    bar_str  = [repmat('#', 1, filled), repmat('-', 1, bar_len - filled)];
    fprintf('\r[%s] %6.2f%% | exp=%d | n=%d | type=%d | run=%d/%d', ...
        bar_str, percent, e, n, t, r, num_runs);
    if done == total; fprintf('\n'); end
end
