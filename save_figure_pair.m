
function save_figure_pair(fig, base_folder, latest_folder, name, formats)
    for i = 1:numel(formats)
        fmt = formats{i};
        switch fmt
            case 'fig'
                savefig(fig, fullfile(base_folder, [name '.fig']));
                savefig(fig, fullfile(latest_folder, [name '.fig']));
            otherwise
                saveas(fig, fullfile(base_folder, [name '.' fmt]));
                saveas(fig, fullfile(latest_folder, [name '.' fmt]));
        end
    end
end
