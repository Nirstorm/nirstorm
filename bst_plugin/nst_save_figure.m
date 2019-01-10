function nst_save_figure(fig_fn, options, hfig)

if nargin < 3
    hfig = gcf;
end

switch options.save_fig_method
    case 'export_fig'
        if ~isfield(options, 'export_fig_dpi')
            options.export_fig_dpi = 90;
        end
        export_fig(fig_fn, '-transparent', sprintf('-r%d', options.export_fig_dpi));
    case 'saveas'
        saveas(hfig, fig_fn);
    otherwise
        error('Unknown figure saving method: %s', options.save_fig_method);
end

end