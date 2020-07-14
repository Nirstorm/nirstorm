function nst_glm_display_model(model,mode)
    n_regressor=size(model.X,2);
    if strcmp(mode,'matrix')
        Tag='DesignModelMatrix';
    else
        Tag='DesignModelTime';
    end
    hFig = findobj(0, 'Type', 'Figure', 'Tag', Tag);
    if isempty(hFig)
        hFig = figure(...
            'MenuBar',     'none', ...
            ... 'Toolbar',     'none', ...
            'Toolbar',     'figure', ...
            'NumberTitle', 'off', ...
            'Name',        sprintf('Design Matrix'), ...
            'Tag',         Tag, ...
            'Units',       'Pixels');
        % Figure already exists: re-use it
    else
        clf(hFig);
        figure(hFig);
    end
    
    % Plot frequency response
    hAxes = axes('Units', 'pixels', 'Parent', hFig, 'Tag', 'AxesImpz');
    
    switch mode
    case 'matrix'
        h = imagesc(1:n_regressor,model.time(end:-1:1),model.X, ...
            'Parent',              hAxes);
        colormap('gray');
        set(gca,'xtick',[]);
        
            % Enable zooming by default
        zoom(hFig, 'on');
        set(hFig, bst_get('ResizeFunction'), @ResizeCallback);
        ResizeCallback(hFig);
    
    case 'timecourse'    
        for i_reg=1:n_regressor
            subplot(n_regressor,1,i_reg,'Parent',hFig)
            plot(model.time,model.X(:,i_reg))
            xlim([model.time(1) model.time(end)])
            title(model.reg_names{i_reg},'interpreter', 'none')
        end
    end  

    function ResizeCallback(hFig, ev)
        % Get figure position
        figpos = get(hFig, 'Position');
        textH = 110;        % Text Height
        marginL = 70;
        marginR = 30;
        marginT = 30;
        marginB = 50;
        axesH = round((figpos(4) - textH) ./ 2);
        % Position axes
        set(hAxes, 'Position', max(1, [marginL, textH + marginB + axesH, figpos(3) - marginL - marginR, axesH - marginB - marginT]));
        %set(hLabel1,  'Position', max(1, [marginL, textH + marginB,         figpos(3) - marginL - marginR, axesH - marginB - marginT]));
        %set(hLabel1,    'Position', max(1, [40,                  1,  round((figpos(3)-40)/2),  textH]));
        %set(hLabel2,    'Position', max(1, [round(figpos(3)/2),  1,  round(figpos(3)/2),       textH]));
    end

end

