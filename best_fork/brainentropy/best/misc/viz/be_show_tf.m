function be_show_tf(TF, TA, FA, wtp, varargin)
    
if (numel(varargin)>0)&(ishandle(varargin{1}))
    ax = varargin{1};
else
    figure; 
    ax = axes;
end

FLIMS       =   FA([1 end]);
if numel(varargin)>1
    FLIMS   =   varargin{2};
end

XLIMS       =   TA([1 end]);
if numel(varargin)>2 && ~isempty(varargin{3})
    XLIMS   =   varargin{3};
end

yoff        =   0;
if numel(varargin)>3 && ~isempty(varargin{4})
    yoff    =   varargin{4};
end

TIT         =   false;
if numel(varargin)>4 && ischar(varargin{5})
    TIT     =   varargin{5};
end

if numel(varargin)>5 && iscell(varargin{6})
    rdgs    =   varargin{6};
end

% TF/ridge plane limits
Ytic  = 2.^(fix( log2(FA(1)) ):fix( log2(FA(end)) ));
temp  = num2str(Ytic','%4.0f');
Yticl = mat2cell( temp, ones(1, numel(Ytic)), size(temp,2) );
minF  = 1;

% Draw plane
switch wtp
    case 'cwt'
        imagesc(TA, log2(FA), TF,'parent', ax );
        if ~TIT
            title(' Average wavelet transform ','FontSize', 16, 'FontWeight', 'bold');
        else
            title(TIT,'FontSize', 16, 'interpreter', 'latex');
        end
    
    case 'rdg'
        imagesc(TA, log2(FA), TF,'parent', ax );
        if ~TIT
            title(' Multivariate ridge plan ','FontSize', 16, 'FontWeight', 'bold');
        else
            title(TIT,'FontSize', 16, 'interpreter', 'latex');
        end
    
        
    case 'mask'
        mask    =   zeros( size(TF) );
        mask([rdgs{:}]) = 1;
        imagesc(TA, log2(FA), mask,'parent', ax );
        hold on
        nfrqs   =   numel( FA );
        for jj  =   1 : numel(rdgs)
            rl  =   rdgs{jj};
            tl  =   TA( fix( mean(rl/nfrqs) ) );
            fl  =   FA( fix( mean( mod(rl, nfrqs) ) ) +3 );
            text(tl, log2(fl), num2str(jj), 'FontSize', 12, 'Color', [1 1 1], 'horizontalAlignment', 'center', 'interpreter', 'latex');
        end
        if TIT
            title('Thresholded ridge plane ','interpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold');
        end
        colormap('gray')
end

% get CLIM
iD      = be_closest(FLIMS, FA);
tiD 	= be_closest(XLIMS, TA);
CLIM    = [min( min(TF(iD+1:end,tiD(1):tiD(end))) ) max( max(TF(iD+1:end,tiD(1):tiD(end))) )];

set(ax, 'XLim', TA([1 end]), 'YGrid','on', 'YLim',log2([minF FA(end)]),...
        'YTick',log2(Ytic(:)), 'YTickLabel',Yticl, 'YDir', 'normal', ...
        'FontName','Helvetica', 'FontSize',16, ...
        'FontWeight', 'bold', 'clim', CLIM, ...
        'Position', get(gca, 'Position')+[0 .02 0 -.02] );

xlabel('time (s)', 'FontSize',16, 'FontWeight', 'demi'); 
ylabel('frequency (Hz)','FontSize',16, 'FontWeight', 'demi');
set(ax, 'ylim', log2(FLIMS), 'xlim', XLIMS )

if yoff
    set(ax,'YtickLabel', '')
    ylabel('')
end


        
return