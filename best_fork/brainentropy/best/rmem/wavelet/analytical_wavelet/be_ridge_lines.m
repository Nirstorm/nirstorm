function [OPTIONS]  =   be_ridge_lines(SCALO, iD, OPTIONS)

% === Get baseline ridges

% time stamps
Tlims       =   be_closest( OPTIONS.optional.BaselineSegment([1 end]), OPTIONS.mandatory.DataTime );
SLICE       =   SCALO(:,Tlims(1):Tlims(2));

% get maxima
nbF         =   numel( OPTIONS.wavelet.freqs_analyzed);
O           =   OPTIONS;
O.ridges.scalo_threshold    =   0;
O.ridges.min_duration       =   3;
[d, GRPS]   =   be_localmaxima(SLICE, O);
GRPS        =   cellfun( @(a) a+( Tlims(1)-1 )*nbF, GRPS, 'uni', 0 );
OPTIONS.automatic.Modality(iD).ridges_baseline  =   GRPS;


% The following is only for inhouse ridge implementation
if strcmp(OPTIONS.ridges.method, 'LillyOlhede2010')
    return
end


% === Get data ridges

% Get time stamps
Tlims       =   be_closest( OPTIONS.optional.TimeSegment([1 end]), OPTIONS.mandatory.DataTime );
SLICE       =   SCALO(:,Tlims(1):Tlims(2));

% get maxima
[d, GRPS]   =   be_localmaxima(SLICE, OPTIONS);
GRPS        =   cellfun(@(a) a + ( Tlims(1)-1 )*nbF, GRPS, 'uni', 0);
OPTIONS.automatic.Modality(iD).ridges_data  =   GRPS;

return

