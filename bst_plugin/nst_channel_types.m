function chan_types = nst_channel_types()
% NST_CHANNEL_TYPES return an enumeration of channel types
%
% CHANNEL_TYPES = NST_CHANNEL_TYPES()
%    CHANNEL_TYPES: struct with numerical fields listing all available
%                   channel types:
%                   - CHANNEL_TYPES.WAVELENGTH
%                   - CHANNEL_TYPES.Hb

chan_types.WAVELENGTH = 1;
chan_types.HB = 2;
end