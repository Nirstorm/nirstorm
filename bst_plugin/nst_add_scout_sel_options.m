function options = nst_add_scout_sel_options(options, label, comment, surf_type, atlas_names, add_hint)

assert(isvarname(label));

if strcmp(surf_type, 'scalp')
    s_isurf = 'iScalp';
elseif strcmp(surf_type, 'cortex')
    s_isurf = 'iCortex';
else
    disp('ERROR: surf_type must either be "scalp" or "cortex"');
    bst_error('surf_type must either be "scalp" or "cortex"');
end

if nargin < 6
   add_hint = 0; 
end

% Parse all subjects in current protocols and make the list of
% all available scalp scouts in atlas "User scouts".
pSubjects = bst_get('ProtocolSubjects');
scout_list = {};
iscout_info = 1;
scout_info = struct([]);
for isubject=1:length(pSubjects.Subject)
    sSubject = pSubjects.Subject(isubject);
    if ~isempty(sSubject.Surface) && ~isempty(sSubject.(s_isurf)) % no surface available
        sSurf = in_tess_bst(sSubject.Surface(sSubject.(s_isurf)).FileName);
        for iatlas=1:length(atlas_names)
            iatlas_selected = strcmp(atlas_names{iatlas}, {sSurf.Atlas.Name});
            if any(iatlas_selected)
                scouts = sSurf.Atlas(iatlas_selected).Scouts;
                for iscout=1:length(scouts)
                    scout_list{end+1} = [sSubject.Name ' > ' atlas_names{iatlas} ' > ' scouts(iscout).Label];
                    % Save info on selected scout to retrieve them after user selection
                    scout_info(iscout_info).sSubject = sSubject;
                    scout_info(iscout_info).isurface = sSubject.(s_isurf);
                    scout_info(iscout_info).iatlas = iatlas_selected;
                    scout_info(iscout_info).iscout = iscout;
                    scout_info(iscout_info).sScout = scouts(iscout);
                    iscout_info = iscout_info + 1;
                end
            end
        end
    end
end

option_name = ['scout_sel_' label];

if add_hint==1
    msg = '<B>If scout does not appear in list, cancel process, use File > Reload process definitions,</B><BR>';
    options.([option_name '_hint1']).Comment = msg;
    options.([option_name '_hint1']).Type = 'label';
    
    msg = '<B>Or run the matlab command: panel_process_select(''ParseProcessFolder'', 1);</B><BR>';
    options.([option_name '_hint2']).Comment = msg;
    options.([option_name '_hint2']).Type = 'label';
    
    msg = '<B>Then run this process again</B><BR>';
    options.([option_name '_hint3']).Comment = msg;
    options.([option_name '_hint3']).Type = 'label';
end

options.(option_name).Comment = comment;
options.(option_name).Type    = 'combobox';
if isempty(scout_list)
    scout_list = {'No scout found'};
end
options.(option_name).Value   = {1, scout_list};

% Keep information on scout selection to retrieve it in after user input
options.([option_name '_info']).Type = 'value';
options.([option_name '_info']).Value = scout_info;
options.([option_name '_info']).Hidden = 1;
end
