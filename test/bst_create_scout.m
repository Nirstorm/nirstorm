function new_scout_sel = bst_create_scout(subject_name, surface_type, scout_name, scout_seed, scout_size, atlas_name)
assert(isvarname(scout_name));

if strcmp(surface_type, 'scalp')
    s_isurf = 'iScalp';
elseif strcmp(surface_type, 'cortex')
    s_isurf = 'iCortex';
else
    disp('ERROR: surf_type must either be "scalp" or "cortex"');
    bst_error('surf_type must either be "scalp" or "cortex"');
end

if nargin < 6
    atlas_name = 'User scouts';
end

[sSubject, iSubject] = bst_get('Subject', subject_name);

sSurf = in_tess_bst(sSubject.Surface(sSubject.(s_isurf)).FileName);
iatlas = strcmp(atlas_name, {sSurf.Atlas.Name});
if isempty(iatlas)
    error(['Atlas not found:' atlas_name]);
end

assert(scout_seed <= size(sSurf.Vertices, 1));
scout_vertices = [scout_seed];
for isize=1:scout_size
    cur_vertices = scout_vertices;
    for ivertex=1:length(cur_vertices)
        vidx = cur_vertices(ivertex);
        scout_vertices = unique([scout_vertices find(sSurf.VertConn(vidx, :))]);
    end
end
scout_idx = size(sSurf.Atlas(iatlas).Scouts,2) + 1;
new_scout = db_template('Scout');
new_scout.Vertices = scout_vertices;
new_scout.Seed = scout_seed;
new_scout.Color = [0,0,0];
new_scout.Label = scout_name;
sSurf.Atlas(iatlas).Scouts(scout_idx) = new_scout;
bst_save(file_fullpath(sSubject.Surface(sSubject.(s_isurf)).FileName), sSurf, 'v7');
db_save();

new_scout_sel.sSubject = bst_get('Subject', subject_name);
new_scout_sel.isurface = sSubject.(s_isurf);
new_scout_sel.iatlas = iatlas;
new_scout_sel.iscout = scout_idx;
new_scout_sel.sScout = new_scout;
end
