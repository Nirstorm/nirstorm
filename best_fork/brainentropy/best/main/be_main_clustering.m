function [OPTIONS, obj] = be_main_clustering(obj, OPTIONS)
% BE_MAIN_CLUSTERING launches the appropriate cortex clustering functions 
% according to the chosen MEM pipeline
%
% Inputs:
% -------
%
%	obj			:	MEM obj structure
%   OPTIONS     :   structure (see bst_sourceimaging.m)
%
%
% Outputs:
% --------
%
%   OPTIONS     :   Updated options fields
%	obj			:	Updated structure
%
%% ==============================================   
% Copyright (C) 2011 - LATIS Team
%
%  Authors: LATIS team, 2011
%
%% ==============================================
% License 
%
% BEst is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    BEst is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with BEst. If not, see <http://www.gnu.org/licenses/>.
% -------------------------------------------------------------------------   
   
global MEMglobal
    
if OPTIONS.optional.groupAnalysis > 0
    [OPTIONS, SCR, CLS, active_probability] = be_groupAnalysis(OPTIONS); %, [hmem 0.05 .8]);
    OPTIONS.automatic.Comment = [OPTIONS.automatic.Comment ' | GrpAnalysis'];
    [dum, studID]   = be_get_id( OPTIONS );
    SCR             = MEMglobal.trMTX{studID} * SCR;
    CLS             = MEMglobal.trMTX{studID} * CLS;
    ALPHA           = MEMglobal.trMTX{studID} * active_probability;
    
elseif ~isfield(OPTIONS.optional.clustering, 'initial_alpha')
    % Sources prescoring - MSP (ref. Mattout et al. 2006) and clustering
    % following different strategies depending of the case cmem, wmem or
    % rmem
    switch OPTIONS.mandatory.pipeline
        case 'cMEM'
            [CLS, SCR, OPTIONS] = be_cmem_clusterize_multim(obj, OPTIONS); %TO DO: CHECK IF DATA SHOULD BE MINUS BSL OR NOT
        case 'wMEM'
            [CLS, SCR, OPTIONS] = be_wfdr_clustering_multim(obj, OPTIONS);
        case 'rMEM'
            [CLS, SCR, OPTIONS] = be_rmem_clusterize_multim(obj, OPTIONS);
    end    
    [ALPHA, CLS, OPTIONS]   = be_scores2alpha(SCR, CLS, OPTIONS);
    
elseif strcmp( OPTIONS.mandatory.pipeline, 'wMEM' )
    ALPHA = OPTIONS.optional.clustering.initial_alpha * ones(1,size(OPTIONS.automatic.Modality(1).selected_jk, 2));
    CLS   = OPTIONS.optional.clustering.clusters * ones(1,size(OPTIONS.automatic.Modality(1).selected_jk, 2));
    SCR   = [];
    
else
    ALPHA = OPTIONS.optional.clustering.initial_alpha * ones(1,size(OPTIONS.automatic.Modality(1).data, 2));
    CLS   = OPTIONS.optional.clustering.clusters * ones(1,size(OPTIONS.automatic.Modality(1).data, 2));
    SCR   = [];
end

% the final scores (SCR), clusters (CLS) and alpha's (ALPHA)
obj.SCR   = SCR;
obj.CLS   = CLS;
obj.ALPHA = ALPHA;







 
