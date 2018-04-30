function [ImageGridAmp, OPTIONS] = be_launch_mem_NIRS(obj, OPTIONS)
% BE_LAUNCH_MEM loops on all time samples and solves the MEM 
%
%   INPUTS:
%       -   obj
%       -   OPTIONS
%
%   OUTPUTS:
%       -   ImageGridAMp    :   MEM solution (matrix NsourcesxNtimesamples)
%       -   OPTIONS
%       -   obj
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

% normalize alpha for each coloum
obj.ALPHA = bsxfun(@rdivide,obj.ALPHA,max(obj.ALPHA,[],1)); % Normalize M

% All samples or a selection?
Data  = [];

if ~isempty(OPTIONS.automatic.selected_samples)        
    Data = [Data;obj.data{1}(:,OPTIONS.automatic.selected_samples(1,:))];                 
elseif ~strcmp(OPTIONS.mandatory.pipeline,'wMEM')
    Data  = obj.data;
end

% time series or wavelet representation ?
nbSmp           = size(Data,2);
ImageSourceAmp  = zeros( length(obj.iModS), 1 );

if strcmp(OPTIONS.mandatory.pipeline,'wMEM')
    obj.time     = OPTIONS.automatic.selected_samples(6,:); % USELESS  ???
    obj.scale    = OPTIONS.automatic.selected_samples(2,:); % USELESS ???
    if ~OPTIONS.wavelet.single_box
        ImageSourceAmp  = sparse(length(obj.iModS), size(obj.data{1},2));
    end
else
    obj.time        = obj.t0+(1:nbSmp)/OPTIONS.automatic.sampling_rate; % USELESS ??
    obj.scale       = zeros(1,nbSmp); % ??????
    ImageSourceAmp  = sparse(length(obj.iModS), nbSmp);  
end

% fixed parameters
entropy_drop    = zeros(1,nbSmp);
final_alpha     = cell(1,nbSmp);
final_sigma     = cell(1,nbSmp);
lambda          = cell(1,nbSmp);

if OPTIONS.solver.parallel_matlab == 1
    warning off
    time_it_starts = tic;
    parfor ii = 1 : nbSmp
        [R, E, A, S, L] = MEM_mainLoop(ii, Data, obj, OPTIONS);
        entropy_drop(ii)    =   E;
        final_alpha{ii}     =   A;
        final_sigma{ii}     =   S;
        lambda{ii}          =   L;
        
        %Store in a matrix
        ImageSourceAmp = ImageSourceAmp + store_solution( R', ii, obj, OPTIONS );
    end
    time_it_ends = toc(time_it_starts);
    if OPTIONS.optional.verbose
        fprintf('%s, Elapsed CPU time is %5.2f seconds.\nBye.\n', OPTIONS.mandatory.pipeline, time_it_ends);
    end
    warning on

else
    if OPTIONS.optional.verbose
        fprintf('%s, MEM at each samples (%d samples, may be done in parallel):', OPTIONS.mandatory.pipeline, nbSmp);
        fprintf('\nMultiresolution sample (j,t): j=0 corresponds to the sampling scale.\n');
    end
    time_it_starts = tic;
    for ii = 1 : nbSmp
        [R, E, A, S, L] = MEM_mainLoop(ii, Data, obj, OPTIONS);
        entropy_drop(ii)	= E;        
        final_alpha{ii}  	= A;
        final_sigma{ii}     = S;
        lambda{ii}          = L;
        % Store in matrix
        ImageSourceAmp      = ImageSourceAmp + store_solution( R', ii, obj, OPTIONS);        
    end
    time_it_ends = toc(time_it_starts);
    if OPTIONS.optional.verbose
        fprintf('%s, Elapsed CPU time is %5.2f seconds.\nBye.\n', OPTIONS.mandatory.pipeline, time_it_ends);
    end
end


if strcmp(OPTIONS.mandatory.pipeline, 'wMEM') && OPTIONS.wavelet.single_box
    ImageGridAmp = [];
    OPTIONS.automatic.wActivation   =   full(ImageSourceAmp);
else
    % store the results where it should:
    ImageGridAmp = zeros( obj.nb_dipoles, size(ImageSourceAmp,2) );

    if  strcmp(OPTIONS.optional.normalization,'adaptive')
        ImageGridAmp(obj.iModS,:) = full(ImageSourceAmp)/OPTIONS.automatic.Modality(1,1).ratioAmp;
    else
        ImageGridAmp(obj.iModS,:) = full(ImageSourceAmp)*OPTIONS.automatic.Modality(1).units_dipoles; %Modified by JSB August 17th 2015
    end
    
    clear ImageSourceAmp
end

OPTIONS.automatic.entropy_drops = entropy_drop;
OPTIONS.automatic.final_alpha   = final_alpha;
OPTIONS.automatic.final_sigma   = final_sigma;
OPTIONS.automatic.lambda        = lambda;


end

% =========================================================================

function [R, E, A, S, L] = MEM_mainLoop(ii, Data, obj, OPTIONS)

    % MEM parameters structure
    obj.ind    = ii;
    %obj.data   = Data(obj.iModall,ii); %.*Dunits;
    obj.data   = Data(:,ii);
    %obj.scores = obj.SCR(:,ii);
    
    if strcmp(OPTIONS.optional.normalization, 'adaptive')       
        obj.Jmne   = OPTIONS.automatic.Modality(1).Jmne(:,ii)./max(abs(OPTIONS.automatic.Modality(1).Jmne(:,ii)));
    end
    
    % check if there's a noise cov for each scale
    if (size(obj.noise_var,3)>1)
        if OPTIONS.optional.verbose
            fprintf('%s, Noise variance at scale %i is selected\n',...
                OPTIONS.mandatory.pipeline,OPTIONS.automatic.selected_samples(2,ii));
        end
        obj.noise_var = squeeze(obj.noise_var(:,:,OPTIONS.automatic.selected_samples(2,ii)) );
    end

    if ~sum(obj.CLS(:,ii))
        if OPTIONS.optional.verbose
            disp(['MEM warning: The distributed dipoles could not be clusterized at sample ' num2str(ii) '. (null solution returned)']);
        end
        % Save empty solution
        J               = zeros( size(obj.iModS) );
        entropy_drop    = NaN;
        act_proba       = NaN;
        act_var         = [];
    else
        
        obj.clusters    = obj.CLS(:,ii);
        obj.active_probability = obj.ALPHA(:,ii);

        % initialize the MEM
        [OPTIONS, mem_structure, act_var] = be_memstruct_NIRS(OPTIONS, obj);
        
        % solve the MEM for the current time
        [J, mem_results_struct] = be_solve_mem(mem_structure);   
            nclus       = max(obj.CLS(:,ii));
            niter       = mem_results_struct.iterations;
            entropy_drop= mem_results_struct.entropy;
            act_proba   = mem_results_struct.active_probability; 
            lambda      = mem_results_struct.lambda;
         if OPTIONS.optional.verbose, fprintf('Sample %3d(%2d,%3.3f):',ii,obj.scale(ii),obj.time(ii)); end;           
        
        % Print output
        if sum(isnan(J))
            fprintf('killed\n');
            %mem_results_struct.intenties = zeros( size(obj.iModS) );
            J(isnan(J)) = 0;
        else
             if OPTIONS.optional.verbose, fprintf('\n\t\t%3d clusters,\n\t\t%3d iter.\n\t\tEntropy drop:%4.1f\n',nclus,niter,entropy_drop); end
        end
        
    end
    R = J;
    E = entropy_drop;
    A = act_proba;
    S = act_var;
    L = lambda;
    
end

function [SOL] = store_solution( vec, ii, obj, OPTIONS )
    
    if strcmp(OPTIONS.mandatory.pipeline, 'wMEM') && OPTIONS.wavelet.single_box
        SOL     =   vec;
    elseif strcmp(OPTIONS.mandatory.pipeline, 'wMEM')
        % Construct wavelet
        nbSmp   =   size(obj.data{1},2);
        scale   =   OPTIONS.automatic.selected_samples(2,ii);
        transl  =   OPTIONS.automatic.selected_samples(3,ii);
        wav     =   zeros( 1, nbSmp );
        wav( nbSmp/2^scale + transl ) = 1;
        wav     =   sparse(be_wavelet_inverse( wav, OPTIONS ));
        SOL     =   sparse(vec) * wav;
        
    else
        SOL     =   sparse( obj.nb_dipoles, size(obj.data, 2) );
        SOL(:,ii) = sparse(vec);
        
    end        
        
end

