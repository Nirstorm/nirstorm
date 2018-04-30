function [ImageGridAmp, OPTIONS] = be_launch_mem2(obj, OPTIONS)


% All samples or a selection?
Data  = obj.data;
if ~isempty(OPTIONS.automatic.selected_samples)
    for ii = 1 : numel(OPTIONS.mandatory.DataTypes)
        Data = obj.data{ii}(:,OPTIONS.automatic.selected_samples(1,:));
    end
end
% time series or wavelet representation ?
nbSmp = size(Data,2);
obj.time = obj.t0+(1:nbSmp)/OPTIONS.automatic.sampling_rate;
obj.scale = zeros(1,nbSmp);

if strcmp(OPTIONS.mandatory.pipeline,'wMEM')
   obj.time = obj.t0+OPTIONS.automatic.selected_samples(5,:);
   obj.scale = OPTIONS.automatic.selected_samples(4,:);
end

% fixed parameters
entropy_drop   = zeros(1,nbSmp);
ImageSourceAmp = zeros(length(obj.iModS), nbSmp);

if OPTIONS.solver.parallel_matlab == 1
    warning off
    tic
    parfor ii = 1 : nbSmp
        [R, E] = MEM_mainLoop(ii, Data, obj, OPTIONS);
        entropy_drop(ii)     = E;
        ImageSourceAmp(:,ii) = R';
    end
    toc
    warning on
else
    fprintf('MEM at each samples (%d samples, may be done in parallel):\nMultiresolution sample (j,t): j=0 corresponds to the sampling scale.\n', nbSmp);
    tic
    for ii = 1 : nbSmp
        [R, E] = MEM_mainLoop(ii, Data, obj, OPTIONS);
        entropy_drop(ii)     = E;
        ImageSourceAmp(:,ii) = R';
    end
    toc
end
% store the results where it should:    
ImageGridAmp = zeros(obj.nb_dipoles, nbSmp);
ImageGridAmp(obj.iModS,:) = ImageSourceAmp;
clear ImageSourceAmp

OPTIONS.automatic.entropy_drops = entropy_drop;
return
end

% =========================================================================

function [R, E] = MEM_mainLoop(ii, Data, obj, OPTIONS)

    % MEM parameters structure
    obj.ind    = ii;
    %obj.data   = Data(obj.iModall,ii); %.*Dunits;
    obj.data   = Data(:,ii);
    obj.scores = obj.SCR(:,ii);
    
    if ~sum(obj.CLS(:,ii))
        disp(['MEM warning: The distributed dipoles could not be clusterized at sample ' num2str(ii) '. (null solution returned)'])
        % Save empty solution
        mem_results_struct.intensities = zeros( size(obj.iModS) );
        entropy_drop = NaN;
    else
         
        obj.clusters = obj.CLS(:,ii);
        obj.active_probability = obj.ALPHA(:,ii);

        % initialize the MEM
        [OPTIONS, mem_structure] = be_memstruct(OPTIONS, obj);

        % solve the MEM for the current time
        [J, mem_results_struct] = be_solve_mem(mem_structure);   
            nclus = max(obj.CLS(:,ii));
            niter = mem_results_struct.iterations;
            entropy_drop = mem_results_struct.entropy;
        fprintf('\t%3d(%2d,%3.3f)\t',ii,obj.scale(ii),obj.time(ii));           
        % Save the solution
        if sum(isnan(J))
            fprintf('killed\n');
            %mem_results_struct.intenties = zeros( size(obj.iModS) );
            mem_results_struct.intensities = J;
            mem_results_struct.intensities(isnan(J)) = 0;
        else
            fprintf('%3d clusters,\t%3d iter.\tEntropy drop:%4.1f\n',nclus,niter,entropy_drop);
            [mem_results_struct.intensities] = deal(J);
        end;
    end
    R = mem_results_struct.intensities;
    E = entropy_drop;
end
