function myCplexProgressCallback(info)

    persistent lastProgressPercent;
    persistent initialTime; % To track elapsed time for overall progress
    
    if isempty(lastProgressPercent)
        lastProgressPercent = 0;
        initialTime = tic; % Start timer when callback is first called for a new solve
        bst_progress('set', 0); % Ensure it starts at 0%
    end
    
    % --- Handle completion status first ---
    if info.cplexstatus == 103 || info.cplexstatus == 104 || info.cplexstatus == 105 || info.cplexstatus == 106 || info.cplexstatus == 107
        % CPX_STAT_OPTIMAL        (103)
        % CPX_STAT_INFEASIBLE     (104)
        % CPX_STAT_UNBOUNDED      (105)
        % CPX_STAT_ABORT_TIME_LIM (107) - Reached time limit
        % Add other terminal statuses you deem complete
        bst_progress('set', 100); % Ensure it shows 100% on completion
        bst_progress('stop');
        lastProgressPercent = 0; % Reset for next potential optimization run
        clear initialTime; % Clear persistent timer
        return; % Exit early as optimization is done
    end
    
    
    % --- Calculate current progress ---
    currentProgressPercent = lastProgressPercent; % Default to no change
    
    % 1. Progression basée sur le temps (toujours utile si le gap est lent au début)
    % Assurez-vous que le totalTimeLimit est accessible ici (e.g., passé comme paramètre à la fonction mère ou global)
    % Ou réglez-le directement ici si fixe.
    totalTimeLimit = 300; % <--- Remplacer par votre vrai timelimit global
    
    elapsed = toc(initialTime);
    timeBasedProgress = min(99, round((elapsed / totalTimeLimit) * 100)); % Cap at 99%
    currentProgressPercent = max(currentProgressPercent, timeBasedProgress);
    
    
    % 2. Progression basée sur le gap MIP (si c'est un MIP)
    if isfield(info, 'mip')
        if isfield(info.mip, 'best_obj') && isfield(info.mip, 'incumbent_obj')
            % Only update if an incumbent solution has been found
            if ~isinf(info.mip.incumbent_obj) % Check if a feasible solution exists
                
                % Avoid division by zero/near-zero
                if abs(info.mip.best_obj) < 1e-9 % If objective is very close to zero
                    if abs(info.mip.incumbent_obj) < 1e-9
                        gap = 0; % Already found optimal close to zero
                    else
                        gap = 1; % Still far from zero
                    end
                else
                    if info.Model.sense == 'maximize'
                        gap = (info.mip.best_obj - info.mip.incumbent_obj) / abs(info.mip.best_obj);
                    else % Minimize
                        gap = (info.mip.incumbent_obj - info.mip.best_obj) / abs(info.mip.incumbent_obj);
                    end
                end
                
                % Convert gap to percentage: lower gap = higher progress
                if ~isnan(gap) && gap >= 0 && gap <= 1
                    gapBasedProgress = min(99, round((1 - gap) * 100)); % Cap at 99%
                    currentProgressPercent = max(currentProgressPercent, gapBasedProgress); % Take the max progress
                end
            end
        end
        
        % 3. Progression basée sur les noeuds explorés (si c'est un MIP et vous avez une estimation maxNodes)
        % maxNodesEstimate = 100000; % <-- Adjust this estimate based on your problem size
        % if isfield(info.mip, 'nodes_processed') && info.mip.nodes_processed > 0
        %     nodeBasedProgress = min(99, round((info.mip.nodes_processed / maxNodesEstimate) * 100));
        %     currentProgressPercent = max(currentProgressPercent, nodeBasedProgress);
        % end
    
    elseif isfield(info, 'lp') % For Linear Programs (LP) or Simplex/Barrier phases
        % Simple progress for LP based on iterations (requires estimate for max iterations)
        % maxIterationsEstimate = 500; % <-- Adjust this estimate
        % if isfield(info.lp, 'itercount') && info.lp.itercount > 0
        %    iterBasedProgress = min(99, round((info.lp.itercount / maxIterationsEstimate) * 100));
        %    currentProgressPercent = max(currentProgressPercent, iterBasedProgress);
        % end
    end
    
    
    % --- Update bst_progress if the value has increased ---
    if currentProgressPercent > lastProgressPercent
        bst_progress('set', currentProgressPercent);
        lastProgressPercent = currentProgressPercent;
        % Optional: Update text with current objective/gap
        % bst_progress('text', sprintf('Optimizing... Gap: %.2f%%', gap * 100));
    end

end