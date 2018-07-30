function varargout = process_nst_compute_group_ttest( varargin )
% process_nst_compute_group_ttest
% Compute a ttest on cB. 
% Mean computiation do the following process : 
% compute mean and covariance of the cB to calculate the following ttest :  
% mean(cB)/ sqrt(var(cB)) with n_subject-1 degree of freedom
%
% @=============================================================================
% This function is part of the Brainstorm software:
% http://neuroimage.usc.edu/brainstorm
% 
% Copyright (c)2000-2017 University of Southern California & McGill University
% This software is distributed under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPLv3
% license can be found at http://www.gnu.org/copyleft/gpl.html.
% 
% FOR RESEARCH PURPOSES ONLY. THE SOFTWARE IS PROVIDED "AS IS," AND THE
% UNIVERSITY OF SOUTHERN CALIFORNIA AND ITS COLLABORATORS DO NOT MAKE ANY
% WARRANTY, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, NOR DO THEY ASSUME ANY
% LIABILITY OR RESPONSIBILITY FOR THE USE OF THIS SOFTWARE.
%
% For more information type "brainstorm license" at command prompt.
% =============================================================================@
%
% Authors: Edouard Delaire, 2018
%  
eval(macro_method);
end



%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'Compute Group Analysis';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'NIRS - wip';
    sProcess.Index       = 1403;
    sProcess.isSeparator = 0;
    sProcess.Description = 'https://github.com/Nirstorm/nirstorm/wiki/%5BWIP%5D-GLM-implementation';
    % todo add a new tutorials
    
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'data'};
    sProcess.OutputTypes = {'data','raw'};
    
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 2;
   
    
%     todo : ask the subjects to include / exclude     
%     sProcess.options.subjects.Comment = 'Inclued Sujbets';
%     sProcess.options.subjects.Type    = 'subjectname';
%     sProcess.options.subjects.Value   = '';

    % === Method for the group analysis
    sProcess.options.method.Comment  = {'Mean',''; ...
                                      'mean', ''};
    sProcess.options.method.Type     = 'radio_linelabel';
    sProcess.options.method.Value    = 'mean';
    
    
    % === TAIL FOR THE TEST STATISTIC
    sProcess.options.tail.Comment  = {'One-tailed (-)', 'Two-tailed', 'One-tailed (+)', ''; ...
                                      'one-', 'two', 'one+', ''};
    sProcess.options.tail.Type     = 'radio_linelabel';
    sProcess.options.tail.Value    = 'two';
    
    
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) 
    Comment = sProcess.Comment; 
end

function OutputFiles = Run(sProcess, sInputs)
    OutputFiles={};
    
    % parse input : 

	for i=1:numel(sInputs)
        s_B(i)=in_bst_data(sInputs(i).FileName);
        B(:,i)=s_B(i).F; % might need to transpose. Have to check how cB is registered 
	end     

    n_chan=size(B,1);
    n_subject=size(B,2);
    
    switch sProcess.options.method.Value 
        case 'mean' 
            [t,df]=compute_mean(B);
    end    

    % Calculate p-values from t-values
    p = ComputePvalues(t, df, 't',   sProcess.options.tail.Value );
   
    % Change folder to inter study. 
    iStudy = sInputs(1).iStudy;

    [tmp, iSubject] = bst_get('Subject', sInputs(1).SubjectName);
    [sStudyIntra, iStudyIntra] = bst_get('AnalysisInterStudy', iSubject);
    
    [ChannelFile] = bst_get('ChannelFileForStudy', iStudy);
    [tmp, iChannelStudy] = bst_get('ChannelForStudy', iStudyIntra);
    db_set_channel(iChannelStudy, ChannelFile, 0, 0);
    
    
    % === OUTPUT STRUCTURE ===
    % Initialize output structure

    sOutput = db_template('statmat');
    sOutput.pmap         = [p;p]';
    sOutput.tmap         = [t;t]';
    sOutput.df           = ones(n_chan,1)*df;
    sOutput.ChannelFlag= ones(1,n_chan);
    sOutput.Correction   = 'no';
    sOutput.Type         = 'data';
    sOutput.Time         = [1];
    sOutput.ColormapType = 'stat2';
    sOutput.DisplayUnits = 't';
    sOutput.Options.SensorTypes = 'NIRS';

    
    % Formating a readable comment such as -Rest +Task
    comment=['Group Stat (' num2str(n_subject) ' subjects), '] ;
    comment=[comment  s_B(1).Comment(1:end-2) ];

    switch sProcess.options.tail.Value 
        case {'one-'}
             comment=[ comment ' < 0 '];      
        case {'two'}
             comment=[ comment ' <> 0 ']; 
        case { 'one+'}
             comment=[ comment ' > 0 ']; 
    end  
    
    sOutput.Comment=comment;
    sOutput = bst_history('add', sOutput, s_B(1).History, '');

    sOutput = bst_history('add', sOutput, 'Group stat', comment);
    OutputFiles{1} = bst_process('GetNewFilename', fileparts(sStudyIntra.FileName), 'pdata_ttest_matrix');
    save(OutputFiles{1}, '-struct', 'sOutput');
    db_add_data(iStudyIntra, OutputFiles{1}, sOutput);

end

function [t,df]=compute_mean(B)
  
    n_chan=size(B,1);
    n_subject=size(B,2);

    mean_B= mean(B,2);
    varB=  var(B,0,2) ;

    df=n_subject-1;
    t=zeros(1,n_chan);
    
    for i = 1:n_chan
        t(i)= mean_B(i) / sqrt(varB(i)) ; 
    end

end



%% ===== COMPUTE P-VALUES ====
% see process_test_parametric2 for more information
function p = ComputePvalues(t, df, TestDistrib, TestTail)
    % Default: two-tailed tests
    if (nargin < 4) || isempty(TestTail)
        TestTail = 'two';
    end
    % Default: F-distribution
    if (nargin < 3) || isempty(TestDistrib)
        TestDistrib = 'f';
    end
    % Nothing to test
    if strcmpi(TestTail, 'no')
        p = zeros(size(t));
        return;
    end
    
    % Different distributions
    switch lower(TestDistrib)
        % === T-TEST ===
        case 't'
            % Calculate p-values from t-values 
            switch (TestTail)
                case 'one-'
                    % Inferior one-tailed t-test:   p = tcdf(t, df);
                    % Equivalent without the statistics toolbox (FieldTrip formula)            
                    p = 0.5 .* ( 1 + sign(t) .* betainc( t.^2 ./ (df + t.^2), 0.5, 0.5.*df ) );
                case 'two'
                    % Two-tailed t-test:     p = 2 * (1 - tcdf(abs(t),df));
                    % Equivalent without the statistics toolbox
                    p = betainc( df ./ (df + t .^ 2), df./2, 0.5);
                    % FieldTrip equivalent: p2 = 1 - betainc( t.^2 ./ (df + t.^2), 0.5, 0.5.*df );
                case 'one+'
                    % Superior one-tailed t-test:    p = 1 - tcdf(t, df);
                    % Equivalent without the statistics toolbox (FieldTrip formula)
                    p = 0.5 .* ( 1 - sign(t) .* betainc( t.^2 ./ (df + t.^2), 0.5, 0.5.*df ) );
            end
            
        % === F-TEST ===
        case 'f'
            v1 = df{1};
            v2 = df{2};
            % Evaluate for which values we can compute something
            k = ((t > 0) & ~isinf(t) & (v1 > 0) & (v2 > 0));
            % Initialize returned p-values
            p = ones(size(t));                    
            % Calculate p-values from F-values 
            switch (TestTail)
                case 'one-'
                    % Inferior one-tailed F-test
                    % p = fcdf(t, v1, v2);
                    p(k) = 1 - betainc(v2(k)./(v2(k) + v1(k).*t(k)), v2(k)./2, v1(k)./2);
                case 'two'
                    % Two tailed F-test
                    % p = 2*min(fcdf(F,df1,df2),fpval(F,df1,df2))
                    p(k) = 2 * min(...
                            1 - betainc(v2(k)./(v2(k) + v1(k).*t(k)), v2(k)./2, v1(k)./2), ...
                            1 - betainc(v1(k)./(v1(k) + v2(k)./t(k)), v1(k)./2, v2(k)./2));
                case 'one+'
                    % Superior one-tailed F-test
                    % p = fpval(t, v1, v2);
                    %   = fcdf(1/t, v2, v1);
                    p(k) = 1 - betainc(v1(k)./(v1(k) + v2(k)./t(k)), v1(k)./2, v2(k)./2);
            end
            
        % === CHI2-TEST ===
        case 'chi2'
            % Calculate p-values from Chi2-values 
            %   chi2cdf(x,n) = gammainc(t/2, n/2)
            switch (TestTail)
                case 'one-'
                    % Inferior one-tailed Chi2-test:    p = gammainc(t./2, df./2);
                    error('Not relevant.');
                case 'two'
                    % Two-tailed Chi2-test
                    error('Not relevant.');
                case 'one+'
                    % Superior one-tailed Chi2-test:    p = 1 - gammainc(t./2, df./2);
                    p = 1 - gammainc(t./2, df./2);
            end
    end
end

