function varargout = process_nst_compute_group_ttest( varargin )
% process_nst_compute_group_ttest
% Compute a ttest on the mean of the B according to a constrast vector using 
% t = cBm / sqrt( c Cov(B)m c^T )  
% with Bm = sum(B)/n and Cov(B)m ! sum(Cov(B))/n^2 
%
% B, cov(B) and the corresponding degrree of freedom are estimated inprocess_nst_compute_glm.
% 
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
    sProcess.InputTypes  = {'matrix'};
    sProcess.OutputTypes = {'data','raw'};
    
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 2;
   
    
%     todo : ask the subjects to include / exclude     
%     sProcess.options.subjects.Comment = 'Inclued Sujbets';
%     sProcess.options.subjects.Type    = 'subjectname';
%     sProcess.options.subjects.Value   = '';
    
    
    sProcess.options.Contrast.Comment = 'Contrast vector';
    sProcess.options.Contrast.Type    = 'text';
    sProcess.options.Contrast.Value   = '[-1 1]';
    
    % === TAIL FOR THE TEST STATISTIC
    sProcess.options.tail.Comment  = {'One-tailed (-)', 'Two-tailed', 'One-tailed (+)', ''; ...
                                      'one-', 'two', 'one+', ''};
    sProcess.options.tail.Type     = 'radio_linelabel';
    sProcess.options.tail.Value    = 'two';
    
    
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) 
    Comment = sProcess.Comment;
    Comment = [ Comment ' C = ' sProcess.options.Contrast.Value ]; 
end

function OutputFiles = Run(sProcess, sInputs)
    OutputFiles={};
    

    % Number of B, and Cov(B) 
    n_beta=0;
    n_covb=0;
    
    % Todo add some verification on the channels, matrix... 
    
    % parse input : 
    for i=1:length(sInputs) 
        name= strsplit(sInputs(i).Comment,' ');
        if( strcmp(name(1), 'covB') == 1)
            n_covb= n_covb +1;

            covB(n_covb)=in_bst_data(sInputs(i).FileName);
            
            name= strsplit(cell2mat(name(end)),'=');
            name= strsplit(cell2mat(name(end)),'_');

            df(n_covb)=str2num(cell2mat(name(1)));
        elseif ( strcmp(name(1), 'B') == 1)
            n_beta= n_beta +1;
            B(n_beta)=in_bst_data(sInputs(i).FileName);
        end
     end   
    if( n_beta~=n_covb || n_beta < 2 ) 
       bst_report('Error', sProcess, sInputs, 'This process require beta and its covariance ');
    end    
    

    [B,covB,df]=compute_mean(B,covB,df);
    
    n_cond=size(B.Value',1);
    n_chan=size(B.Value',2);
    
    % exctract the constrast vector. 
    if( strcmp( sProcess.options.Contrast.Value(1),'[') && strcmp( sProcess.options.Contrast.Value(end),']') )
        % The constrast vector is in a SPM-format : 
        % sProcess.options.Contrast.Value = '[1,0,-1]'
        C=strsplit( sProcess.options.Contrast.Value(2:end-1),',');
        C=str2num(cell2mat(C));
    else 
        C=[];
        % Parse the input 
        % Accepted form are : 'X event1 Y event2' 
        % where X,Y is a sign(+,-), a signed number(-1,+1) or a decimal(+0.5,-0.5)
        % event12 are the name of regressor present in the design matrix
        % A valid input can be '- rest + task' ( spaces are important)
        
        expression='[-+]((\d+.\d+)|((\d)+)|)\s+\w+';
        [startIndex,endIndex] = regexp( sProcess.options.Contrast.Value , expression ); 

        for(i=1:length(startIndex))
            word=sProcess.options.Contrast.Value(startIndex(i):endIndex(i)); % can be '-rest','+rest'..
            
            [evt_ind_start,evt_ind_end]=regexp(word, '\s+\w+' );            
            evt_name=word(evt_ind_start+1:evt_ind_end);

            
            % Find the weight of the regressor 
            if strcmp(word(1:evt_ind_start-1),'+')
                evt_coef=1;
            elseif strcmp(word(1:evt_ind_start-1),'-')
                evt_coef=-1;
            else
                evt_coef=str2double(word(1:evt_ind_start));
            end
            
            
            %Find the position of the regressor            
            ind=find(strcmp(B.Description,evt_name))
            if( isempty(ind) || ~isnumeric(evt_coef) )
               bst_report('Error', sProcess, sInputs, [ 'Event ' evt_name ' has not been found']);
               return;
            end
            
            C(ind)=evt_coef;
        end

        if isempty(C)
            bst_report('Error', sProcess, sInputs, 'The format of the constrast vector (eg [-1 1] ) is not recognized');
            return
        end
   end
    
    % Add zero padding for the trend regressor 
    if length(C) < n_cond
       C= [C zeros(1, n_cond - length(C)) ]; 
    end    
     
    B.Value=C*B.Value';
    t=zeros(1,n_chan);
    
    for i = 1:n_chan
        t(i)= B.Value(i) / sqrt( C*covB.Value(:,:,i)*transpose(C) ) ; 
    end
    
    % Calculate p-values from t-values
    p = ComputePvalues(t, df, 't',   sProcess.options.tail.Value );
    df=ones(n_chan,1)*df;
    
    
    % Change folder to inter study. 
    iStudy = sInputs(1).iStudy;

    [tmp, iSubject] = bst_get('Subject', sInputs.SubjectName);
    [sStudyIntra, iStudyIntra] = bst_get('AnalysisInterStudy', iSubject);
    
    [ChannelFile] = bst_get('ChannelFileForStudy', iStudy);
    [tmp, iChannelStudy] = bst_get('ChannelForStudy', iStudyIntra);
    db_set_channel(iChannelStudy, ChannelFile, 0, 0);
    
    
    % === OUTPUT STRUCTURE ===
    % Initialize output structure
    
    
    sOutput = db_template('statmat');
    sOutput.pmap         = [p;p]';
    sOutput.tmap         = [t;t]';
    sOutput.df           = df;
    sOutput.ChannelFlag= ones(1,n_chan);
    sOutput.Correction   = 'no';
    sOutput.Type         = 'data';
    sOutput.Time         = [1];
    sOutput.ColormapType = 'stat2';
    sOutput.DisplayUnits = 't';
    sOutput.Options.SensorTypes = 'NIRS';

    
    % Formating a readable comment such as -Rest +Task
    comment=['Group Stat (' num2str(n_beta) ' subjects), '] ;
    for i=1:n_cond
        if ( C(i) < 0)
            if( C(i) == -1 )
                comment=[  comment  ' - ' cell2mat(B.Description(i)) ' '];
            else
                comment=[  comment num2str(C(i)) ' ' cell2mat(B.Description(i)) ' '];
            end
        elseif ( C(i) > 0 )
            if( C(i) == 1)
                comment=[  comment  ' + ' cell2mat(B.Description(i)) ' '];  
            else
                comment=[  comment  ' + ' num2str(C(i)) ' ' cell2mat(B.Description(i)) ' '];
            end 
        end     
    end
    
    switch sProcess.options.tail.Value 
        case {'one-'}
             comment=[ comment ' < 0 '];      
        case {'two'}
             comment=[ comment ' <> 0 ']; 
        case { 'one+'}
             comment=[ comment ' > 0 ']; 
    end  
    
    sOutput.Comment=comment;
    sOutput = bst_history('add', sOutput, B.History, '');

    sOutput = bst_history('add', sOutput, 'ttest computation', comment);
    OutputFiles{1} = bst_process('GetNewFilename', fileparts(sStudyIntra.FileName), 'pdata_ttest_matrix');
    save(OutputFiles{1}, '-struct', 'sOutput');
    db_add_data(iStudyIntra, OutputFiles{1}, sOutput);

end

function [meam_B,mean_covB,mean_df]=compute_mean(B,covB,df)
   
    % Initialize the structure
    meam_B=B(1);
    mean_covB=covB(1);
    mean_df=df(1);
    
    n=numel(B);
    
    % Calculate the mean
    for i =2:n
        meam_B.Value=meam_B.Value+B(i).Value;
        mean_covB.Value=mean_covB.Value+covB(i).Value;
        mean_df=mean_df+df(i);
    end
    
    
    meam_B.F=meam_B.F/n;
    mean_covB.F=mean_covB.F/n;
    mean_df=mean_df/n;
   
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

