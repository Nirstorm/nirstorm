function varargout = process_nst_compute_ttest( varargin )
% process_nst_compute_ttest
% Compute a ttest according to a spm-like constrast vector using 
% t = cB / sqrt( c Cov(B) c^T )  
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
    sProcess.Comment     = 'Compute Subject Analysis';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'NIRS - wip';
    sProcess.Index       = 1402;
    sProcess.isSeparator = 0;
    sProcess.Description = 'https://github.com/Nirstorm/nirstorm/wiki/%5BWIP%5D-GLM-implementation';
    % todo add a new tutorials
    
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'matrix'};
    sProcess.OutputTypes = {'data','raw'};
    
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
   
    
    sProcess.options.Contrast.Comment = 'Contrast vector';
    sProcess.options.Contrast.Type    = 'text';
    sProcess.options.Contrast.Value   = '[-1 1]';
    
    % === tail for the ttest 
    sProcess.options.tail.Comment  = {'One-tailed (-)', 'Two-tailed', 'One-tailed (+)', ''; ...
                                      'one-', 'two', 'one+', ''};
    sProcess.options.tail.Type     = 'radio_linelabel';
    sProcess.options.tail.Value    = 'one+';

    
    
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) 
    Comment = sProcess.Comment;
    Comment = [ Comment ' C = ' sProcess.options.Contrast.Value ]; 
end

function OutputFiles = Run(sProcess, sInputs)
    OutputFiles={};
    
    % check if Beta and covV are in sInputs
    check_beta=0;
    check_covb=0;
    
    % parse input : 
    for i=1:length(sInputs) 
        name= strsplit(sInputs(i).Comment,' ');
        if( strcmp(name(1), 'covB') == 1)
            covB=in_bst_data(sInputs(i).FileName);
            
            name= strsplit(cell2mat(name(end)),'=');
            name= strsplit(cell2mat(name(end)),'_');

            df=str2num(cell2mat(name(1)));
            check_covb=1;
        elseif ( strcmp(name(1), 'B') == 1)
            B=in_bst_data(sInputs(i).FileName);
            check_beta=1;
        end
     end   
    if( ~check_beta || ~check_covb ) 
       bst_report('Error', sProcess, sInputs, 'This process require beta and its covariance ');
    end    
    
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
    
    p = ComputePvalues(t, df, 't',   sProcess.options.tail.Value );
    df=ones(n_chan,1)*df;
    
    
    % Saving the output.
    iStudy = sInputs.iStudy;


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
    comment='T-test : ';
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
    
    sOutput.Comment=comment;
    sOutput = bst_history('add', sOutput, B.History, '');

    sOutput = bst_history('add', sOutput, 'ttest computation', comment);
    OutputFiles{1} = bst_process('GetNewFilename', fileparts(sInputs(1).FileName), 'pdata_ttest_matrix');
    save(OutputFiles{1}, '-struct', 'sOutput');
    db_add_data(iStudy, OutputFiles{1}, sOutput);

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

