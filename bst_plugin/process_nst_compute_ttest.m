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
    sProcess.Comment     = 'Compute ttest';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'wip - GLM';
    sProcess.Index       = 1402;
    sProcess.isSeparator = 0;
    sProcess.Description = 'http://neuroimage.usc.edu/brainstorm/Tutorials/NIRSFingerTapping#Movement_correction';
    % todo add a new tutorials
    
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'matrix'};
    sProcess.OutputTypes = {'data','raw'};
    
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
   
    
    sProcess.options.Contrast.Comment = 'Contrast vector';
    sProcess.options.Contrast.Type    = 'text';
    sProcess.options.Contrast.Value   = '[-1 1]';
    
    sProcess.options.Student.Comment={' One-tailed or two-tailed hypothesis? <br /> One tail','two tail'};
    sProcess.options.Student.Type ='radio';
    sProcess.options.Student.Value=1;
    sProcess.options.Student.Hidden  = 0;

    
    
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) 
    Comment = sProcess.Comment;
    Comment = [ Comment ' C = ' sProcess.options.Contrast.Value ]; 
end

function OutputFiles = Run(sProcess, sInputs)
    OutputFiles={};
    
    % parse input : 
    for i=1:length(sInputs) 
        name= strsplit(sInputs(i).Comment,' ');
        if( strcmp(name(1), 'covB') == 1)
            covB=in_bst_data(sInputs(i).FileName);
            
            name= strsplit(cell2mat(name(end)),'=');
            name= strsplit(cell2mat(name(end)),'_');

            df=str2num(cell2mat(name(1)));
        elseif ( strcmp(name(1), 'B') == 1)
            B=in_bst_data(sInputs(i).FileName);
        else
           bst_report('Error', sProcess, sInputs, [ 'The file ' sInputs(i).FileName ' is not recognized. Please only input the B and covB matrix' ]);
        end    
    end
    n_chan=size(B.Value',2);
    C=[1 -1 0]; 
    disp( int2str(df));
    
    B=C*B.Value';
    t=zeros(1,n_chan);
    
    for i = 1:n_chan
        t(i)= B(i) / sqrt( C*covB.Value(:,:,i)*transpose(C) ) ; 
    end
    
    if(  sProcess.options.Student.Value == 1)
        p=tcdf(-abs(t), df);
    else
        p=2*tcdf(-abs(t), df);
    end    
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
    sOutput.Time         = 1;
    sOutput.ColormapType = 'stat2';
    sOutput.DisplayUnits = 'F';
%     sOutput.nComponents  = StatA.nComponents;
%     sOutput.GridAtlas    = StatA.GridAtlas;
%     sOutput.Freqs        = StatA.Freqs;
%     sOutput.TFmask       = StatA.TFmask;

    sOutput = bst_history('add', sOutput, 'ttest computation', FormatComment(sProcess));
    OutputFiles{1} = bst_process('GetNewFilename', fileparts(sInputs(1).FileName), 'pdata_ttest_matrix');
    save(OutputFiles{1}, '-struct', 'sOutput');
    db_add_data(iStudy, OutputFiles{1}, sOutput);
    
    
    
    
    
end
