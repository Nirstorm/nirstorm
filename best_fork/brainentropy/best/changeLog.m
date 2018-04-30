%% BEST CHANGE LOG %%

% 05 - 12 - 2014
- be_install 	            :   automatic update added

% 27 - 10 - 2014 	    :	Changes requested by Francois
- panel_brainentropy 	    :	remove use of caller stack
- process_inverse_mem	    : 	change position in menu liste from 325 to 328 (after the separator)
				change the process name to "Compute sources: BEst"

% 10 - 10 - 2014 	    :   New version of the packge 	
- be_pipelineoptions        :	permute DEF and OPT in the call to be_struct_copy_fields 		
- be_main		    :	be_main reserved for plugin, replaced by be_main_call
	
% 18 - 09 - 2014
- be_check_caller           :   bug correction - now brainstorm mode is set if any caller is a fcn within BrainstormHomeDir
- be_check_timedef          :   bug correction - OPTIONS.optional.minW replaced by minW
- be_look_for_brainentropy  :   added a line to move the MEM process to the right spot in brainstorm

% 03 - 09 - 2014            
- process_inverse_mem       :   creation
- be_main                   :   detect if called from installation or not  
                                (avoid infinite loop)
                                
% 11 - 08 - 2014
- be_ridgefilter            :   fix bug when called in standalone - select the right channels

% 09 - 08 - 2014
- panel_brainentropy  	    : 	fix bug when multiple files are selected - doesnt look for default 					baseline
                            :	fix bug when the field OPTIONS.mandatory.version is missing

% 06 - 08 - 2014
- be_cwavelet               :   fix bug when no signal padding is necessary - unset variable
- be_cwsparse               :   fix bug when no signal padding is necessary - unset variable

% 29 - 07 - 2014 
- be_main                   :	fix bug on stand alone version: field MSP_min_window had no default value

% 28 - 05 - 2014 
- panel_brainentropy        :   limit panel size at first instance
- panel_brainentropy        :   avoid deprecated java calls warnings
- be_cwavelet               :   replace built-in padarray.m
- be_cwsparse               :   replace built-in padarray.m
- be_bpfilter               :   check if sig proc toolbox is installed. If not, filter with bst 				functions

% 26 - 05 - 2014
- panel_brainentropy        :   fixed bug with rMEM selection
- be_look_for_wavelab       :   no more automatic download - print instructions for manual download
- be_check_timedef          :   add the FLAG output argument

% 24 - 05 - 2014
- panel_brainentropy        :   restrict selection range to [4 5] for wMEM noisecov method
- panel_brainentropy        :   restrict selection range to [1 4] for cMEM noisecov method
- panel_brainentropy        :   restrict selection range to [1 5] for rMEM noisecov method
- panel_brainentropy        :   emptyroom option only available when expert mode is on
- panel_brainentropy        :   fix bug with frequency range definition in rMEM
- panel_brainentropy        :   show references at first panel instance
- be_wMEM_PIPELINEOPTIONS   :   changed default noise cov method for wMEM to 5
- be_main                   :   new label in TF files indicating to use be_vizr.m to display
- be_main                   :   set default display to 0
- be_vizr                   :   text now appears in white
- be_cmem_solver            :   comment of the result file now has prefix cMEM
- be_rmem_solver            :   comment of the result file now has prefix rMEM


% 20 - 03 -2014
- be_look_for_wavelab       :   corrected small bug, did not add the correct 
                                folder to matlab paths (line 40)
- be_wmem_pipelineoptions   :   ajout de l'option par defaut DEF.solver.NoiseCov_method
                                afin d'initialiser la matrice de covariance 
                                des donnees avec l'echelle 1 des coefficients
                                temps-echelle
- panel_brainentropy        :   correction du champ scales du pipeline wMEM
                                ne perd plus seon contenu lorsqu'on switch en 
                                mode expert

% 13 - 03 - 2014
- be_selected_coeffs        :   remplacer ''(OPTIONS.automatic.selected_samples(6,:)-1)/fs'' 
                                par ''(OPTIONS.automatic.selected_samples		 		 			(2,:)-1)fs/2'' (lignes 117 et 118)'' 
                                J-M please double-check
- panel_brainentropy        :   ajout de l'option emptyroom

% 07 - 03 - 2014
- be_wavelet_inverse.m 		:	remplacer "iwt_po" par "IWT_PO" (ligne 48)
- be_discrete_wavelet_transform :	remplacer "fwt_po" par "FWT_PO" (ligne 44)
- be_struct_copy_fields		: 	ajouter le champ "varargin" dans l'appel � 
                                la fonction "be_struct_copy_fields". Ceci   
                                permet de sp�cifier si les options choisies 
                                dans l'appel � la fonction doivent "override" 
                                les options par d�faut (ce qui n'�tait pas 
                                le cas avant)
- be_pipelineoptions		: 	changer l'appel "be_struct_copyfields()" en 
                                ajoutant le param�tre 1 � la fin (override)     
                                (ligne 23)

 
