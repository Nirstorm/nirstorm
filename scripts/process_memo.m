function varargout = process_memo( varargin )
eval(macro_method);
end


function sProcess = GetDescription() %#ok<DEFNU>

sProcess.options.text.Comment = 'Text: ';
sProcess.options.text.Type    = 'text';
sProcess.options.text.Value   = 'the text';

sProcess.options.flag.Comment = 'Check box';
sProcess.options.flag.Type    = 'checkbox';
sProcess.options.flag.Value   =  0;

sProcess.options.float.Comment = 'Value with unit';
sProcess.options.float.Type    = 'value';
sProcess.options.float.Value   = {50, 'unit', 0}; % last if nb of subunit digits, use 0 for integer

sProcess.options.choice.Comment = 'Choice';
sProcess.options.choice.Type    = 'combobox';
choices = choice_enum();
sProcess.options.choice.Value   = {choices.Choice1,...
                                   fieldnames(choices)};

end

%% ===== RUN =====
function OutputFile = Run(sProcess, sInputs) %#ok<DEFNU>
   OutputFile = {};
    
   %% Retrieve options
   fag = sProcess.options.flag.Value;
   
   float_val = sProcess.options.float.Value{1};
   
   choices = choice_enum();
   choice_int = sProcess.options.choice.Value{1};
   choice_str = sProcess.options.choice.Value{2}{sProcess.options.choice.Value{1}};

   switch choice_int
       case choices.Choice1
       case choices.Choice2
       otherwise
           error('Invalid choice value');
   end
end

function enum = choice_enum()
enum.Choice1 = 1;
enum.Choice2 = 2;
end
