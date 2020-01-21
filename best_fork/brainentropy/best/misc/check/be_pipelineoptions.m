function OPT = be_pipelineoptions(OPT, pipeline)

if nargin==1
    pipeline = OPT.mandatory.pipeline;
end

switch pipeline
    
    case 'cMEM'

        DEF = be_cmem_pipelineoptions();   
        
    case 'wMEM'
        
        DEF = be_wmem_pipelineoptions();
        
    case 'rMEM'

        DEF = be_rmem_pipelineoptions();
        
end

OPT = be_struct_copy_fields(OPT, DEF, [], 0);
OPT.mandatory.pipeline = pipeline;

return

