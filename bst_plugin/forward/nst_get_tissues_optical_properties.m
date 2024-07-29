function  prop = nst_get_tissues_optical_properties(method,wavelength)
%==========================================================================
%                       INITIALISATION: SET DEFAULTS
%==========================================================================
p = inputParser;
addRequired(p,'method',@isscalar);
addRequired(p,'wavelength',@ischar);
parse(p,method,wavelength);


%==========================================================================
%                       GET OPTICAL PROPERTIES
%==========================================================================
switch method
    case 1
        % Strangman_2003: mus was calculated from mus', conversion in mm-1
        % NB: another useful references : Yamada 2009 JBO 14(6) has a set of absorption and scattering coeff at multiple wavelength for the five layers

        
        if strcmp(wavelength,'685')
            % 690 nm : [mua, mus, g, n]
            prop=[0       0       1    1      % medium 0: the environment
                0.0159  10      0.92 1.37   % medium 1: skin
                0.0101  12.5    0.92 1.37   % medium 2: skull
                0.0004  0.125   0.92 1.37   % medium 3: CSF
                0.0178  15.625  0.92 1.37   % medium 4: gray matter
                0.0178  15.625  0.92 1.37]; % medium 5: white matter
            
       elseif strcmp(wavelength,'690')
            % 690 nm : [mua, mus, g, n]
            prop=[0       0       1    1      % medium 0: the environment
                0.0159  10      0.92 1.37   % medium 1: skin
                0.0101  12.5    0.92 1.37   % medium 2: skull
                0.0004  0.125   0.92 1.37   % medium 3: CSF
                0.0178  15.625  0.92 1.37   % medium 4: gray matter
                0.0178  15.625  0.92 1.37]; % medium 5: white matter
            
        elseif strcmp(wavelength,'830')
            % 830 nm : [mua, mus, g, n]
            prop=[0      0      1     1       % medium 0: the environment
                0.0191 8.25   0.92  1.37    % medium 1: skin
                0.0136 10.75  0.92  1.37    % medium 2: skull
                0.0026 0.125   0.92 1.37    % medium 3: CSF
                0.0186 13.875 0.92  1.37    % medium 4: gray matter
                0.0186 13.875 0.92  1.37]; % medium 5: white matter
        end      

end