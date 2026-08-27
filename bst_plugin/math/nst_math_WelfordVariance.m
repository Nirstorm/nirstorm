classdef nst_math_WelfordVariance < handle
    properties (SetAccess = private)
        Count = 0      % Total number of matrix samples processed
        N     = 0      % Expected total number of sample
        Mean           % Matrix of running means [M x N]
    end
    
    properties (Access = private)
        M2             % Matrix of running sum of squares [M x N]
    end
    
    properties (Dependent)
        Variance       % Matrix of sample variances [M x N]
        StdDev         % Matrix of standard deviations [M x N]
    end
    
    methods

        % Constructor: initialize matrix
        function obj = nst_math_WelfordVariance(N, matrixSize)

            obj.N = N;
            
            % Pre-allocate state matrices
            obj.Mean = zeros(matrixSize);
            obj.M2 = zeros(matrixSize);
            
        end


        % Update the estimator with exactly ONE new 2D matrix sample
        function update(obj, X)
            % Ensure input is a valid 2D sample
            if isempty(X) || ~ismatrix(X)
                error('Input must be a single non-empty 2D matrix.');
            end
            
            obj.Count = obj.Count + 1;
            
            % Element-wise Welford's algorithm (fully vectorized across the 2D grid)
            delta       = X - obj.Mean;
            obj.Mean    = obj.Mean + delta / obj.Count;
            delta2      = X - obj.Mean;
            
            % Update M2
            obj.M2 = obj.M2 + delta .* delta2;
        end
        
        % Getter for Variance matrix
        function varVal = get.Variance(obj)
            if obj.Count > 1
                varVal = obj.M2 / (obj.Count - 1);
            else
                varVal = zeros(size(obj.Mean));
            end
        end
        
        % Getter for Standard Deviation matrix
        function stdVal = get.StdDev(obj)
            stdVal = sqrt(obj.Variance);
        end
    end
end
