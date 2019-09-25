classdef ArFitTest < matlab.unittest.TestCase
    
    properties
        sigma
        N
    end
    
    methods(TestMethodSetup)
        function setup(testCase)
            testCase.sigma =  sqrt(0.2);
            testCase.N =  1024;
        end
    end    
    
        methods(Test)
            
            function test_AR1(testCase)
                
                ar_coeffs=linspace(0.02,0.99,20);
                
                for coeff=ar_coeffs
                    A=[1 coeff];
                    x= testCase.sigma*randn(1,testCase.N);
                    y = filter(1,A,x); 

                    W=nst_math_fit_AR(y,1);
                    assert( all( abs(W-A) < 0.1))
                end
                
            end
             
            function test_ARP(testCase)
                                
                for order = 1:7
                    for ntrials=1:100
                        l=cumsum(rand(1,order+1)) ;
                        l=l/l(end);

                        A= l(end:-1:1);

                        x= testCase.sigma*randn(1,testCase.N*5);
                        y = filter(1,A,x); 

                        W=nst_math_fit_AR(y,order);
                        assert( all( abs(W-A) < 0.1),'Failed to estimate order %s',num2str(order))
                    end
                end
                
            end

         end
            
end