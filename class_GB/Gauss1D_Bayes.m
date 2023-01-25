classdef Gauss1D_Bayes < Atom
    %#ok<*PROPLC>
    properties
       para
       count_samp
       list_mus
       list_ramdas
       list_nus
       list_sigmas
    end
    
    methods (Access = public)
        
        function gb = learn(gb,in)
            num_samp = numel(in);
            for i_samp  = 1 : num_samp
                gb.para = gb.fcn_learn(gb.para,in(i_samp));
                gb      = gb.savePara;
            end
        end
        
        function gb = savePara(gb)
            gb.list_mus     = [gb.list_mus,gb.para.muS];
            gb.list_ramdas = [gb.list_ramdas,gb.para.ramdaS];
            gb.list_nus    = [gb.list_nus,gb.para.nuS];
            gb.count_samp  = gb.count_samp + numel(gb.para.nuS);
            
        end
       
        function gb = calc_list_sigmas(gb)
            gb.list_ramdas(gb.list_ramdas<=0) = 0.1;
            gb.list_sigmas =  sqrt(1./gb.list_ramdas);
        end
    end
    
    methods (Access = public)
        
        function gb = Gauss1D_Bayes
            %% ParameterInitialization
            gb = gb.initPara;            
        end
        
        function gb = initPara(gb)
            %initPara
            m       = gb.randomv(1);
            beta    = gb.randomv(-1);
            a       = gb.randomv(-1);
            b       = gb.randomv(1);
            mH      = gb.randomv(1);
            betaH   = gb.randomv(1);
            aH      = gb.randomv(1);
            bH      = gb.randomv(1);
            
            muS         = m;
            ramdaS      = (beta*a)/((1+beta)*b);
            nuS         = 2*a;
                        
            list_var = who;
            list_var = gb.rmlist(list_var,'gb');            
            gb.para  = gb.generateStruct(list_var,2);  
            
            gb.count_samp  = 0;            
            gb.list_mus    = [];
            gb.list_ramdas = [];
            gb.list_nus    = [];
        end
        
    end
    
    methods (Access = public) %% post learning
        
        function pdf1    = getPDF(gb,i_samp,x)
            if nargin < 3
                x    = gb.getX;
            end
            list_mus    = gb.list_mus;
            list_nus    = gb.list_nus;
            gb          = gb.calc_list_sigmas;
            list_sigmas = gb.list_sigmas;
            samp_end    = numel(list_sigmas);
            if nargin < 2
                i_samp = samp_end;
            elseif i_samp > gb.count_samp
                i_samp = samp_end;
            end
            pd = makedist('tLocationScale',...
                'mu',list_mus(i_samp),...
                'sigma',list_sigmas(i_samp),...
                'nu',max(list_nus(i_samp),1));
            pdf1 = pdf(pd,x);
        end
        
        function [IAF,prob] = decideIAF(gb,pdf1,x)
            if nargin < 3
                x       = gb.getX;
            end
            [m,I]   = max(pdf1);
            prob    = m;
            IAF     = x(I);
        end
        
        function data_out = getIAF(gb)
            num_samp = gb.count_samp;
            data_out = zeros(num_samp,2);
            for i_samp = 1 : num_samp
                pdf1 = gb.getPDF(i_samp);
                [data_out(i_samp,1),data_out(i_samp,2)] = ...
                    decideIAF(gb,pdf1);
            end
        end
    end
    methods (Static)
        function x = getX
            x = [5:20];
        end        
        
        function out = randomv(in)
            if in < 0
                out = abs(randn(abs(in)));
            else
                out = randn(in);
            end
        end
        
        function para = fcn_learn(para,in)
            N = numel(in);
            
            para.betaH = N + para.beta;
            para.mH    = 1/(para.betaH) * (sum(in)+para.beta*para.m);
            para.aH    = N/2 + para.a;
            para.bH    = 1/2*(sum(in.^2)+para.beta*para.m^2-para.betaH*para.mH^2)+para.b;
            
            para.muS    = para.m;
            para.ramdaS = (para.beta*para.a)/((1+para.beta)*para.b);
            para.nuS    = 2 * para.a;
            
            %%% update
            para.beta = para.betaH;
            para.m    = para.mH;
            para.a    = para.aH;
            para.b    = para.bH;
        end
        
    end
end