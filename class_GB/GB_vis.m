classdef GB_vis
    methods (Static)
        %%
        function plotDist(gb,k)
            if nargin < 2
                k = 3;
            end
            %% dist
            num_samp= gb.count_samp;
            vi      = visualize_data;
            collist = vi.genGrad(vi.para_col.col3(:,1),num_samp);
            
            list_mus    = gb.list_mus;
            list_nus    = gb.list_nus;
            list_sigmas = gb.list_sigmas;
            list_leg    = [];
            vi.figure;
            hold on;
            for i_samp = 1 : k : num_samp
                pd = makedist('tLocationScale',...
                    'mu',list_mus(i_samp),...
                    'sigma',list_sigmas(i_samp),...
                    'nu',max(list_nus(i_samp),1));
                x    = [5:20];
                pdf1 = pdf(pd,x);
                plot(x,pdf1,'Color',collist(:,i_samp),'LineWidth',1.5);
                list_leg = [list_leg,{sprintf('%d samples',i_samp)}];
            end
            vi.setFig(-4,10)
            vi.setLabel('Frequency [Hz]','Probability');
            drawnow; 
        end
        %%
        function plotRes(gb)
            vi = visualize_data;
            vi.fig;
            vi.sp(3,1,1);
            plot(gb.list_mus,'LineWidth',1.5,'Color',vi.para_col.col4(:,1));
            title('Average');
            ylabel('Frequency [Hz]');
            fcnRoutine;
            vi.sp(3,1,2);
            plot(gb.list_sigmas,'LineWidth',1.5,'Color',vi.para_col.col4(:,2));
            title('sigma');
            ylabel('Frequency [Hz]');
            fcnRoutine;
            vi.sp(3,1,3);
            plot(gb.list_nus,'LineWidth',1.5,'Color',vi.para_col.col4(:,3));
            title('nu');
            xlabel('Number of samples')
            fcnRoutine;
            
            function fcnRoutine
                vi.setFig(-4,8);
                yl = ylim;
                ylim([yl(1) yl(2) + 1]);
            end
            drawnow limitrate; 
        end
        %%
        function plotProb(gb)
            vi          = visualize_data;
            data_prob   = gb.getIAF;
            vi.figure;
            hold on;
            vi.sp(2,1,1);
            plot(data_prob(:,1),'Color',vi.para_col.col(:,1),'LineWidth',1.5);
            vi.setFig(-4,10)
            vi.setLabel('Number of samples','Probability')
            vi.sp(2,1,2);
            plot(data_prob(:,2),'Color',vi.para_col.col(:,1),'LineWidth',1.5);
            vi.setFig(-4,10)
            yl = ylim;
            ylim([yl(1) yl(2) + 1]);
            vi.setLabel('Number of samples','Peak Freq.')
            drawnow limitrate; 
        end
    end
end