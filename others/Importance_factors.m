            stds=Su;
            if HYPLAS_flag
            else
%                 designRV=RandVarProp(1,:)==1;
%                 componts={rv_data_comp{designRV}};
%                 stds(designRV)=Su(designRV).*arq.x0([componts{:}]);
            end

            en = v_out0{4};

            den = dv_out0{4};
            dMen = dM_out{4};
            dSen = dS_out{4};
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % STRESS IMPORTANCE FACTOR COMPUTATION
            % ImpF = (dF/dv*Sv)/F
            i_out=2;
            [VM,imax] = max(abs(M_out{i_out}));
            en = v_out0{i_out}(imax);
            Men = M_out{i_out}(imax);
            Sen = S_out{i_out}(imax);
            max(abs(M_out{i_out})+3*S_out{i_out})
            
            den = dv_out0{i_out}(:,imax);
            dMen = dM_out{i_out}(:,imax);
            dSen = dS_out{i_out}(:,imax);
            v0_Mv_Sv=[en,Men,Sen]
            variacao=[den.*stds',dMen.*stds',dSen.*stds']
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            DeltEnergy=[dMen'; dSen'].*([1/Men;1/Sen]*stds);
            DeltEnergy=[den';dMen'; dSen'].*([1;1;1]*stds);
            figure,bar(abs(DeltEnergy'))
            legend('Deterministic','Mean','S.D.')
            title(['R.V. importance for \it{x} = ' num2str(xopt,'%0.5g ')])
            xlabel('Random Variables')
            ylabel('Output Variation (%)')
            disp('Test of Random Var Importance done!....')
