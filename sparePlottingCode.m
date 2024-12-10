%% Spare plotting code

%         figure(1)
%         clf
%         hold on
%         for i = [1 2 3]%[1 3]
%             plot(x_c,v(i:eqn_ma:N*eqn_ma),'-o')
%         end
%         hold off
%         text(0,1,sprintf('t=%f',t))
% %         xlim([0.75*l_e l_e+l_sep+0.25*l_e])
%         xlabel('x')
%         ylabel('c')
%         legend('R','O')
%         drawnow
% %             
%         figure(2)
%         clf
%         hold on
%         for i = [1 2 3]%[1 3]
%             plot(x_c,v2(i:eqn_mi:N*eqn_mi),'-o')
%         end    
% %         for i = eqn_ma+1:eqn_mi
% %             plot(x_c,v2(i:eqn_mi:N*eqn_mi),'-o')
% %         end    
% %         plot(x_c,v2(index_HX:eqn_mi:N*eqn_mi)+v2(index_X:eqn_mi:N*eqn_mi),'-o')
% %         plot(x_c,v2(index_HY:eqn_mi:N*eqn_mi)+v2(index_Y:eqn_mi:N*eqn_mi),'-o')
%         hold off
%         xlabel('x')
%         ylabel('c_{mi}')
%         legend('R','Cl','O')
%         drawnow
% %     
%         figure(3)
%         clf
%         hold on
%         plot(x_c,v(eqn_ma:eqn_ma:N*eqn_ma),'-o')
%         plot(x_c,v2(eqn_ma:eqn_mi:N*eqn_mi)*F/V_T/C_S,'-o')
%         plot(x_c,phi_el-v(eqn_ma:eqn_ma:N*eqn_ma)-v2(eqn_ma:eqn_mi:N*eqn_mi)*F/V_T/C_S,'-o')
%         plot(x_c,phi_el,'-o')
%         hold off
%         xlabel('x')
%         ylabel('\phi')
%         legend('\phi','\phi_{St}','\phi_{Donnan}','\phi_{electrode}')
%         drawnow
% %     
%         figure(4)
%         clf
%         hold on
%         plot(x_c,v2(s_mo+1:eqn_mi:N*eqn_mi),'-o')
%         hold off
%         xlabel('x')
%         ylabel('\sigma_{ele}')
%         drawnow
    % 
%         figure(5)
%         clf
%         hold on
%         plot(x_f(2:end-1)',flux_1_E,'-o')
%         plot(x_f(2:end-1)',flux_1_D,'-o')
%         hold off
%         legend('ELE','Diff')
%         xlabel('x')
%         ylabel('Flux')
%         drawnow
    %     
%         figure(6)
%         clf
%         hold on
%         for i = index_H
%             plot(x_c,-log10(v(i:eqn_ma:N*eqn_ma)/1e3),'-o')
%         end
%         plot(x_c,-log10(v(index_H:eqn_ma:N*eqn_ma)/1e3) -log10(v(index_OH:eqn_ma:N*eqn_ma)/1e3),'-o')
%         hold off
%         text(3e-3,14,sprintf('t=%.2f',t))
%     %     xlim([0 L-2*l_res])
%     %     ylim([12 16])
%         xlabel('x')
%         ylabel('pH_{mA},(pH+pOH)_{mA}')
%         drawnow
%     
%         figure(7)
%         clf
%         hold on
%         for i = index_H
%             plot(x_c,-log10(v2(i:eqn_mi:N*eqn_mi)/1e3),'-o')
%         end
%         plot(x_c,-log10(v2(index_H:eqn_mi:N*eqn_mi)/1e3) -log10(v2(index_OH:eqn_mi:N*eqn_mi)/1e3),'-o')
% %         plot(x_c,-log10(v2(index_H:eqn_mi:N*eqn_mi).*v2(index_X:eqn_mi:N*eqn_mi)./v2(index_HX:eqn_mi:N*eqn_mi)/1e3),'-o')
% %         plot(x_c,-log10(v2(index_H:eqn_mi:N*eqn_mi).*v2(index_Y:eqn_mi:N*eqn_mi)./v2(index_HY:eqn_mi:N*eqn_mi)/1e3),'-o')
%         plot(x_c,-log10(v2(index_H:eqn_mi:N*eqn_mi).*v2(index_COO:eqn_mi:N*eqn_mi)./v2(index_COOH:eqn_mi:N*eqn_mi)/1e3),'-o')
% 
%         hold off
%         text(3e-3,14,sprintf('t=%.2f',t))
%         legend('pH','pH+pOH','pKa5','pKa9','pKa9.5')
%     %     ylim([12 16])
%         xlabel('x')
%         ylabel('pH & pKa')
%         drawnow

%         figure(8)
%         clf
%         hold on
%         plot(t_plotting_array,c_eff_array,'g')
%         xlabel('t')
%         ylabel('c_eff')
% 
%         figure(9)
%         clf
%         hold on
% 
%         plot(x_c,v2(index_COO:eqn_mi:end),'-o')
%         plot(x_c,v2(index_COO:eqn_mi:end)+v2(index_COOH:eqn_mi:end),'-o') % Total
% 
%         hold off
%         text(3e-3,14,sprintf('t=%.2f',t))
%         xlim([0 2*l_e+l_sep])
%         legend('COO','COOH+COO')
% %         ylim([12 16])
%         xlabel('x')
%         ylabel('Concentration[mM]')
%         drawnow








%         figure(1)
%         clf
%         hold on
%         summation = 0;
%         for i = 1:s_mo
%             summation = summation + z_mo_input(i).*v(i:eqn_ma:N*eqn_ma);
%         %             plot(x_c,v(i:eqn_ma:N*eqn_ma),'-o')
%         end
%         plot(x_c,summation,'-o')
%         hold off
%         text(3e-3,10,sprintf('t=%f',t))
%         %         xlim([0.75*l_e l_e+l_sep+0.25*l_e])
%         xlabel('x')
%         ylabel('c')
%         drawnow
%             
%         figure(2)
%         clf
%         hold on
%         for i = 1:s_mo
%             plot(x_c,v2(i:eqn_mi:N*eqn_mi),'-o')
%         end    
% %         for i = eqn_ma+1:eqn_mi
% %             plot(x_c,v2(i:eqn_mi:N*eqn_mi),'-o')
% %         end    
% %         plot(x_c,v2(index_HX:eqn_mi:N*eqn_mi)+v2(index_X:eqn_mi:N*eqn_mi),'-o')
% %         plot(x_c,v2(index_HY:eqn_mi:N*eqn_mi)+v2(index_Y:eqn_mi:N*eqn_mi),'-o')
% 
%         hold off
%         xlabel('x')
%         ylabel('c_{mi}')
%         drawnow
% % % 
%         figure(3)
%         clf
%         hold on
%         plot(x_c,v(eqn_ma:eqn_ma:N*eqn_ma),'-o')
%         plot(x_c,v2(eqn_ma:eqn_mi:N*eqn_mi)*F/V_T/C_S,'-o')
%         plot(x_c,phi_el-v(eqn_ma:eqn_ma:N*eqn_ma)-v2(eqn_ma:eqn_mi:N*eqn_mi)*F/V_T/C_S,'-o')
%         plot(x_c,phi_el,'-o')
%         hold off
%         xlabel('x')
%         ylabel('\phi')
%         legend('\phi','\phi_{St}','\phi_{Donnan}','\phi_{electrode}')
%         drawnow
% % 
%     figure(4)
%     clf
%     hold on
%     plot(x_c,v2(eqn_ma:eqn_mi:N*eqn_mi),'-o')
%     hold off
%     xlabel('x')
%     ylabel('\sigma_{ele}')
%     drawnow

%     figure(7)
%     clf
%     hold on
%     for i = 3
%         plot(x_c,-log10(v2(i:eqn_mi:N*eqn_mi)/1e3),'-o')
%     end
%     plot(x_c,-log10(v2(3:eqn_mi:N*eqn_mi)/1e3) -log10(v2(4:eqn_mi:N*eqn_mi)/1e3),'-o')
%     plot(x_c,-log10(v2(3:eqn_mi:N*eqn_mi).*v2(7:eqn_mi:N*eqn_mi)./v2(6:eqn_mi:N*eqn_mi)/1e3),'-o')
%     plot(x_c,-log10(v2(3:eqn_mi:N*eqn_mi).*v2(9:eqn_mi:N*eqn_mi)./v2(8:eqn_mi:N*eqn_mi)/1e3),'-o')
% 
%     hold off
%     text(3e-3,14,sprintf('t=%.2f',t))
%     legend('pH','pH+pOH','pKa5','pKa9')
% %     ylim([12 16])
%     xlabel('x')
%     ylabel('pH & pKa')
%         drawnow