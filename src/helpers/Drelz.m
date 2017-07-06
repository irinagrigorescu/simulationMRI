function D = Drelz(t, T1, M0)
% % % % IRINA GRIGORESCU
% % % % 
% % % % INPUT:
% % % %     t = time (ms)
% % % %     T1 = longitudinal relaxation (ms)
% % % %     M0 = proton density
% % % % OUTPUT:
% % % %     D  = [      0        ;
% % % %                 0        ;
% % % %            M0 * (1-E1) ] ;


E1 = exp(-t/T1);

D  = [      0        ;
            0        ;
       M0 * (1-E1) ] ;
             
end