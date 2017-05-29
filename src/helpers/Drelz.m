function D = Drelz(t, T1, M0)

E1 = exp(-t/T1);

D  = [      0        ;
            0        ;
       M0 * (1-E1) ] ;
             
end