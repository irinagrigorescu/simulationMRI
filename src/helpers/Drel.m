function D = Drel(t, T1, T2)

E1 = exp(-t/T1);
E2 = exp(-t/T2);

D  = [  E2    0    0   ;
         0   E2    0   ;
         0    0    E1 ] ;
             
end