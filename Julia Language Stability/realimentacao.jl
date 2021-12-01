function realimentacao(A,B,C,D)
    ngen=3;
    npss=3;
    C = [0 0 0 1 0 0 0 0 0 0 0 0;0 0 0 0 1 0 0 0 0 0 0 0;0 0 0 0 0 1 0 0 0 0 0 0]; 
    
    input = [1,2,3]
    output = [1,2,3]

    TW  = [3,3,3];
    NB  = [3,3,3];
    
    MAT_KP = Diagonal([14.17; 34.72; 31.19]);
    MAT_T1 = Diagonal([1.05; 0.23; 0.15]);
    MAT_T2 = Diagonal([0.59; 0.03; 0.02]);
    sys=ss(A,B,C,D);
    sys=tf(sys);
    for in01 = 1:ngen
        out01 = in01;
        wh_in01 = input[in01];
        wh_out01 = output[out01];
        tw=TW[in01];
        nl=NB[in01];
        kp=MAT_KP[in01,out01];
        t1=MAT_T1[in01,out01];
        t2=MAT_T2[in01,out01];

        num = [tw,0];
        den = [tw,1];

        num1 = [t1,1];
        den1 = [t2,1];
        
        
        

        for kk= 1:nl
            num=conv(num1,num);
            den=conv(den1,den);
        end

        numh = num*kp;
        denh = den;
        
        sys_h = tf(numh,denh);

        sys1=sys[wh_out01,wh_in01];
    
        f=feedback(sys1,sys_h);
        step(f);       
              
            
    
           

           
     end
    
    
end