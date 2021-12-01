function decomposicao(mp)
    #global j1=mp[1:4*ngen,1:4*ngen];
    #j2=mp[1:4*ngen,4*ngen+1:4*ngen+2*nbus];
    #j3=mp[ngen*4+1:ngen*4+2*nbus,1:4*ngen];
    #j4=mp[ngen*4+1:ngen*4+2*nbus,4*ngen+1:4*ngen+2*nbus];
    #b1=mp[1:4*ngen,4*ngen+2*nbus+1:ngen*4+nbus*2+ngen*2];
    #b2=mp[ngen*4+1:ngen*4+2*nbus,4*ngen+2*nbus+1:ngen*4+nbus*2+ngen*2];
    #A=j1-j2*inv(j4)*j3;
    #B=b1-j2*inv(j4)*b2;
    ng=ngen;
    nb=nbus;
    ns=ngen;
    wp1 = 4*ng;
    wp2 = wp1 + 2*nb;
    wp3 = wp2 #+ ns;

    wp4 = 4*ng;
    wp5 = wp4 + 2*nb;
    wp6 = wp5 + 2*ng;

    aux_11 = 1:wp1;
    aux_21 = (wp1+1):(wp2)
    aux_31 = (wp2+1):(wp3)

    aux_12 = 1:wp4
    aux_13 = (wp4+1):(wp5)
    aux_14 = (wp5+1):(wp6)
    
    j1 = mp[aux_11,aux_12]
    j2 = mp[aux_11,aux_13]
    b1 = mp[aux_11,aux_14]
    j3 = mp[aux_21,aux_12]
    j4 = mp[aux_21,aux_13]
    b2 = mp[aux_21,aux_14]
    c1 = mp[aux_31,aux_12]
    c2= mp[aux_31,aux_13]
    c3 = mp[aux_31,aux_14]
    A = j1 - j2*inv(j4)*j3;
    B = b1 - j2*inv(j4)*b2;
    C = c1 - c2*inv(j4)*j3;
    D = c3 - c2*inv(j4)*b2;
    return A,B,C,D;    
end