function degrau(A,B)
    D = [0 0 0 0 0 0;0 0 0 0 0 0; 0 0 0 0 0 0];
    C = [0 0 0 1 0 0 0 0 0 0 0 0;0 0 0 0 1 0 0 0 0 0 0 0;0 0 0 0 0 1 0 0 0 0 0 0]; 

    sys=ss(A,B,C,D);
    sys=minreal(tf(sys));
    feedback(sys[1,2],sys2)
    y,t,x=step(sys);
    return y,t
end