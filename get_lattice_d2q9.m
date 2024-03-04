function d2q9=get_lattice_d2q9()
    d2q9.cs2 = 1/3;
    w0 = 16/36;
    w1 = 4/36;
    w2 = 1/36;
    d2q9.w = [w0 w1 w1 w1 w1 w2 w2 w2 w2];
    d2q9.vx = [0 1 0 -1 0 1 -1 -1 1];
    d2q9.vy = [0 0 1 0 -1 1 1 -1 -1];
    d2q9.opp = [0 3 4 1 2 7 8 5 6]+1;
end