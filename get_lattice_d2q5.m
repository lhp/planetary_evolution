function d2q5=get_lattice_d2q5()
    w0 = 1/3;
    wn = 1/6;
    d2q5.w = [w0 wn wn wn wn];
    d2q5.vx = [0 1 0 -1 0];
    d2q5.vy = [0 0 1 0 -1];
    d2q5.cs2 = 1/3;
end