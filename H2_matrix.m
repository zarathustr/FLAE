function H = H2_matrix(h)
    Hx2 = h(1);  Hy2 = h(2);  Hz2 = h(3);
    H = [
        Hy2,  Hz2,    0, -Hx2;
        Hz2, -Hy2,  Hx2,    0;
          0,  Hx2,  Hy2,  Hz2;
       -Hx2,    0,  Hz2, -Hy2];
end