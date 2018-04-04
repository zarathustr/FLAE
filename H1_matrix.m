function H = H1_matrix(h)
    Hx1 = h(1);  Hy1 = h(2);  Hz1 = h(3);
    H = [
        Hx1,   0, -Hz1, Hy1;
          0, Hx1,  Hy1, Hz1;
       -Hz1, Hy1, -Hx1,   0;
        Hy1, Hz1,    0,-Hx1];
end