function H = H3_matrix(h)
    Hx3 = h(1);  Hy3 = h(2);  Hz3 = h(3);
    H = [
        Hz3, -Hy3,  Hx3,    0;
       -Hy3, -Hz3,    0,  Hx3;
        Hx3,    0, -Hz3,  Hy3;
          0,  Hx3,  Hy3,  Hz3];
end