% Fast Linear Quaternion Attitude Estimator from Vector Observations
% Published in IEEE Transactions on Automation Science and Engineering
% Authors: Jin Wu, Zebo Zhou et al.
% Copytight (c) 2016

function [Q, iter] = FLAE_newton( r_base, b_base, weights )

    MM = zeros(3, 3);
   
    for i = 1 : length(b_base(1, :))
        
        MM = MM + weights(i) * b_base(:, i) * r_base(:, i)';
    end
    
    Hx1 = MM(1, 1);    Hx2 = MM(1, 2);    Hx3 = MM(1, 3);
    Hy1 = MM(2, 1);    Hy2 = MM(2, 2);    Hy3 = MM(2, 3);
    Hz1 = MM(3, 1);    Hz2 = MM(3, 2);    Hz3 = MM(3, 3);
    
    W = [Hx1 + Hy2 + Hz3, -Hy3 + Hz2, -Hz1 + Hx3, -Hx2 + Hy1;
        -Hy3 + Hz2, Hx1 - Hy2 - Hz3, Hx2 + Hy1, Hx3 + Hz1;
        -Hz1 + Hx3, Hx2 + Hy1, Hy2 - Hx1 - Hz3, Hy3 + Hz2;
        -Hx2 + Hy1, Hx3 + Hz1, Hy3 + Hz2, Hz3 - Hy2 - Hx1];
   

   c = det(W);
   b = 8 * Hx3 * Hy2 * Hz1 - 8 * Hx2 * Hy3 * Hz1 - ...
       8 * Hx3 * Hy1 * Hz2 + 8 * Hx1 * Hy3 * Hz2 + ...
       8 * Hx2 * Hy1 * Hz3 - 8 * Hx1 * Hy2 * Hz3;
   a = -2 * Hx1 * Hx1 - 2 * Hx2 * Hx2 - 2 * Hx3 * Hx3 - ...
        2 * Hy1 * Hy1 - 2 * Hy2 * Hy2 - 2 * Hy3 * Hy3 - ...
        2 * Hz1 * Hz1 - 2 * Hz2 * Hz2 - 2 * Hz3 * Hz3;

   lambda = 1.0;
   old_lambda = 1.0;
   iter = 0;
 
   while(abs(old_lambda - lambda) > 1e-8 && iter <= 50)
       old_lambda = lambda;
       lambda = lambda - ((lambda^2 + a / 2)^2 + b * lambda + c - a^2 / 4) / ...
                          (4 * lambda * (lambda^2 + a / 2) + b);
       iter = iter + 1;
   end
   
   G = W - lambda * eye(4);
   
   
   
   pivot = G(1, 1);  
   G(1, :) = G(1, :) / pivot;
   G(2, :) = G(2, :) - G(2, 1) * G(1, :);
   G(3, :) = G(3, :) - G(3, 1) * G(1, :);
   G(4, :) = G(4, :) - G(4, 1) * G(1, :);

   
   pivot = G(2, 2);
   G(2, :) = G(2, :) / pivot;
   G(1, :) = G(1, :) - G(1, 2) * G(2, :);
   G(3, :) = G(3, :) - G(3, 2) * G(2, :);
   G(4, :) = G(4, :) - G(4, 2) * G(2, :);
   
   pivot = G(3, 3);
   G(3, :) = G(3, :) / pivot;
   G(1, :) = G(1, :) - G(1, 3) * G(3,:);
   G(2, :) = G(2, :) - G(2, 3) * G(3,:);
   G(4, :) = G(4, :) - G(4, 3) * G(3,:);
   
   q = [G(1, 4); G(2, 4); G(3, 4); -1];
 
   Q = q ./ norm(q);

end

