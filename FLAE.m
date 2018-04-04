% Fast Linear Quaternion Attitude Estimator from Vector Observations
% Published in IEEE Transactions on Automation Science and Engineering
% Authors: Jin Wu, Zebo Zhou et al.
% Copytight (c) 2016

function Q = FLAE(r_base, b_base, weights, method)

    MM = zeros(3, 3);
   
    for i = 1 : length(b_base(1, :))
        MM = MM + weights(i) * b_base(:, i) * r_base(:, i)';
    end

    hx = MM(:, 1)';
    hy = MM(:, 2)';
    hz = MM(:, 3)';
    WW = H1_matrix(hx) + H2_matrix(hy) + H3_matrix(hz);

    if(strcmp(method, 'symbolic'))
        
        c = det(WW);
        b = - 8 * det(MM');
        a = - 2 * trace(MM * MM');

        T0 = 2 * a^3 + 27 * b^2 - 72 * a * c;
        T1 = (T0 + sqrt(- 4 * (a^2 + 12 * c)^3 + T0^2))^(1 / 3);
        T2 = sqrt(- 4 * a + 2^(4 / 3) * (a^2 + 12 * c) / T1 + 2^(2 / 3) * T1);
   
        lambda1 =   0.20412414523193150818310700622549 * (T2 - sqrt(-T2^2 - 12 * a - 12 * 2.4494897427831780981972840747059 * b / T2));
        lambda2 =   0.20412414523193150818310700622549 * (T2 + sqrt(-T2^2 - 12 * a - 12 * 2.4494897427831780981972840747059 * b / T2));
        lambda3 = - 0.20412414523193150818310700622549 * (T2 + sqrt(-T2^2 - 12 * a + 12 * 2.4494897427831780981972840747059 * b / T2));
        lambda4 = - 0.20412414523193150818310700622549 * (T2 - sqrt(-T2^2 - 12 * a + 12 * 2.4494897427831780981972840747059 * b / T2));
    
        lambda = max(real([lambda1 lambda2 lambda3 lambda4]));
    
        G = WW - lambda * eye(4);

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
    elseif(strcmp(method, 'eig'))
        [V, D] = eig(WW);
        q = V(:, 4);
        Q = q ./ norm(q);
    elseif(strcmp(method, 'newton'))
        
        c = det(WW);
        b = - 8 * det(MM');
        a = - 2 * trace(MM * MM');
        
        lambda = 1.0;
        old_lambda = 1.0;
        iter = 0;
 
        while(abs(old_lambda - lambda) > 1e-8 && iter <= 50)
            old_lambda = lambda;
            lambda = lambda - ((lambda^2 + a / 2)^2 + b * lambda + c - a^2 / 4) / ...
                          (4 * lambda * (lambda^2 + a / 2) + b);
            iter = iter + 1;
        end
        
        G = WW - lambda * eye(4);

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
end

