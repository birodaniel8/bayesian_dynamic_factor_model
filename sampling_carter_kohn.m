function x = sampling_carter_kohn(F, P, G, Q, j)
%SAMPLING_CARTER_KOHN - Carter & Kohn (1994) sampling algorithm for sampling Kalman Filtered states 
%
% Inputs:
%    F: (T x k) state estimations from Kalman Filter (updated values)
%    P: (k x k x T) state covariance from Kalman Filter (updated values)
%    G: (k x k) state transition matrix
%    Q: (k x k) process noise matrix
%    j: (int) size of the block for which the Q matrix is positive definite
%       default: size(Q,1) ie. whole matrix
%
% Outputs:
%    x: (T x k) sampled states

%------------- BEGIN CODE --------------

if nargin == 4, j = size(Q, 1); end

if size(G, 3) > 1 || size(Q, 3) > 1
    x = sampling_carter_kohn_with_time_varying_parameters(F, P, G, Q, j);
else
    F = F';
    x = [];
    Q_star = Q(1:j, 1:j);
    G_star = G(1:j, :);
    for s = size(F, 2):-1:1
        if s == size(F, 2)
            x(:, s) = mvnrnd(F(:, s), (P(:, :, s) + P(:, :, s).')/2, 1)';  % the covariance matrix is calculated as (A+A')/2 to ensure positive definiteness
        else
            F_star = F(:, s) + P(:, :, s) * G_star' * inv(G_star*P(:, :, s)*G_star'+Q_star) * (x(1:j, s+1) - G_star * F(:, s));
            P_star = P(:, :, s) - P(:, :, s) * G_star' * inv(G_star*P(:, :, s)*G_star'+Q_star) * G_star * P(:, :, s);
            P_star = (P_star + P_star.') / 2;
            x(:, s) = mvnrnd(F_star, P_star, 1)';
        end
    end
    x = x';
end

%------------- END OF CODE --------------
end