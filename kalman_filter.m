function [x, P] = kalman_filter(z, H, R, G, Q, x0, P0)
%KALMAN_FILTER - Perform the Kalman filter estimation
%
% Inputs:
%    z: (T x m) matrix of observations
%    H: (m x k x (T)) observation model matrix
%    R: (m x m x (T)) masurement noise covariance
%    G: (k x k x (T)) state transition matrix
%    Q: (k x k x (T)) process noise matrix
%    x0: (k x 1) initial state values 
%    P0: (k x k) initial state covariance matrix
%
% Outputs:
%    x: (T x k) estimated states
%    P: (k x k x T) estimated state covariance matrixes
% 

% Note:
%    In any of the input matrixes (H, R, G, Q) is more then 3D, then a time varying parameter Kalman Filter is called

%------------- BEGIN CODE --------------

if size(H, 3) > 1 || size(R, 3) > 1 || size(G, 3) > 1 || size(Q, 3) > 1
    [x, P] = kalman_filter_with_time_varying_parameters(z, H, R, G, Q, x0, P0);
else
    z = z';
    N = size(z, 2);
    x = zeros(size(H, 2), N);
    P = zeros(size(H, 2), size(H, 2), N);
    
    % Loop through and perform the Kalman filter equations recursively
    for i = 1:N
        if i == 1
            x(:, i) = G * x0;  % Predict the state vector from the initial values
            P(:, :, i) = G * P0 * G' + Q;  % Predict the covariance from the initial values
        else          
            x(:, i) = G * x(:, i-1);  % Predict the state vector
            P(:, :, i) = G * P(:, :, i-1) * G' + Q;  % Predict the covariance
        end
        if all(isnan(z(:, i)))
            % if all of them is missing at time 'i', there is not update and the updated values are the predicted ones
        elseif any(isnan(z(:, i)))
            keep = ~isnan(z(:, i));
            K = P(:, :, i) * H(keep, :)' / (H(keep, :) * P(:, :, i) * H(keep, :)' + R(keep,keep));  % Calculate the Kalman gain matrix only with the observed data
            x(:, i) = x(:, i) + K * (z(keep, i) - H(keep, :) * x(:, i));  % Update the state vector only with the observed data
            P(:, :, i) = (eye(size(P, 1)) - K * H(keep, :)) * P(:, :, i);  % Update the covariance only with the observed data
        else
            K = P(:, :, i) * H' / (H * P(:, :, i) * H' + R);  % Calculate the Kalman gain matrix
            x(:, i) = x(:, i) + K * (z(:, i) - H * x(:, i));  % Update the state vector
            P(:, :, i) = (eye(size(P, 1)) - K * H) * P(:, :, i);  % Update the covariance
        end
    end
    x = x';
end

%------------- END OF CODE --------------
end