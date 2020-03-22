function [x, P] = kalman_filter_with_time_varying_parameters(z, H, R, G, Q, x0, P0)
%KALMAN_FILTER - Perform the Kalman filter estimation with time varying parameters
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

%------------- BEGIN CODE --------------

z = z';
N = size(z, 2);
x = zeros(size(H, 2), N);
P = zeros(size(H, 2), size(H, 2), N);

% Convert all input matrixes into 3D objects
if size(H, 3) == 1, H = repmat(H, 1, 1, N); end
if size(R, 3) == 1, R = repmat(R, 1, 1, N); end
if size(G, 3) == 1, G = repmat(G, 1, 1, N); end
if size(Q, 3) == 1, Q = repmat(Q, 1, 1, N); end

% Loop through and perform the Kalman filter equations recursively
for i = 1:N
    if i == 1
        x(:, i) = G(:, :, i) * x0;  % Predict the state vector from the initial values
        P(:, :, i) = G(:, :, i) * P0 * G(:, :, i)' + Q(:, :, i);  % Predict the covariance from the initial values
    else
        x(:, i) = G(:, :, i) * x(:, i-1);  % Predict the state vector
        P(:, :, i) = G(:, :, i) * P(:, :, i-1) * G(:, :, i)' + Q(:, :, i);  % Predict the covariance
    end
    if all(isnan(z(:, i)))
        % if all of them is missing at time 'i', there is not update and the updated values are the predicted ones
    elseif any(isnan(z(:, i)))
        keep = ~isnan(z(:, i));
        K = P(:, :, i) * H(keep, :, i)' / (H(keep, :, i) * P(:, :, i) * H(keep, :, i)' + R(keep, keep, i));  % Calculate the Kalman gain matrix only with the observed data
        x(:, i) = x(:, i) + K * (z(keep, i) - H(keep, :, i) * x(:, i));  % Update the state vector only with the observed data
        P(:, :, i) = (eye(size(P, 1)) - K * H(keep, :, i)) * P(:, :, i);  % Update the covariance only with the observed data
    else
        K = P(:, :, i) * H(:, :, i)' / (H(:, :, i) * P(:, :, i) * H(:, :, i)' + R(:, :, i));  % Calculate the Kalman gain matrix
        x(:, i) = x(:, i) + K * (z(:, i) - H(:, :, i) * x(:, i));  % Update the state vector
        P(:, :, i) = (eye(size(P, 1)) - K * H(:, :, i)) * P(:, :, i);  % Update the covariance
    end
end
x = x';

%------------- END OF CODE --------------
end