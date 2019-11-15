function dy = cdiff(y, h)
%CDIFF Calculates numerical derivative
%
%   cdiff(y, h) calculates numerical derivative of $y = f(x)$ using
%   the two-point formula
%   \[
%      \frac{df}{dx} = \frac{f(x + h) - f(x - h)}{2h}.
%   \]
%   where $h$ is a small, positive number.

n = length(y);
dy = zeros(1, n);
for i = 2:n-1
   dy(i) = (y(i+1) - y(i-1)) / (2 * h);
end
dy(end) = [];
dy(1) = [];
end
