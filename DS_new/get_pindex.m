function index = get_pindex(tuning)

theta = linspace(0, 2*pi, length(tuning) + 1);
theta = theta(1:end-1);
[X, Y] = pol2cart(theta, tuning);
[angle, ~] = cart2pol(sum(X), sum(Y));
[~, index] = max(cos(theta - angle));

end
