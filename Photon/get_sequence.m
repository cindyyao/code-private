function seq = get_sequence(steps)

if mod(steps, 2) == 1
    a(1, :) = ceil(steps/2):steps;
    a(2, :) = ceil(steps/2):-1:1;
    seq = reshape(a, 1, steps+1);
    seq = seq(2:end);
else
    a(1, :) = steps/2+1:steps;
    a(2, :) = steps/2:-1:1;
    seq = reshape(a, 1, steps);
end
