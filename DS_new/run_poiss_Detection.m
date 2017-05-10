function precision = run_poiss_Detection(c50, alpha1, phi, alpha2, r_base, ymax, c, repeat_n)

x(1) = 0;
x(2) = log10(c);
THETA = linspace(0, 2*pi, 50); % direction
for dir = 1:length(THETA)
    theta = THETA(dir);
    input = ymax.* x(2).^alpha1./(x(2).^alpha1 + c50.^alpha1);
    gain = (0.5 + 0.5 * cos(theta + phi)).^alpha2;
    modelOutput(dir,:) = input.*gain + r_base;
end

correct = zeros(repeat_n, length(THETA));

for dir = 1:length(THETA)
    for repeat = 1:repeat_n
        for j = 1:2
            theta = THETA(dir);
            input = ymax.* x(j).^alpha1./(x(j).^alpha1 + c50.^alpha1);
            gain = (0.5 + 0.5 * cos(theta + phi)).^alpha2;
            output = input.*gain + r_base;
            for i = 1:length(output)
                output(i) = random('poiss', output(i));
            end

            output_r = repmat(output, length(THETA), 1);
            p = r_base.^output.*exp(-r_base)./factorial(output);
            p_null = prod(p, 2);
            p = modelOutput.^output_r.*exp(-modelOutput)./factorial(output_r);
            p_stim = max(prod(p, 2));
            l(j) = p_stim/p_null;
        end
        if l(1) < l(2)
            correct(repeat, dir) = 1;
        elseif l(1) == l(2)
            correct(repeat, dir) = 0.5;
        end
    end
end
precision = sum(correct(:))/repeat_n/length(THETA);