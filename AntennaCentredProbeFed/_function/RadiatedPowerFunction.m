function Pr_norm = RadiatedPowerFunction(n, a3e0mn, a3o0mn, b3e0mn, b3o0mn)
    double precision; format long;
    %m = 0;

    Pr_norm = ((1:n) .* ((1:n) + 1) ./ (2 * (1:n) + 1)) * (a3e0mn .* conj(a3e0mn) + a3o0mn .* conj(a3o0mn) + b3e0mn .* conj(b3e0mn) + b3o0mn .* conj(b3o0mn));

end
