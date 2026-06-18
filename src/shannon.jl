# Q12.1 — Entropie de Shannon H(X) = -Σ P(x) log₂(P(x))
function shannon_entropy(data, n_bins, lo, hi)
    width = (hi - lo) / n_bins
    counts = zeros(Int, n_bins)
    for x in data
        i = floor(Int, (x - lo) / width) + 1
        i = clamp(i, 1, n_bins)
        counts[i] += 1
    end
    n = length(data)
    H = 0.0
    for c in counts
        c == 0 && continue
        p = c / n
        H -= p * log2(p)
    end
    return H
end

# Q12.1 — H(X,Y,Z,<v²>) = H(X) + H(Y) + H(Z) + H(<v²>)
function total_entropy(molecules, domain)
    xs  = [m.pos[1] for m in molecules]
    ys  = [m.pos[2] for m in molecules]
    zs  = [m.pos[3] for m in molecules]
    v2s = [sum(m.velocity .^ 2) for m in molecules]

    hx = -domain.Lx / 2; hX = domain.Lx / 2
    hy = -domain.Ly / 2; hY = domain.Ly / 2
    hz = -domain.Lz / 2; hZ = domain.Lz / 2

    H  = shannon_entropy(xs,  10, hx, hX)
    H += shannon_entropy(ys,  10, hy, hY)
    H += shannon_entropy(zs,  10, hz, hZ)
    H += shannon_entropy(v2s, 200, 0.0, 200_000.0)
    return H
end
