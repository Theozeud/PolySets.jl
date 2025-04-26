function Plots.plot(ps::PolysSet; 
                    xlim=(-1.0, 1.0), 
                    npoints=200, 
                    labels=:auto, 
                    kwargs...)
    x = range(xlim[1], xlim[2]; length=npoints)
    y = evaluate(ps, x)
    plt = plot(; kwargs...)
    for i in 1:npolys(ps)
        label = labels === :auto ? "poly $i" : labels[i]
        @views yi = y[i,:]
        plot!(plt, x, yi, label=label; kwargs...)
    end
    return plt
end


function Plots.plot!(plt, ps::PolysSet; xlim=(-1.0, 1.0), npoints=200, labels=:auto, kwargs...)
    x = range(xlim[1], xlim[2]; length=npoints)
    y = evaluate(ps, x)
    for i in 1:npolys(ps)
        label = labels === :auto ? "poly $i" : labels[i]
        @views yi = y[i,:]
        plot!(plt, x, yi , label=label; kwargs...)
    end
    return plt
end