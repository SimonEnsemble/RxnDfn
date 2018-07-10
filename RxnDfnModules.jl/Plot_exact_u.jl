#Creates a series of plots to be made into a gif
using Plots

function gif_maker(t, x, u, exact_u)
    Nₜ = ceil(Int, st.tf / sample_time)
    @gif for i = 1:Nₜ
        scatter(x,u[:,i], ylim=(1.05*minimum(u), 1.05*maximum(u)), label="numerical", xlabel="x", ylabel="u", title=@sprintf("t = %.2f", t[i]))
        plot!(x, exact_u.(x,t[i]), label="exact")
    end
end
