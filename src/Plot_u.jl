#Creates a series of plots to be made into a gif
using Plots
using Printf

function gif_maker(t, x, u, st, sample_time)
    Nₜ = ceil(Int, st.tf / sample_time)
    @gif for i = 1:Nₜ
        plot(x,u[:,i], ylim=(1.05*minimum(u), 1.05*maximum(u)), legend=false, xlabel="x", ylabel="u", title=@sprintf("t = %.2f", t[i]))
    end
end
