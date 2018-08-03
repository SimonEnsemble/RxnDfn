function sample_variables(st::SpaceTime, discretization::Discretization, sample_time::Float64)
    # initialize u_sample and t_sample to be filled below with specific time steps for graphing purposes
    num_samples = ceil(Int, st.tf / sample_time) # determines size of u_sample and t_sample
    sample_freq = floor(Int, discretization.Nₜ / (num_samples - 1)) # don't count the 0 sample
    u_sample = zeros(Float64, discretization.Nₓ, num_samples) # nb_unknowns = size of u, num_samples, how many u's we will store
    t_sample = zeros(Float64, num_samples)
    sample_counter = 1 # for if loop below to identify column of u_sample, increments below after u_sample is imput

    return sample_freq, u_sample, t_sample, sample_counter
end
