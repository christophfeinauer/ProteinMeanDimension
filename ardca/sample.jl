using ArDCA
using HDF5
using JLD
using ProgressMeter
using Base.Threads: @threads, threadid


function add_pseudocount!(arnet, q=21, pc=0.01)
    arnet.p0[:] = (1-pc) * arnet.p0[:] .+ pc*(1/q);
    return arnet
end


function load_ardca(ardca_file)
    arnet = load(ardca_file)["arnet"]
    return arnet
end



function create_samples(ardca_file::String, out_file::String, nsamples::Int; pc::Float64=0.01, overwrite=false)

    if isfile(out_file) && !overwrite
        error("file $out_file exists and overwrite is false")
    end

    samples, logp = create_samples(ardca_file, nsamples; pc=pc)

    if isfile(out_file)
        rm(out_file)
    end

    # transpose samples for python
    samples = permutedims(samples, (3, 2, 1))
    h5write(out_file, "samples", convert.(Int8, samples .- 1))
    h5write(out_file, "logp", logp)

end


function create_samples_folder(folder_path::String; nsamples::Int = 1000, pc = 0.01, overwrite=false)

    for model_path in sort(filter(s->endswith(s, ".jld"), readdir(folder_path, join=true)), by=s->filesize(s))
        out_file = model_path*"_samples.h5"
        if isfile(out_file) && !overwrite
            continue
        end
        println("creating samples for $model_path")
        create_samples(model_path, out_file, nsamples; pc=pc, overwrite=overwrite)
    end
end

function create_samples(ardca_file::String, nsamples::Int; pc::Float64=0.01)

    arnet = load_ardca(ardca_file)
    N = length(arnet.idxperm)
    q = length(arnet.H[1])
    add_pseudocount!(arnet, q, pc)

    samples = zeros(Int64, nsamples * N, q, N)
    logp = zeros(nsamples * N* q)

    @showprogress for sample_id in 1:nsamples
        @threads for i in 1:N
            k = (sample_id-1)*N + i
            seq = rand(1:q, N)
            for a in 1:q
                seq[i] = a
                samples[k, a, :] .= seq[:]
                logp[(k-1)*q+a] = sum(log.(arnet(seq))) 
            end
        end
    end

    return samples, logp

end
