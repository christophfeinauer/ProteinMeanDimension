using ArDCA
using JLD
using DCAUtils
using StatsBase


function train_folder(datafolder=joinpath(@__DIR__, "../alignments"), outfolder=joinpath(@__DIR__, "models"); lambdaH=0.0001, lambdaJ=0.01, theta=0.2)


    for alignmentpath in sort(filter(s->endswith(s, "a2m"), readdir(datafolder, join=true)), by=f->filesize(f))
        outpath = joinpath(outfolder, basename(alignmentpath)*".ardca.jld")
        if isfile(outpath) && filesize(outpath) > 0
            continue
        end
        arnet, _ = ardca(alignmentpath; lambdaJ=lambdaJ, lambdaH=lambdaH, theta=theta, permorder=:NATURAL)
        save(outpath, "arnet", arnet)
    end

end

function train_folder_loglambdaJsweep(datafolder=joinpath(@__DIR__, "../data"), outfolder=joinpath(@__DIR__, "models"); lambdaH=0.0001, lambdaJstart=-6, lambdaJstop=-1, npoints=60, theta=0.2)

    for alignmentpath in sort(filter(s->endswith(s, "a2m"), readdir(datafolder, join=true)), by=f->filesize(f))
        for lambdaJ in exp10.(range(lambdaJstart, stop=lambdaJstop, length=npoints))
            outpath = joinpath(outfolder, basename(alignmentpath)*".ardca."*"lambdaJ_$lambdaJ"*".jld")
            if isfile(outpath) && filesize(outpath) > 0
                continue
            end
            arnet, _ = ardca(alignmentpath; lambdaJ=lambdaJ, lambdaH=lambdaH, theta=theta, permorder=:NATURAL)
            save(outpath, "arnet", arnet)
        end
    end

end

function calculate_sr(pc=0.01, q=21)
    proteins = ["BRCA1", "GAL4", "UBC9", "SUMO1"]

    for protein in proteins
        # get experimental values
        exp_file = filter(f->startswith(basename(f), protein) && endswith(f, "exp"), readdir(joinpath(@__DIR__, "../data"), join=true))
        if length(exp_file) != 1
            error("found several or no exp files for $protein")
        end
        exp_file = exp_file[1]
        ranks = map(line->parse(Float64, split(line, "/")[end]), filter(f->startswith(f, ">"), readlines(exp_file)))
        Zexp = read_fasta_alignment(exp_file, 1.0)
        for model_file in filter(f->startswith(basename(f), protein) && endswith(f, "jld") && contains(f, "lambdaJ"),  readdir(joinpath(@__DIR__, "models"), join=true))
            out_file = model_file*".sr"
            (isfile(out_file) && filesize(out_file) > 0) && continue
            model = load(model_file)["arnet"]
            model.p0[:] = (1-pc)*model.p0[:] .+ pc*q
            preds = [sum(log.(model(Zexp[:, m]))) for m in 1:size(Zexp, 2)]
            sr = corspearman(preds, ranks)

            @show model_file, sr

            open(out_file, "w") do fid
                println(fid, "$model_file $sr")
            end

        end
    end
end
