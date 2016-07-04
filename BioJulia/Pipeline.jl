using Bio

macro output_timings(t,name)
    quote
        @printf "BioJulia,%s,%f\n" $name minimum($t)
    end
end

macro timeit(ex,name)
    quote
        t = Float64[]
        for i in 1:6
            e = 1000*(@elapsed $(esc(ex)));
            if i > 1
                # warm up on first iteration
                push!(t, e)
            end
        end
        @output_timings t $name
    end
end
