ENV["JULIA_PROJECT"] = @__DIR__

function safe_read(cmd)
    try
        return read(cmd, String)
    catch e
        return "[ERROR] $(cmd): $(e)"
    end
end

function run_in(dir, cmd)
    cd(dir) do
        println(safe_read(cmd))
    end
end

function r_cmd()
    if Sys.which("Rscript") !== nothing
        return `Rscript Pipeline.R`
    end

    conda = Sys.which("conda")
    if conda !== nothing
        return Cmd([conda, "run", "-n", "mitos-benchmarks", "Rscript", "Pipeline.R"])
    end

    return Cmd(["/bin/false"]) # will be caught by safe_read
end

println("""
MIToS pipeline (Julia)
======================
""")

run_in("MIToS", `julia Pipeline.jl`)

println("""
ProDy pipeline (Python + C)
===========================
""")

run_in("ProDy", `python Pipeline.py`)

println("""
Bio3D pipeline (R)
==================
""")

run_in("Bio3D", r_cmd())

println("""
Biopython pipeline (Python)
===========================
""")

run_in("BioPython", `python Pipeline.py`)
