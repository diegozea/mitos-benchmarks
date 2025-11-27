println("""
MIToS (Julia)
=============
""")

cd("MIToS")

println(read(`julia ReadMSA.jl`, String))
println(read(`julia Entropy.jl`, String))
println(read(`julia PID.jl`, String))
println(read(`julia MIp.jl`, String))

cd("..")

println("""
ProDy (Python + C)
==================
""")

function safe_read(cmd)
    try
        return read(cmd, String)
    catch e
        return "[ERROR] $(cmd): $(e)"
    end
end

cd("ProDy")

println(safe_read(`python ReadMSA.py`))
println(safe_read(`python Entropy.py`))
println(safe_read(`python PID.py`))
println(safe_read(`python MIp.py`))

cd("..")

println("""
Bio3D (R)
=========
""")

cd("Bio3D")

println(safe_read(`Rscript ReadMSA.R`))
println(safe_read(`Rscript Entropy.R`))
println(safe_read(`Rscript PID.R`))

cd("..")

println("""
Biopython (Python)
==================
""")

cd("BioPython")

println(read(`python ReadMSA_Stockholm.py`, String))
println(read(`python ReadMSA_Fasta.py`, String))
pid_cmd = setenv(Cmd(["python", "PID.py"]), ENV)
println(read(pid_cmd, String))
println(read(`python Entropy_MI.py`, String))
println(read(`python StockholmAnnotations.py`, String))

cd("..")
