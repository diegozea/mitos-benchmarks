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

cd("ProDy")

println(read(`python ReadMSA.py`, String))
println(read(`python Entropy.py`, String))
println(read(`python PID.py`, String))
println(read(`python MIp.py`, String))

cd("..")

println("""
Bio3D (R)
=========
""")

cd("Bio3D")

println(read(`Rscript ReadMSA.R`, String))
println(read(`Rscript Entropy.R`, String))
println(read(`Rscript PID.R`, String))

cd("..")
