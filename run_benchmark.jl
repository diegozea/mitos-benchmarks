println("""
MIToS (Julia)
=============
""")

cd("MIToS")

println(readstring(`julia ReadMSA.jl`))
println(readstring(`julia Entropy.jl`))
println(readstring(`julia PID.jl`))
println(readstring(`julia MIp.jl`))

cd("..")

println("""
ProDy (Python + C)
==================
""")

cd("ProDy")

println(readstring(`python ReadMSA.py`))
println(readstring(`python Entropy.py`))
println(readstring(`python PID.py`))
println(readstring(`python MIp.py`))

cd("..")

println("""
Bio3D (R)
=========
""")

cd("Bio3D")

println(readstring(`Rscript ReadMSA.R`))
println(readstring(`Rscript Entropy.R`))
println(readstring(`Rscript PID.R`))

cd("..")
