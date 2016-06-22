println("""
MIToS (Julia)
=============
""")

cd("MIToS")

println(readall(`julia ReadMSA.jl`))
println(readall(`julia Entropy.jl`))
println(readall(`julia PID.jl`))
println(readall(`julia MIp.jl`))

cd("..")

println("""
ProDy (Python + C)
==================
""")

cd("ProDy")

println(readall(`python ReadMSA.py`))
println(readall(`python Entropy.py`))
println(readall(`python PID.py`))
println(readall(`python MIp.py`))

cd("..")

println("""
Bio3D (R)
=========
""")

cd("Bio3D")

println(readall(`Rscript ReadMSA.R`))
println(readall(`Rscript Entropy.R`))
println(readall(`Rscript PID.R`))

cd("..")
