using BenchmarkTools
using JLD

res = JLD.load("result_SmallBenchmark.jld")

for (t,d,l) in [("#### MSA module","MIToS.MSA.","small_bench_msa"),
              ("#### PDB module","MIToS.PDB.","small_bench_mitos_pdb"),
              ("#### Information module","MIToS.Information.","small_bench_information")]
    println("")
    println(t)
    for (k,v) in res["small_result"][l]
        println("")
        println("|", k, "| |")
        println("|-|-|")
        for (kk,vv) in minimum(v)
            println("|", replace(kk,d,"") , "|", BenchmarkTools.prettytime(time(vv)), "|")
        end
    end
end

println("#### Pfam/pipeline\n")
println("|-|-|")
for (k,v) in minimum(res["small_result"]["small_bench_pfam"])
    println("|", replace(k,"MIToS.Pfam.","") , "|", BenchmarkTools.prettytime(time(v)), "|")
end

