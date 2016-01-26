module ParSpMatVec

const spmatveclib  = abspath(joinpath(splitdir(Base.source_path())[1],"..","deps","builds","ParSpMatVec"))

include("A_mul_B.jl")
include("Ac_mul_B.jl")

end  # module ParSpMatVec

