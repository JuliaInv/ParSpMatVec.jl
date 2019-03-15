module ParSpMatVec
using SparseArrays
using LinearAlgebra
const spmatveclib  = abspath(joinpath(splitdir(Base.source_path())[1],"..","deps","builds","ParSpMatVec"))

include("A_mul_B.jl")
include("Ac_mul_B.jl")

export isBuilt
function isBuilt()
	A = sprandn(10,10,0.25);
	x = rand(10);
	y = rand(10);
	try 
		ParSpMatVec.Ac_mul_B!( 1.0, A, x, 0.0, y, 1);
		return true; 
	catch
		warn("ParSpMatVec is not built");
		return false;
	end
end


end # module
