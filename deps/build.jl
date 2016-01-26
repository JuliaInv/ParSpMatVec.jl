using BinDeps

@BinDeps.setup

# set to true to support intel fortran compiler
useIntelFortran = false

# construct absolute path
depsdir = splitdir(Base.source_path())[1]
builddir = joinpath(depsdir,"builds")
srcdir  = joinpath(depsdir,"src")

if !isdir(builddir)
	mkdir(builddir)
end

@unix_only begin
	oldDir = pwd()
	cd(srcdir)
	@build_steps begin
		if useIntelFortran
			run(`ifort --O3 -xHost -fPIC -fpp -openmp -integer-size 64 -diag-disable=7841 -shared Ac_mul_B.f90  A_mul_B.f90 -o ../builds/ParSpMatVec`)
		else
			run(`gfortran -O3 -fPIC -cpp -fopenmp -fdefault-integer-8 -shared  Ac_mul_B.f90  A_mul_B.f90 -o ../builds/ParSpMatVec`)
		end
	end
	cd(oldDir)
end

