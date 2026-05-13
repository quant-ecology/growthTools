
# Works using R 4.3.3

# fails using R 4.4.0 ("In R 4.4.0, there were subtle internal changes in R objects that RCall.jl relies on.
# DiffEqR wraps R functions in Julia’s RFunction type.
# The current DiffEqR version you have is likely not compatible with R 4.4.0, so RFunction fails when trying to access its r field.")

# trying 4.5.3
# fails, same error. 

# ====== R Script to Install a Specific Version of DiffEqR on Mac ======

# 1. Install required R packages
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
if (!requireNamespace("JuliaCall", quietly = TRUE)) install.packages("JuliaCall")

# installed/updated all dependencies recommended

library(remotes)
library(JuliaCall)

# 2. Setup Julia
# If you already installed Julia manually, point to its path. Example path for Julia 1.9:
julia_path <- "/Applications/Julia-1.10.app/Contents/Resources/julia/bin" 
julia_setup(JULIA_HOME = julia_path, installJulia = FALSE)

# Loaded Julia in bash, ran:
# import Pkg
# Pkg.add("DiffEqBase")
# Pkg.add("OrdinaryDiffEq")  # often required by DiffEqR
# Pkg.add("DiffEqOperators") # sometimes needed

julia_command("import DiffEqBase")
 

# 3. Specify the DiffEqR version you want
#diff_eqr_version <- "v1.2.0"  # <-- change this to the version you want
diff_eqr_version <- "v2.0.0"  # <-- change this to the version you want

# 4. Install that specific version from GitHub
remotes::install_github(paste0("SciML/DiffEqR@", diff_eqr_version))

# 5. Load DiffEqR and test
library(diffeqr)

# Verify the installed version
cat("Installed DiffEqR version:", as.character(packageVersion("diffeqr")), "\n")
# Installed DiffEqR version: 2.0.0


# ====== Simple test ======

# this used to fail...
de<-diffeqr::diffeq_setup()

f <- function(u,p,t) {
  return(1.01*u)
}

u0 <- 1/2
tspan <- c(0., 1.)

prob = de$ODEProblem(f, u0, tspan)
sol = de$solve(prob)

sol$.(0.2)

plot(sol$t,sol$u,"l")

f <- function(u,p,t) {
  du1 = p[1]*(u[2]-u[1])
  du2 = u[1]*(p[2]-u[3]) - u[2]
  du3 = u[1]*u[2] - p[3]*u[3]
  return(c(du1,du2,du3))
}

u0 <- c(1.0,0.0,0.0)
tspan <- c(0.0,100.0)
p <- c(10.0,28.0,8/3)
prob <- de$ODEProblem(f, u0, tspan, p)
sol <- de$solve(prob)

mat <- sapply(sol$u,identity)

udf <- as.data.frame(t(mat))

matplot(sol$t,udf,"l",col=1:3)


# Alternative: work directly in Julia code?
# but then how best to pass parameters...

# Define the ODE function in Julia
julia_command('using DifferentialEquations; f(u,p,t) = -u; u0 = 1.0; tspan = (0.0, 1.0); prob = ODEProblem(f,u0,tspan); sol = solve(prob);')


# Pull the solution back into R
sol <- julia_eval("sol")
print(sol)


cmd<-paste()

