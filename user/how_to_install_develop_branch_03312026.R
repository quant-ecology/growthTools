## I recommend getting the development version of growthTools, if you aren't using
# it already, it has lots of new features. Here's how to get it for now:

#### growthTools depends on mleTools, which you can install via:
#pak::pak("ctkremer/mleTools")

#### Install development version of growthTools:

# remove any old versions of the package before installing
remove.packages('growthTools')

# updated code for installing package:
install.packages("pak")
library(pak)
pak::pak("quant-ecology/growthTools@develop")
