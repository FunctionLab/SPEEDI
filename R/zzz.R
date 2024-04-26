# Add our new version of the addGeneIntegrationMatrix function to the ArchR namespace
# This makes it so we don't have to define the namespace of all of the ArchR functions
# inside the function
environment(addGeneIntegrationMatrix_SPEEDI) <- asNamespace('ArchR')
environment(log_warning) <- asNamespace('logr')

# We reset the RNG to its default to counteract ArchR automatically loading an
# alternative RNG when using more than 1 thread in an interactive session
.onLoad <- function(...) {
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
}
