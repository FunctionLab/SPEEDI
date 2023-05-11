# Add our new version of the addGeneIntegrationMatrix function to the ArchR namespace
# This makes it so we don't have to define the namespace of all of the ArchR functions
# inside the function
environment(addGeneIntegrationMatrix_SPEEDI) <- asNamespace('ArchR')
