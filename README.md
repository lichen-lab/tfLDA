# tfLDA
Application of topic models to a compendium of ChIP-Seq datasets uncovers recurrent transcriptional regulatory modules


# Maintainer
Li Chen <li.chen@auburn.edu>


# Install circMeta
```r
install.packages("devtools")
library(devtools)
install_github("lichen-lab/tfLDA")
```


# Descriptions for circJuncDE

## Usage
tfLDA(obj, topic = c(5, 10, 20), seed = 1234, iterations = 500, 
    burnin = NULL, alpha = 50, eta = 0.1, alphaByTopic = TRUE)

## Arguments
*  obj: A tfLDA object, which can be created by createtfLDAobj function
*  topic: a series of preset number of topics
*  seed: seed for performing MCMC
*  iterations: iterations to run MCMC. Default is 500
*  bunrin: burnin for MCMC. Default is NULL
*  alpha: Dirichlet hyperparameter for topic distribution for each document
*  eta: Dirichlet hyperparameter for motif distribution for each topic
*  alphaByTopic: if TRUE, alpha=alpha/number of topic

                 
                
## Output values
* a tfLDA object contains the same output values of lda R package

# Examples of applying tfLDA for K562 and HeLaS3 




