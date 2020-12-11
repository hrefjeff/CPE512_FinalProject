# Amazon Web Services

# Problem Introduction

This paper compares and contrasts execution of parallelization of a program on the DMC environment and AWS environment. The Synchronous 2D Heat Transfer program is used and the results of the program running on both environments are analyzed. Given the accessible nature of AWS I wanted to know how much faster or slower the performance of a program ran on its environment compared to the DMC environment.

# Parallelization Approach Employed

I deployed the 2-D temperature program with paired locally blocking MPI_Send and MPI_Recv routines to both the DMC and AWS computing environments. The programming language chosen for this experinemnt is C. AWS computing environments come with the OpenMPI C wrapper compiler mpicc by default and installing mpic++ was not explained in detail. Upon trying to install it I ran into some issues where the nodes that were spun up during the process did not install correctly. I ultimately had to revert back to the C version.

# Empirical Method and Results

Uploading the 2D locally blocking with paired MPI send/recv program to each environment. Running the code and collecting execution time information.

**Excel charts of both will be uploaded here**

# General Conclusions

The general consensus is that the AWS environment out performed the DMC environment when using the most expensive option of computing available. The c5x.large did better.

## How effective were your parallel implementations?

## What were sources of bottlenecks?
## What are the limits of your approach?


## Do you have other ideas that may provide better results in the future?
## Can your parallel method be applied to other problems?
# Appendix
Source Code for your problem, example program runs


## Helpful Links

1. https://docs.aws.amazon.com/parallelcluster/latest/ug/tutorials_03_batch_mpi.html <----- THE GOAT DOC
1. https://docs.aws.amazon.com/parallelcluster/latest/ug/tutorials_01_hello_world.html
1. https://docs.aws.amazon.com/parallelcluster/latest/ug/networking.html#awsbatch-networking
1. https://slurm.schedmd.com/squeue.html#lbAJ
1. https://console.aws.amazon.com/cloudwatch/home?region=us-east-1#dashboards:name=parallelcluster-finalproject-us-east-1 

## Topics

1. AWS ParallelCluster
1. AWS Batch
1. Getting your code onto the 
1. Script explanation
1. Cloudwatch logging
1. What are processes?

## Open questions

1. How do I run an MPI cpp job in AWS slurm
1. How do I monitor my runs
1. How do I collect billing/cost information associated with runs
1. 









# Random resources

1. https://aws.amazon.com/premiumsupport/knowledge-center/batch-job-stuck-runnable-status/
1. https://docs.aws.amazon.com/parallelcluster/latest/ug/install.html
1. https://gist.github.com/sean-smith/d29e83b64d7ac3f2195e349447d2dbfe 
1. 
