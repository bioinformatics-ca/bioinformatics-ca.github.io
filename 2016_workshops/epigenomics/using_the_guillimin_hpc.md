# Using the Guillimin HPC

This document provides information on how to get started with using the Guillimin HPC.


## Connecting to the Guillimin HPC scheduler:
On Linux/Mac, use the following command from the command line:
<code>
ssh username@guillimin.clumeq.ca
</code>
The username and password are the ones that we provided on a separate piece of paper.

On Windows, connect to the ```guillimin.clumeq.ca``` server using the Putty tool that you previously installed.
  
Note: Your home folder is limited to 10GB of storage. The tutorials we will do shouldn't make you go over this quota, but to be aware of the remaining space you have, you can use the command ```myquota```. 



## Loading modules
All available software are organized as modules that you can load. For example, to start using Java 7, you need to load the Java 7 module:

**1- Import the CVMFS list of modules**
module use /cvmfs/soft.mugqic/CentOS6/modulefiles

**2- Get the list of available modules**
module avail

**3- Load your module of choice**
module load mugqic/java/openjdk-jdk1.7.0_60

**4- You can now run the module you loaded directly**
```
[class99@lg-1r17-n04 ~]$ javac -version
javac 1.7.0_60-ea
```





## Launching a job on the scheduler
**Note:** Complete documentation on the Guillimin scheduler is available here: http://www.hpc.mcgill.ca/index.php/starthere/81-doc-pages/91-guillimin-job-submit 

In order to launch a script as a compute job directly from the command line, you can specify it by creating a script that will be executed, like the ```testjob.sh``` script in this example:
```qsub -l nodes=1:ppn=1 -A bem-651-ae testjob.sh```

where the content of ```testjob.sh``` could be:
```
#!/bin/bash
#PBS -l nodes=1:ppn=2
#PBS -A bws-221-af
#PBS -o convert_genome.output
#PBS -e convert_genome.error
#PBS -N convert_genome

module use /cvmfs/soft.mugqic/CentOS6/modulefiles/
module load mugqic/bowtie2/2.2.4
perl bismark_genome_preparation --bowtie2 --verbose .
```

If you do not want to put your command in a script, you can write it directly in the prompt using the ```echo``` command like this:
```echo "touch test.txt" | qsub -l nodes=1:ppn=1 -A bem-651-ae -d .```

Important parameters for the qsub command in this workshop are the following:
* **nodes** : total number of machines to use for this job
* **ppn** : total number of cores to use for this job
* **-A bem-651-ae** is the name of the allocation to use, in this case 'bem-651-ae'
* **-d workdir** The working directory, in other words the path in which this command will be executed.

Please keep in mind that you need always need to load the modules you want to use in your scripts.
 
