# Using the Guillimin HPC

This document provides information on how to get started with using the Guillimin HPC.


## Connecting to the Guillimin HPC scheduler:

**Host:** username@guillimin.clumeq.ca
**Login / Password:** Provided to you on a separate piece of paper at the beginning of the workshop.

**Note:** Your home folder is limited to 10GB of storage. The tutorials we will do shouldn't make you go over this quota, but to be aware of the remaining space you have, you can use the command ```myquota```.


#### Linux/Mac
Use the following command from the command line:
```
ssh username@guillimin.clumeq.ca
```
The username and password are the ones that we provided on a separate piece of paper.

#### Windows

Connect to the ```guillimin.clumeq.ca``` server using the **Putty** tool that you installed for this workshop.
  
 



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
 
Jobs are submitted to the scheduler using the qsub command. The important qsub parameters for this workshop are the following:
* **nodes** : total number of machines to use for this job
* **ppn** : total number of cores to use for this job
* **-A bem-651-ae** : is the name of the allocation to use, in this case 'bem-651-ae'
* **-d workdir** : The working directory, in other words the path in which this command will be executed 

Please keep in mind that you need always need to load the modules you want to use in your scripts.
 
 

### Using a bash script

You can specify a bash script to launch as a compute job directly from the command line, like the ```testjob.sh``` script in this example:
```
qsub -l nodes=1:ppn=1 -A bem-651-ae testjob.sh
```

where the content of ```testjob.sh``` could be:
```
#!/bin/bash
#PBS -l nodes=1:ppn=2
#PBS -A bem-651-ae

module use /cvmfs/soft.mugqic/CentOS6/modulefiles/
module load mugqic/bowtie2/2.2.4
perl bismark_genome_preparation --bowtie2 --verbose .
```


### Specifying the whole command from the prompt

If you do not want to put your command in a script, you can write it directly in the prompt using the ```echo``` command like this:
```echo "touch test.txt" | qsub -l nodes=1:ppn=1 -A bem-651-ae -d .```


## Downloading a file from Guillimin to your computer

#### Linux/Mac
Use the following from the command line:
```
cd /path/to/destination/folder
scp user99@guillimin.clumeq.ca:/home/user99/myfile.txt .
```

To download a full directory, add the **-R** option
```
cd /path/to/destination/folder
scp -R user99@guillimin.clumeq.ca:/home/user99/myfile.txt .
```

#### Windows
Use the WinSCP software that you installed for this workshop, connecting to ```guillimin.clumeq.ca```.
