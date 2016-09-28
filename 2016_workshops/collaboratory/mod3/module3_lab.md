---
layout: post2
permalink: /bioinformatics_on_big_data_module3_lab/
title: Bioinformatics on Big Data 2016 Module 3 Lab
header1: Bioinformatics on Big Data 2016
header2: Module 3 Lab
image: CBW_bigdata_icon.jpg
---

# Bioinformatics on Big Data Module 3 Lab

This lab was created by George Mihaiescu

## Introduction 

### Description of the lab

Welcome to the lab for Working Reproducibly in the Cloud! This lab will take you through the steps to setup and configure your virtual machine with Docker packages and will show you how to access data in the Cloud.

After this lab, you will be able to:

* Setup and launch a virtual machine
* Setup and launch a Docker container
* Configure a VM with Docker packages 
* Learn how to access protected data in the Cloud as well as non-protected/public data

Things to know before you start:

The lab may take between 1-2 hours, depending on your familiarity with Cloud Computing and alignment tasks.
   
### Requirements

* Laptop connected to OICR's wifi network  
* Web browser
* Collaboratory credentials 

**Note:** The Collaboratory credentials you are given for the workshop will only work during the workshop.

## Log In to the Collaboratory

In your browser, go to <https://console.cancercollaboratory.org>.  Log in using your provided credentials.

![image_a](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/collaboratory/mod3/mod3_a.png?raw=true)

Once logged in, the first page open will be the "Overview Page" that shows how many resources the project you are part of has access to, as well as the current usage.

![image_a](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/collaboratory/mod3/mod3_b.png?raw=true)

## Before Launching a Virtual Machine

### Create a SSH Key-pair

#### In the Collaboratory

This is the method we will be using for this workshop.

In the bar on the left of the page, click on "Access and Security."  At the top of the page, under "Access and Security", select the "Key Pairs" tab.  Click on the "Create Key Pair" button.  Name your key-pair and click on the "Create Key Pair" button.

![image_aa](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/collaboratory/mod3/mod3_c.png?raw=true)

#### On Windows with PuTTY

If you are using a Windows computer and do not have access to Collaboratory the [PuTTY Key Generator](https://the.earth.li/~sgtatham/putty/latest/x86/puttygen.exe) can be used to create SSH Key-pairs.  Once the program is installed, open it and click on "Generate."

<img src="https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/collaboratory/mod3/mod3_aa.png?raw=true" class="center">

<img src="https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/collaboratory/mod3/mod3_bb.png?raw=true" class="center">

Once the key is generated, you can add a Key comment (ie your name and the date the key was created) and key passphrases.  The key can be saved as a public key or a private key.

<img src="https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/collaboratory/mod3/mod3_cc.png?raw=true" class="center">

#### On Mac/Linux

If you are using a Unix system and do not have access to Collaboratory in Terminal, the command `ssh-keygen` can be used to create SSH key-pairs.  You will be prompted to enter a file name for the key and a passphrase.  Leaving the file name blank will save the key /home/user/.ssh/id_rsa file.

```
ssh-keygen
Generating public/private rsa key pair.
Enter file in which to save the key (/home/user/.ssh/id_rsa):
Enter passphrase (empty for no passphrase): 
Enter same passphrase again: 
Your identification has been saved in testkey.txt.
Your public key has been saved in testkey.txt.pub.
The key fingerprint is:
SHA256:OnHW+VAa/n05qiS2V0HqHbaWAln9J0njgy9zbwYFIa8
```

You will need to change the permissions on the key.  Use `ls` to view the files in your .ssh directory and `chmod` to change the permissions.  Use `cat` to view the contents of the key file.

```
ls ~/.ssh
id_rsa
id_rsa.pub
chmod 400 ~/.ssh/id_rsa
cat ~/.ssh/id_rsa.pub
```

### Import an Existing SSH Key-pair

Alternatively, you can import a key-pair by hitting the "Import Key Pair" button.  You will need to name your key-pair and paste the public key into the text box.  Click on "Import Key Pair" when done.

![image_a](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/collaboratory/mod3/mod3_d.png?raw=true)

If you already have a SSH key-pair on your laptop that you want to use, feel free to upload the public key part by clicking the "Import Key Pair" button.
Otherwise, choose the "Create Key Pair" option and download the private part of the SSH key that will be created.

### Customize Your Security Groups

You will need to know your IP address for this.  To find you IP address, open a new tab or window and go to Google and search for "what is my ip".

![image_a](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/collaboratory/mod3/mod3_e.png?raw=true)

Return to the Collaboratory page.  Select the "Security Groups" tab and click on the "Create Security Group" button.  Name your security group (ie ssh_yourname) and write a description.  Click on "Create Security Group".

![image_a](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/collaboratory/mod3/mod3_f.png?raw=true)

You will need to allow SSH access from your IP address.  Beside the name for the security group you just created, click on "Manage Rules".  Click on the "Add Rule" button.  

In the dropdown menus and boxes, select or enter:  
* Custom TCP Rule   
* Ingress  
* Port  
* 22  
* CIDR  
* your IP address  

Repeat this step and add a second rule with allowing TCP port 80:
* Custom TCP Rule   
* Ingress  
* Port  
* 80  
* CIDR  
* your IP address  

![image_a](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/collaboratory/mod3/mod3_g.png?raw=true)


### Choose Your Flavor

In the menu on the left, select "Instances."  Click on the "Launch Instance" button.

![image_a](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/collaboratory/mod3/mod3_h.png?raw=true)

Make sure you are in the "Details" tab.  

In the dropdown menus and boxes, select or enter:

* Nova  
* Name you instance  
* c1.micro  
* 1  
* Boot from image  
* Ubuntu 16.04 - 2016.04.25 (297.9 MB)  

![image_a](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/collaboratory/mod3/mod3_i.png?raw=true)

Select the "Access and Security" tab.  Select the key pair you previously created and check the box beside "ssh".  

![image_a](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/collaboratory/mod3/mod3_j.png?raw=true)

Select the "Networking" tab.  Choose the appropriate network.

![image_a](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/collaboratory/mod3/mod3_k.png?raw=true)

Launch the instance by hitting the "Launch" button.

![image_a](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/collaboratory/mod3/mod3_l.png?raw=true)

It will take a few minutes for the instance to start.

To view your instances, in the left hand menu, click on "Instances".

![image_a](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/collaboratory/mod3/mod3_m.png?raw=true)

## Log Into Your Instance

### Mac/Linux Instructions

You will need to change the file permissions for your private SSH key.

```bash
 chmod 400 path_to_private_key
```

Where `path_to_private_key` is the path to where you have saved your key.


Then, to log in:

```bash
 ssh -i path_to_private_key ubuntu@142.1.177.XXX
```

XXX is the last octet from the floating IP address you assigned to the instance.

### Windows Instructions

To configure Putty, start Putty and do the following:

* Fill in the "Host name" field with 142.1.177.XXX.

XXX is the last octet from the floating IP address you assign to the instance.

 
<img src="../../../resources/Putty_Basic_Options.png" alt="Basic Putty Options" class="center">

* In the left hand categories,under the Connection category choose Data.  In the auto-login username field write ***ubuntu***.

<img src="../../../resources/Putty_Data_Options.png" alt="Putty Data Options" class="center"> 

* In the left hand categories, in the Connection category next to SSH click on the **+**. Click on Auth. In the private-key file for authentication field, hit browse and find your private key.

<img src="../../../resources/Putty_Auth_Options.png" alt="Putty Auth Options" class="center">

* In the left hand categories, click on Session.  In the Saved Sessions field write **Collaboratory** and click save.

**Now that Putty is configured**, all you have to do is start PuTTY and double-click on "Collaboratory" to login.

## Customize Your Virtual Machine

You will need to upgrade your package index and existing packages by running:

```
sudo apt-get update && sudo apt-get upgrade
```


## Docker Installation

Run the following commands to install the Docker engine software and required dependencies:

```
sudo apt-get install -y apt-transport-https ca-certificates unzip
sudo apt-key adv --keyserver hkp://p80.pool.sks-keyservers.net:80 --recv-keys 58118E89F3A912897C070ADBF76221572C52609D

echo 'deb https://apt.dockerproject.org/repo ubuntu-xenial main' | sudo tee /etc/apt/sources.list.d/docker.list

sudo apt-get update && sudo apt-get -y install linux-image-extra-$(uname -r) linux-image-extra-virtual

sudo apt-get -y install docker-engine

sudo service docker start

sudo docker run hello-world
```

This will install Docker and run hello-world to test the installation.

## Run a Bioinformatics Tool in Docker

We will first need a data file.  To get the file from the ftp server, we will use the `wget` command.

```
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA12878/alignment/NA12878.chrom20.ILLUMINA.bwa.CEU.low_coverage.20121211.bam
```

We will need the Docker container with the "bamstats" tool.  We will use the Docker `pull` command to retrieve the container.

```
sudo docker pull quay.io/briandoconnor/dockstore-tool-bamstats
```

Run the container, mounting the sample file inside.

```
sudo docker run -it -v `pwd`/NA12878.chrom20.ILLUMINA.bwa.CEU.low_coverage.20121211.bam:/NA12878.chrom20.ILLUMINA.bwa.CEU.low_coverage.20121211.bam -v /tmp:/home/ubuntu quay.io/briandoconnor/dockstore-tool-bamstats
```

In this command, `-it` means run interactively and `-v` maps a file or directory from the VM inside the Docker container.  Whatever is created inside /home/ubuntu (inside the container) will be in the host tmp directory.  This will allow the files that are create to survive after the container is terminated.

The OS inside the Docker container is different than that of the host VM.  You can check this with:

```
cat /etc/lsb-release
```

Recall that we chose Ubuntu 16.04 for our VM.

Inside the Docker container, execute the bamstats binary against the sample file.

```
cd && /usr/local/bin/bamstats 4 /NA12878.chrom20.ILLUMINA.bwa.CEU.low_coverage.20121211.bam
```

Exit the docker container by typing "exit" and go to "/tmp" where the report was created.

```
exit
cd /tmp 
unzip bamstats_report.zip 
sudo python3 -m http.server
```

Visit the page to see the statistics for that sample BAM:
	<http://142.1.177.XXX/bamstats_report.html>
	
XXX is the last octet from the floating IP address you assigned to the instance.


## Access Data in the Cloud

**Note that you need DACO approval to access ICGC data stored in the Cloud**

**The following exercises will be done by the instructor only**

### Create a Configuration File

First, determine how many cores and how much memory you can allocate to the download.

```
cat /proc/cpuinfo | grep -c processor
free -g
```

On a VM with 8 cores and 55 GB of RAM, you can allocate 7 cores and 49 GB of RAM (7 GB per thread).  Create a text file `application.properties` that contains the access token and the number of cores and memory per core.  Use `cat` to view the file contents.  

```
cat application.properties
accessToken=XXX
transport.parallel=7
transport.memory=7
```

### Download Protected Data

Pull the Docker container containing the storage client.

```
 docker pull icgc/icgc-storage-client
```

Initiate the download, mounting the application.properties file as well as the destination directory for the download (-v /tmp/:/data),  so it survives the termination of the Docker container:

```
sudo docker run -v /tmp/:/data -v /home/ubuntu/application.properties:/icgc/icgc-storage-client/conf/application.properties --privileged icgc/icgc-storage-client bin/icgc-storage-client --profile collab download --object-id 6329334b-dcd5-53c8-98fd-9812ac386d30 --output-dir /data
```

It takes around 30 min for a 120 GB file to be downloaded using a VM with 8 cores and 56 GB of RAM, including the time needed by the storage client to perform an automated checksum to verify downloaded data integrity.

The download speed depends on the disk IO which is shared with other VMs running on the same physical server, as well as other shared resources (network, storage cluster).

### Important Notes

ICGC data stored in AWS S3 is only available from EC2 instances.

Download of ICGC protected data from S3 is only available within the “us-east-1” EC2 region of AWS.

 ICGC data stored in Cancer Genome Collaboratory is only available from VMs inside Collaboratory, or from external IP addresses listed on the user enrollment form.

It is the responsibility of the users to protect access to the EC2 or Collaboratory VMs, as well as the restricted data they have access to.







