---
layout: post2
permalink: /bioinformatics_on_big_data_module4_lab/
title: Bioinformatics on Big Data 2016 Module 4 Lab
header1: Bioinformatics on Big Data 2016
header2: Module 4 Lab
image: CBW_bigdata_icon.jpg
---

# Bioinformatics on Big Data Module 4 Lab

This lab was created by George Mihaiescu

## Introduction 

### Description of the lab

*write description here something like:*

Welcome to the lab for Sharing and Scaling a VM! 

After this lab, you will be able to:

* Snapshot a VM and share it with other cloud tenants.  
* Launch a new VM from a shared image.  
* Scale out your VMs to meet the needs of this task.  

Things to know before you start:

The lab may take between 1-2 hours, depending on your familiarity with Cloud Computing and VMs. Don't worry if you don't complete the lab! It will be made available for you to complete later.   
There are a few thought-provoking Questions or Notes pertaining to sections of the lab. These are optional, and may take more time, but are meant to help you better understand the the tasks you are performing. These questions will be denoted by boxes, as follows:

   Question(s):

   * Thought-provoking question goes here
   
### Requirements

* Laptop connected to OICR's wifi network  
* Web browser  
* The VM from Module 3 Lab



## Customize Your VM

Install some packages.

```
apt-get –y install vim atop
``

## Prepare to Take a Snapshot

SSH into the VM configured in Module 3 Lab

### For Mac/Linux Users

To log back in, open terminal and type:   

```bash
 ssh -i path_to_private_key ubuntu@142.1.177.XXX
```
XXX is the last octet from the floating IP address you assigned to the instance.

### For Windows Users

To log back in, open PuTTY and double-click on "Collaboratory" to login.

### Clean Up Its History

```
rm -i ~/.bash_history
```
**Note:** the `-i` option in `rm` will prompt you to answer yes or no to remove a file.  It's a good idea to use `-i` with `rm` so that you don't inadvertently delete the wrong thing.

Remove your public key from the default user’s “authorized_keys” file, so they don’t have free access to virtual machines based on the snapshot.

```
rm –i ~/.ssh/authorized_keys
```

**Note:** after removing your public key you will not be able to log in again into this virtual machine, and you will have to start a new one from the snapshot 

## Take the Snapshot

In your browser, go to <https://console.cancercollaboratory.org>.  Log in using your provided credentials.

On the lefthand side, expand the "Compute" menu and then click on "Instances."  In the table, click on "Create Snapshot" for your instance.

Enter a descriptive name for your snapshot, preferably "your_name_snap" and click "Create Snapshot."

![image_a](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/collaboratory/mod4/mod4_a.png?raw=true)

## Terminate the Old Instance

To terminate your old instance, place a check mark in the box beside the name of your instance and click on the red "Terminate Instances" button.  In the dialogue box, verify the instance name and click on "Terminate Instance."

![image_a](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/collaboratory/mod4/mod4_b.png?raw=true)

## Start a New Image from the Snapshot

On the lefthand side, from the expanded "Compute" menu, select "Images."  In the table, you should see the snapshot that you just took.  Click on "Launch Instance" for your snapshot.

Choose the same or a larger flavor when starting an instance from a snapshot. Assign a floating IP, then SSH into it and verify the packages you installed before are present.

![image_a](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/collaboratory/mod4/mod4_c.png?raw=true)

## Install the CLI Tool to Share Your Snapshot

In order to share the snapshot you took with other Openstack projects, you will need to use a CLI tool, python-glanceclient.

```
sudo apt-get -y install python-glanceclient
```

## Download the Environment Variables

In your browser on the Collaboratory website, in the lefthand side in the expanded "Compute" menu, select "Access and Security."

Click on the "API Access" tab.

??? what are they downloading

![image_a](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/collaboratory/mod4/mod4_d.png?raw=true)

## Source the Credentials

Create a file and paste the contents of the download, then source them.

![image_a](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/collaboratory/mod4/mod4_e.png?raw=true)

??? Are students creating a file in the VM?  How are students creating the file and sourcing it?

## Share the Snapshot

Execute the following commands to share the snapshot:

```
glance image-list  -> find the ID of the snapshot
```

Share it with the project ID a295a12382ef451b9bf2353d624455bb

```
glance member-create b42a556a-9525-4575-9828-9f01833a939d a295a12382ef451b9bf2353d624455bb
```

## Using `cloud-init` to Bootstrap a VM

Launch a new Ubuntu 16.04 instance using flavor c1.micro and paste the following in the “Post-creation” tab, “Direct Input” window:

```
#!/bin/sh
# This script will install Docker with all its dependencies, download a Docker container and a sample file
# Run bamstats on the sample file and then serve the report over HTTP

sudo apt-get install -y apt-transport-https ca-certificates unzip
sudo apt-key adv --keyserver hkp://p80.pool.sks-keyservers.net:80 --recv-keys 58118E89F3A912897C070ADBF76221572C52609D
echo 'deb https://apt.dockerproject.org/repo ubuntu-xenial main' | sudo tee /etc/apt/sources.list.d/docker.list
sudo apt-get update && sudo apt-get -y install linux-image-extra-$(uname -r) linux-image-extra-virtual
sudo apt-get -y install docker-engine
sudo service docker start
sudo docker run hello-world
sudo wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA12878/alignment/NA12878.chrom20.ILLUMINA.bwa.CEU.low_coverage.20121211.bam -O /NA12878.chrom20.ILLUMINA.bwa.CEU.low_coverage.20121211.bam
sudo docker pull quay.io/briandoconnor/dockstore-tool-bamstats
sudo docker run -v `pwd`/NA12878.chrom20.ILLUMINA.bwa.CEU.low_coverage.20121211.bam:/NA12878.chrom20.ILLUMINA.bwa.CEU.low_coverage.20121211.bam -v /tmp:/home/ubuntu quay.io/briandoconnor/dockstore-tool-bamstats  /bin/bash -c "cd; /usr/local/bin/bamstats 4 /NA12878.chrom20.ILLUMINA.bwa.CEU.low_coverage.20121211.bam"
cd /tmp && unzip /tmp/bamstats_report.zip
cd /tmp && sudo python3 -m http.server 80 &
```

![image_a](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/collaboratory/mod4/mod4_f.png?raw=true)

Monitor the bootstrap progress by checking the boot console output in “Instance details” under the “Log” tab.

![image_a](https://github.com/bioinformatics-ca/bioinformatics-ca.github.io/blob/master/2016_workshops/collaboratory/mod4/mod4_g.png?raw=true)

