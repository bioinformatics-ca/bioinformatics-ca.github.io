# Laptop Programs to Install

1) Install latest version of R which can be downloaded from http://probability.ca/cran/.  If you need help on installing R see below for assistance.

1b) Download and install the most recent version of R Studio desktop from http://www.rstudio.com/.

2) BioConductor core packages. To do this, open R and type at the > prompt, then wait for prompt and type second command:

source("http://bioconductor.org/biocLite.R");

biocLite();

3) A robust text editor. For Windows/PC - notepad++ (http://notepad-plus-plus.org/). For Linux - gEdit (http://projects.gnome.org/gedit/). For Mac – TextWrangler (http://www.barebones.com/products/textwrangler/download.html)

4) A file decompression tool. For Windows/PC – 7zip (http://www.7-zip.org/). For Linux – gzip (http://www.gzip.org). For Mac – already there.

5) fastqc – This tool is available for Windows/Mac/Linux here: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/

6) A robust internet browser such as Firefox or Safari (Internet Explorer and Chrome are not recommended because of Java issues).

7) Java -The visualization program that we will be using (IGV) requires Java. Check if you have Java installed: https://www.java.com/verify/ and download Java if you do not have it installed.

8) Integrative Genomics Viewer 2.3 (IGV) - Once java is installed, go to http://www.broadinstitute.org/igv/ and register in order to get access to the downloads page. Once you have gained access to the downloads page, click on the appropriate launch button that matches the amount of memory available on your laptop (if you have space, 1.2GB is good, more is better). Chrome: Chrome does not launch "java webstart" files by default. Instead, the launch buttons below will download a "jnlp" file. This should appear in the lower left corner of the browser. Double-click the downloaded file to run. Windows users: To run with more than 1.2 GB you must install 64-bit Java. This is often not installed by default even with the latest Windows 7 machines with many GB of memory. In general trying to launch with more memory than your OS/Java combination supports will result in the obscure error "could not create virtual machine".

9) SSH client - Mac and Linux users already have a command line ssh program that can be run from the terminal. For Windows users, please download PuTTY:http://www.chiark.greenend.org.uk/~sgtatham/putty/download.html

10) SCP/SFTP client - We will be moving data from the servers to the student laptops for visualization. Mac and Linux users already have a command line scp and sftp program. For Windows users, please install WinSCP: http://winscp.net/eng/download.php

11) A PDF viewer (Adobe Acrobat or equivalent).
