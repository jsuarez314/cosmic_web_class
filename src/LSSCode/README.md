# Hacking Xiao-Dong Li's LSS Code

```
  ____       _           _____ _        _      _              
 |  _ \     | |         / ____| |      | |    | |             
 | |_) | ___| |_ __ _  | (___ | | _____| | ___| |_ ___  _ __  
 |  _ < / _ \ __/ _` |  \___ \| |/ / _ \ |/ _ \ __/ _ \| '_ \ 
 | |_) |  __/ || (_| |  ____) |   <  __/ |  __/ || (_) | | | |
 |____/ \___|\__\__,_| |_____/|_|\_\___|_|\___|\__\___/|_| |_|
                                                              
```

This is a hacked version of the repo at
[https://github.com/xiaodongli1986/LSS_Code](https://github.com/xiaodongli1986/LSS_Code)
to calculate only the Beta-Skeleton of a given set of points.

Original version requires the "ifort" compiler (Intel - Fortran compiler).
This version can compile using "gfortran".

If the number of points to calculate the Beta Skeleton graph is big enough,
the program does a non-exact calculation of the graph.

You may be interested in the NGL library, it takes ~100 more time to do
the exact calculations.

##########################################################################

# Install:

Requirements: gfortran, Unix-like OS.
Clone the repo, go to the source directory (src/). Then just type "make".

      cd src
      make

It will create an excecutable file located in the binary folder (bin/)

It is recomendable to export the path of the binary file to the hidden 
file ".bashrc" in the home directory.
If HackingLSSCode is in the Home directory, it will works like:

   emacs ~/.bashrc

add the line

    export PATH=$PATH:$HOME/HackingLSSCode/bin

save, exit and load the new path by typing:

      source ~/.bashrc

      
# How to run:

Usage:

	./BSK_calc -input intpufilename -output outputfilname -beta beta -printinfo printinfo -numNNB numNNB 

Example:

	./bin/LSS_BSK_calc -input FC_N20000_REAL.cat -output test01 -beta 1.0 -printinfo True -numNNB 300

	input:	Mandatory. The input file to read and do calculations.
	output:	Mandatory. The output file to save the graph.
	beta:	Mandatory. real number >= 1.0
	printinfo: Optional. True or False
	numNNB: Mandatory. The number of neighbours to perform the search of the Beta-Neighbours.

Xiao Dong Li suggets:	

|beta   | numNNB |  Observations                           |
|------:|-------:|-----------------------------------------|
|1.0	|   100	 |I found 300 nummNNB gives best accuracy. |
|3.0	|    25  |                                         |
|5.0	|    20  |                                         |
|10.0	|    15  |                                         |