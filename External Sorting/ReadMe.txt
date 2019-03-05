ReadMe File

- For MySort:
	- navigate to folder ~/cs553-pa2a/
	- for cleaning all the class files, execute:
		make clean
	- For compilation and execution of the programs, just run below command:
		make mysort
		- When you execute above command, it will first compile java files, then submit slurm jobs for 2gb and 20gb input files respectively.
	- output will be stored in respective log files which are as follows:
		- mysort2GB.log : 2GB input file logs.
		- mysort20GB.log : 20GB input file logs.
		
- For running linux sort:
	- navigate to folder ~/cs553-pa2a/ and execute below command.
		- make linsort
	- Output:
		- output will be stored in files linsort2GB.log and linsort20GB.log respectively.