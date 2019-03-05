JFLAGS = -g
JC = javac
J = jar
SB1 = sbatch hadoopsort8GB.slurm
SB2 = sbatch hadoopsort20GB.slurm
SB3 = sbatch hadoopsort80GB.slurm

SB4 = sbatch sparksort8GB.slurm
SB5 = sbatch sparksort20GB.slurm
SB6 = sbatch sparksort80GB.slurm

.SUFFIXES: .java .class
.java.class:
	$(JC) $*.java
hadoopsort:	
	$(SB1)
	$(SB2)
	$(SB3)

sparksort:	
	$(SB4)
	$(SB5)
	$(SB6)

CLASSES = \
        *.java 

default: classes

classes: $(CLASSES:.java=.class)

clean:
	$(RM) *.class
	$(RM) *.jar


