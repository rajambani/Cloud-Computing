JFLAGS = -g
JC = javac
SB1 = sbatch mysort2GB.slurm
SB2 = sbatch mysort20GB.slurm
SB3 = sbatch linsort2GB.slurm
SB4 = sbatch linsort20GB.slurm

.SUFFIXES: .java .class
.java.class:
	$(JC) $*.java
mysort:
	$(JC) *.java
	$(SB1)
	$(SB2)

linsort:
	$(SB3)
	$(SB4)

CLASSES = \
        MySort.java 

default: classes

classes: $(CLASSES:.java=.class)

clean:
	$(RM) *.class

