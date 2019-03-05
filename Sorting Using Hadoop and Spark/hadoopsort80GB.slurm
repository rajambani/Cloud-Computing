#!/bin/bash

#SBATCH --job-name="HadoopSort80GB"
#SBATCH --output="HadoopSort80GB.log"
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --partition=compute

javac -classpath $(hadoop classpath) HadoopSort.java
jar cvf HadoopSort.jar *.class
hadoop jar HadoopSort.jar HadoopSort /input/data-80GB /user/rambani/output-hadoop -D mapred.reduce.tasks=1

START_TIME=$SECONDS
hadoop jar /opt/hadoop-2.9.0/share/hadoop/mapreduce/hadoop-mapreduce-examples-2.9.0.jar teravalidate /user/rambani/output-hadoop /user/rambani/report-hadoop
ELAPSED_TIME=$(($SECONDS - $START_TIME))

hadoop fs -get /user/rambani/report-hadoop/part-r-00000

echo "Total time for TeraValidate: $ELAPSED_TIME seconds"

