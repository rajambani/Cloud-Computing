#!/bin/bash

#SBATCH --job-name="SparkSort20GB"
#SBATCH --output="SparkSort20GB.log"
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --partition=compute

rm *.class
javac -classpath /opt/spark-2.3.0-bin-hadoop2.7/jars/spark-core_2.11-2.3.0.jar:/opt/spark-2.3.0-bin-hadoop2.7/jars/spark-sql_2.11-2.3.0.jar:/opt/spark-2.3.0-bin-hadoop2.7/jars/scala-compiler-2.11.8.jar:/opt/spark-2.3.0-bin-hadoop2.7/jars/scala-library-2.11.8.jar SparkSort.java
jar cvf SparkSort.jar *.class

spark-submit --class SparkSort --master yarn --deploy-mode client --driver-memory 1g --executor-memory 1g --executor-cores 4 --num-executors 4 SparkSort.jar /input/data-20GB /user/rambani/output-spark

START_TIME=$SECONDS
hadoop jar /opt/hadoop-2.9.0/share/hadoop/mapreduce/hadoop-mapreduce-examples-2.9.0.jar teravalidate /user/rambani/output-spark /user/rambani/report-spark
ELAPSED_TIME=$(($SECONDS - $START_TIME))

hadoop fs -get /user/rambani/report-spark/part-r-00000

echo "Total time for TeraValidate: $ELAPSED_TIME seconds"

