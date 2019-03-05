import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import scala.Tuple2;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.api.java.function.FlatMapFunction;
import org.apache.spark.api.java.function.PairFunction;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/*
 * @author: Raj Ambani - A20396925
 */

public class SparkSort 
{
	public static void main(String[] args) 
	{
		double startTime = System.currentTimeMillis();
		//Initialize Java Spark context.
		JavaSparkContext context = new JavaSparkContext("local[*]", "SparkSort");
		//JavaSparkContext context = new JavaSparkContext("local[*]", "SparkSort");

		//Read input file by using command line arguments.
		JavaRDD<String> fileLine = context.textFile(args[0]);
		//System.out.println("Test Here file read");

		//map records to key value pair.
		JavaPairRDD<String, String> recordPair = fileLine.mapToPair(new KeyValuePair());

		//Sort records based on their keys.
		boolean sortingType = true;
		JavaPairRDD<String, String> sortedPairs = recordPair.sortByKey(sortingType);

		//Creating anonymous class to sort the records based on their keys.
		//Use of FlatMap to combine key and values and then save to output file
		JavaRDD<String> finalOutput = sortedPairs.flatMap(new FlatMapFunction<Tuple2<String, String>, String>() 
		{
			private static final long serialVersionUID = -123479245454L;
			public Iterator<String> call(Tuple2<String, String> tuple) throws Exception 
			{
				List<String> recordList = new ArrayList<String>();
				//System.out.println("Test Here 1");
				String record = tuple._1() + "  " + tuple._2().trim() + "\t";
				//recordList.add(record);
				recordList.add(record);
				return recordList.iterator();
			}
		});

		//System.out.println("-----------------------------Start of saving output file----------------------------");
		
		//finalOutput.reduce();
		
		
		//Save output to text file.
		finalOutput.saveAsTextFile(args[1]);
		
		double endTime = System.currentTimeMillis();
		double timeTaken = (endTime - startTime)/1000;

		System.out.println("File line size: "+ fileLine.count());
		System.out.println("----------------------------------------Total time taken by my Spark Sort implementation: "+ timeTaken);
		
		//Closing SparkContext object.
		context.close();
	}
}

//This class is mainly for creating key value pairs
class KeyValuePair implements PairFunction<String, String, String> {

	private static final long serialVersionUID = -1L;

	@Override
	public Tuple2<String, String> call(String record) throws Exception 
	{
		// TODO Auto-generated method stub
		//System.out.println("Test Here inside key value pair class");
		String key = record.substring(0, 10);
		String value = record.substring(10);
		return new Tuple2<String, String>(key, value);
	}


}