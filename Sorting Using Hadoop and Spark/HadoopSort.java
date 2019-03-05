import java.io.IOException;
import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.Path;
import org.apache.hadoop.io.LongWritable;
import org.apache.hadoop.io.Text;
import org.apache.hadoop.mapreduce.Job;
import org.apache.hadoop.mapreduce.Mapper;
import org.apache.hadoop.mapreduce.Reducer;
import org.apache.hadoop.mapreduce.lib.input.*;
import org.apache.hadoop.mapreduce.lib.output.*;

/*
 * @author: Raj Ambani - A20396925
 */
public class HadoopSort {

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		
		long startTime = System.currentTimeMillis();
		
		try 
		{
			Configuration conf = new Configuration();
			Job job = Job.getInstance(conf, "Hadoop TeraSort");
			
			//set input file path
			FileInputFormat.addInputPath(job, new Path(args[0]));
			//Set output file path
			FileOutputFormat.setOutputPath(job, new Path(args[1]));
			
			//Set main class
			job.setJarByClass(HadoopSort.class);
			//set mapper class
			job.setMapperClass(HadoopSortMapper.class);
			//set reducer class
			job.setReducerClass(HadoopSortReducer.class);
			
			//set class for key string
			job.setOutputKeyClass(Text.class);
			//set class for output string
			job.setOutputValueClass(Text.class);
			
			//Set input format.
			job.setInputFormatClass(TextInputFormat.class);
			//set output format
			job.setOutputFormatClass(TextOutputFormat.class);
			
			//set reducer task
			//job.setNumReduceTasks(100);
			
			//This will ensure that all the jobs have completed execution before the end time is calculated.
			if(job.waitForCompletion(true))
			{
				long totalTime = System.currentTimeMillis() - startTime;
				double timeTaken = ( (double) totalTime ) / 1000;
				System.out.println("Total time taken by HadoopSort to sort " + args[0] + " = " + timeTaken);
			}
			
		} 
		catch (IOException e) 
		{
			// TODO Auto-generated catch block
			e.printStackTrace();
		} 
		catch (ClassNotFoundException e) 
		{
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		catch (InterruptedException e) 
		{
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}
	
}

class HadoopSortMapper extends Mapper<LongWritable, Text, Text, Text>
{
	@Override
	protected void map(LongWritable key, Text value, Mapper<LongWritable, Text, Text, Text>.Context context)
			throws IOException, InterruptedException 
	{
		// TODO Auto-generated method stub
		
		//super.map(key, value, context);
		Text myKey = new Text();
		//Set first 10 bytes as key for sorting purpose.
		myKey.set(value.toString().substring(0, 10));
		
		Text myValue = new Text();
		//set bytes from 11 to 100 as value for that key.
		myValue.set(value.toString().substring(10));
		
		//Write key value pair
		context.write(myKey, myValue);
	}
}

class HadoopSortReducer extends Reducer<Text, Text, Text, Text>
{
	protected void reduce(Text key, Iterable<Text> values, Context con)
			throws IOException, InterruptedException 
	{
		// TODO Auto-generated method stub
		//super.reduce(arg0, arg1, arg2);
		
		Text myValue = new Text();
		for(Text val:values)
		{
			myValue = val;
		}
		//Reducing all different values to one specific key.
		con.write(key, myValue);
		
	}
}