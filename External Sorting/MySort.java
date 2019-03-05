//package TeraSort;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.TreeMap;
import java.util.TreeSet;

public class MySort 
{
	static String inputFile = "input/input1gb", outputFile = "output";
	static String tmpPath = "/tmp/temp/"; 
	static String outputPath = "/tmp/output/";
	static int fileCount = 0;
	static long inputFileSize;
	//static String outputPath = "C:\\Users\\Raj Ambani\\Desktop\\CS 553 Cloud Computing\\Assignment 2\\gensort-win-1.5\\";
	static long memorySize = 250000000; //1GB
	static int numOfThreads = 4;
	static long chunkSize;
	public static void main(String[] args) 
	{
		//For cluster
		if(args.length > 0)
			inputFile = args[0];

		/*System.out.println("Cleaning tmp directory...");
		File dir = new File(MySort.tmpPath);
		if(dir != null && dir.listFiles() !=null)
		{
			for(File file: dir.listFiles()) 
				if (!file.isDirectory()) 
					file.delete();
		}*/

		//Create directories.
		new File(tmpPath).mkdir();
		new File(outputPath).mkdir();

		long startTime = System.currentTimeMillis();
		createChunks();
		double totalTime = (double) System.currentTimeMillis() - startTime;

		System.out.println("Total time for external sort: "+ (double)(totalTime/1000) + " seconds");
		//Move output file to output folder.
		/*File f = new File("/tmp/output/output");
		if(f.renameTo(new File("~/cs553-pa2a/output/output")))
			System.out.println("Output file moved successfully...");
		else
			System.out.println("Output file movement failed");*/

	}

	static void createChunks()
	{
		File f = new File(inputFile);
		//int totalBytes = 0;
		//This will get chunk size.
		inputFileSize = f.length();
		chunkSize = getChunkSize(inputFileSize);

		System.out.println("Splitting and sorting input file "+ inputFile +" into chunks");
		FileReader fr;
		try 
		{
			fr = new FileReader(f);
			BufferedReader br = new BufferedReader(fr);
			String line = null;
			//int fileCount = 0;
			ArrayList<String> buff;
			//Read first line
			line = br.readLine();
			//System.out.println("line length: "+ line.length());
			double start = System.currentTimeMillis();
			while(line != null)
			{
				buff = new ArrayList<String>();
				long currentChunkSize = 0;
				while(line != null && ((currentChunkSize + line.length()) <= chunkSize))
				{
					//System.out.println("Line: "+ line);
					//totalBytes += line.length();
					buff.add(line);
					currentChunkSize += line.length() + 2;
					line = br.readLine();

				}
				//System.out.println("currentChunkSize: "+currentChunkSize);
				//System.out.println("------------Buff: "+ buff);
				writeToFile(buff, fileCount);	// Write current chunk to tmp file
				buff = null;	// freeing current buff
				fileCount += 1; // increementing tmpfile counter
			}
			double end = (System.currentTimeMillis() - start) / 1000;
			System.out.println("Total time to sort and create chunks: "+ end);

			br.close();
			fr.close();
			System.out.println("Number of temp files created: "+fileCount);
			//System.out.println("Chunks created with size: "+ chunkSize);
			//Start of sorting individual chunks.
			//startSorting();

			//Start of merging sorted files.
			//System.out.println("Starting merging of sorted files...");
			Thread t[] = new Thread[MySort.numOfThreads];
			MySort.numOfThreads = 1;
			for(int i=0; i<MySort.numOfThreads; i++)
			{
				System.out.println("Starting merging of sorted files...");
				t[i] = new Thread(new MyMergerLinux());
				t[i].start();
			}
			for(int i=0; i<MySort.numOfThreads; i++)
			{
				try 
				{
					t[i].join();
				} 
				catch (InterruptedException e)
				{
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
			System.out.println("Files are merged");
		} 
		catch (FileNotFoundException e)
		{
			// TODO Auto-generated catch block
			e.printStackTrace();
		} 
		catch (IOException e)
		{
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}
	//Returns chunk size
	static long getChunkSize(long fileSize)
	{
		//For cleaning any unused memory.
		System.gc();
		long availaleFreeMemory = Runtime.getRuntime().freeMemory();
		//if()
		//For cluster
		//memorySize = availaleFreeMemory/numOfThreads/2;
		//MySort.memorySize = availaleFreeMemory;
		System.out.println("File Size:"+ fileSize + " MemorySize: "+memorySize +" Concurrency: "+ numOfThreads);
		long numOfChunks = fileSize / memorySize;
		System.out.println("Available Free Memory: "+ availaleFreeMemory);
		long chunkSize;
		chunkSize = memorySize;
		System.out.println("Number of Chunks: "+ numOfChunks +" Chunk Size: "+chunkSize);
		return chunkSize;
	}

	//created files and write data to that file.
	static void writeToFile(ArrayList<String> buff, int fileCount)
	{
		//Creating temp files to store intermediate data.
		String fileName = MySort.tmpPath +"temp"+ fileCount;
		//System.out.println("Current file written: "+ fileName +":" + fileCount);

		//sort here
		MyMergeSortLinux sort = new MyMergeSortLinux();
		buff = sort.sortData(buff);

		File f = new File(fileName);
		try 
		{
			FileWriter fw = new FileWriter(f);
			BufferedWriter bw = new BufferedWriter(fw);
			for(String s:buff)
			{
				//System.out.println("s: "+s);
				bw.write(s);
				bw.newLine();
				/*if(s.contains(" "))
				{
					fw.write("\n");
				}*/
			}
			bw.flush();
			bw.close();
			//fw.flush();
			fw.close();
		}
		catch (IOException e) 
		{
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}

} // end of class

//MySort class is defined for performing sorting using multiple threads.
class MyMergeSortLinux extends Thread
{
	//File tmpFiles[];
	long start, end;
	//String tmpPath = "/tmp/";

	//Constructor.
	public MyMergeSortLinux() 
	{
		// TODO Auto-generated constructor stub
	}
	MyMergeSortLinux(long start, long end)
	{
		//this.tmpFiles = tmpFiles;
		this.start=start;
		this.end=end;
	}

	@Override
	public void run()
	{

	}



	public ArrayList<String> sortData(ArrayList<String> buff)
	{
		ArrayList<String> sortedBuff = new ArrayList<>();
		//Map<String, String> sortedMap = new TreeMap<>();
		TreeSet<String> sortedSet = new TreeSet<>();

		for(String s:buff)
		{
			if(s.length() > 10)
				sortedSet.add(s);
			//sortedMap.put(s.substring(0, 10), s.substring(10));
		}

		sortedBuff.addAll(sortedSet);
		/*Iterator<Map.Entry<String, String>> it = sortedMap.entrySet().iterator();
		while(it.hasNext())
		{
			Map.Entry<String, String> currentEntry = it.next();
			String key = currentEntry.getKey();
			sortedBuff.add(key+ "" + currentEntry.getValue());
		}*/
		return sortedBuff;
	}

}// end of class

class MyMergerLinux extends Thread
{
	//Constructor
	MyMergerLinux()
	{

	}

	@Override
	public void run()
	{
		// TODO Auto-generated method stub
		/*File f = new File(MySort.tmpPath);
		File tmpFiles[] = f.listFiles();*/
		//For cluster -6
		//int numOfFiles = tmpFiles.length;
		int numOfFiles = MySort.fileCount;
		BufferedReader br[] = new BufferedReader[numOfFiles];
		//long totalBytes = 0;
		try
		{
			//String lineOfEachFile[] = new String[numOfFiles];
			TreeMap<String, Integer> lineMap = new TreeMap<String, Integer>();
			for(int x=0; x<numOfFiles; x++)
			{
				//reading each sorted temp files in our code.
				FileReader fr = new FileReader(MySort.tmpPath + "temp" +x);
				br[x] = new BufferedReader(fr);
				//lineOfEachFile[x] = br[x].readLine();
				lineMap.put(br[x].readLine(), x);
				//br[i].close();
				//fr.close();

			}


			FileWriter fw = new FileWriter(MySort.outputPath + MySort.outputFile);
			BufferedWriter bw = new BufferedWriter(fw);

			//long loopSize = MySort.chunkSize / MySort.memorySize;
			long loopSize = MySort.inputFileSize/100;
			int i=0, minIndex=0;
			String min;
			//int temp=0;
			while(i < loopSize)
			{
				/*if(i%1000000 == 0)
					System.out.println("merging i: "+ i);*/
				min = (String)lineMap.firstKey();
				minIndex = (int)lineMap.get(min);

				//if(min != null)
				//{
				bw.write(min + "\r\n");
				//lineOfEachFile[minIndex] = br[minIndex].readLine();
				lineMap.remove(min);
				if((min = br[minIndex].readLine()) != null)
					lineMap.put(min, minIndex);
				/*}
				else
					lineOfEachFile[minIndex] = br[minIndex].readLine();*/
				++i;
			}
			bw.close();
			/*for(long i=0; i<loopSize; i++)
			{
				//Assuming first record is the smallest one.
				String initial = lineOfEachFile[0];
				long minIndex = 0;
				//Iterate over all the files.
				for(long j=0; j<numOfFiles; j++)
				{
					if(initial != null && lineOfEachFile[(int) j] != null)
					{
						//This loop will determine the smallest record out of the list present.
						if((initial.substring(0, 10).compareTo(lineOfEachFile[(int)j].substring(0, 10))) > 0)
						{
							initial = lineOfEachFile[(int) j];
							minIndex = j;
						}
					}
				}
				if(initial != null)
				{
					//Write the smallest record to output file.
					totalBytes += 1;
					bw.write(initial);
					bw.write("\r\n");
					//change the pointer to next line in current file.
					//bw.newLine();
				}

				lineOfEachFile[(int) minIndex] = br[(int) minIndex].readLine();
			}

			//Close all the read write buffers.
			for(long j=0; j<numOfFiles; j++)
			{
				//CLosing input buffer.
				br[(int) j].close();
			}
			//Closing output buffer.
			//bw.flush();
			bw.close();
			//closing output file.
			//fw.flush();
			fw.close();
			System.out.println("Total bytes merged: "+ totalBytes);*/
		}
		catch (Exception e) 
		{
			// TODO: handle exception
			e.printStackTrace();
		}
	}
}
