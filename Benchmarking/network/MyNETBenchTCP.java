//package NetworkBenchmark;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.net.ServerSocket;
import java.net.Socket;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

public class MyNETBenchTCP 
{
	static int numberOfThreads;
	static String hostName;
	static double packetSize;
	//Fixed buffer size of 1GB
	final static double memorySize = 1073741824; //1073741824
	static String Protocol;
	static double Concurrency;
	static double BlockSize;
	static double MyNetBenchValue;
	static double TheoreticalValue;
	static double MyNETBenchEfficiency;

	public static void main(String[] args) 
	{

		try 
		{
			FileReader fr = new FileReader(new File(args[0]));
			BufferedReader br = new BufferedReader(fr);

			FileReader fr1 = new FileReader(new File("server.txt"));
			BufferedReader br1 = new BufferedReader(fr1);
			hostName = br1.readLine().trim();
			br1.close();
			fr1.close();

			int mode = Integer.parseInt(args[1]);
			
			List<String> lines = new ArrayList<String>();

			String line = br.readLine();
			while (line != null) 
			{
				lines.add(line);
				//System.out.println(line);
				line = br.readLine();
			}
			System.out.println(lines.toString() + mode);

			switch(mode)
			{
			case 1:
				server(lines);
				break;
			case 2:
				client(lines);
			}
		}
		catch (FileNotFoundException e) 
		{
			e.printStackTrace();
		} 
		catch (IOException e) 
		{
			e.printStackTrace();
		}
	}

	//TCP Server code
	static void server(List<String> lines)
	{

		if(lines.get(0).equalsIgnoreCase("TCP"))
		{
			//Initialize TCP Server Socket
			ServerSocket tcpServer = null;
			try
			{
				//Start server at port 4091
				tcpServer = new ServerSocket(6061);
				System.out.println("Server Started...");
				while(true)
				{
					//Start accepting connections.
					new TCPServer(tcpServer.accept()).start();
				}
			}
			catch(Exception ex)
			{
				ex.printStackTrace();
			}
			finally
			{
				try 
				{
					//Close TCP server connection.
					tcpServer.close();
				} 
				catch (IOException e) 
				{
					e.printStackTrace();
				}
			}
		}
	}

	//TCP Client code
	static void client(List<String> lines)
	{
		//Fetch number of threads.
		numberOfThreads = Integer.parseInt(lines.get(2));
		//packet size
		packetSize = Integer.parseInt(lines.get(1));

		double iterations =  memorySize/packetSize/numberOfThreads;
		Thread thread[] = new Thread[numberOfThreads];

		//Create threads and send data in parallel.
		for(int i=0; i<numberOfThreads; i++)
		{
			thread[i] = new Thread(new TCPClient(i));
			thread[i].start();
		}

		try 
		{
			for(int j=0; j<numberOfThreads; j++)
				thread[j].join();
		} 
		catch (InterruptedException e) 
		{
			e.printStackTrace();
		}

		//This part is for calculating total time for all the threads.
		TCPClient tcpClient = new TCPClient();
		double timerPerThread[] = tcpClient.getTimePerThread();
		double totalTime = 0;;
		for(int i=0; i<timerPerThread.length; i++)
		{
			totalTime += timerPerThread[i];
		}

		System.out.println("Threads: "+numberOfThreads);
		double actualTime = totalTime / numberOfThreads;
		//double throughput = ((memorySize*100) / numberOfThreads / (1024*1024)) / (actualTime/1000);
		double throughput = (double) 8 *( ((double)(memorySize*100) / (1024*1024)) / (double)(totalTime/1000));
		//memorySize*100
		double latency = (totalTime / (memorySize*100));
		if(MyNETBenchTCP.packetSize == 1)
		{
			latency = (totalTime / (memorySize/1000));
		}

		System.out.println("Time taken: " + (totalTime/1000) +"sec " + "Throughput: "+throughput + "MbPS");
		System.out.println("Latency: "+latency + " ms");
		
		//set variables to be written to file
		Protocol = "TCP";
		Concurrency = numberOfThreads;
		BlockSize = packetSize;
		MyNetBenchValue = throughput;
		TheoreticalValue = 56000;
		MyNETBenchEfficiency = (double)((MyNetBenchValue/TheoreticalValue) * 100);

		if(packetSize == 1)
		{
			MyNetBenchValue = latency;
			TheoreticalValue = 0.0007;
			MyNETBenchEfficiency = (double)((MyNetBenchValue/TheoreticalValue) * 100);	
		}
		
		writeToFile();
	}
	
	//Method to write output to file.
	static void writeToFile()
	{
		try 
		{
			String basePath = new File("").getAbsolutePath();
		    System.out.println(basePath);
			FileWriter fw = new FileWriter(new File("network-TCP-1-1thread.out.dat"), true);
			BufferedWriter write = new BufferedWriter(fw);
			write.write("\n"+ Protocol +"\t\t\t"+ Concurrency +"\t\t\t\t"+ BlockSize +"\t\t"+ MyNetBenchValue +"\t\t"+ TheoreticalValue 
					+"\t\t\t"+ MyNETBenchEfficiency +"%");
			write.close();
			fw.close();
		} 
		catch (IOException e) 
		{
			e.printStackTrace();
		}
	}
}

class TCPServer extends Thread
{
	Socket tcpClient = null;
	String readLine;

	public TCPServer(Socket tcpClient) 
	{
		this.tcpClient = tcpClient;
	}

	@Override
	public void run() 
	{
		System.out.println("Connection Established");
		try 
		{
			PrintWriter pw = new PrintWriter(tcpClient.getOutputStream(), true);
			BufferedReader br = new BufferedReader(new InputStreamReader(tcpClient.getInputStream()));
			if(MyNETBenchTCP.packetSize == 1)
			{
				while((readLine = br.readLine()) != null)
				{
					//write data to client
					pw.println(readLine);
					//pw.flush();
				}
			}
			else
			{
				while((readLine = br.readLine()) != null)
				{
					//Receiving data from client.
					//pw.println(readLine);
					br.readLine();
					//pw.flush();
				}
			}
			System.out.println("End of data receiving");
		}
		catch (IOException e)
		{
			e.printStackTrace();
		}
		finally
		{
			try 
			{
				//End client connection.
				tcpClient.close();
			}
			catch (IOException e) 
			{
				e.printStackTrace();
			}
		}
	}
}

class TCPClient implements Runnable
{

	static double timePerThread[] = new double[MyNETBenchTCP.numberOfThreads];
	int currentThread;
	Socket tcpSocket;

	//Default constructor.
	public TCPClient() 
	{

	}
	public TCPClient(int currentThread) 
	{
		this.currentThread = currentThread;
	}

	@Override
	synchronized public void run() 
	{
		try 
		{
			//Establish connection to server at specified port.
			tcpSocket = new Socket(MyNETBenchTCP.hostName, 6061);
			//tcpSocket = new Socket("localhost", 6061);

			PrintWriter pw = new PrintWriter(tcpSocket.getOutputStream(), true);
			BufferedReader br = new BufferedReader(new InputStreamReader(tcpSocket.getInputStream()));
			byte buffer[]=new byte[(int) MyNETBenchTCP.packetSize];

			//Loop for creating a packet of specified size.
			new Random().nextBytes(buffer);
			/*for(int i=0; i<MyNETBenchTCP.packetSize; i++)
			{
				buffer[i] = (byte) (Math.random() * MyNETBenchTCP.packetSize);
			}*/

			double loopSize = MyNETBenchTCP.memorySize/(MyNETBenchTCP.packetSize*MyNETBenchTCP.numberOfThreads);
			double startTime, endTime;

			//For Ping pong 1-byte
			if(MyNETBenchTCP.packetSize == 1)
			{
				startTime = System.currentTimeMillis();
				//Send and Receive data.
				long  i;
				for(i=0; i<loopSize/1000; i++)
				{
					//Send data
					pw.println(buffer);
					//receive data
					//br.readLine();
					//System.out.println("i: "+i + " loopSize:"+loopSize);
					//System.out.println("i:"+i);
				}
				//System.out.println("i:"+i);
				endTime = System.currentTimeMillis();
				double totalTime = endTime - startTime;
				timePerThread[currentThread] = totalTime;
			}
			
			else
			{
				startTime = System.currentTimeMillis();
				//Send and Receive data.
				long  i;
				for(int j=0; j<100; j++) // Loop for 100x operations.
				{
					for( i=0; i<loopSize; i++)
					{
						//Send data
						pw.println(buffer);
						//pw.flush();
						//receive data
						//br.readLine();
						//System.out.println("i: "+i + " loopSize:"+loopSize);
					}
					//System.out.println("Before readline client");
					//Receive ACK.
					//br.readLine();
					//System.out.println("after readline client");
				}
				endTime = System.currentTimeMillis();
				double totalTime = endTime - startTime;
				timePerThread[currentThread] = totalTime;
			}
		
	} 
	catch (IOException e)
	{
		e.printStackTrace();
	}
	finally
	{
		try 
		{
			//Close TCP Connection here.
			tcpSocket.close();
		} 
		catch (IOException e) 
		{
			e.printStackTrace();
		}
	}
}

	public double[] getTimePerThread()
	{
		return TCPClient.timePerThread;
	}
	

}

