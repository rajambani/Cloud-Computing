
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.net.DatagramPacket;
import java.net.DatagramSocket;
import java.net.InetAddress;
import java.net.SocketException;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;


public class MyNETBenchUDP 
{
	static int numberOfThreads;
	static String hostName;
	static int port = 5013;
	static double packetSize;
	//Fixed buffer size of 1GB
	final static double memorySize = 1000000000;
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
			int mode = Integer.parseInt(args[1]);

			FileReader fr1 = new FileReader(new File("server.txt"));
			BufferedReader br1 = new BufferedReader(fr1);
			hostName = br1.readLine().trim();
			br1.close();
			fr1.close();

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
				//hostName = args[2];
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

	//UDP Server code
	static void server(List<String> lines)
	{
		if(lines.get(0).equalsIgnoreCase("UDP"))
		{
			//Start Server.
			new UDPServer().start();
			//System.out.println("Server Started...");
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
			thread[i] = new Thread(new UDPClient(i));
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
		UDPClient udpClient = new UDPClient();
		double timerPerThread[] = udpClient.getTimePerThread();
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

		System.out.println("Time taken: " + (totalTime/1000) +" sec " + "Throughput: "+throughput + " MbPS");
		System.out.println("Latency: "+latency + " ms");
		
		//set variables to be written to file
		Protocol = "UDP";
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
			write.write("\n"+ Protocol +"\t\t\t"+ Concurrency +"\t\t\t\t"+ BlockSize +"\t\t\t"+ MyNetBenchValue +"\t\t"+ TheoreticalValue 
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

//Creating multiple threads for the server.
class UDPServer extends Thread
{
	DatagramSocket udpSocket = null;

	//Default Constructor for creating objects.
	public UDPServer() 
	{

	}

	@Override
	public void run() 
	{
		//MyNETBenchUDP.port = 5001;

		try
		{
			udpSocket = new DatagramSocket(MyNETBenchUDP.port);
			System.out.println("UDP Server started...");
			byte packetSize[] = new byte[(int) MyNETBenchUDP.packetSize];
			//Start receiving packets from client
			DatagramPacket packet = new DatagramPacket(packetSize, packetSize.length);
			udpSocket.receive(packet);
			InetAddress hostAddress = null;
			byte[] data = new byte[packet.getLength()];
			new Random().nextBytes(data);
			while(data != null)
			{
				data = packet.getData();
			}
			
			hostAddress = packet.getAddress();
			MyNETBenchUDP.port = packet.getPort();
			// Create new Packet to be sent
			packet = new DatagramPacket(data, data.length, hostAddress, MyNETBenchUDP.port);
			udpSocket.send(packet);

		}
		catch (SocketException e) 
		{
			e.printStackTrace();
		} 
		catch (IOException e) 
		{
			e.printStackTrace();
		} 
		finally 
		{
			udpSocket.close();
		}

	}
}

//Create multiple threads for UDP Client conenction.
class UDPClient extends Thread
{

	static double timePerThread[] = new double[MyNETBenchUDP.numberOfThreads];
	int currentThread;
	DatagramSocket udpSocket;

	//Default constructor.
	public UDPClient() 
	{

	}
	public UDPClient(int currentThread) 
	{
		this.currentThread = currentThread;
	}

	@Override
	synchronized public void run() 
	{
		try 
		{
			//Establish connection to server at specified port.
			//Create Socket
			udpSocket = new DatagramSocket();
			//Set host address.
			InetAddress hostAddress = InetAddress.getByName(MyNETBenchUDP.hostName);
			//InetAddress hostAddress = InetAddress.getByName("localhost");
			//MyNETBenchUDP.port = 5001;

			byte buffer[] = new byte[(int) MyNETBenchUDP.packetSize];
			new Random().nextBytes(buffer);
			/*//Loop for creating a packet of specified size.
			for(int i=0; i<MyNETBenchUDP.packetSize; i++)
			{
				buffer[i] = (byte) (Math.random() % MyNETBenchUDP.packetSize);
			}*/

			//Create UDP Packet
			DatagramPacket packet = new DatagramPacket(buffer, buffer.length,hostAddress,MyNETBenchUDP.port);
			

			//Calculate loop size as per the nuber of threads and the packetsize.
			double loopSize = MyNETBenchUDP.memorySize/(MyNETBenchUDP.packetSize*MyNETBenchUDP.numberOfThreads);
			double startTime, endTime;

			//For 1-byte ping pong
			if(MyNETBenchUDP.packetSize == 1)
			{
				startTime = System.currentTimeMillis();
				//Send and Receive data.
				for(double i=0; i<loopSize/1000; i++)
				{
					//Send data to server
					udpSocket.send(packet);
					//Receive data from server.
					//udpSocket.receive(packet);
					//System.out.println("i:"+i);
				}
				endTime = System.currentTimeMillis();
			}
			//For other conditions
			else
			{
				startTime = System.currentTimeMillis();
				//Send and Receive data.
				for(int j=0; j<100; j++) // Loop for 100x operations.
				{
					for(double i=0; i<loopSize; i++)
					{
						//Send data to server
						udpSocket.send(packet);
						//Receive data from server.
						//udpSocket.receive(packet);
						//System.out.println("i:"+i);
					}
				}
				endTime = System.currentTimeMillis();
			}

			double totalTime = endTime - startTime;
			timePerThread[currentThread] = totalTime;
		} 
		catch (IOException e)
		{
			e.printStackTrace();
		}
	}

	public double[] getTimePerThread()
	{
		return UDPClient.timePerThread;
	}

}
