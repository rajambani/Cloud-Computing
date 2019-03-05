#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <pthread.h>
#include <math.h>
#include <sys/time.h>

char *file = NULL;
char *outputFileName = NULL;
char type[5];
long int memorySize = 1000000000; //1GB of data
long int concurrency,BlockSize;
double MyRAMBenchValue, TheoreticalValue = 68.256, theoreticalLatency=0.01406, Latency;
//Workload Concurrency BlockSize MyRAMBenchValue TheoreticalValue MyRAMBenchEfficiency


//Function Prototypes are defined here.
//Functions for Read-Write random access
void performRWR(char line[10][10]);
void * randomOperation(void *threadDetails);

//Functions for Read-Write in sequential manner.
void performRWS(char line[10][10]);
void * sequentialOperation(void *threadDetails);

//Function to write output to file.
void writeToFile();

/*structure for maintaining current thread details*/
typedef struct Thread_Details 
{
	int startOffset,endOffset;
	int threadsize;
	long int blocksize;
} Thread_Details;


//FILL ME IN
int main(int argc, char *argv[]) 
{

    	//Read input file
	FILE *fptr;
	file = argv[1];
	outputFileName = argv[2];
	char ch;
	char line[10][10];// = (char *)malloc(sizeof(char) * 10);
	printf("Input file name: %s", argv[1]);
	
	/*  opening the file for reading */
	fptr = fopen(file, "r");
	if (fptr == NULL)
	{
		printf("Cannot open file \n");
		exit(0);
	}
	ch = fgetc(fptr);
	int lineCount=0, wordCount;
	char *word;// = (char *)malloc(sizeof(char) * 10);
	while (ch != EOF)
	{
		word = (char *)malloc(sizeof(char) * 10);
		wordCount=0;
		while(ch != EOF && ch != '\n')
		{
			word[wordCount] = ch;
        		wordCount++;
			ch = fgetc(fptr);
			//strcat(word,ch);
		}
		strcpy(line[lineCount],word);
		//printf("\n%s %s %d %d", line[lineCount], word, lineCount, wordCount);
		//free(word);
		lineCount++;
		ch = fgetc(fptr);
	}
	fclose(fptr);
	//printf("\n%s", line[2]);
	strcpy(type,line[0]);

	//Perform Read-Write random
	if(strcmp(line[0], "RWR") == 0)
	{
		performRWR(line);	
	}
	//Perform Read-Write Sequential
	else
	{
		performRWS(line);
	}

	return 0;
}

//Function for calculating Read-Write random access.
void performRWR(char line[10][10])
{
	//printf("\n%s", line[2]);
	int numOfThreads = atoi(line[2]);
	concurrency = numOfThreads;
	long int blockSize = atoi(line[1]);
	//parameter to be written in file.
	BlockSize = blockSize;
	long int loopSize = memorySize/blockSize;
	
	//Creating thread array.
	pthread_t threads[numOfThreads];
	//thread structure variable used to manintain thread related information like offset and other details.
	Thread_Details *threadDetails = (Thread_Details*) malloc(numOfThreads * sizeof (Thread_Details));

	//timer variables.
	struct timeval startTime, endTime;
	double totalTime, throughput, latency;

	printf("\n%d %d ", numOfThreads, blockSize);
	int i;
	//int   tid;

	//Start system timer.
	gettimeofday(&startTime, NULL);

	for(i=0; i<numOfThreads; i++)
	{
		//pthread_id_np_t   tid;
		pthread_t tid = pthread_self();
		printf("\nThread id: %d", tid);
		//current thread starting offset address.
		threadDetails[i].startOffset = i * (loopSize / numOfThreads);
		//current thread ending offset address.
		threadDetails[i].endOffset = (i+1) * (loopSize / numOfThreads) - 1;
		threadDetails[i].blocksize = blockSize;
		threadDetails[i].threadsize = numOfThreads;
		//Creating threads.
		pthread_create(&threads[i], NULL, randomOperation, (void*)&threadDetails[i]);
	}
	for(i=0; i<numOfThreads; i++)
	{
		//Thread join ensures that all the thread execution has completed and can move to next processing.
		pthread_join(threads[i], NULL);
	}
	//End System timer.
	gettimeofday(&endTime, NULL);
	//Total time in seconds.
	totalTime = (double) (endTime.tv_sec - startTime.tv_sec) + (double)(endTime.tv_usec - startTime.tv_usec) / 1000000 ;
	//Calculate throughput
	//memorySize *= 100;
	throughput = (((double)( memorySize*100)/1073741824) / (totalTime )); // *100 1048576Mbps 1073741824Gbps
	//Calculate Latency in microseconds
	latency = ((totalTime)/((memorySize*100)/blockSize))*1000; //*100 microseconds
	MyRAMBenchValue = throughput;

	// For i-byte latency
	if(blockSize == 1)
	{
		//printf("\n Inside if");
		//Total time in microseconds.
		totalTime = (double)(endTime.tv_sec - startTime.tv_sec)*1000000 + (double) (endTime.tv_usec - startTime.tv_usec) ;
		//Calculate throughput
		//memorySize *= 100;
		throughput = (((double)( memorySize/10)/1073741824) / (double) (totalTime )); // *100 1048576Mbps 1073741824Gbps
		//Calculate Latency in microseconds
		latency = (double)((totalTime)/(double)((memorySize/10)/blockSize)); // /10 for bytes
		MyRAMBenchValue = throughput;	
		Latency = latency;
	}
	printf("\n%i\t\t%d\t%s\t%f\t%f\t%f\n", numOfThreads,blockSize, "RWR", totalTime,throughput,latency);
	writeToFile();
	free(threadDetails);	
}

//Function for performing random read and write operations using single/multiple threads.
void * randomOperation(void *threadDetails)
{
	Thread_Details *threadInfo = (Thread_Details *) threadDetails;
	long int blockSize = threadInfo -> blocksize;
	long int loopSize = memorySize/blockSize;
	long int blocksPerThread = loopSize/(threadInfo->threadsize);
	long int randomOffset,i;
	srand((unsigned)time(NULL));
	int start, end, j;
	
	// For 1-byte latency
	if(blockSize == 1)
	{	
		char *source = (char *) malloc(memorySize/concurrency);
		char *destination =  (char *) malloc(memorySize/concurrency);
		start = threadInfo->startOffset;
		end = threadInfo->endOffset;
		memset(source, '.', blockSize);
		for(i=start; i<=end/10; i++)
		{
			//Do random offset selection.
			randomOffset = rand()%blocksPerThread + (start);
			//Read and Write memory at some random offset.
			memcpy(destination+(randomOffset*blockSize), source, blockSize);
		}
		free(source);
		free(destination);
		return;
	}
	
	/*Used to set the starting value(seed), so time(NULL) will make use of 
	 *computer's internal clock to control the choice of the seed.
	 */
	for(j=0; j<100; j++)	
	{
		char *source = (char *) malloc(memorySize/concurrency);
		char *destination =  (char *) malloc(memorySize/concurrency);
		start = threadInfo->startOffset;
		end = threadInfo->endOffset;
		memset(source, '.', blockSize);
		for(i=start; i<=end; i++)
		{
			//Do random offset selection.
			randomOffset = rand()%blocksPerThread + (start);
			//Read and Write memory at some random offset.
			memcpy(destination+(randomOffset*blockSize), source+(randomOffset*blockSize), blockSize);
		}
		free(source);
		free(destination);
	}
	//free(source);
	//free(destination);
}

//Function for calculating Read-Write sequential access.
void performRWS(char line[10][10])
{
	//printf("\n%s", line[2]);
	int numOfThreads = atoi(line[2]);
	concurrency = numOfThreads;
	long int blockSize = atoi(line[1]);
	//parameter to be written in file.
	BlockSize = blockSize;
	long int loopSize = memorySize/blockSize;
	
	//Creating thread array.
	pthread_t threads[numOfThreads];
	//thread variable used to manintain thread related information like offset and other details.
	Thread_Details *threadDetails = (Thread_Details*) malloc(numOfThreads * sizeof (Thread_Details));

	//timer variables.
	struct timeval startTime, endTime;
	double totalTime, throughput, latency;
	//Start system timer.
	gettimeofday(&startTime, NULL);

	printf("\n %d %d ", numOfThreads, blockSize);
	int i;
	for(i=0; i<numOfThreads; i++)
	{
		//current thread starting offset address.
		threadDetails[i].startOffset = i * (loopSize / numOfThreads);
		//current thread ending offset address.
		threadDetails[i].endOffset = (i+1) * (loopSize / numOfThreads) - 1;
		threadDetails[i].blocksize = blockSize;
		threadDetails[i].threadsize = numOfThreads;
		//Creating threads.
		pthread_create(&threads[i], NULL, sequentialOperation, (void*)&threadDetails[i]);
	}
	for(i=0; i<numOfThreads; i++)
	{
		//Thread join ensures that all the thread execution has completed and can move to next processing.
		pthread_join(threads[i], NULL);
	}
	//End System timer.
	gettimeofday(&endTime, NULL);
	
	//For 1-byte latency
	if(blockSize == 1)
	{
		//printf("\n Inside if RWS");
		//Total time in seconds.
		totalTime = (double) ((endTime.tv_sec - startTime.tv_sec)*1000000) + (double) (endTime.tv_usec - startTime.tv_usec) ;
		//Calculate throughput
		//memorySize *= 100;
		throughput = (((double)( memorySize/10)/1073741824) /(double) (totalTime )); // *100 1048576Mbps 1073741824Gbps
		//Calculate Latency in microseconds
		latency = (double)((totalTime)/(double)((memorySize/10)/blockSize)); // /10 for bytes
		MyRAMBenchValue = throughput;	
		Latency = latency;
	}
	else
	{
		//Total time in seconds.
		totalTime =(double)(endTime.tv_sec - startTime.tv_sec) + (double)(endTime.tv_usec - startTime.tv_usec) / 1000000 ;
		//Calculate throughput
		throughput = (double)((( memorySize*100)/1073741824) /(double) (totalTime )); //1073741824 GBPS
		//Calculate Latency in microseconds
		latency = (double)((totalTime)/(double)((memorySize*100)/blockSize))*1000; 
		MyRAMBenchValue = throughput;
	}

	printf("\n%i\t\t%d\t%s\t%f\t%f\t%f\n", numOfThreads,blockSize, "RWS", totalTime,throughput,latency);
	writeToFile();
	free(threadDetails);	
}

//Function for performing sequential read and write operations using single/multiple threads.
void * sequentialOperation(void *threadDetails)
{
	Thread_Details *threadInfo = (Thread_Details *) threadDetails;
	long int blockSize = threadInfo -> blocksize;	
	int start, end, j, i;
	
	// For 1-byte latency
	if(blockSize == 1)
	{	
		//printf("\n Inside if");
		char *source = (char *) malloc(memorySize/concurrency);
		char *destination =  (char *) malloc(memorySize/concurrency);
		start = threadInfo->startOffset;
		end = threadInfo->endOffset;
		memset(source, '.', blockSize);
		for(i=start; i<=(end/10); i++)
		{
			//Read and Write memory at sequential offset.
			memcpy(destination+(i*blockSize), source+(i*blockSize), blockSize);
		}
		free(source);
		free(destination);
		return;
	}

	char *source = (char *) malloc(memorySize);
	char *destination =  (char *) malloc(memorySize);
	memset(source, 'a', blockSize);
	for(j=0; j<100; j++)	
	{
		
		start = threadInfo->startOffset;
		end = threadInfo->endOffset;
		for(i=start; i<=end; i++)
		{
			//Read and Write memory at sequential offset.
			memcpy(destination+(i*blockSize), source+(i*blockSize), blockSize);
		}
	}
	free(source);
	free(destination);
}

void writeToFile()
{
	FILE *f;
	char cwd[256];
	if (getcwd(cwd, sizeof(cwd)) == NULL)
	      perror("\ngetcwd() error");
	else
	      printf("\ncurrent working directory is: %s\n", cwd);

	strcat(cwd,"/output");
	chdir(cwd);
	getcwd(cwd, sizeof(cwd));
	printf("\ncurrent working directory is: %s\n", cwd);
    	f = fopen("memory-RWR-1-1thread.out.dat", "a");
	printf("\n %s\n", outputFileName);
	double MyRAMBenchEfficiency = MyRAMBenchValue / TheoreticalValue;

	if(BlockSize == 1)
	{
		MyRAMBenchValue = Latency;
		TheoreticalValue = theoreticalLatency;
		MyRAMBenchEfficiency = (theoreticalLatency-Latency) / theoreticalLatency;
	}
	
	//printf(f, "%ld, %ld, %d, %f, %lf, %f\n", memorySize, concurrency, BlockSize, MyRAMBenchValue,TheoreticalValue, MyRAMBenchEfficiency);
	fprintf(f, "%s\t\t%d\t\t%ld\t\t%lf\t\t%lf\t\t%lf\n", type, concurrency, BlockSize, MyRAMBenchValue,TheoreticalValue, (MyRAMBenchEfficiency*100));
	fclose(f);
}

