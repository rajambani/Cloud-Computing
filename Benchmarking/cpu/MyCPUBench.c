#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <pthread.h>
#include <sys/time.h>

const long long size = (1000000000000/6); //12 zeroes

int concurrency=0;
double myCPUBenchValue=0;
double theoreticalValue=1, efficiency=0;
char *file = NULL;
char type[10];
char *outputFileName = NULL;


//Function Prototypes
//Quarter Precision functions
void performQPOps(char *type);
void *charOps(void *arg);

//Half Precision functions
void performHPOps(char *type);
void *shortOps(void *arg);

//Single Precision functions
void performSPOps(char *type);
void *intOps(void *arg);

//Double Precision functions
void performDPOps(char *type);
void *doubleOps(void *arg);

//Writing results to file
void writeToFile();


int main(int argc, char *argv[]) 
{

	//Read input file
	FILE *fptr;
	file = argv[1];
	outputFileName = argv[2];
	char ch;
	char line[3];
	printf("Here it is: %s \n", argv[1]);
	outputFileName = argv[2];
	theoreticalValue = 2*1*2.3;
	
	/*  opening the file for reading */
	fptr = fopen(file, "r");
	if (fptr == NULL)
	{
		printf("Cannot open file \n");
		exit(0);
	}
	ch = fgetc(fptr);
	int cnt=0;
	while (ch != EOF)
	{
		printf ("%c", ch);
		//if(ch != ' ')		
		line[cnt] = ch;
		cnt++;
		ch = fgetc(fptr);
	}
	fclose(fptr);
	//printf("%c %c %c", line[0], line[1], line[3]);

	if(line[0] == 'Q' || line[0] == 'q')
	{
		performQPOps(line);
		theoreticalValue *= 128;
	}
	else if(line[0] == 'H' || line[0] == 'h')
	{
		performHPOps(&line[0]);
		theoreticalValue *= 64;
	}
	else if(line[0] == 'S' || line[0] == 's')
	{
		performSPOps(&line[0]);
		theoreticalValue *= 32;
	}
	else if(line[0] == 'D' || line[0] == 'd')
	{
		performDPOps(&line[0]);
		theoreticalValue *= 16;
	}

	//char type[0] = line[0], type[1]=line[1];
	concurrency = line[3]-48;
	type[0] = line[0];
	type[1] = line[1];
	//strcpy(type,);
	//strcat(type,line[1]);
	writeToFile();
	
    return 0;
}

//Perform Quarter Precision Operations.
void performQPOps(char *type)
{
	printf("%c %c %c", type[0], type[1], type[3]);
	struct timeval startTime, endTime;
	double totalTime,gigaOps;
	//Converting num of threads from char to int.
	int numOfThreads = type[3] - 48;
	//printf("\nTHreads:%d ", numOfThreads);
	int i;
	//Number of iterations
	double numOfOps= size ;//* 2ul;
	pthread_t *threads = (pthread_t*) malloc(numOfThreads * sizeof (pthread_t));
	
	//Fetch End time
	gettimeofday(&startTime, NULL);
	
	//Create multiple threads as enterd in the input file
	for(i=0; i<numOfThreads;i++)
	{
		//pthread_create(&threads[i],NULL,charOps,(void*)(int)numOfThreads);
		pthread_create(&threads[i], NULL, charOps, (void*)numOfThreads);
	}
	//This will ensure execution of all thread is completed before calculating timings.
	for(i=0; i<numOfThreads; i++)
	{
		pthread_join(threads[i], NULL);
	}
	
	//Fetch End time
	gettimeofday(&endTime, NULL);

	//total time in seconds
	totalTime = (double) (endTime.tv_sec - startTime.tv_sec) + (double) ((endTime.tv_usec - startTime.tv_usec) / 1000000);
	//number of giga flops
	gigaOps = (double) ((6*size) / (double) (totalTime * 1e9));
	myCPUBenchValue = gigaOps;
	
	printf("\n Size:%ld Threads:%d TotalTime:%.6fsec GigaOps:%f \n", (size*2), numOfThreads, totalTime, gigaOps);
	free(threads);

}
//This function will execute character operations for 1 trillion times.
void * charOps(void *arg)
{
	printf("\n executing QP");
	int nofOfThreads = (int)arg;
	long long i;
	char x1 = 'a';
	char x2 = 'z';
	char res;
	double loop = (double)(size/nofOfThreads);
	/* strong scaling in executing instructions*/
	for (i = 0; i < loop; i++) // 1 ops
	{
		x2-x1; // 1 ops
	}
	//printf("n i:%f ",i);
}

//Perform Half Precision Operations.
void performHPOps(char *type)
{
	printf("%c %c %c", type[0], type[1], type[3]);
	struct timeval startTime, endTime;
	double totalTime,gigaOps;
	//Converting num of threads from char to int.
	int numOfThreads = type[3] - 48;
	//printf("\nTHreads:%d ", numOfThreads);
	int i;
	//Number of iterations
	double numOfOps= size * 2ul;
	pthread_t *threads = (pthread_t*) malloc(numOfThreads * sizeof (pthread_t));
	
	//Fetch End time
	gettimeofday(&startTime, NULL);
	
	//Create multiple threads as enterd in the input file
	for(i=0; i<numOfThreads;i++)
	{
		//pthread_create(&threads[i],NULL,charOps,(void*)(int)numOfThreads);
		pthread_create(&threads[i], NULL, shortOps, (void*)numOfThreads);
	}
	//This will ensure execution of all thread is completed before calculating timings.
	for(i=0; i<numOfThreads; i++)
	{
		pthread_join(threads[i], NULL);
	}
	
	//Fetch End time
	gettimeofday(&endTime, NULL);

	//total time in seconds
	totalTime = (double) (endTime.tv_sec - startTime.tv_sec) + (double) ((endTime.tv_usec - startTime.tv_usec) / 1000000);
	//number of giga flops
	gigaOps = (double) ((6*size) / (double) (totalTime * 1e9));
	myCPUBenchValue = gigaOps;
	
	printf("\n Size:%ld Threads:%d TotalTime:%.6fsec GigaOps:%f \n", (size*2), numOfThreads, totalTime, gigaOps);
	free(threads);
}
//This function will execute short type operations for 1 trillion times.
void * shortOps(void *arg)
{
	int nofOfThreads = (int)arg;
	long int i;
	short x1 = 10;
	short x2 = 20;
	short res;
	double loop = (double)(size/nofOfThreads);
	/* strong scaling in executing instructions*/
	for (i = 0; i < loop; i++) 
	{
		x2-x1;
	}
	//printf("n i:%f ",i);
}

//Perform Single Precision Operations.
void performSPOps(char *type)
{
	printf("%c %c %c", type[0], type[1], type[3]);
	struct timeval startTime, endTime;
	double totalTime,gigaOps;
	//Converting num of threads from char to int.
	int numOfThreads = type[3] - 48;
	//printf("\nTHreads:%d ", numOfThreads);
	int i;
	//Number of iterations
	double numOfOps= size * 2ul;
	pthread_t *threads = (pthread_t*) malloc(numOfThreads * sizeof (pthread_t));
	
	//Fetch End time
	gettimeofday(&startTime, NULL);
	
	//Create multiple threads as enterd in the input file
	for(i=0; i<numOfThreads;i++)
	{
		//pthread_create(&threads[i],NULL,charOps,(void*)(int)numOfThreads);
		pthread_create(&threads[i], NULL, intOps, (void*)numOfThreads);
	}
	//This will ensure execution of all thread is completed before calculating timings.
	for(i=0; i<numOfThreads; i++)
	{
		pthread_join(threads[i], NULL);
	}
	
	//Fetch End time
	gettimeofday(&endTime, NULL);

	//total time in seconds
	totalTime = (double) (endTime.tv_sec - startTime.tv_sec) + (double) ((endTime.tv_usec - startTime.tv_usec) / 1000000);
	//number of giga flops
	gigaOps = (double) ((6*size) / (double) (totalTime * 1e9));
	myCPUBenchValue = gigaOps;
	
	printf("\n Size:%ld Threads:%d TotalTime:%.6fsec GigaOps:%f \n", (size*2), numOfThreads, totalTime, gigaOps);
	free(threads);

}
//This function will execute integer type operations for 1 trillion times.
void * intOps(void *arg)
{
	int nofOfThreads = (int)arg;
	double i;
	int x1 = 10;
	int x2 = 20;
	int res;
	double loop = (double)(size/nofOfThreads);
	/* strong scaling in executing instructions*/
	for (i = 0; i < loop; i++) 
	{
		x2-x1;
	}
	//printf("n i:%f ",i);
}

//Perform Double Precision Operations.
void performDPOps(char *type)
{
	printf("%c %c %c", type[0], type[1], type[3]);
	struct timeval startTime, endTime;
	double totalTime,gigaOps;
	//Converting num of threads from char to int.
	int numOfThreads = type[3] - 48;
	//printf("\nTHreads:%d ", numOfThreads);
	int i;
	//Number of iterations
	double numOfOps= size * 2ul;
	pthread_t *threads = (pthread_t*) malloc(numOfThreads * sizeof (pthread_t));
	
	//Fetch End time
	gettimeofday(&startTime, NULL);
	
	//Create multiple threads as enterd in the input file
	for(i=0; i<numOfThreads;i++)
	{
		//pthread_create(&threads[i],NULL,charOps,(void*)(int)numOfThreads);
		pthread_create(&threads[i], NULL, doubleOps, (void*)numOfThreads);
	}
	//This will ensure execution of all thread is completed before calculating timings.
	for(i=0; i<numOfThreads; i++)
	{
		pthread_join(threads[i], NULL);
	}
	
	//Fetch End time
	gettimeofday(&endTime, NULL);

	//total time in seconds
	totalTime = (double) (endTime.tv_sec - startTime.tv_sec) + (double) ((endTime.tv_usec - startTime.tv_usec) / 1000000);
	//number of giga flops
	gigaOps = (double) ((6*size) / (double) (totalTime * 1e9));
	myCPUBenchValue = gigaOps;
	
	printf("\n Size:%ld Threads:%d TotalTime:%.6fsec GigaOps:%f \n", (size*2), numOfThreads, totalTime, gigaOps);
	free(threads);

}
//This function will execute double type operations for 1 trillion times.
void * doubleOps(void *arg)
{
	int nofOfThreads = (int)arg;
	double i;
	double x1 = 10.010;
	double x2 = 20.015;
	double res;
	double loop = (size/nofOfThreads);
	/* strong scaling in executing instructions*/
	for (i = 0; i < loop; i++) 
	{
		x2-x1;
	}
	//printf("n i:%f ",i);
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
    	f = fopen("cpu_SP_1thread.out.dat", "a");
	//printf("\n %s", outputFileName);
	efficiency = myCPUBenchValue / theoreticalValue;
	fprintf(f, "%s\t\t%d\t\t%s\t\t%fGigaOps\t\t%fGigaOps\t\t%f %\n", type, concurrency,type, myCPUBenchValue, theoreticalValue, (efficiency*100));
	fclose(f);
}
