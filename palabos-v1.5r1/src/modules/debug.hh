#ifndef DEBUG_HH
#define DEBUG_HH

// C++ INCLUDES
#include <execinfo.h> // For complete stacktrace on exception see GNU libc manual
#include <unistd.h> // For gethostname
#include <sys/types.h> // For getpid
#include <iostream>


namespace plb{
	void printTrace (void){			// Obtain a stacktrace and print it to stdout.
	if(global::mpi().isMainProcessor()){
		try{
			void *Array[10];
			int size;
			char **strings;
			int i;
			size = backtrace (Array, 10);
			strings = backtrace_symbols (Array, size);
			std::cout << "Obtained "<< size << "stack frames" << std::endl;
			for (i = 0; i < size; i++){printf ("%s\n", strings[i]);}
			free (strings);
		}
		catch(std::exception& e){ std::cout << "Exception Caught: " << e.what() << "\n"; throw; }
		}
	}

	void waitGDB(void){
		int i = 0;
		char hostname[256];
		gethostname(hostname, sizeof(hostname));
		printf("PID %d on %s ready for attach\n", getpid(), hostname);
		fflush(stdout);
		unsigned int seconds;
		seconds = 10;
		while (0 == i){ usleep(seconds); }
	}
}
#endif //DEBUG_HH
