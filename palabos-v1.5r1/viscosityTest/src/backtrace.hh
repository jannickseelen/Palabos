#ifndef BACKTRACE_HH
#define BACKTRACE_HH

// C++ INCLUDES
#include <execinfo.h> // For complete stacktrace on exception see GNU libc manual
#include <iostream>

namespace plb{
	void printTrace (void){			// Obtain a backtrace and print it to stdout.
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
}
#endif //BACKTRACE_HH
