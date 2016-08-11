#ifndef DEBUG_HH
#define DEBUG_HH

// C++ INCLUDES
#include <execinfo.h> // For complete stacktrace on exception see GNU libc manual
#include <unistd.h> // For gethostname
#include <sys/types.h> // For getpid
#include <iostream>
// Palabos INCLUDES
#include <core/plbLogFiles.h>

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
		catch(std::exception& e){
			std::string ex = e.what();
			std::string line = std::to_string(__LINE__);
			global::log().entry("[STACKTRACE]: "+ex+" [FILE:"+__FILE__+",LINE:"+line+"]");
			throw e; }
	}
}
#ifdef PLB_MPI_PARALLEL
void waitGDB(void){
	int i = 0;
	char hostname[256];
	gethostname(hostname, sizeof(hostname));
	const int pid = getpid();
	std::string mesg = "PID "+std::to_string(pid)+" on "+hostname+" ready for GDB attach";
	std::cout << mesg << std::endl;
	global::log().entry("[GDB]: "+mesg+ " [FILE:"+__FILE__+",LINE:"+std::to_string(__LINE__)+"]");
	fflush(stdout);
	unsigned int seconds;
	seconds = 10;
	while (0 == i){ }
}
#endif

template<typename T>
std::string adr_string(const T& var){
	const void * address = static_cast<const void*>(var);
	std::stringstream ss;
	ss << address;
	return ss.str();
}
}
#endif //DEBUG_HH
