#ifndef DEBUG_HH
#define DEBUG_HH

// C++ INCLUDES
#include <execinfo.h> // For complete stacktrace on exception see GNU libc manual
#include <unistd.h> // For gethostname
#include <sys/types.h> // For getpid
#include <iostream>
#include <sstream>
#include <signal.h>
#include <string>
#include <sys/time.h>
#include <sys/resource.h>



// Palabos INCLUDES
#include <core/plbLogFiles.h>

namespace plb{

void printTrace (void){			// Obtain a stacktrace and print it to stdout.
try{
	int length = 100;
	void *array[length];
	size_t size;

	//char **strings;
	//int i;
	size = backtrace(array, length);

	fprintf(stderr, "[ERROR]: EXCEPTION STACKTRACE \n");
	backtrace_symbols_fd(array, size, STDERR_FILENO);

	global::log("[ERROR]: EXCEPTION STACKTRACE \n");

	void* const* buffer = nullptr;
	char** trace;
	trace = backtrace_symbols(buffer, size);
	for(int i = 0; i<length; i++){
		global::log(std::to_string(**trace));
	}
	free(trace);
}
catch(const std::exception& e){
	std::string ex = e.what();
	std::string line = std::to_string(__LINE__);
	global::log().entry("[STACKTRACE]: "+ex+" [FILE:"+__FILE__+",LINE:"+line+"]");
	throw e; }
}

void exHandler(const std::exception& e, const std::string& file, const std::string& function, const int& line){
	std::string location = "[FILE: "+file+", FUNC: "+function+", LINE: "+std::to_string(line)+"]";
	std::string ex = e.what();
	std::string mesg = "[EXCEPTION] "+ ex + location;
	global::log(mesg);
	std::cerr << mesg << std::endl;
	printTrace();
	throw e;
}

void memUsage(){
	try{
		struct rusage usage;
		int answer = getrusage(RUSAGE_SELF, &usage);
		if(answer == 0){
			std::string mesg;
			mesg = std::to_string(usage.ru_utime.tv_sec) + " user CPU time used in seconds";
			global::log(mesg);
			pcout << mesg << std::endl;
			mesg = std::to_string(usage.ru_stime.tv_sec) + " system CPU time used in seconds";
			global::log(mesg);
			pcout << mesg << std::endl;
			mesg = std::to_string(usage.ru_maxrss) + " maximum resident set size";
			global::log(mesg);
			pcout << mesg << std::endl;
			mesg = std::to_string(usage.ru_ixrss) + " integral shared memory size";
			global::log(mesg);
			pcout << mesg << std::endl;
			mesg = std::to_string(usage.ru_idrss) + " integral unshared data size";
			global::log(mesg);
			pcout << mesg << std::endl;
			mesg = std::to_string(usage.ru_isrss) + " integral unshared stack size";
			global::log(mesg);
			pcout << mesg << std::endl;
			mesg = std::to_string(usage.ru_minflt) + " page reclaims (soft page faults)";
			global::log(mesg);
			pcout << mesg << std::endl;
			mesg = std::to_string(usage.ru_majflt) + " page faults (hard page faults)";
			global::log(mesg);
			pcout << mesg << std::endl;
			mesg = std::to_string(usage.ru_nswap) + " swaps";
			global::log(mesg);
			pcout << mesg << std::endl;
			mesg = std::to_string(usage.ru_inblock) + " block input operations";
			global::log(mesg);
			pcout << mesg << std::endl;
			mesg = std::to_string(usage.ru_oublock) + " block output operations";
			global::log(mesg);
			pcout << mesg << std::endl;
			mesg = std::to_string(usage.ru_msgsnd) + "IPC messages sent";
			global::log(mesg);
			pcout << mesg << std::endl;
			mesg = std::to_string(usage.ru_msgrcv) + "IPC messages received";
			global::log(mesg);
			pcout << mesg << std::endl;
			mesg = std::to_string(usage.ru_nsignals) + " signals received";
			global::log(mesg);
			pcout << mesg << std::endl;
			mesg = std::to_string(usage.ru_nvcsw) + " voluntary context switches";
			global::log(mesg);
			pcout << mesg << std::endl;
			mesg = std::to_string(usage.ru_nivcsw) + " involuntary context switches";
			global::log(mesg);
			pcout << mesg << std::endl;
		}
	}
	catch(const std::exception& e){exHandler(e,__FILE__,__FUNCTION__,__LINE__);}
}

void sigHandler(int sig) {
	try{
		int length = 100;
		void *array[length];
		size_t size;

		// get void*'s for all entries on the stack
		size = backtrace(array, length);

		std::string mesg = "";
		if(sig > 0 && sig < 32){
		switch(sig){
			case 1: mesg = "[ERROR]: (Signal Hangup) Report that user's terminal is disconnected.\n Signal used to report the termination of the controlling process. \n";
			case 2: mesg = "[ERROR]: (Signal Interrupt) Program interrupted \n";
			case 3: mesg = "[ERROR]: (Signal Quit) Terminate process and generate core dump \n";
			case 4: mesg = "[ERROR]: (Signal Illegal Instruction) Generally indicates that the executable file \n is corrupted or use of data where a pointer to a function was expected. \n";
			case 5: mesg = "[ERROR]: (Signal Trace Trap) \n";
			case 6: mesg = "[ERROR]: (Signal Abort) Process detects error and reports by calling abort. \n";
			case 7: mesg = "[ERROR]: (Signal BUS error) Indicates an access to an invalid address. \n";
			case 8: mesg = "[ERROR]: (Signal Floating-Point Exception) Erroneous arithmetic operation, \n such as zero divide or an operation resulting in overflow (not necessarily with a floating-point operation). \n";
			case 9: mesg = "[ERROR]: (Signal Kill) Cause immediate termination \n";
			case 10: mesg = "[ERROR]: (Signal Userdefined 1) \n";
			case 11: mesg = "[ERROR]: (Signal Segmentation Violation) Invalid access to storage: \n When a program tries to read or write outside the memory it has allocated. \n";
			case 12: mesg = "[ERROR]: (Signal Userdefined 2) \n";
			case 13: mesg = "[ERROR]: (Signal Broken Pipe) Error condition like trying to write to a socket which is not connected. \n";
			case 14: mesg = "[ERROR]: (Signal Alarm Clock) Alarm clock (POSIX) Indicates expiration of a timer. Used by the alarm() function. \n";
			case 15: mesg = "[ERROR]: (Signal Termination) This signal can be blocked, handled, and ignored. Generated by kill command. \n";
			case 16:  mesg = "[ERROR]: (Signal Stack fault) \n";
			case 17: mesg = "[ERROR]: (Signal Child status has changed) Signal sent to parent process \n whenever one of its child processes terminates or stops. \n";
			case 18: mesg = "[ERROR]: (Signal	Continue) Signal sent to process to make it continue. \n";
			case 19: mesg = "[ERROR]: (Signal 	Stop, unblockable) Stop a process. This signal cannot be handled, ignored, or blocked. \n";
			case 20: mesg = "[ERROR]: (Signal Keyboard stop) Interactive stop signal. This signal can be handled and ignored. (ctrl-z) \n";
			case 21: mesg = "[ERROR]: (Signal Background read from tty) \n";
			case 22: mesg = "[ERROR]: (Signal Background write to tty) \n";
			case 23: mesg = "[ERROR]: (Signal Urgent condition on socket) Signal sent when urgent or out-of-band data arrives on a socket. \n";
			case 24: mesg = "[ERROR]: (Signal CPU limit exceeded) \n";
			case 25: mesg = "[ERROR]: (Signal size limit exceeded) \n";
			case 26: mesg = "[ERROR]: (Signal 	Virtual Time Alarm) Indicates expiration of a timer. \n";
			case 27: mesg = "[ERROR]: (Signal 	Profiling alarm clock) Indicates expiration of a timer. Use for code profiling facilities. \n";
			case 28: mesg = "[ERROR]: (Signal 	Window size change) \n";
			case 29: mesg = "[ERROR]: (Signal I/O now possible) Pollable event occurred (System V) \n Signal sent when file descriptor is ready to perform I/O (generated by sockets) \n";
			case 30: mesg = "[ERROR]: (Signal 	Power failure restart) \n";
			case 31: mesg = "[ERROR]: (Signal 	Bad system call) \n";
		}
		}
		else{  mesg = "[ERROR]: (Signal "+std::to_string(sig)+") \n";}
		// print out all the frames to stderr
		global::log(mesg);
		fprintf(stderr, mesg.c_str(), sig);

		memUsage();

		backtrace_symbols_fd(array, size, STDERR_FILENO);

		void* const* buffer = nullptr;
		char** trace;
		trace = backtrace_symbols(buffer, size);
		for(int i = 0; i<length; i++){
			global::log(std::to_string(**trace));
		}
		free(trace);
		exit(1);
	}
	catch(const std::exception& e){exHandler(e,__FILE__,__FUNCTION__,__LINE__);}
}

void installSigHandler(){
	try{
		#ifdef __APPLE__
			signal(SIGVTALRM, sigHandler);
			signal(SIGTSTP, sigHandler);
		#elif __linux__
			signal(SIGSTKFLT, sigHandler);
			signal(SIGPWR, sigHandler);
			signal(SIGVTALRM, sigHandler);
			signal(SIGTSTP, sigHandler);
		#else
			signal(SIGSTP, sigHandler);
			signal(SIGVALRM, sigHandler);
		#endif
		signal(SIGHUP, sigHandler);
		signal(SIGINT, sigHandler);
		signal(SIGQUIT, sigHandler);
		signal(SIGILL, sigHandler);
		signal(SIGTRAP, sigHandler);
		signal(SIGABRT, sigHandler);
		signal(SIGFPE, sigHandler);
		signal(SIGKILL, sigHandler);
		signal(SIGUSR1, sigHandler);
		signal(SIGSEGV, sigHandler);
		signal(SIGUSR2, sigHandler);
		signal(SIGPIPE, sigHandler);
		signal(SIGALRM, sigHandler);
		signal(SIGTERM, sigHandler);
		signal(SIGCHLD, sigHandler);
		signal(SIGCONT, sigHandler);
		signal(SIGSTOP, sigHandler);
		signal(SIGTTIN, sigHandler);
		signal(SIGTTOU, sigHandler);
		signal(SIGURG, sigHandler);
		signal(SIGXCPU, sigHandler);
		signal(SIGXFSZ, sigHandler);
		signal(SIGPROF, sigHandler);
		signal(SIGWINCH, sigHandler);
		signal(SIGIO, sigHandler);
		signal(SIGSYS, sigHandler);
		signal(SIGBUS, sigHandler);
		signal(SIGTRAP, sigHandler);
	}
	catch(const std::exception& e){exHandler(e,__FILE__,__FUNCTION__,__LINE__);}
}

#ifdef PLB_MPI_PARALLEL
void waitGDB(void){
	try{
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
	catch(const std::exception& e){exHandler(e,__FILE__,__FUNCTION__,__LINE__);}
}
#endif

template<typename T>
std::string adr_string(const T& var){
	std::string ans ="";
	try{
		const void * address = static_cast<const void*>(var);
		std::stringstream ss;
		ss << address;
		ans = ss.str();
	}
	catch(const std::exception& e){exHandler(e,__FILE__,__FUNCTION__,__LINE__);}
	return ans;
}

std::string safe_string(const T& val){
	std::string ans ="";
	try{
		if(std::isnan(val)){ ans ="NaN";}
		else if(val>std::numeric_limits<T>::max()){ans = "Inf";}
		else if(val<std::numeric_limits<T>::min()){ ans = "-Inf"; }
		else if(!val){	ans = "NULL";}
		else{ ans = std::to_string(val);}
	}
	catch(const std::exception& e){exHandler(e,__FILE__,__FUNCTION__,__LINE__);}
	return ans;
}

std::string array_string(const Array<T,3>& array){
	std::string ans = "";
	try{
		std::string a1 = safe_string(array[0]);
		std::string a2 = safe_string(array[1]);
		std::string a3 = safe_string(array[2]);
		ans = "["+a1+", "+a2+", "+a3+"]";
	}
	catch(const std::exception& e){exHandler(e,__FILE__,__FUNCTION__,__LINE__);}
	return ans;
}

std::string box_string(const Box3D& box){
	std::string ans = "";
	try{
		std::string a1 = safe_string(box.x0);
		std::string a2 = safe_string(box.x1);
		std::string a3 = safe_string(box.y0);
		std::string a4 = safe_string(box.y1);
		std::string a5 = safe_string(box.z0);
		std::string a6 = safe_string(box.z1);
		ans = "["+a1+", "+a2+"]["+a3+", "+a4+" ]["+a5+", "+a6+"]";
	}
	catch(const std::exception& e){exHandler(e,__FILE__,__FUNCTION__,__LINE__);}
	return ans;
}


}
#endif //DEBUG_HH
