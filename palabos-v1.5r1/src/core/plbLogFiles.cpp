/* This file is part of the Palabos library.
 *
 * Copyright (C) 2011-2015 FlowKit Sarl
 * Route d'Oron 2
 * 1010 Lausanne, Switzerland
 * E-mail contact: contact@flowkit.com
 *
 * The most recent release of Palabos can be downloaded at 
 * <http://www.palabos.org/>
 *
 * The library Palabos is free software: you can redistribute it and/or
 * modify it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * The library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "parallelism/mpiManager.h"
#include "core/plbLogFiles.h"
#include "core/util.h"
#include <sys/stat.h>
#include <stdio.h>
#include <utility>
#include <fstream>
#include <exception>
#include <string>
#include <ctime>

namespace plb {

namespace global {

PlbLogFile::PlbLogFile(PlbLogFile const& rhs)
{ }

PlbLogFile& PlbLogFile::operator=(PlbLogFile const& rhs)
{
    return *this;
}

inline bool exists(const std::string& name) {
	struct stat buffer;
	return (stat (name.c_str(), &buffer) == 0);
}

inline std::string time_string(){
	time_t rawtime = time(0);
	struct tm* now  = localtime(&rawtime);
	const std::string day = std::to_string(now->tm_mday);
	const std::string month = std::to_string(now->tm_mon+1);
	const std::string hour = std::to_string(now->tm_hour);
	const std::string minute = std::to_string(now->tm_min);
	const std::string second = std::to_string(now->tm_sec);
	const std::string time = "["+day+"/"+month+" "+hour+":"+minute+":"+second+"]";
	return time;
}

void PlbLogFile::init(bool parallel_){
	this->parallel = parallel_;
	this->indentation=0;
	this->indentSpaces="";
	if (global::mpi().isMainProcessor()){this->fName=global::directories().getLogOutDir()+"Main.Log";}
	else{this->fName = global::directories().getLogOutDir()+"Rank"+util::val2str(global::mpi().getRank())+".log";}
	if(exists(this->fName)){ remove(fName.c_str());}
}

void PlbLogFile::push(std::string sectionName)
{
    entry("\nSECTION " + sectionName);
    indentation += 4;
    indentSpaces = std::string(indentation, ' ');
}

void PlbLogFile::pop()
{
    indentation -= 4;
    indentSpaces = std::string(indentation, ' ');
}

void PlbLogFile::entry(std::string entryText) {
	std::ofstream ofile;
	if(!ofile.is_open()){
		ofile.open(this->fName.c_str(),std::ofstream::out | std::ofstream::app);
	}
	ofile << indentSpaces;
	ofile << time_string();
	ofile << entryText;
	ofile << "\n";
}

void PlbLogFile::flushEntry(std::string entryText) {
	std::ofstream ofile;
	if(!ofile.is_open()){
		ofile.open(this->fName.c_str(),std::ofstream::out | std::ofstream::app);
	}
	ofile << indentSpaces;
	ofile << time_string();
	ofile << entryText;
	ofile << std::endl;
	ofile.flush();
	ofile.close();
}

LogFileCollection::LogFileCollection(bool parallel_)
    : parallel(parallel_)
{ }

LogFileCollection::~LogFileCollection() {
    std::map<std::string, PlbLogFile*>::iterator it=collection.begin();
    for (; it != collection.end(); ++it) {
        delete it->second;
    }
}

LogFileCollection::LogFileCollection(LogFileCollection const& rhs)
{ }

LogFileCollection& LogFileCollection::operator=(LogFileCollection const& rhs) {
    return *this;
}

PlbLogFile& LogFileCollection::get(std::string nameOfLogFile) {
    std::map<std::string, PlbLogFile*>::iterator it=collection.find(nameOfLogFile);
    if (it == collection.end()) {
        PlbLogFile* logfile = &global::log();
		logfile->init(parallel);
        collection.insert(make_pair(nameOfLogFile, logfile));
        return *logfile;
    }
    else {
        return *it->second;
    }
}

PlbLogFile& logfile(std::string nameOfLogFile) {
    static LogFileCollection collection(true);
    return collection.get(nameOfLogFile);
}

PlbLogFile& logfile_nonparallel(std::string nameOfLogFile) {
    static LogFileCollection collection(false);
    return collection.get(nameOfLogFile);
}

}  // namespace global

}  // namespace plb
