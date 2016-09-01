#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <palabos3D.h>
#include "myheaders3D.h"
#include <memory>
#include <string>

namespace plb{

template<typename T>
class Constants{
private:
public:
	static int objCount;

	Constants();

	~Constants();

	void initialize(const std::string& fileName);

// Properties
	static std::string parameterXmlFileName, obstacle_file, obstacle_mat, wall_file, wall_mat;
	static plint testIter, testRe, testTime, maxRe, minRe, maxGridLevel, referenceResolution, margin,
		borderWidth, extraLayer, blockSize, envelopeWidth;
	static T initialTemperature, gravitationalAcceleration, u0lb, epsilon, maxT, imageSave;
	static bool dynamicObstacle, dynamicWall, test;
	static Array<std::string,2> wall, obstacle;
	static Array<T,2> wall_data;
	static Array<T,4> obstacle_data;
	static Precision precision;
	static std::unique_ptr<Constants<T> > c;
private:
	static bool master;
};

template<typename T>
int Constants<T>::objCount= 0;

template<typename T>
std::string Constants<T>::parameterXmlFileName= "";

template<typename T>
std::string Constants<T>::obstacle_file= "";

template<typename T>
std::string Constants<T>::obstacle_mat= "";

template<typename T>
std::string Constants<T>::wall_file= "";

template<typename T>
std::string Constants<T>::wall_mat= "";

template<typename T>
plint Constants<T>::extraLayer= 0;

template<typename T>
plint Constants<T>::borderWidth= 0;

template<typename T>
plint Constants<T>::maxGridLevel= 0;

template<typename T>
plint Constants<T>::envelopeWidth= 0;

template<typename T>
plint Constants<T>::referenceResolution= 0;

template<typename T>
plint Constants<T>::maxRe= 0;

template<typename T>
plint Constants<T>::minRe= 0;

template<typename T>
plint Constants<T>::testRe= 0;

template<typename T>
plint Constants<T>::margin= 0;

template<typename T>
plint Constants<T>::testIter= 0;

template<typename T>
plint Constants<T>::testTime= 0;

template<typename T>
plint Constants<T>::blockSize= 0;

template<typename T>
T Constants<T>::maxT= 0;

template<typename T>
T Constants<T>::initialTemperature= 0;

template<typename T>
T Constants<T>::gravitationalAcceleration= 0;

template<typename T>
T Constants<T>::u0lb= 0;

template<typename T>
T Constants<T>::epsilon= 0;

template<typename T>
T Constants<T>::imageSave= 0;

template<typename T>
bool Constants<T>::dynamicObstacle= false;

template<typename T>
bool Constants<T>::dynamicWall= false;

template<typename T>
bool Constants<T>::test= false;

template<typename T>
bool Constants<T>::master= false;

template<typename T>
Array<T,2> Constants<T>::wall_data = Array<T,2>();

template<typename T>
Array<T,4> Constants<T>::obstacle_data = Array<T,4>();

template<typename T>
Precision Constants<T>::precision = FLT;

template<typename T>
std::unique_ptr<Constants<T> > Constants<T>::c(nullptr);

}//namespace plb

#endif //CONSTANTS_H
