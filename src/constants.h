#ifndef CONSTANTS_H
#define CONSTANTS_H

#include <palabos3D.h>
#include "myheaders3D.h"
#include <memory>
#include <string>

namespace plb{

template<typename T>
struct Object{
	std::string fileName;
	Array<T,3> start;
	Array<T,3> dim;
	int referenceDirection;
	double density;
	std::string material;
	double initialTemperature;
	bool dynamicMesh;
};

template<typename T>
struct Param{
	T u;
	T length;
	T resolution;
	T lx;
	T ly;
	T lz;
	T dx;
};

template<typename T>
class Constants{
private:
public:
	static int objCount;

	Constants();

	~Constants();

	void initialize(const std::string& fileName);

	void createMinimalSettings();

// Properties
	static Object<T> obstacle, wall;
	static Param<T> physical, lb;
	static std::string parameterXmlFileName;
	static plint testIter, ibIter, testRe, testTime, maxRe, minRe, maxGridLevel, margin,
		borderWidth, extraLayer, blockSize, envelopeWidth;
	static T initialTemperature, gravitationalAcceleration, epsilon, maxT, imageSave;
	static bool test;
	static Precision precision;
	static std::unique_ptr<Constants<T> > c;
private:
	static bool master;
};

template<typename T>
int Constants<T>::objCount= 0;

template<typename T>
Object<T> Constants<T>::obstacle;

template<typename T>
Object<T> Constants<T>::wall;

template<typename T>
Param<T> Constants<T>::physical;

template<typename T>
Param<T> Constants<T>::lb;

template<typename T>
std::string Constants<T>::parameterXmlFileName= "";

template<typename T>
plint Constants<T>::extraLayer= 0;

template<typename T>
plint Constants<T>::borderWidth= 0;

template<typename T>
plint Constants<T>::maxGridLevel= 0;

template<typename T>
plint Constants<T>::envelopeWidth= 0;

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
plint Constants<T>::ibIter= 0;

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
T Constants<T>::epsilon= 0;

template<typename T>
T Constants<T>::imageSave= 0;

template<typename T>
bool Constants<T>::test= false;

template<typename T>
bool Constants<T>::master= false;

template<typename T>
Precision Constants<T>::precision = FLT;

template<typename T>
std::unique_ptr<Constants<T> > Constants<T>::c(nullptr);

}//namespace plb

#endif //CONSTANTS_H
