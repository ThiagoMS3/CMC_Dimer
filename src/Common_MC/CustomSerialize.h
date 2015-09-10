#ifndef CUSTOMSERIALIZATION_H_
#define CUSTOMSERIALIZATION_H_

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>

#include <iostream>
#include <fstream>

// Definitions
using namespace std;
using namespace boost;

template<typename T>
inline void ExportObject(T& object, const char filename[])
{
	ofstream out;

	out.open(filename);
	boost::archive::text_oarchive out_obj(out);
	out_obj << object;
	out.close();
}

template<typename T>
inline void ImportObject(T& object, const char filename[])
{
	ifstream in;

	in.open(filename);
	boost::archive::text_iarchive in_obj(in);
	in_obj >> object;
	in.close();
}

template<typename T>
inline void ExportObjectText(T& object, const char filename[])
{
	ofstream out;

	out.open(filename);
	boost::archive::text_oarchive out_obj(out);
	out_obj << object;
	out.close();
}

template<typename T>
inline void ImportObjectText(T& object, const char filename[])
{
	ifstream in;

	in.open(filename);
	boost::archive::text_iarchive in_obj(in);
	in_obj >> object;
	in.close();
}

#endif
