#ifndef FORALLOPERATIONS_H
#define FORALLOPERATIONS_H
#include <iostream>
#include <vector>

using namespace std;

template<class T>
class forAllOperations{

	public:
		forAllOperations<T>();
		virtual ~forAllOperations<T>();

		vector<T> temp1dvector;
		vector<vector<T>> tempvector;

};

#endif


#define forAll(tempvector) \
for(unsigned int i=0;i<tempvector.size();i++)  \
for(unsigned int j=0;j<tempvector[i].size();j++)

#define forAllInternal(tempvector) \
for(unsigned int i=1;i<tempvector.size()-1;i++)  \
for(unsigned int j=1;j<tempvector[i].size()-1;j++)

#define forAllInternalUCVs(tempvector) \
for(unsigned int i=1;i<tempvector.size()-2;i++)  \
for(unsigned int j=1;j<tempvector[i].size()-1;j++)


#define forAllInternalVCVs(tempvector) \
for(unsigned int i=1;i<tempvector.size()-1;i++)  \
for(unsigned int j=1;j<tempvector[i].size()-2;j++)
