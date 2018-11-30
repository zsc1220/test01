#pragma once

#define ERRCODE(x) (errcode = ((errcode>100) ? (errcode) : (x))) 

typedef struct CheckingData
{
	char* ID;
	float value;
	std::string type;
};

//float calmatrix(CheckingData *, int , int , int , float *);
void calmatrix(CheckingData *, int, int, float *,float *);
