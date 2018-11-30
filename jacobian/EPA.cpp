#include <iostream>
#include <string>
#include <math.h>
#include <Eigen/Dense>  

#include "calibration.h"

using namespace Eigen;

extern "C" {
	#include "epanet2.H"
}

int main()
{
	int a;
	int  errcode = 0;
	char *f1 = new char, *f2 = new char, *f3 = new char;
	printf("输入文件:\n");//D:\work\操作\供水\epaneth\1abc.inp
	printf("D:\\work\\操作\\供水\\epaneth\\1abc.inp\n");
	f1 = const_cast<char*>("D:\\work\\操作\\供水\\epaneth\\NET.inp");
	//scanf("%s", f1);
	printf("输出文件:\n");
	printf("D:/work/操作/供水/epaneth/1.txt\n");
	f2 = const_cast<char*>("D:\\work\\操作\\供水\\epaneth\\1.txt");
	//scanf("%s", f2);
	//char _f3[2] = "0";
	//f3 = &_f3[2];
	f3 = const_cast<char*>("\0");
	printf("calculating...\n");

	//errcode = ENepanet(f1, f2, f3, NULL);

	errcode = ENopen(f1, f2, f3);
	ERRCODE(ENsolveH());

	CheckingData CD[10];

	char *node = const_cast<char*>("2");
	CD[4].ID = node;
	CD[4].value = 25.96;
	CD[4].type = "presure";
	char *node2 = const_cast<char*>("4");
	CD[2].ID = node2;
	CD[2].value = 28.99;
	CD[2].type = "presure";
	char *node3 = const_cast<char*>("5");
	CD[5].ID = node3;
	CD[5].value = 33.79;
	CD[5].type = "presure";
	char *pipe = const_cast<char*>("1");
	CD[0].ID = pipe;
	CD[0].value = -240.63;
	CD[0].type = "flow";
	char *pipe2 = const_cast<char*>("4");
	CD[1].ID = pipe2;
	CD[1].value = 140.06;
	CD[1].type = "flow";
	char *pipe3 = const_cast<char*>("5");
	CD[3].ID = pipe3;
	CD[3].value = 120.06;
	CD[3].type = "flow";

	float *checkc = new float;
	*checkc = 10;
	float *checkq = new float;
	*checkq = 10;
	int i = 0;
	while (*checkc >= 0.01 || *checkq >= 0.01)
	{
		calmatrix(CD, 3, 1, checkc, checkq);
		i++;
	}
	printf("%d\n",i);

	delete checkc, checkq;
	ERRCODE(ENsolveQ());
	//calmatrix(CD, 2, 1, 1);
	ERRCODE(ENreport());
	ENclose();
	//return(errcode);

	if (errcode == 0) { printf("校核！\n"); }
	else { printf("出错！\n"); }

	system("pause");
}