#include <iostream>
#include <string>
#include <math.h>
#include <Eigen/Dense>  

#include "calibration.h"

using namespace Eigen;

extern "C" {
#include "epanet2.H"
}
/*
typedef struct CheckingData
{
	char* ID;
	float value;
	std::string type;
};*/

void calmatrix(CheckingData *CD, int length, int count , float *checkc,float *checkq)
/*
length:监测值个数
count:压力监测点数量

*checkc:摩阻系数校核精度
*checkq:用水量校核精度
*/
{
	int i = 0, m = 0, n = 0, a[6] = { 0 };
	
	a[0] = ENnetSize(0);//节点个数
	a[1] = ENnetSize(1);//水池个数
	a[2] = ENnetSize(2);//连接个数
	a[3] = ENnetSize(3);//管段个数
	a[4] = ENnetSize(4);//水泵个数
	a[5] = ENnetSize(5);//连接点个数
	MatrixXd
		A(a[0], a[3] + a[4]),
		B(a[3] + a[4], a[3] + a[4]),
		S(a[3] + a[4], a[3] + a[4]),
		D(a[3] + a[4], a[3] + a[4]),
		JhC(a[0], a[3]),
		JqC(a[3] + a[4], a[3]),
		JhQ(a[0], a[0]),
		JqQ(a[3] + a[4], a[0]);
	A.fill(0); B.fill(0); S.fill(0); D.fill(0); JhC.fill(0); JqC.fill(0);
	n = a[3];

	float _d = 0, _q = 0, _h = 0, _C = 0,
		*d = &_d, *q = &_q, *h = &_h, *C = &_C;
	int _N1 = 0, _N2 = 0,
		*N1 = &_N1, *N2 = &_N2,
		_code = 0, *code = &_code;
	
	for (i = 1; i <= a[2]; i++)
	{
		ENgetlinktype(i, code);
		if (*code <= 2 && *code != -1)
		{
			ENgetlinknodes(i, N1, N2);
			ENgetnodetype(*N1, code);
			if (*code == 0) A(*N1 - 1, i - 1) = -1;
			ENgetnodetype(*N2, code);
			if (*code == 0) A(*N2 - 1, i - 1) = 1;
			ENgetlinkvalue(i, 8, q);
			ENgetlinkvalue(i, 0, d);
			if (*d != 0)
			{
				ENgetlinkvalue(i, 10, h);
				ENgetlinkvalue(i, 2, C);
				if (*h != 0) {B(m, m) = abs(*q / (1.852 * (*h))); }
				S(m, m) = *q / (*C);
				D(m, m) = 4.871*(*q) / (1.852*(*d));
				m++;
			}
			else
			{
				ENgetpumpvalue(i, h);
				ENgetlinkvalue(i, 12, C);
				if (*h != 0 && *q != 0 ) { B(n, n) = 1 / (2 * (*h) * abs(*q)); }
				S(n, n) = 0;
				D(n, n) = 0;
				n++;
			}
		}
	}

	/*
	printf("矩阵A：\n");
	std::cout << A << "\n" <<std::endl;
	printf("矩阵B：\n");
	std::cout << B << "\n" << std::endl;
	printf("矩阵S：\n");
	std::cout << S << "\n" << std::endl;
	printf("矩阵D：\n");
	std::cout << D << "\n" << std::endl;
	*/

	//节点水压对管道阻力系数雅克比矩阵
	JhC = (((A * B*(A.transpose())).inverse()) * (A*S)).block(0, 0, a[0], a[3]);
	//管道流量对管道阻力系数的雅克比矩阵
	JqC = (S - ((B*(A.transpose()))*((A*B*(A.transpose()))).inverse()*(A*S))).block(0, 0, a[3] + a[4], a[3]);

	//节点水压对节点流量的雅克比矩阵
	JhQ = -(A*B*(A.transpose())).inverse();
	//管道流量对节点流量的雅克比矩阵
	JqQ = B * (A.transpose())*((A*B*(A.transpose())).inverse());
	
	float _pattern = 0, *pattern = &_pattern;
	float *patterns = new float[a[0]];
	memset(patterns,0,a[0]*sizeof(patterns));
	int jgroup = 0;
	for (i = 1; i <= a[5]; i++)
	{
		ENgetnodetype(i, code);
		if (*code == 0) 
		{
			*pattern = 0;
			ENgetnodevalue(i, 2, pattern);
			int c = *pattern;
			if (jgroup <= c) jgroup = c;
			ENgetnodevalue(i, 9, q);
			patterns[c] += *q;
		}
	}

	/*
	printf("矩阵JhC：\n");
	std::cout << JhC << "\n" << std::endl;
	printf("矩阵JqC：\n");
	std::cout << JqC << "\n" << std::endl;
	
	printf("矩阵JhQ：\n");
	std::cout << JhQ << "\n" << std::endl;
	printf("矩阵JqQ：\n");
	std::cout << JqQ << "\n" << std::endl;
	*/

	//摩阻系数校核矩阵dertC()//用户用水量校核矩阵dertH(), 
	MatrixXd	dertC(length, 1), dertQ(jgroup, 1),
		Jhqi(length + a[3], a[3]), dertHQ(length, jgroup),
		dertHQC(length + a[3], 1), w(length + a[3], length + a[3]);
	dertHQC.fill(0); dertHQ.fill(0);
	w.setIdentity(length + a[3], length + a[3]);
	Jhqi.fill(0),w = w * 0.01;
	m = 0, n = count;
	int _nodeIndex = 0, *nodeIndex = &_nodeIndex,
		_cIndex = 0, *cIndex = &_cIndex, linkindex = 0;
	float _presure = -1, *presure = &_presure;

	for (int j = 0; j < length; j++)
	{
		if (CD[j].type == "presure") //Jh（C1）
		{
			ENgetnodeindex(CD[j].ID, cIndex);
			if (*cIndex != 0)
			{
				linkindex = 1;
				for (i = 1; i < *cIndex; i++)
				{
					ENgetnodetype(i, code);
					if (*code != 0) linkindex++;
				}
				ENgetnodevalue(*cIndex, 11, presure);
				w(m, m) = pow(0.3, -2);
				Jhqi.row(m) = JhC.row(*cIndex - linkindex);
				dertHQC(m, 0) = CD[j].value - *presure;

				int c1 = 0, c2 = count;
				for (i = 1; i <= a[5]; i++)
				{
					ENgetnodetype(i, code);
					if (*code == 0)
					{
						*pattern = 0, *q = 0;
						ENgetnodevalue(i, 2, pattern);
						int c = *pattern;
						ENgetnodevalue(i, 9, q);
						dertHQ(m, c - 1) += (*q / patterns[c]) * JhQ(*cIndex - linkindex, c1);
						c1++;
					}
				}
				m++;
			}
			else
			{
				//压力监测点错误
			}
		}
		else//Jq(C1)
		{
			ENgetlinkindex(CD[j].ID, cIndex);
			if (*cIndex != 0)
			{
				linkindex = 1;
				for (i = 1; i < *cIndex; i++)
				{
					ENgetlinktype(i, code);
					if (*code > 1) linkindex++;
				}
				float _flow = 0, *flow = &_flow;
				ENgetlinkvalue(*cIndex, 8, flow);
				w(n, n) = pow(2, -2);
				Jhqi.row(n) = JqC.row(*cIndex - linkindex);
				dertHQC(n, 0) = CD[j].value - *flow;

				int c1 = 0, c2 = count;
				for (i = 1; i <= a[5]; i++)
				{
					ENgetnodetype(i, code);
					if (*code == 0)
					{
						*pattern = 0, *q = 0;
						ENgetnodevalue(i, 2, pattern);
						int c = *pattern;
						ENgetnodevalue(i, 9, q);
						dertHQ(n, c - 1) += (*q / patterns[c]) * JqQ(*cIndex - linkindex, c1);
						c1++;
					}
				}
				n++;
			}
			else
			{
				//流量监测点错误
			}
		}
	}

	
	for (i = 0; i < a[3]; i++)
	{
		Jhqi(length + i, i) = 1;
	}
	/*
	printf("矩阵Jhqi：\n");
	std::cout << Jhqi << "\n" << std::endl;
	
	printf("矩阵w：\n");
	std::cout << w << "\n" << std::endl;
	
	printf("矩阵dertHQC：\n");
	std::cout << dertHQC << "\n" << std::endl;
	printf("矩阵dertHQ：\n");
	std::cout << dertHQ << "\n" << std::endl;
	*/

	//管道摩阻系数改正值矩阵
	dertC = (Jhqi.transpose()*w*Jhqi).inverse()*Jhqi.transpose()*w*dertHQC;
	
	//节点需水量改正值矩阵
	dertQ = (dertHQ.transpose()*(w.block(0, 0, length, length))*dertHQ).inverse()
		*dertHQ.transpose()*(w.block(0, 0, length, length))*(dertHQC.block(0, 0, length, 1));
	
	/*
	printf("矩阵dertC：\n");
	std::cout << dertC << "\n" << std::endl;
	printf("矩阵dertQ：\n");
	std::cout << dertQ << "\n" << std::endl;
	*/
	float dertc2 = 0, dertq2 = 0;
	for (i = 0; i < a[3]; i++) 
	{
		dertc2 += (dertC(i, 0)* dertC(i, 0));
	}
	for (i = 0; i < jgroup; i++)
	{
		dertq2 += (dertQ(i, 0)*dertQ(i, 0));
	}

	dertc2 = sqrt(dertc2);
	dertq2 = sqrt(dertq2);

	*checkc = (float)dertc2;
	*checkq = (float)dertq2;
	m = 0;
	for (i = 1; i <= a[2]; i++) 
	{
		ENgetlinkvalue(i, 0, d);
		if (*d != 0)
		{
			ENgetlinkvalue(i, 2, C);
			ENsetlinkvalue(i, 2, *C + dertC(m, 0));
			m++;
		}
	}
	for (i = 1; i <= a[5]; i++)
	{
		ENgetnodetype(i, code);
		if (*code == 0)
		{
			*pattern = 0;
			ENgetnodevalue(i, 2, pattern);
			int c = *pattern;
			ENgetnodevalue(i, 9, q);
			if(dertQ(c - 1, 0)!=NULL)ENsetnodevalue(i, 1, *q + (*q / patterns[c])*dertQ(c - 1, 0));
		}
	}

	delete[] patterns;
	patterns = NULL;

	int errcode = 0;
	if (dertc2 > 0.01 || dertq2 > 0.01) { ERRCODE(ENsolveH());}
	else 
	{
		for (i = 1; i <= a[2]; i++)
		{
			ENgetlinkvalue(i, 0, d);
			if (*d != 0)
			{
				ENgetlinkvalue(i, 2, C);
				std::cout << *C << "\n" << std::endl;
			}
		}
	}
}