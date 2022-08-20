// Homework4a.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//
#include<iostream>
#include<cmath>
#include<ctime>
#include<fstream>
#define PI 3.14159265358979323846
using namespace std;

double intrucation(int N)
{
	double axis = 0, plane = 0, volume = 0;
	double q2 = 0.5, q = pow(q2, 0.5);
	double x2 = 0, s2 = 0, r2 = 0;
	double fc = 4 * PI * N - 2 * PI * q * log((N + q) / (N - q)) + 1 / q2;//积分项以及求和n=0的项
	for (int i = 1; i <= N; i++)
	{
		x2 = i * i;
		axis += 1 / (x2 - q2);//计算轴上的值
		for (int j = 1; j <= pow(N * N - x2, 0.5); j++)
		{
			s2 = x2 + j * j;
			plane+= 1 / (s2 - q2);//计算面上的值
			for (int k = 1; k <= pow(N * N - s2, 0.5); k++)
			{
				r2 = s2 + k * k;
				volume+= 1 / (r2 - q2);//计算体上的值
			}
		}
	}
	return 6 * axis + 12 * plane + 8 * volume - fc;
}
/*
int main()
{
	cout << intrucation(100);
}
*/

int main()
{
	double r=0;
	ofstream fout("Homework4aResult.txt", ios::app);
	for(int i=500;i<=6000;i=i+500)
	{
		time_t start=time(0);
		r=intrucation(i);
		cout<<"When Lambda="<<i<<" , the result is "<<r<<" , time used: "<<time(0)-start<<endl;
		fout<<"When Lambda="<<i<<" , the result is "<<r<<endl;
	}
	return 0;
}

/*
int main()
{
	cout << intrucation(3000);
	return 0;
}
*/
// 运行程序: Ctrl + F5 或调试 >“开始执行(不调试)”菜单
// 调试程序: F5 或调试 >“开始调试”菜单

// 入门使用技巧: 
//   1. 使用解决方案资源管理器窗口添加/管理文件
//   2. 使用团队资源管理器窗口连接到源代码管理
//   3. 使用输出窗口查看生成输出和其他消息
//   4. 使用错误列表窗口查看错误
//   5. 转到“项目”>“添加新项”以创建新的代码文件，或转到“项目”>“添加现有项”以将现有代码文件添加到项目
//   6. 将来，若要再次打开此项目，请转到“文件”>“打开”>“项目”并选择 .sln 文件
