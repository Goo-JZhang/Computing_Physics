#include <iostream>
#include<fstream>
#include <stdlib.h>
using namespace std;

int main()
{
    int i,j,k;
    double A[20][6][6];
    double a, b;
    a = -1;
    b = 1;
    srand(0);
    for(i = 0; i < 20; i++)
        for(j = 0; j < 6; j++)
            for(k = 0; k < 6; k++)
            {
                A[i][j][k] = (double) rand()/RAND_MAX*(b-a)+a;
            }
    ofstream outfile;
    outfile.open("random_matrix_C_Windows.txt");
    //outfile.open("random_matrix_C_Mac.txt");
    for(j = 0; j < 6; j++)
    {
        for(k = 0; k < 6; k++)
            {
                outfile<<A[0][j][k]<<" ";
            }
        outfile<<endl;
    }
    outfile.close();
    return(0);
}
