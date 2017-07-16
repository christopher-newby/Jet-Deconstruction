#include <iostream>
#include <iomanip>
#include <math.h>
#include <fstream>
#include "JD3header.h"

using namespace std;

//Code for running through a file of events
void FileRunS4(); //Looks at a file with 4 jets signal
void FileRunB4(); //Looks at a file with 4 jets background
void FileRunS5(); //Looks at a file with 5 jets signal
void FileRunB5(); //Looks at a file with 5 jets background
void FileRunB6(); //Looks at a file with 6 jets background

int main()
{
    ofstream fileout("C:\\Users\\use\\Desktop\\Output2.dat");
    fileout<<"File of Chi outputs:\n--------------------------------\n\n";
    fileout.close();

    FileRunS4();
    FileRunB4();
    FileRunS5();
    FileRunB5();
    FileRunB6();

    return 0;
}

void FileRunS4() //Looks at a file with formatting t,x,y,z,btag
{
    ifstream filein("C:\\Users\\use\\Desktop\\S4j.dat", ios::in);

    ofstream fileout("C:\\Users\\use\\Desktop\\Output2.dat", ios::app);

    int n=4;//number of jets in event
    double t[n],x[n],y[n],z[n];
    int btag[n],i=1;

    fileout<<"Signal 4 jets:\n---------------------\n";

    if(filein.is_open())
    {
        while(!filein.eof())
        {
            filein>>t[i%n]>>x[i%n]>>y[i%n]>>z[i%n]>>btag[i%n];

            if(i%n==0)
            {
                fourvector v1(t[0],x[0],y[0],z[0]);
                fourvector v2(t[1],x[1],y[1],z[1]);
                fourvector v3(t[2],x[2],y[2],z[2]);
                fourvector v4(t[3],x[3],y[3],z[3]);

                microjet mJ[maxsize];

                microjet m1(v1,btag[0]);
                microjet m2(v2,btag[1]);
                microjet m3(v3,btag[2]);
                microjet m4(v4,btag[3]);

                mJ[0]=m1;
                mJ[1]=m2;
                mJ[2]=m3;
                mJ[3]=m4;

                fileout<<"Chi = "<<Chi(mJ,n)<<endl;
            }

            i++;
        }
    }

    fileout<<endl;

    filein.close();
    fileout.close();
}

void FileRunB4() //Looks at a file with formatting t,x,y,z,btag
{
    ifstream filein("C:\\Users\\use\\Desktop\\B4j.dat", ios::in);

    ofstream fileout("C:\\Users\\use\\Desktop\\Output2.dat", ios::app);

    int n=4;//number of jets in event
    double t[n],x[n],y[n],z[n];
    int btag[n],i=1;

    fileout<<"Background 4 jets:\n---------------------\n";

    if(filein.is_open())
    {
        while(!filein.eof())
        {
            filein>>t[i%n]>>x[i%n]>>y[i%n]>>z[i%n]>>btag[i%n];

            if(i%n==0)
            {
                fourvector v1(t[0],x[0],y[0],z[0]);
                fourvector v2(t[1],x[1],y[1],z[1]);
                fourvector v3(t[2],x[2],y[2],z[2]);
                fourvector v4(t[3],x[3],y[3],z[3]);

                microjet mJ[maxsize];

                microjet m1(v1,btag[0]);
                microjet m2(v2,btag[1]);
                microjet m3(v3,btag[2]);
                microjet m4(v4,btag[3]);

                mJ[0]=m1;
                mJ[1]=m2;
                mJ[2]=m3;
                mJ[3]=m4;

                fileout<<"Chi = "<<Chi(mJ,n)<<endl;
            }

            i++;
        }
    }

    fileout<<endl;

    filein.close();
    fileout.close();
}

void FileRunS5() //Looks at a file with formatting t,x,y,z,btag
{
    ifstream filein("C:\\Users\\use\\Desktop\\S5j.dat", ios::in);

    ofstream fileout("C:\\Users\\use\\Desktop\\Output2.dat", ios::app);

    int n=5;//number of jets in event
    double t[n],x[n],y[n],z[n];
    int btag[n],i=1;

    fileout<<"Signal 5 jets:\n---------------------\n";

    if(filein.is_open())
    {
        while(!filein.eof())
        {
            filein>>t[i%n]>>x[i%n]>>y[i%n]>>z[i%n]>>btag[i%n];

            if(i%n==0)
            {
                fourvector v1(t[0],x[0],y[0],z[0]);
                fourvector v2(t[1],x[1],y[1],z[1]);
                fourvector v3(t[2],x[2],y[2],z[2]);
                fourvector v4(t[3],x[3],y[3],z[3]);
                fourvector v5(t[4],x[4],y[4],z[4]);

                microjet mJ[maxsize];

                microjet m1(v1,btag[0]);
                microjet m2(v2,btag[1]);
                microjet m3(v3,btag[2]);
                microjet m4(v4,btag[3]);
                microjet m5(v5,btag[4]);

                mJ[0]=m1;
                mJ[1]=m2;
                mJ[2]=m3;
                mJ[3]=m4;
                mJ[4]=m5;

                fileout<<"Chi = "<<Chi(mJ,n)<<endl;
            }

            i++;
        }
    }

    fileout<<endl;

    filein.close();
    fileout.close();
}

void FileRunB5() //Looks at a file with formatting t,x,y,z,btag
{
    ifstream filein("C:\\Users\\use\\Desktop\\B5j.dat", ios::in);

    ofstream fileout("C:\\Users\\use\\Desktop\\Output2.dat", ios::app);

    int n=5;//number of jets in event
    double t[n],x[n],y[n],z[n];
    int btag[n],i=1;

    fileout<<"Background 5 jets:\n---------------------\n";

    if(filein.is_open())
    {
        while(!filein.eof())
        {
            filein>>t[i%n]>>x[i%n]>>y[i%n]>>z[i%n]>>btag[i%n];

            if(i%n==0)
            {
                fourvector v1(t[0],x[0],y[0],z[0]);
                fourvector v2(t[1],x[1],y[1],z[1]);
                fourvector v3(t[2],x[2],y[2],z[2]);
                fourvector v4(t[3],x[3],y[3],z[3]);
                fourvector v5(t[4],x[4],y[4],z[4]);

                microjet mJ[maxsize];

                microjet m1(v1,btag[0]);
                microjet m2(v2,btag[1]);
                microjet m3(v3,btag[2]);
                microjet m4(v4,btag[3]);
                microjet m5(v5,btag[4]);

                mJ[0]=m1;
                mJ[1]=m2;
                mJ[2]=m3;
                mJ[3]=m4;
                mJ[4]=m5;

                fileout<<"Chi = "<<Chi(mJ,n)<<endl;
            }

            i++;
        }
    }

    fileout<<endl;

    filein.close();
    fileout.close();
}

void FileRunB6() //Looks at a file with formatting t,x,y,z,btag
{
    ifstream filein("C:\\Users\\use\\Desktop\\B6j.dat", ios::in);

    ofstream fileout("C:\\Users\\use\\Desktop\\Output2.dat", ios::app);

    int n=6;//number of jets in event
    double t[n],x[n],y[n],z[n];
    int btag[n],i=1;

    fileout<<"Background 6 jets:\n---------------------\n";

    if(filein.is_open())
    {
        while(!filein.eof())
        {
            filein>>t[i%n]>>x[i%n]>>y[i%n]>>z[i%n]>>btag[i%n];

            if(i%n==0)
            {
                fourvector v1(t[0],x[0],y[0],z[0]);
                fourvector v2(t[1],x[1],y[1],z[1]);
                fourvector v3(t[2],x[2],y[2],z[2]);
                fourvector v4(t[3],x[3],y[3],z[3]);
                fourvector v5(t[4],x[4],y[4],z[4]);
                fourvector v6(t[5],x[5],y[5],z[5]);

                microjet mJ[maxsize];

                microjet m1(v1,btag[0]);
                microjet m2(v2,btag[1]);
                microjet m3(v3,btag[2]);
                microjet m4(v4,btag[3]);
                microjet m5(v5,btag[4]);
                microjet m6(v6,btag[5]);

                mJ[0]=m1;
                mJ[1]=m2;
                mJ[2]=m3;
                mJ[3]=m4;
                mJ[4]=m5;
                mJ[5]=m6;

                fileout<<"Chi = "<<Chi(mJ,n)<<endl;
            }

            i++;
        }
    }

    fileout<<endl;

    filein.close();
    fileout.close();
}
