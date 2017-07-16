#include "JD3header.h"
#include <bitset>
#include <cmath>
#include <iostream>
#include <fstream>


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////LRsets///////////////////////////////////////////////////////////////////////////////////////////////////
LRsets::LRsets()        // Default constructor with set {1,2,3}
{
    intset[0]=3;
    len=3;
    count=0;
    countend=false;

    for(int i=1;i<=maxsize;i++)
    {
        if(i<=len)
            intset[i]=i;
        else
            intset[i]=0;
    }
    //ctor
}

LRsets::LRsets(int (&set)[maxsize+1])       // Constructor.  Use set passed to the constructor.
{
    len = set[0];
    count=0;
    countend=false;

    for(int i=1;i<=maxsize;i++)
    {
        if(i<=len)
            intset[i]=set[i];
        else
            intset[i]=0;
    }
}//end of constructor

void LRsets::reset(int (&set)[maxsize+1])       // Resets data according to set, as in constructor.
{
    len = set[0];
    count=0;
    countend=false;

    for(int i=1;i<=maxsize;i++)
    {
        if(i<=len)
            intset[i]=set[i];
        else
            intset[i]=0;
    }
}//end of reset

// Returns next division of set into setA and setB.
// First call of s.next(setA,setB) gives setA = set, setB = {}.
// Further calls give each subset of "set" for setA, the complement for setB.
// Call 2^N is the last split into and returns setA = {}, setB = set.
// Further calls are in error and produce N = -1.
bool LRsets::next(int (&setA)[maxsize+1],int (&setB)[maxsize+1])
{

    if(countend)
    {
        setA[0]=0;
        setB[0]=0;
    }
    else
    {
        int lenL,lenR;
        bitset<maxsize> b(count);       //turns the count into a binary number

        lenL=1;
        lenR=1;

        for(int i=0;i<len;i++)
        {
            if(b[i]==0)
            {
                setA[lenL]=intset[i+1];
                lenL++;
            }
            else
            {
                setB[lenR]=intset[i+1];
                lenR++;
            }
        }

        setA[0]=lenL-1;
        setB[0]=lenR-1;

    //  This is a test to make sure the count does not exceed the total number of combinations for splitting up the list
    //  If we are at the end of the splitting subset, then the code just returns a set of length 0
    //              (as long as the user realizes that the 0th element of the returned array is the length)
        unsigned long max=int(pow(2,len));
        count++;

        if(count==max)
            countend = 1;
        else
            countend = 0;

    }//end of long else

    return countend;

}//end of next
//////////////////////////////////////End of LRSets////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////Partitions///////////////////////////////////////////////////////////////////////////////////////////////
Partitions::Partitions(int (&set)[maxsize+1])   // Constructor which uses set passed to it
{
    for(int i=0;i<maxsize+1;i++)
    {
        if(i<=set[0])
            intset[i]=set[i];
        else
            intset[i]=0;
    }

    for(int i=0;i<maxsize+1;i++)
    {
        count[i]=0;
        maxcount[i]=0;
    }

    maxcount[0]=int(pow(2,intset[0]-1))-1;
    c=0;

    countend=false;

    MaxBell=BellNumber(intset[0]);

    //ctor
}

void Partitions::reset(int (&set)[maxsize+1])   // Resets the class with the input set
{
    for(int i=0;i<maxsize+1;i++)
    {
        if(i<set[0])
            intset[i]=set[i];
        else
            intset[i]=0;
    }

    for(int i=0;i<maxsize+1;i++)
    {
        count[i]=0;
        maxcount[i]=0;
    }

    maxcount[0]=int(pow(2,intset[0]-1))-1;
    c=0;

    countend=false;
}//end of reset

// cycles through the partitions of input set
// Returns a 2D array of the specific partition set you are on
// governed by the rules laid out at beginning of class.
// Returns a boolian which if false until you reach the end of the
// total number of possible partitions, then it returns true.
// If you try to call after all partitions, then it returns empty sets.
bool Partitions::next(int (&set)[maxsize+1][maxsize+1])
{
    c++;

    nestLevel=0;
    maxcount[nestLevel]=int(pow(2,intset[0]-1))-1;
    nestwhile(set,intset);
    set[0][0]=nestLevel+1;

    for(int i=nestLevel;i<maxsize;i++)
        count[i]=0;

    count[nestLevel-1]++;

    if(c>=MaxBell)
        countend=true;
    else
        countend=false;

    return countend;
}//end of next

// Function which generates the elements of the partitions
// calls itself recursively
void Partitions::nestwhile(int (&TotSet)[maxsize+1][maxsize+1],int Rset[maxsize+1])
{
        int setA[maxsize+1],setB[maxsize+1];
        int lenA=1,lenB=1;
        int countinteger=count[nestLevel]+maxcount[nestLevel]+1;
        bitset<maxsize> b(countinteger);  //turns the c1 into a binary number

        for(int i=0;i<Rset[0];i++)
        {
            if(b[i]==1)
            {
                setA[lenA]=Rset[i+1];
                lenA++;
            }
            else
            {
                setB[lenB]=Rset[i+1];
                lenB++;
            }
        }
        setA[0]=lenA-1;
        setB[0]=lenB-1;

        for(int i=0;i<=setA[0];i++)
            TotSet[nestLevel+1][i]=setA[i];

        if(setB[0]==1)
        {
            ++nestLevel;
            maxcount[nestLevel]=0;
            TotSet[nestLevel+1][0]=setB[0];
            TotSet[nestLevel+1][1]=setB[1];
        }
        else if(setB[0]!=0)
        {
            ++nestLevel;
            maxcount[nestLevel]=int(pow(2,setB[0]-1))-1;
            nestwhile(TotSet,setB);
        }
        else if(setB[0]==0)
            return;
}//end of nestwhile


// This returns the binomial coefficient of a choose b
int Partitions::BinomialCoeff(int a,int b)
{
    return factorial(a)/(factorial(b)*factorial(a-b));
}

// This returns a!
int Partitions::factorial(int a)
{
    int hold=1;

    for(int i=1;i<=a;i++)
    {
        hold=hold*i;
    }

    return hold;
}

// This subroutine calculates the total number of partitions
int Partitions::BellNumber(int a)
{
    int b[maxsize+1];
    for(int i=0;i<maxsize+1;i++)
        b[i]=0;

    b[0]=1;

    for(int i=0;i<maxsize;i++)
    {
        for(int k=0;k<=i;k++)
            b[i+1]+=BinomialCoeff(i,k)*b[k];
    }

    return b[a];
}

//////////////////////////////////////End of Partitions////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////fourvector//////////////////////////////////////////////////////////////////////////////////////////////
fourvector::fourvector()    // Default constructor, sets all values to 0
{
    for(int i=0;i<4;i++)
        vector[i]=0;
    // default constructor
}

fourvector::fourvector(double invec[4])
{
    for(int i=0;i<4;i++)
        vector[i]=invec[i];
    //Constructor which sets the internal vector v(i)=invec[i]
}

fourvector::fourvector(double t, double x, double y, double z)
{
    vector[0]=t;
    vector[1]=x;
    vector[2]=y;
    vector[3]=z;
    //Constructor which sets the internal vector v(i)=(t,x,y,z)
}

// v.in(mu) for input, v.in(mu) = r
double& fourvector::in(const int& mu)
{
    return vector[mu];
}//end of in

// v(mu) for output,  r = v(mu)
const double fourvector::operator()(const int& mu) const
{
    return vector[mu];
}// end of operaor()

// Sum of two vectors: v3 = v1 + v2
fourvector fourvector::operator+ (const fourvector& v2) const
{
    double hold[4];

    for(int i=0;i<4;i++)
        hold[i]=vector[i]+v2(i);

    return fourvector(hold);
}//end of operator+

// Difference of two vectors: v3 = v1 - v2
fourvector fourvector::operator- (const fourvector& v2) const
{
    double hold[4];

    for(int i=0;i<4;i++)
        hold[i]=vector[i]-v2(i);

    return fourvector(hold);
}//end of operator-

// Vector times scalar: v2 = v1*c
fourvector fourvector::operator* (const double& a) const
{
    double hold[4];

    for(int i=0;i<4;i++)
        hold[i]=vector[i]*a;

    return fourvector(hold);
}//end of operator*

// -1 times vector: v2 = - v1
fourvector fourvector::operator- () const
{
    double hold[4];

    for(int i=0;i<4;i++)
        hold[i]=-1.0*vector[i];

    return fourvector(hold);
}//end of operator-

// Scalar times vector: v2 = c*v1
fourvector operator*(const double& a, const fourvector& v)
{
    double hold[4];

    for(int i=0;i<4;i++)
        hold[i]=a*v(i);

    fourvector h(hold);

    return h;
}//end of operator*

// Dot product of two 4-vectors: dot(v1,v2)
double dot(const fourvector& v1, const fourvector& v2)
{
    double hold=0;

    hold+=v1(0)*v2(0);

    for(int i=1;i<4;i++)
        hold-=v1(i)*v2(i);

    return hold;
}//end of dot

// Square of one 4-vector: square(v) = dot(v,v)
double square(const fourvector& v)
{
    double hold=0;

    hold+=v(0)*v(0);

    for(int i=1;i<4;i++)
        hold-=v(i)*v(i);

    return hold;
}//end of square
//////////////////////////////////////End of fourvector////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////microjet/////////////////////////////////////////////////////////////////////////////////////////////////
microjet::microjet()                                                   // Default constructor which sets p to (0,0,0,0) and btag to 0
{
    btag=0;
}//end of default ctor

microjet::microjet(fourvector pinit,int btinit)                        // Constructor which sets btag to btinit and p to pinit
{
    p=pinit;
    btag=btinit;
}//end of fourvector microjet ctor

microjet::microjet(double t, double x, double y, double z, int btinit) // Constructor which sets btag to btinit and p to (t,x,y,z)
{
    fourvector hold(t,x,y,z);

    p=hold;
    btag=btinit;
}//end of 4 doubles microjet ctor

microjet::microjet(double invec[4], int btinit)                        // Constructor which sets btag to btinit and p to invec (with invec[0]=t,...)
{
    fourvector hold(invec);

    p=hold;
    btag=btinit;
}//end of array microjet
//////////////////////////////////////End of microjet//////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////jetdata//////////////////////////////////////////////////////////////////////////////////////////////////
jetdata::jetdata()              // Default constructor which sets all values, except for momenta, to the very big number from parameters
{
    set[0]=0;

    parameters para;

    MsqK=para.verybig;
    kTK=para.verybig;

//    calculatingParams();        // Calculates the parameters k,y,mu,phi

    y=para.verybig;
    phi=para.verybig;
    ksq=para.verybig;
    k=para.verybig;
    musq=para.verybig;

    yccL=para.verybig;
    phiccL=para.verybig;
    yccR=para.verybig;
    phiccR=para.verybig;
    //Default ctor
}

jetdata::jetdata(microjet m1[maxsize], int length)  // Constructor for first call of jetdata
{
    for(int i=0;i<length;i++)
    {
        mList[i]=m1[i];
        p=p+mList[i].p;
    }

    set[0]=length;
    for(int i=1;i<=set[0];i++)
        set[i]=i;

    parameters para;

    MsqK=para.verybig;
    kTK=para.verybig;

    calculatingParams();        // Calculates the parameters k,y,mu,phi

    yccL=para.verybig;
    phiccL=para.verybig;
    yccR=para.verybig;
    phiccR=para.verybig;

    //end of ctor
}

jetdata::jetdata(jetdata jK, int s1[maxsize+1]) // Initializes all the momentum parts of jetdata
{
    set[0]=s1[0];

    for(int i=1;i<=s1[0];i++)
    {
        set[i]=i;
        mList[i-1]=jK.mList[s1[i]-1];
        p=p+mList[i-1].p;
    }

    MsqK=jK.musq;
    kTK=jK.k;

    calculatingParams();        // Calculates the parameters k,y,mu,phi

    parameters para;

    yccL=para.verybig;
    phiccL=para.verybig;
    yccR=para.verybig;
    phiccR=para.verybig;

    //end of ctor
}

void jetdata::calculatingParams()   //calculates k, musq, y, and phi
{
    phi=atan2(p(1),p(2));
    k=sqrt(pow(p(1),2.0)+pow(p(2),2.0));
    ksq=pow(k,2.0);
    musq=pow(p(0),2.0)-pow(p(3),2.0)-ksq;
    y=1.0/2.0*log((p(0)+p(3))/(p(0)-p(3)));
}//end of calculatingParams
//////////////////////////////////////End of jetdata///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////parameters///////////////////////////////////////////////////////////////////////////////////////////////
parameters::parameters()
{
    pi=3.1415926535898;
    MassZ=91.188;
    asMZ=0.118;
    cA=3.0;
    cF=4.0/3.0;
    TR=1.0/2.0;
    Nf=5;
    b0=(33.0 - double(2*Nf))/(12.0*pi);
    verybig=10000;
    Rfatjet=1.2;
    Rfatjetsq=pow(Rfatjet,2.0);
    pTmin=200.0;
    npdfg=2.0;
    npdfH=2.0;
    kappapsq=4.0;
    kappanpsq=4.0;
    cR=2.0;
    nR=1.0;
    cnp=1.0;
    nnp=1.5;
    massH=120.0;
    DeltaMH=10.0;
    pTtagcut=15.0;
    PTb=0.6;
    PTnotb=0.02;
    //  end of ctor
}

double alphas(double Qsq) // Alpha_s as a function of the scale.
{
    parameters para;
    double hold,MZsq;
    MZsq=pow(para.MassZ,2.0);
    hold=para.asMZ/(1.0+para.asMZ*para.b0*log(Qsq/MZsq));

    return hold;
}//end of alphas

//////////////////////////////////////End of parameters////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////Probability Functions////////////////////////////////////////////////////////////////////////////////////
double Sudakov(const jetdata& jJ)             // Sudakov exponent
{
    flavor f;
    double kTsq, massjetsq, massmothersq, kTmother, thetaLsq, thetaRsq;
    parameters para;

    f=jJ.f;
    kTsq=jJ.ksq;
    massjetsq=jJ.musq;
    massmothersq=jJ.MsqK;
    kTmother=jJ.kTK;


    if(f==gluon)  //If for setting the left and right color connected thetas for a gluon parent
    {
        if(jJ.yccL==para.verybig)
            thetaLsq=para.Rfatjetsq;
        else
            thetaLsq=pow(jJ.y-jJ.yccL,2.0)+pow(Dphi(jJ.phi,jJ.phiccL),2.0);
        if(jJ.yccR==para.verybig)
            thetaRsq=para.Rfatjetsq;
        else
            thetaRsq=pow(jJ.y-jJ.yccR,2.0)+pow(Dphi(jJ.phi,jJ.phiccR),2.0);
    }
    else if(f==quark||f==bottom)    //If for setting the left and right color connected thetas for a quark or bottom parent
    {
        if(jJ.yccR==para.verybig)
            thetaRsq=para.Rfatjetsq;
        else
            thetaRsq=pow(jJ.y-jJ.yccR,2.0)+pow(Dphi(jJ.phi,jJ.phiccR),2.0);
    }
    else if(f==qbar||f==bbar)   //If for setting the left and right color connected thetas for a anti-quark or anti-bottom parent
    {
        if(jJ.yccL==para.verybig)
            thetaLsq=para.Rfatjetsq;
        else
            thetaLsq=pow(jJ.y-jJ.yccL,2.0)+pow(Dphi(jJ.phi,jJ.phiccL),2.0);
    }

    double hold=1.0;
    double Sggg,Sgqq,Sqgq,S=0;
    double h1,h2;

    // If/else for the parent kinematic part of the sudakov factor
    if(massmothersq==para.verybig&&kTmother==para.verybig)
        h1=kTsq;
    else
        h1=sqrt(kTsq)*massmothersq/(2.0*kTmother);

    // If/else for the theta part of the sudakov factor, depending upon parent flavor
    if(f==gluon)
        h2=sqrt(thetaLsq*thetaRsq);
    else if(f==quark||f==bottom)
        h2=thetaRsq;
    else if(f==qbar||f==bbar)
        h2=thetaLsq;

    if(f==gluon)// Sudakov for gluon
    {
        Sggg=para.cA/(para.pi*pow(para.b0,2.0))*(log(alphas(massjetsq)/(alphas(h1)))*(1/alphas(h2*kTsq)-11*para.b0/12)+1/alphas(massjetsq)-1/alphas(h1));
        if(Sggg<0)
            Sggg=0;
        Sgqq=para.TR/(3.0*para.pi*para.b0)*log(alphas(massjetsq)/alphas(h1));

        S=Sggg+para.Nf*Sgqq;
    }// End of Sudakov for gluon
    else// Sudakov for quark or qbar
    {
        Sqgq=para.cF/(para.pi*pow(para.b0,2.0))*(log(alphas(massjetsq)/alphas(h1))*(1/alphas(h2*kTsq)-3*para.b0/4)+1/alphas(massjetsq)-1/alphas(h1));
        if(Sqgq<0)
            Sqgq=0;
        S=Sqgq;
    }// End of Sudakov for quark or qbar

    hold=exp(-S);

    return hold;
}//end of Sudakov

double SplittingProb(const jetdata& jJ, const jetdata& jL, const jetdata& jR)    // Probability factor at vertex
{
    double H;
    double kfac,angfac;
    parameters para;

//  This section checks to see if there is an angle factor needed
    if((jR.f==quark||jR.f==bottom)&&(jL.f==qbar||jL.f==bbar))
        angfac=1.0;
    else
        angfac=anglefactor(jL, jR, jJ);

//  Calculating the momenta dependent factor of Hxxx
    if(jJ.f==quark||jJ.f==qbar||jJ.f==bottom||jJ.f==bbar)
    {
        if(jL.f==quark||jL.f==bottom)
            kfac=8.0*para.pi*para.cF*alphas(jJ.musq)/jJ.musq*jJ.k/jR.k*(1.0+pow(jL.k/jJ.k,2.0));
        else
            kfac=8.0*para.pi*para.cF*alphas(jJ.musq)/jJ.musq*jJ.k/jL.k*(1.0+pow(jR.k/jJ.k,2.0));
    }
    else if(jJ.f==gluon&&jL.f==gluon&&jR.f==gluon)
        kfac=8.0*para.pi*para.cA*(alphas(jJ.musq)/jJ.musq)*(jJ.ksq/(jR.k*jL.k))*pow((1.0-jR.k*jL.k/jJ.ksq),2.0);
    else
        kfac=8.0*para.pi*para.TR*alphas(jJ.musq)/jJ.musq*(jL.ksq+jR.ksq)/jJ.ksq;

    H=kfac*angfac;

    return H;
}//end of SplittingProb

double anglefactor(const jetdata& jL, const jetdata& jR, const jetdata& jK)          // Angle factor g(y,phi) with jets.
{
   double ys,yh,yk,phis,phih,phik;
   double thetaSHsq,thetaSKsq,thetaHKsq;
   double hold=1.0;
   bool test=false;     // value for whether there doesn't exist a color connected partner
                        // true = no color connection
   parameters para;

   if(jK.f==gluon)// Angle factor for gluon parent
   {
       if(jL.k>jR.k)//test for which jet is harder
       {
           if(jK.yccR!=para.verybig)
           {
               ys=jR.y;
               yh=jL.y;
               yk=jK.yccR;
               phis=jR.phi;
               phih=jL.phi;
               phik=jK.phiccR;
           }
           else
               test=true;
        }
       else
       {
           if(jK.yccL!=para.verybig)
           {
               ys=jL.y;
               yh=jR.y;
               yk=jK.yccL;
               phis=jL.phi;
               phih=jR.phi;
               phik=jK.phiccL;
           }
           else
               test=true;
       }
   }//end of angle factor parameters for gluon parent
   if(jK.f==quark||jK.f==qbar||jK.f==bottom||jK.f==bbar)// Angle factor for quark parent
   {
       if(jK.f==quark||jK.f==bottom)
       {
            ys=jR.y;
            yh=jL.y;
            phis=jR.phi;
            phih=jL.phi;

            if(jK.yccR!=para.verybig)
            {
                yk=jK.yccR;
                phik=jK.phiccR;
            }
            else
                test=true;
       }
       else
       {
            ys=jL.y;
            yh=jR.y;
            phis=jL.phi;
            phih=jR.phi;

            if(jK.yccL!=para.verybig)
            {
                yk=jK.yccL;
                phik=jK.phiccL;
            }
            else
                test=true;
       }
   }//end of angle factors for quark or antiquark parent

   if(!test)
   {
       thetaSHsq=pow(ys-yh,2.0)+pow(Dphi(phis,phih),2.0);
       thetaSKsq=pow(ys-yk,2.0)+pow(Dphi(phis,phik),2.0);
       thetaHKsq=pow(yh-yk,2.0)+pow(Dphi(phih,phik),2.0);

       hold=thetaHKsq/(thetaSHsq+thetaSKsq);
   }
   else
        hold=1.0;

   return hold;
}//end of anglefactor

double Dphi(const double& phi1, const double& phi2)            //calculates the difference between two phi's and returns a value between -pi and pi
{
    double hold;
    parameters para;

    hold=phi1-phi2;

    while(hold>para.pi||hold<-para.pi)
    {
        if(hold>para.pi)
            hold-=2*para.pi;
        else
            hold+=2*para.pi;
    }

    return hold;
}//end of Dphi
//////////////////////////////////////End of Probability Functions/////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////Splitting Subroutines////////////////////////////////////////////////////////////////////////////////////
double Chi(microjet inMjet[maxsize],int length)  //Returns the observable for an initial list of microjets
{
    double chireturn,Greturn,Hreturn;
    jetdata jInit(inMjet,length);

    Greturn=Gstart(jInit);
    Hreturn=Hstart(jInit);

    chireturn=Hreturn/Greturn;

    return chireturn;

}//end of Chi

double SplitJet(const jetdata& jJ)                           //Splits the jet given to it and returns the probability
{
    double prob=0;
    double pL=0,pR=0;
    int s1[maxsize+1],sL[maxsize+1],sR[maxsize+1];
    bool EndofSets;
    parameters para;
    double hprob;

    unsigned long max=int(pow(2,jJ.set[0]));

    s1[0]=jJ.set[0];
    for(int i=1;i<=maxsize;i++) //initializes the set which will be split into L and R sets
    {
        if(i<=s1[0])
            s1[i]=i;
        else
            s1[i]=0;
    }

    if(s1[0]==1&&((2*jJ.musq/jJ.k<jJ.MsqK/jJ.kTK&&jJ.MsqK!=para.verybig)||(jJ.MsqK==para.verybig&&jJ.musq<jJ.ksq)/*Kinematic Check*/))
    {
        /*Final State Stuff (single element in jetdata jJ)*/
        if(jJ.f!=Higgs)
            prob=FinalState(jJ);
    }
    else if(((2*jJ.musq/jJ.k<jJ.MsqK/jJ.kTK&&jJ.MsqK!=para.verybig)||(jJ.MsqK==para.verybig&&jJ.musq<jJ.ksq))/*Kinematic Check*/)
    {
        LRsets LR1(s1);
        unsigned int count=0;
        EndofSets=false;

        //While for looping through all different splittings of the input microjets
        // extra check is to reduce calculation time since, if there is no parent, then the probability for
        //  a splitting with a certain set of jets right is the same as for those same set of jets left
        //  as long the other jets remain the same
        while(!EndofSets&&!(jJ.MsqK==para.verybig&&count==int(max/2.0)))
        {
            EndofSets=LR1.next(sL,sR);
            count++;

            if((sL[0]!=0&&sR[0]!=0)/*Split Check*/)
            {
                jetdata jL(jJ,sL);
                jetdata jR(jJ,sR);

                /*Checking parent flavors*/
                if(jJ.f==gluon)
                {
                    pL=0;
                    pR=0;
                    for(int i=0;i<3;i++)
                    {
                        switch(i)//switch for daughter flavors
                        {
                            case 0:     //g->g + g
                                jL.f=gluon;
                                jR.f=gluon;

                                jL.yccL=jJ.yccL;
                                jL.yccR=jR.y;
                                jL.phiccL=jJ.phiccL;
                                jL.phiccR=jR.phi;

                                jR.yccL=jL.y;
                                jR.yccR=jJ.yccR;
                                jR.phiccL=jL.phi;
                                jR.phiccR=jJ.phiccR;
                                break;//end of case g->g + g
                            case 1:     //g->qb + q
                                jL.f=qbar;
                                jR.f=quark;

                                jL.yccL=jJ.yccL;
                                jL.yccR=para.verybig;
                                jL.phiccL=jJ.phiccL;
                                jL.phiccR=para.verybig;

                                jR.yccL=para.verybig;
                                jR.yccR=jJ.yccR;
                                jR.phiccL=para.verybig;
                                jR.phiccR=jJ.phiccR;
                                break;//end of case g->qb + q
                            case 2:     //g->bb + b
                                jL.f=bbar;
                                jR.f=bottom;

                                jL.yccL=jJ.yccL;
                                jL.yccR=para.verybig;
                                jL.phiccL=jJ.phiccL;
                                jL.phiccR=para.verybig;

                                jR.yccL=para.verybig;
                                jR.yccR=jJ.yccR;
                                jR.phiccL=para.verybig;
                                jR.phiccR=jJ.phiccR;
                                break;//end of case g->qb + q
                        }

                        pL=SplitJet(jL);
                        pR=SplitJet(jR);

                        if(jL.f==qbar)//For the 4 quark flavors
                            prob=prob+(para.Nf-1)*pL*pR*SplittingProb(jJ,jL,jR)*Sudakov(jJ);
                        else
                            prob=prob+pL*pR*SplittingProb(jJ,jL,jR)*Sudakov(jJ);
                    }//end of for-loop which cycles through flavors for gluon parent
                }//end of if for gluon parent
                else if(jJ.f==quark)
                {
                    pL=0;
                    pR=0;

                    //q->q + g
                    jL.f=quark;
                    jR.f=gluon;

                    jL.yccL=para.verybig;
                    jL.yccR=jR.y;
                    jL.phiccL=para.verybig;
                    jL.phiccR=jR.phi;

                    jR.yccL=jL.y;
                    jR.yccR=jJ.yccR;
                    jR.phiccL=jL.phi;
                    jR.phiccR=jJ.phiccR;

                    //end of setting parameters for q->q + g

                    pL=SplitJet(jL);
                    pR=SplitJet(jR);

                    prob=prob+pL*pR*SplittingProb(jJ,jL,jR)*Sudakov(jJ);
                }//end of if for quark parent
                else if(jJ.f==qbar)
                {
                    pL=0;
                    pR=0;

                    //qb->g + qb
                    jL.f=gluon;
                    jR.f=qbar;

                    jL.yccL=jJ.yccL;
                    jL.yccR=jR.y;
                    jL.phiccL=jJ.phiccL;
                    jL.phiccR=jR.phi;

                    jR.yccL=jL.y;
                    jR.yccR=para.verybig;
                    jR.phiccL=jL.phi;
                    jR.phiccR=para.verybig;

                    //end of setting parameters for qb->g + qb

                    pL=SplitJet(jL);
                    pR=SplitJet(jR);

                    prob=prob+pL*pR*SplittingProb(jJ,jL,jR)*Sudakov(jJ);
                }//end of if for anti-quark parent
                else if(jJ.f==bottom)
                {
                    pL=0;
                    pR=0;

                    //q->q + g
                    jL.f=bottom;
                    jR.f=gluon;

                    jL.yccL=para.verybig;
                    jL.yccR=jR.y;
                    jL.phiccL=para.verybig;
                    jL.phiccR=jR.phi;

                    jR.yccL=jL.y;
                    jR.yccR=jJ.yccR;
                    jR.phiccL=jL.phi;
                    jR.phiccR=jJ.phiccR;

                    //end of setting parameters for q->q + g

                    pL=SplitJet(jL);
                    pR=SplitJet(jR);

                    prob=prob+pL*pR*SplittingProb(jJ,jL,jR)*Sudakov(jJ);
                }//end of if for bottom parent
                else if(jJ.f==bbar)
                {
                    pL=0;
                    pR=0;

                    //qb->g + qb
                    jL.f=gluon;
                    jR.f=bbar;

                    jL.yccL=jJ.yccL;
                    jL.yccR=jR.y;
                    jL.phiccL=jJ.phiccL;
                    jL.phiccR=jR.phi;

                    jR.yccL=jL.y;
                    jR.yccR=para.verybig;
                    jR.phiccL=jL.phi;
                    jR.phiccR=para.verybig;

                    //end of setting parameters for qb->g + qb

                    pL=SplitJet(jL);
                    pR=SplitJet(jR);

                    prob=prob+pL*pR*SplittingProb(jJ,jL,jR)*Sudakov(jJ);
                }//end of if for anti-quark parent
                else if(jJ.f==Higgs&&(fabs(sqrt(square(jL.p+jR.p))-para.massH)<para.DeltaMH))
                {
                    pL=0;
                    pR=0;

                    //H->q + qb
                    jL.f=bottom;
                    jR.f=bbar;

                    jL.yccL=para.verybig;
                    jL.yccR=jR.y;
                    jL.phiccL=para.verybig;
                    jL.phiccR=jR.phi;

                    jR.yccL=jL.y;
                    jR.yccR=para.verybig;
                    jR.phiccL=jL.phi;
                    jR.phiccR=para.verybig;

                    //end of setting parameters for H->q + qb

                    pL=SplitJet(jL);
                    pR=SplitJet(jR);

                    hprob=pL*pR*16*pow(para.pi,2.0)/(4*para.massH*para.DeltaMH);

                    prob=prob+hprob;
                }//end of if for Higgs parent

            }//end of if for set length of 0

        }//end of long while for splitting initial jetdata's momenta
    }//end of long else for jetdata length >1

    if(jJ.MsqK==para.verybig&&jJ.set[0]!=1)//To account for not doing half of the splittings
        prob=prob*2.0;

    return prob;

}//end of SplitJet

double Gstart(const jetdata& jJ)                             //Retruns prob for main jet starting with G
{
    int s1[maxsize+1],sL[maxsize+1],sR[maxsize+1];
    double prob=0.0,pL=0,pR=0,Qsq,cfactor=0;
    bool EndofSplit;
    parameters para;

    Qsq=jJ.ksq+jJ.musq;

    s1[0]=jJ.set[0];

    for(int i=1;i<=s1[0];i++)
        s1[i]=i;

    EndofSplit=false;
    LRsets LR1(s1);

    //while for splitting the inital set of jets into a hard and initial state sets
    while(!EndofSplit)
    {
        EndofSplit=LR1.next(sL,sR);//sL is the initial state and sR is the hard jet

        if(sR[0]!=0)
        {
            if(sL[0]!=0)
            {
                jetdata jL(jJ,sL);
                jetdata jR(jJ,sR);

                jR.f=gluon;
                jL.f=gluon;

                jR.MsqK=para.verybig;
                jR.kTK=para.verybig;
                jL.MsqK=para.verybig;
                jL.kTK=para.verybig;

                if(jL.ksq<Qsq/4.0)//Kinematic check
                {
                    pL=InitStateRad(jL,Qsq);
                    pR=SplitJet(jR);
                    cfactor=para.npdfg*pow(pow(para.pTmin,2.0)/jR.ksq,para.npdfg)*1.0/jR.ksq;
                }
                else
                {
                    pL=0;
                    pR=0;
                    cfactor=0;
                }
            }//end of IF for no initial state radiation
            else
            {
                jetdata jR(jJ,sR);
                jR.f=gluon;
                jR.MsqK=para.verybig;
                jR.kTK=para.verybig;

                pR=SplitJet(jR);
                pL=1.0;

                cfactor=para.npdfg*pow(pow(para.pTmin,2.0)/jR.ksq,para.npdfg)*1.0/jR.ksq;
            }

            prob=prob+pL*pR*cfactor;
        }//end of IF for nothing in hard jet
    }//end of while for looping through possible combinations of hard and init-state radiation

    return prob;
}//end of Gstart

double Hstart(const jetdata& jJ)                             //Retruns prob for main jet starting with H
{
    int s1[maxsize+1],sL[maxsize+1],sR[maxsize+1];
    double prob=0.0,pL=0,pR=0,cfactor=0;
    double mHsq,Qsq;
    bool EndofSplit;
    parameters para;

    Qsq=jJ.ksq+jJ.musq;
    mHsq=pow(para.massH,2.0);

    s1[0]=jJ.set[0];

    for(int i=1;i<=s1[0];i++)
        s1[i]=i;

    EndofSplit=false;
    LRsets LR1(s1);

    //while for splitting the inital set of jets into a hard and initial state sets
    while(!EndofSplit)
    {
        EndofSplit=LR1.next(sL,sR);//sL is the initial state and sR is the hard jet

        if(sR[0]!=0)
        {
            if(sL[0]!=0)
            {
                jetdata jL(jJ,sL);
                jetdata jR(jJ,sR);

                jR.f=Higgs;
                jL.f=gluon;

                jR.MsqK=para.verybig;
                jR.kTK=para.verybig;
                jL.MsqK=para.verybig;
                jL.kTK=para.verybig;

                if(jL.ksq<Qsq/4.0)//Kinematic check
                {
                    pL=InitStateRad(jL,Qsq);
                    pR=SplitJet(jR);
                    cfactor=para.npdfH*pow((pow(para.pTmin,2.0)+mHsq)/(jR.ksq+mHsq),para.npdfH)*1.0/(jR.ksq+mHsq);
                }
                else
                {
                    pL=0;
                    pR=0;
                    cfactor=0;
                }
            }//end of IF for no initial state radiation
            else
            {
                jetdata jR(jJ,sR);
                jR.f=Higgs;
                jR.MsqK=para.verybig;
                jR.kTK=para.verybig;

                pR=SplitJet(jR);
                pL=1.0;

                cfactor=para.npdfH*pow((pow(para.pTmin,2.0)+mHsq)/(jR.ksq+mHsq),para.npdfH)*1.0/(jR.ksq+mHsq);
            }

            prob=prob+pL*pR*cfactor;
        }//end of IF for nothing in hard jet
    }//end of while loop over possible combinations of hard and init-state radiation

    return prob;
}//end of Hstart

double InitStateRad(const jetdata& jJ, double& Qsq)                       //Returns prob for initial state radiation
{
    double prob=0.0,hold=1.0,His=0.0,H1,H2,Q;
    int s1[maxsize+1],sPart[maxsize+1][maxsize+1];
    bool EndofParts=false;
    parameters para;

    Q=sqrt(Qsq);

    s1[0]=jJ.set[0];
    for(int i=1;i<=s1[0];i++)
        s1[i]=i;

    Partitions P1(s1);

    //while loop for generating all possible combinations of partitions of the initial state microjets
    while(!EndofParts)
    {
        EndofParts=P1.next(sPart);

        hold=1.0;
        for(int i=1;i<=sPart[0][0];i++)     //for loop for generating all the jetdatas' momenta
        {
            jetdata jP(jJ,sPart[i]);
            jP.f=gluon;
            jP.MsqK=para.verybig;
            jP.kTK=para.verybig;

            H1=8.0*para.pi*para.cA*alphas(jP.ksq+para.kappapsq)/(jP.ksq+para.kappapsq)*1.0/pow(1+para.cR*jP.k/Q,para.nR);
            H2=16.0*para.pi*para.cnp*pow(para.kappanpsq,para.nnp-1.0)/pow(jP.ksq+para.kappanpsq,para.nnp);

            His=H1+H2;

            hold=hold*His*SplitJet(jP);
        }

        prob=prob+hold;
    }//end of while

    return prob;

}//end of InitStateRad



double FinalState(const jetdata& jF)                          //Returns prob for final state
{
    double pS,pb=1.0;
    parameters para;

    pS=Sudakov(jF);//final state sudakov factor

    // B-tag probablilities section
    switch (jF.mList[0].btag)
    {
        case -1:
            if(jF.f==bottom||jF.f==bbar)
                pb=1.0-para.PTb;
            else
                pb=1.0-para.PTnotb;
            break;
        case 0:
            pb=1.0;
            break;
        case 1:
            if(jF.f==bottom||jF.f==bbar)
                pb=para.PTb;
            else
                pb=para.PTnotb;
            break;
    }

    return pS*pb;
}//end of FinalState
//////////////////////////////////////End of Splitting Subroutines/////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
