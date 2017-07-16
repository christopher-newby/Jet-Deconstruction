// Specification file for classes etc.
// used for shower deconstruction code
// D.E. Soper, C. Newby, C. Jackson
// 6 November 2011

using namespace std;
const int maxsize = 12;

class LRsets
{
// Uses arrays of integers of size maxsize+1 to represent sets of integers.
// set[0]=N is the size of the set,
// {set[1],...,set[N]} are the elements.
// Class LRsets returns subsets setA and setB of set, where setA U setB = set.
public:
  LRsets();                                                // Default constructor with set {1,2,3}
  LRsets(int (&set)[maxsize+1]);                           // Constructor.  Use set passed to the constructor.
  void reset(int (&set)[maxsize+1]);                       // Resets data according to set, as in constructor.
  bool next(int (&setA)[maxsize+1],int (&setB)[maxsize+1]);// Returns next division of set into setA and setB; value = 0.
                                                           // First call of s.next(setA,setB) gives setA = set, setB = {}.
                                                           // Further calls give each subset of "set" for setA, the complement for setB.
                                                           // Call 2^N is the last split into and returns setA = {}, setB = set, value = 1.
                                                           // Further calls are in error and return setA = {}, setB = {}, value 1.
private:
  int len;                //length of the whole set
  int intset[maxsize];    //interal set which is either given or set by the default compiler
  unsigned long count;    //a count of which step in the splitting the class is at
  bool countend;          //boolian for when the count is done
};

class Partitions
{
// Uses arrays of integers of size maxsize+1 to represent different sets of integers.
// set[0]=N is the size of the set on input,
// {set[1],...,set[N]} are the elements.
// The Class Partitions returns (through next) a 2D array which contains the info of how many parts it contains,
//          the length of each part, and the specific elements in that part
//      ->Notation for returned array is:
//                                          set[0][0]=number of partitions
//                                          set[i][0]=length of partition i
//                                          set[i][j]=jth element of partition i
    public:
        Partitions(int (&set)[maxsize+1]);              // Constructor which uses set passed to it
        void reset(int (&set)[maxsize+1]);              // Resets the class with the input set
        bool next(int (&set)[maxsize+1][maxsize+1]);    // cycles through the partitions of input set
                                                        // Returns a 2D array of the specific partition set you are on
                                                        // governed by the rules laid out at beginning of class.
                                                        // Returns a boolian which if false until you reach the end of the
                                                        // total number of possible partitions, then it returns true.
                                                        // If you try to call after all partitions, then it returns empty sets.
    private:
        int intset[maxsize+1];                          // Internal input set
        int count[maxsize+1];                           // Array of counts, location depends upon which level of the splitting you are in
        int maxcount[maxsize+1];                        // Array of maximum number of counts for each splitting level
        void nestwhile(int (&TotSet)[maxsize+1][maxsize+1],int Rset[maxsize+1]);    // Function which generates the elements of the partitions
                                                                                    // calls itself recursively
        int nestLevel;                                  // The level of splitting (or how many times you have recursively called nestwhile)
        bool countend;                                  // Is the return value for next()
        int BellNumber(int a);                          // This subroutine calculates the total number of partitions
        int BinomialCoeff(int a,int b);                 // This returns the binomial coefficient of a choose b
        int factorial(int a);                           // This returns a!
        int MaxBell;                                    // Is the maximum bell number for the input set length
        int c;                                          // Internal count which counts up to the Bell number
};

class fourvector
{
public:
  fourvector();				                                // Default constructor, sets all values to 0
  fourvector(double invec[4]);                              //Constructor which sets the internal vector v(i)=invec[i]
  fourvector(double t, double x, double y, double z);       //Constructor which sets the internal vector v(i)=(t,x,y,z)
  double& in(const int& mu);						        // v.in(mu) for input, v.in(mu) = r
  const double operator()(const int& mu) const;	        // v(mu) for output,  r = v(mu)
  fourvector operator+ (const fourvector& v2) const;        // Sum of two vectors: v3 = v1 + v2
  fourvector operator- (const fourvector& v2) const;        // Difference of two vectors: v3 = v1 - v2
  fourvector operator* (const double& a) const;		        // Vector times scalar: v2 = v1*c
  fourvector operator- () const;		                    // -1 times vector: v2 = - v1
private:
  double vector[4];				                            // The vector, v[0], v[1], v[2], v[3]
};

// Scalar times vector: v2 = c*v1
fourvector operator*(const double& a, const fourvector& v);
// Dot product of two 4-vectors: dot(v1,v2)
double dot(const fourvector& v1, const fourvector& v2);
// Square of one 4-vector: square(v) = dot(v,v)
double square(const fourvector& v);

enum flavor {gluon,quark,qbar,bottom,bbar,Higgs};           //Enumberated data type for the flavor of the jet

class microjet
{
public:
  microjet();                                                   // Default constructor which sets p to (0,0,0,0) and btag to 0
  microjet(fourvector pinit,int btinit);                        // Constructor which sets btag to btinit and p to pinit
  microjet(double t, double x, double y, double z, int btinit); // Constructor which sets btag to btinit and p to (t,x,y,z)
  microjet(double invec[4], int btinit);                        // Constructor which sets btag to btinit and p to invec (with invec[0]=t,...)
  fourvector p;     //internal momentum of the microjet
  int btag;         //whether the microjet was labeled as a b or not
                    //      btag=0 means not even attempted to tag
                    //      btag=-1 is not tagged as b
                    //      btag=1 means was tagged as b
};

class jetdata
{
public:
  jetdata();              // Default constructor which sets all values, except for momenta, to the very big number from parameters
  jetdata(microjet mIn[maxsize], int length);// Constructor for first call of jetdata
  jetdata(jetdata jK, int s1[maxsize+1]);    // Initializes all the momentum parts of jetdata
  int set[maxsize+1];     // List of constituents of the jet.
  microjet mList[maxsize];// List of all microjets in jet (starting from i=0 to maxsize)
  fourvector p;           // Momentum of the jet.
  flavor f;               // Flavor of the jet.
  double yccL;            // rapidity of the left color connected partner of the jet.
  double phiccL;          // aximuthal angle of the left color connected partner of the jet.
  double yccR;            // rapidity of the right color connected partner of the jet.
  double phiccR;          // aximuthal angle of the right color connected partner of the jet.
  double MsqK;            // squared mass of the mother of the jet.
  double kTK;             // absolute value of the transverse momentum of the mother of the jet.

/////////////////////Commonly used parameters
  double k,ksq;             // Absolute value of transverse momentum and squared
  double musq;              // Virtuality
  double y;                 // rapidity
  double phi;               // Azimuthal angle
private:
  void calculatingParams();     //calculates k, musq, y, and phi
};

class parameters
{
public:
  parameters();
  double pi;              // 3.1415926535898.
  double MassZ;           // Mass of Z-boson, 91.188 GeV.
  double asMZ;            // Alpha_s at MZ = 0.118;
  double cA;              // SU(3) group theory parameter, = 3.
  double cF;              // SU(3) group theory parameter, = 4/3.
  double TR;              // SU(3) group theory parameter, = 1/2.
  int Nf;                 // Number of quark flavors, normally 5.
  double b0;              // Coefficient in beta function, (33 - 2*Nf)/12/pi.
  double verybig;         // A large number.
  double Rfatjet;         // Radius of the fat jet, default 1.2.
  double Rfatjetsq;       // Square of Rfatjet.
  double pTmin;           // Smallest pT of Z-boson allowed, default 200.0.
  double npdfg;           // Power from pdfs for hard scattering, background process; default 2.0.
  double npdfH;           // Power from pdfs for hard scattering, signal process; default 2.0.
  double kappapsq;        // Cutoff parameter for perturbative IS radiaton; default 4.0;
  double kappanpsq;       // Cutoff parameter for non-erturbative IS radiaton; default 4.0;
  double cR;              // Scale parameter for perturbative IS radiation, default 2.0.
  double nR;              // Denominator power for perturbative IS radiation, default 1.0.
  double cnp;             // Coefficient for non-perturbative IS radiation, default 1.0.
  double nnp;             // Denominator power for non-perturbative IS radiation, default 1.5.
  double massH;           // Mass of Higgs boson, default 120.0.
  double DeltaMH;         // Allowed variation of Higgs mass, default 10.0.
  double pTtagcut;        // Minimum pT of microjets that can carry b-tags, default 15 GeV
  double PTb;             // Probability for positive b-tag of a b-quark; default 0.6.
  double PTnotb;          // Probability for positive b-tag of parton not a b-quark; default 0.02;
};

double alphas(double Qsq); // Alpha_s as a function of the scale.

double Sudakov(const jetdata& jJ);          // Sudakov exponent

double SplittingProb(const jetdata& jJ, const jetdata& jL, const jetdata& jR);       // Probability factor at vertex

double anglefactor(const jetdata& jL, const jetdata& jR, const jetdata& jK);          // Angle factor g(y,phi) with jets.

double Dphi(const double& phi1, const double& phi2);            //calculates the difference between two phi's and returns a value between -pi and pi


double /*Observable*/Chi(microjet inMjet[maxsize], int length);  //Returns the observable for an initial list of microjets
double SplitJet(const jetdata& jJ);                              //Splits the jet given to it and returns the probability
double Gstart(const jetdata& jJ);                                //Retruns prob for main jet starting with G
double Hstart(const jetdata& jJ);                                //Retruns prob for main jet starting with H
double InitStateRad(const jetdata& jJ, double& Qsq);             //Returns prob for initial state radiation
double FinalState(const jetdata& jF);                             //Returns prob for final state (should only send a jet with one momenta from initial set)
