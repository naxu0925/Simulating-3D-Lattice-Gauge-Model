
// This program calculates the bin average for the original dataset

#include <iostream>
#include <fstream>
#include <math.h>
#include <cmath>
#include <iomanip>
#include <cassert>
#include <string>
#include <strstream>
using namespace std; 

const int NBIT = 64;
const int nsample = 1000;

int dd, ll, nbin, nline,nn;
int qtime, qpnt, tot_pnt;
int init, istp, mstp;
double jj1=-1.0e0, ti,tf;
string filename, confname, outputname;

int *tm_table;
double *tt_table;
double **data;   
double *enrg, *enrg2; 
double *maga, *mag2;
double *data_cv, *data_chi; 
double *tx, *ty, *tz, *t_order;
// variables for Random Numer Generator
int iir, jjr, kkr, nnr, seed[4];

void readinput(void);
void readdata(void);
void bootstrapping(double vec[], int n, int nboot, double& ave, double& err);
void initran(void);
double ran(void);

void aveanderr(double [],int ,double &, double &);
void aveanderr(long double [],int ,long double &, long double &);

   
int main () { 
  int i, t;
  double beta; 
  double av1, av2, av3, av4, av5, av6, av7, av8, av9;
  double er1, er2, er3, er4, er5, er6, er7, er8, er9;
  ofstream datafl;

  initran();

  readinput();
  readdata();
  outputname=filename;
  outputname.erase(0,4);
  outputname.insert(0,"res");
  datafl.open(outputname.c_str(), ios::app);
  for (t = 0; t < tot_pnt; ++t) {
    beta = 1.0e0/tt_table[t];

    for (i = 0; i < nbin; ++i) {
      enrg[i] = data[ t + i * tot_pnt ][0];
      enrg2[i] = data[ t + i * tot_pnt ][1];
      data_cv[i] = ( enrg2[i] - pow(enrg[i],2) );
      maga[i] = data[ t + i * tot_pnt ][2];
      mag2[i] = data[ t + i * tot_pnt ][3];
      data_chi[i] = ( mag2[i] - pow(maga[i],2) );
        tx[i] = abs(data [ t + i * tot_pnt][4]);
        ty[i] = data [ t + i * tot_pnt][5];
        tz[i] = data [ t + i * tot_pnt][6];
        t_order[i] = data [ t + i * tot_pnt][7];
    }

    bootstrapping(enrg, nbin, nsample, av1, er1);
    bootstrapping(data_cv, nbin, nsample, av2, er2);
    av2 *= 3.0*double(nn)*pow(beta,2);
    er2 *= 3.0*double(nn)*pow(beta,2);

    bootstrapping(mag2, nbin, nsample, av3, er3);
    bootstrapping(data_chi, nbin, nsample, av4, er4);
    av4 *= 3.0*double(nn)*beta;
    er4 *= 3.0*double(nn)*beta;
   // aveanderr(maga, nbin, av5, er5);
    bootstrapping(tx, nbin, nsample, av6, er6);
   // bootstrapping(ty, nbin, nsample, av7, er7);
   // bootstrapping(tz, nbin, nsample, av8, er8);
    bootstrapping(t_order, nbin, nsample, av9, er9);
      
      datafl<< setw(18) << setprecision(14) << tt_table[t]
	   << setw(20) << setprecision(14) << av1
	   << setw(20) << setprecision(14) << er1 
	   << setw(22) << setprecision(14) << av2
	   << setw(20) << setprecision(14) << er2
	   << setw(20) << setprecision(14) << av3
	   << setw(20) << setprecision(14) << er3
	   << setw(22) << setprecision(14) << av4
	   << setw(20) << setprecision(14) << er4
       << setw(20) << setprecision(14) << av9
       << setw(20) << setprecision(14) << er9
       << setw(20) << setprecision(14) << av6
       << setw(20) << setprecision(14) << er6
       << endl;

  } // end of t-loop      
  datafl.close();

  for (i = 0; i < nline; ++i) delete [] data[i]; delete [] data;
  delete [] enrg;
  delete [] enrg2;
  delete [] data_cv;
  delete [] maga;
  delete [] mag2;
  delete [] data_chi;
  delete [] tt_table;
  delete [] tx;
  delete [] ty;
  delete [] tz;
  delete [] t_order;
    return 0;
}//end of 'main' 




/*----------------------------------*/
void readdata(void){
/*----------------------------------*/
 fstream datafile; 
 int i;
 double tmp1,tmp2, tmp3, tmp4, tmp5, tmp6,tmp7,tmp8,tmp9,tmp10;
    
 datafile.open(filename.c_str()); 
 i = 0;
 if ( datafile.is_open() ) {
   while ( datafile >> tmp2 >> tmp3 >> tmp4 >> tmp5 >> tmp6 >> tmp7 >> tmp8 >> tmp9 >> tmp10 ) ++i;
 } else {
   cout << "can't open file data.txt" << endl; 
   exit (1);
 }
 datafile.close();  

 nline = i;
 if (nline%tot_pnt != 0) {
   cout << "# of data is WRONG! \n";
   exit (1);
 }

 nbin = nline/tot_pnt;

 data = new double *[nline];
 for (i = 0; i < nline; ++i) data[i] = new double [8];

 tt_table = new double [tot_pnt];
 enrg = new double [nbin];
 enrg2 = new double [nbin];
 data_cv = new double [nbin];
 maga = new double [nbin];
 mag2 = new double [nbin];
 data_chi = new double [nbin];
 tx= new double [nbin];
 ty= new double [nbin];
 tz= new double [nbin];
 t_order= new double [nbin];

 datafile.open(filename.c_str());
 for ( i = 0; i < nline; ++i ) {
   datafile >> tmp2 >> data[i][0] >> data[i][1] >> data[i][2] >> data[i][3]
            >>data[i][4]>>data[i][5]>>data[i][6]>>data[i][7];
 }
 datafile.close();

 datafile.open(filename.c_str());
 for ( i = 0; i < tot_pnt; ++i ) { 
   datafile >>tt_table[i] >> tmp3 >> tmp4 >> tmp5 >> tmp6 >> tmp7 >> tmp8
            >> tmp9 >> tmp10;
 }
 datafile.close();

}//end of function 'readdata'



/*----------------------------------*/
void readinput(void){
/*----------------------------------*/
//  Read file input.txt,
//    which has the format
//       ll dd
//      init istp mstp nbin
//      qtime  qpnt  
//     filename confname
//     ti tf
/*----------------------------------*/
    ifstream inputfile("input.txt");
    
    if(inputfile.is_open())
    {
        while(!inputfile.eof())
        {
            inputfile >> ll >> dd;
            inputfile >> init >> istp >> mstp >> nbin;
            inputfile >> qtime >> qpnt;
            inputfile >> filename >>confname;
            inputfile >> ti >> tf;
        }
    }
    else{
        cout<<"Can not open input.txt !" << endl;
        exit(1); }
    inputfile.close();

    assert((qtime >= qpnt) && (qtime%qpnt == 0));

  nn = int(pow(1.0*ll,1.0*dd));
  tot_pnt = qpnt + 1;

}//end of function 'readinput' 



/*------------------------------------*/
void bootstrapping(double vec[], int n, int nboot, double& ave, double& err)
/*------------------------------------*/
{
  int i, j;
  double sum;

  ave = 0.0e0;
  for (i = 0; i < n; ++i) ave += vec[i];
  ave /= double(n);

  err = 0.0e0;
  for ( j = 0; j < nboot; ++j) {
    sum = 0.0e0;
    for ( i = 0; i < n; ++i) sum += vec[ int(ran()*n) ];
    sum /= double(n);
    err += (sum - ave)*(sum - ave);
  }
  err = sqrt(err/double(nboot));

}



/*------------------------------------*/
void aveanderr(double vec[], int n, double& ave, double& err){
/*------------------------------------*/

 int i;
 
 ave = 0.0e0; err = 0.0e0;
 for ( i = 0; i < n; ++i ) {
       ave += vec[i];
       err += vec[i]*vec[i];
 }

 ave/=double(n);
 err/=double(n);
 err = sqrt( (err-ave*ave)/double(n-1) ) ; 
}//end of 'aveanderr'



/*------------------------------------*/
void aveanderr(long double vec[], int n,long double& ave,long double& err){
/*------------------------------------*/

 int i;
 
 ave=0.0e0; err=0.0e0;
 for ( i=0; i<n; i++ ) {
       ave+=vec[i];
       err+=vec[i]*vec[i];
 }

 ave/=double(n);
 err/=double(n);
 err = sqrt( fabs(err-ave*ave)/double(n-1) ) ; 
}//end of 'aveanderr'



void initran(void) {
  ifstream seedfile ("seed.txt");
  int i;

  if ( seedfile.is_open() ) {
    while ( ! seedfile.eof() ){
      for ( i=0; i<4; i++ ) seedfile >> seed[i];
    }
  } else {
    cout << "Unable to open the seed file" << endl;
    exit (1);
  } //end of 'if seedfile.is_open()'
  seedfile.close();

  iir=1+int(fabs(seed[0]));
  jjr=1+int(fabs(seed[1]));
  kkr=1+int(fabs(seed[2]));
  nnr=seed[3];
}// end of member function 'Cadbt_qmc::initran'


double ran(void) {
  int mzran;
  double r;

  mzran=iir-kkr;
  if(mzran<0) mzran+=2147483579;
  iir=jjr;
  jjr=kkr;
  kkr=mzran;
  nnr=69069*nnr+1013904243;
  mzran=mzran+nnr;
  r=0.5e0+double((0.23283064e-9)*mzran);

  return r;
}// end of member function 'Cadbt_qmc::ran'
