//File: to.cxx
// 1. This is a 3D toric code model with multispin coding and Metropolis algorithm.
// J1<0, ferromagnetic
// spins are in variables spn[i][3],i~[0-nn]
// By naxu 12/26/2013, PTL!

//g++ to_quench.cxx to.cxx

#include <cstdlib>
#include <stdlib.h>
#include <iostream>
#include <fstream> // Provides fstream
#include <math.h>  // Provides pow(), fabs(), sqrt(), etc.
#include <time.h>  // Provides time()
#include <cassert> // Provides assert()
#include <iomanip> // Provides setw(), setprecision()
#include <string>
#include <strstream>
#include "corr.h" 
using namespace std;

namespace naxu_to_quench_namespace{
   
to_quench::~to_quench(void)
    {
        int i;
        for(i=0;i<3;i++)
            delete []spn[i];
    /*    delete [] enrj;
        delete [] enrj2;
        delete [] maga;
        delete [] mag2;
        delete [] tm_table;
        delete [] tx;
        delete [] ty;
        delete [] tz;
        delete [] t_order2;
     */
        delete [] corr;
        
    }
    
void to_quench:: writedata(void)
    {
        ofstream outfile;
        outfile.open(filename.c_str(), ios::app);
        int t;
        if(outfile.is_open())
        {
            for (t = 0; t < tot_pnt; ++t) {
             /*   enrj[t] /= double(mstp);
                enrj2[t] /= double(mstp);
                maga[t] /= double(mstp);
                mag2[t] /= double(mstp);
                tx[t]/=double(mstp);
                ty[t]/=double(mstp);
                tz[t]/=double(mstp);
                t_order2[t]/=double(mstp);
                outfile<< setw(20) << setprecision(12) << double(1.0)/double(qtime)
                <<setw(20) << setprecision(12) << tt_table[t]
                << setw(20) << setprecision(12) << enrj[t]
                << setw(22) << setprecision(12) << enrj2[t]
                << setw(20) << setprecision(12) << maga[t]
                << setw(22) << setprecision(12) << mag2[t]
                << setw(20) << setprecision(12) << tx[t]
                << setw(20) << setprecision(12) << ty[t]
                << setw(20) << setprecision(12) << tz[t]
                << setw(22) << setprecision(12) << t_order2[t]<< endl;
              */
                outfile<<setw(20)<<setprecision(12)<<double(1.0)/double(qtime)
                       <<setw(20)<<setprecision(12)<<tt_table[t]
                       <<setw(24)<<setprecision(14)<<corr[t]/double(mstp)
                <<endl;
            }
        } else {
            cout << "can't open file" << filename << endl;
            exit (1);
        }
        outfile.close();
        
    }


    
    
    
void to_quench:: clean(void)
    {
        int t;
        for(t=0;t<tot_pnt;t++)
        {
           // enrj[t]=enrj2[t]=maga[t]=mag2[t]=tx[t]=ty[t]=tz[t]=t_order2[t]=0.0e0;
            corr[t]=0.0e0;
           
        }
        
    }
 

    void to_quench:: measure(int t)
    {
        int i,j,k,mm, mm_new, too[3], too2[3];
        double tty=0.0e0;
        
        
    /*    ej1=energy();
        enrj[t]+=double(ej1)/double(3*nn);
        enrj2[t]+=pow(double(ej1)/double(3*nn),2.0e0);
        m1=0;
        for(i=0;i<3;i++) to[i]=0.0e0;
        
       for(k=0;k<ll;k++)
        {
            
            for(j=0;j<ll;j++)
            {
                too[0]=1; too[1]=1; too[2]=1;
                for(i=0;i<ll;i++)
                {
                    mm=num(i,j,k);
                    too[0]*=spn[mm][0];
                    too[1]*=spn[mm][1];
                    too[2]*=spn[mm][2];
                    m1+=spn[mm][0]+spn[mm][1]+spn[mm][2];
                }
               to[0]+=double(too[0])/double(ll*ll);
               
            }
            
            for(i=0;i<ll;i++)
            {
                too[0]=1;too[1]=1; too[2]=1;
                for(j=0;j<ll;j++)
                {
                    mm=num(i,j,k);
                    too[0]*=spn[mm][0];
                    too[1]*=spn[mm][1];
                    too[2]*=spn[mm][2];
                }
                to[1]+=double(too[1])/double(ll*ll);
            }

        }
        
        for(i=0;i<ll;i++)
        {
            for(j=0;j<ll;j++)
            {
              too[0]=1;too[1]=1; too[2]=1;
                for(k=0;k<ll;k++)
                {
                    mm=num(i,j,k);
                    too[0]*=spn[mm][0];
                    too[1]*=spn[mm][1];
                    too[2]*=spn[mm][2];
                }
                to[2]+=double(too[2])/double(ll*ll);
            }
        }
     
    maga[t]+=double(m1)/double(nn*3);
    mag2[t]+=pow(double(m1)/double(nn*3),2.0e0);
        tx[t]+=to[0];
        ty[t]+=to[1];
        tz[t]+=to[2];
        t_order2[t]+=(to[0]*to[0]+to[1]*to[1]+to[2]*to[2])/3.0e0; */
      
        for(k=0;k<ll;k++)
        {
            for(i=0;i<ll;i++)
            {
                too[1]=1; too2[1]=1;
                for(j=0;j<ll;j++)
                {
                    mm=num(i,j,k);
                    too[1]*=spn[mm][1];
                    mm_new=num((i+ll/2)%ll,j,(k+ll/2)%ll);
                    too2[1]*=spn[mm_new][1];
                }
                tty=tty+double(too[1]*too2[1])/double(ll*ll);
            }
        }
        
        corr[t]=corr[t]+tty;
        
        
    }
    
void to_quench::newrun(void)
    {
     int j, steps;
        
        reset_temp();
        
        steps = 32;
        readconf();
        for (j = 0; j < steps; ++j) mcsweep();
        writeconf();
        
    } // end of member function 'newrun'
    
    
  
void to_quench:: decrease_temp(void)
    {
        tem=tem-dt;
    }
    
void to_quench::reset_temp(void)
    {
        tem=t_i;
    }
 
 void to_quench:: mcsweep(void)
    
    {
        int i,j,k;
        int d_e;
        double r;
        for(k=0; k<3*nn; k++)
        {
            i=int(nn*ran());
            j=int(3*ran());
            d_e=e_change(i,j);
          if(d_e!=0)
          {
            r=ran();
            if(r<prob(d_e))
                spn[i][j]=-1*spn[i][j];
          }
            else
            {
              if(ran()>0.5)
                  spn[i][j]=-1*spn[i][j];
            }
        }
        
    }

    
  int to_quench:: e_change(int mm, int dim) // Calculate the change in energy when flip one spn[mm][0~2]
  {
    int delta_e=0;
    int i,j,k,i1,j1,k1,i2,j2,k2;
    k=mm/(ll*ll);
    j=mm%(ll*ll)/ll;
    i=mm%ll;
    i1=(i+1)%ll;
    j1=(j+1)%ll;
    k1=(k+1)%ll;
    i2=(i+ll-1)%ll;
    j2=(j+ll-1)%ll;
    k2=(k+ll-1)%ll;
    if(dim==0)
    {
        delta_e=2*spn[mm][0]*(spn[num(i,j1,k)][0]*spn[mm][1]*spn[num(i1,j,k)][1]
                              +spn[num(i,j2,k)][1]*spn[num(i,j2,k)][0]*spn[num(i1,j2,k)][1]
                              +spn[mm][2]*spn[num(i1,j,k)][2]*spn[num(i,j,k1)][0]
                              +spn[num(i,j,k2)][2]*spn[num(i,j,k2)][0]*spn[num(i1,j,k2)][2] );

        
    }
    else if(dim==1)
    {
        delta_e=2*spn[mm][1]*(spn[num(i1,j,k)][1]*spn[mm][0]*spn[num(i,j1,k)][0]
                              +spn[num(i2,j,k)][0]*spn[num(i2,j,k)][1]*spn[num(i2,j1,k)][0]
                              +spn[mm][2]*spn[num(i,j,k1)][1]*spn[num(i,j1,k)][2]
                              +spn[num(i,j,k2)][2]*spn[num(i,j,k2)][1]*spn[num(i,j1,k2)][2] );
    }
    else if(dim==2)
    {
        delta_e=2*spn[mm][2]*(spn[mm][0]*spn[num(i1,j,k)][2]*spn[num(i,j,k1)][0]
                              +spn[num(i2,j,k)][0]*spn[num(i2,j,k)][2]*spn[num(i2,j,k1)][0]
                              +spn[mm][1]*spn[num(i,j1,k)][2]*spn[num(i,j,k1)][1]
                              +spn[num(i,j2,k)][1]*spn[num(i,j2,k)][2]*spn[num(i,j2,k1)][1]);
        
    }
      
      
      
    return delta_e;
  }

  
  double to_quench::prob(int d_e)
  {
    double probability;
      probability=exp(-1.0*double(d_e)/double(tem));
     return probability;
  }



 int to_quench :: num(int i1, int j1, int k1)
  {
    int mm;
    mm=i1+j1*ll+k1*ll*ll;
    return mm;
  }


  double to_quench :: energy(void)  // Calculate the energy of the system(j0=-1 included);
  {
    int i,j,k;
    int ii,jj,kk;
    int m;
    double ee;
    
      ee=0.e0;
    for(k=0; k<ll; k++)
      {
      for(j=0; j<ll; j++)
	{
        for(i=0; i<ll; i++)
          {
            m=i+j*ll+k*ll*ll;
            ii=(i+1)%ll;
	        jj=(j+1)%ll;
	        kk=(k+1)%ll;
            ee+=-1.0*(spn[m][0]*spn[m][1]*spn[num(ii,j,k)][1]*spn[num(i,jj,k)][0]
	      +spn[m][1]*spn[m][2]*spn[num(i,jj,k)][2]*spn[num(i,j,kk)][1]
                     +spn[m][0]*spn[m][2]*spn[num(ii,j,k)][2]*spn[num(i,j,kk)][0]);
             	  }
	}
      }
     
      return ee;
        }

    void to_quench :: gaugeupdate(int mm)
    {
        int i,j,k,i1,j1,k1,i2,j2,k2;
        int dim;
        k=mm/(ll*ll);
        j=mm%(ll*ll)/ll;
        i=mm%ll;
        i2=(i+ll-1)%ll;
        j2=(j+ll-1)%ll;
        k2=(k+ll-1)%ll;
        for(dim=0;dim<3;dim++)
            spn[mm][dim]=-1*spn[mm][dim];
        spn[num(i2,j,k)][0]=-1*spn[num(i2,j,k)][0];
        spn[num(i,j2,k)][1]=-1*spn[num(i,j2,k)][1];
        spn[num(i,j,k2)][2]=-1*spn[num(i,j,k2)][2];
        
        
    }
    

  
void to_quench::initialize(void)
  {
    initran(); modify_seed();//initialize random numbers;
    //-------------------------------------------------//
      readinput(); // readinput
    nn=int(pow(ll,dd));
    tot_pnt=qpnt+1;
    dt=(t_i-t_f)/qtime;
    //------------------------------------//
    allocate_arrays(); // allocate arrays
    //--------------------------------------//
    t_table();
    tem=t_i;   //generate temperature table
    //-----------------------------------------//
    if(init==0)  //Initial configuration
      initconf();
      else
	readconf();
  } // End of initialization


void to_quench::writeconf(void) const
  {
    ofstream conf;
	conf.open(confname.c_str(), ios::trunc);
    int i;

    if ( conf.is_open() ) {
      for (i = 0; i < nn; ++i) conf << spn[i][0] <<"\n"<< spn[i][1]<<"\n" << spn[i][2]<<"\n";
    } else {
      cout << "In writeconf, can't open conf.txt" << endl;
      exit (1);
    }
    conf.close();
  } // end of member function 'writeconf'



  void to_quench::readconf(void)
  {
    ifstream conf;
	conf.open(confname.c_str());
    int i;
    
    if ( conf.is_open() ) {
      for (i = 0; i < nn; ++i) conf >> spn[i][0] >> spn[i][1] >> spn[i][2];
    } else {
      cout << "In readconf, can't open conf.txt" << endl;
      exit (1);
    }
    conf.close();
  } // end of member function 'readconf'



 
  void to_quench::initconf(void)
  {
    int i,j;
    double r;
    
    for ( i=0; i<nn; i++)
      {
        for ( j=0; j<3; j++)
	 {
         if(ran()<0.5)
             spn[i][j]=1;
         else
             spn[i][j]=-1;
	
	 }
        }
      
  }

  
void to_quench:: t_table(void)
  {
    int i,n=0;
    tem=t_i;
    for(i=0;i<=qtime;i++){
      if(i%(qtime/qpnt) == 0){
        tm_table[n] = i;
        tt_table[n] = tem;
        n++;
	  }
      tem = tem - dt;
    } 
    assert(n == qpnt+1);
  }


  void to_quench:: allocate_arrays(void)
  {
    int i;
    
    spn= new int *[nn];
    for(i=0;i<nn; i++) spn[i]=new int [3];
  /*  enrj = new double [tot_pnt];
    enrj2= new double [tot_pnt];
    maga = new double [tot_pnt];
    mag2 = new double [tot_pnt];
    tx=new double[tot_pnt];
    ty=new double[tot_pnt];
    tz=new double[tot_pnt];
    t_order2=new double[tot_pnt]; */
    tm_table= new int [tot_pnt];
    tt_table= new double[tot_pnt];
      corr= new double [tot_pnt];
  }
  

  
void to_quench:: readinput(void)
{
  // Read input file input.txt, which has the following format:
  // ll dd
  // init istp mstp nbin
  // qtime qpnt
  // filename confname
  // t_i t_f
  j0=-1.0e0;
  ifstream inputfile("input.txt");
  if(inputfile.is_open())
    {
      while(!inputfile.eof())
	{
	  inputfile >> ll >> dd;
          inputfile >> init >> istp >> mstp >> nbin;
          inputfile >> qtime >> qpnt;
          inputfile >> filename >>confname;
          inputfile >> t_i >> t_f;
	}
    }
  else{
    cout<<"Can not open input.txt !" << endl;
    exit(1); }
  inputfile.close();
}

void to_quench::modify_seed(void)
  {
    ofstream newseed ("seed.txt", ios::trunc);
    time_t systime;
    int i;

    systime = time(NULL);
    if ( newseed.is_open() ){
      for ( i=0; i<4; i++ ) newseed << int(systime * ran() * sqrt(2.0) ) << endl;
    } else {
      cout << "Can't open seed, can't modify" << endl;
      exit (1);
    }
    newseed.close();
  } // end of member function 'modify_seed'


double to_quench::ran(void)
  {
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
  } // end of member function 'ran'




  void to_quench::initran(void) 
  {
    ifstream seedfile("seed.txt");
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
  } // end of member function 'initran'


    void to_quench:: checkspn(void)
    {
        int i;
        for(i=0;i<nn;i++)
            cout<<spn[i][0]<<"/"<<spn[i][1]<<"/"<<spn[i][2]<<endl;
    }

}
