//File: corr.h
//Header file for the class to_quench
//By Na Xu, 4/23/2014


#include <strstream>
#include <string>
using namespace std;
namespace naxu_to_quench_namespace{

  class to_quench{

  public:
    // variables for Random Number Generator
    int iir, jjr, kkr, nnr, seed [4];
    
    // variables for Monte Carlo simulations
    int mstp, nbin; //monte carlo steps, bin numbers
    int qtime; // quench time
    int qpnt;  // # of measurement
    int istp;  // equilibration MC steps
    int init;  // if init=0 means stars from scratch, else start from saved conf
    
    void initran (void);
    void modify_seed(void);
    double ran(void);
    
    void initialize(void);
    void readinput(void);
    void allocate_arrays(void);
    void lattice(void);
    void t_table(void);
    void initconf(void);
    void readconf(void);
    void writeconf(void) const;
    void newrun(void);

    void mcsweep(void);
    void decrease_temp(void);
    void reset_temp(void);
    void clean(void);
    void checkstop(void) const;
    void measure(int t);
    void gaugeupdate(int mm);
    void writedata(void);
    double energy(void);
    int  e_change(int mm, int dim);
    int  num(int i1, int j1,int k1);
    double prob(int d_e);
    void checkspn(void);
    ~to_quench();
     
  private:
    int dd;        
    int ll;
    int nn;
    int nb;
    int tot_pnt;
    double j0;
    double t_i;
    double t_f;
    double tem;
    double dt;
    std:: string filename, confname;
    int **spn;
    

    double ej1;
    double m1;
    double to[3];
    int *tm_table;
    double *tt_table;
    double *enrj;
    double *enrj2;   
    double *maga;
    double *mag2;
  /*  double *tx, *ty, *tz;
    double *tx_2, *ty_2, *tz_2;
    double *t_order, *t_order2;
   */
      double *corr;


  } ;// end of of class "to_quench"

} // end of namespace "naxu_to_quench_namespace"


#endif
