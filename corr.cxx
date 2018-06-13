//
//  to.cxx
//  
//
//  Created by Na Xu on 1/3/14， PTL!
//
//

#include "corr.h"
#include <iostream>
#include <math.h>
using namespace std;
using namespace naxu_to_quench_namespace;

int main()
{
    to_quench quench;
    int interval, qtime, istp, mstp, nbin;
    int t,n,m,j,k;
    quench.initialize();
    interval=quench.qtime/quench.qpnt;
    istp=quench.istp;
    mstp=quench.mstp;
    nbin=quench.nbin;
    qtime=quench.qtime;
    for(j=0; j< istp; j++)
    {
        quench.mcsweep();
    }
    if(istp !=0)
        quench.writeconf();
    
    for(k=0; k<nbin; k++)
    {
        quench.clean();
        for(j=0; j<mstp; j++)
        {
            quench.newrun();
            n=0;
            for(t=0; t<=qtime; t++)
            {
                quench.mcsweep();
                if((t%interval)==0)
                {
                    quench.measure(n);
                    ++n;
                }
                quench.decrease_temp();
            }
        }
        quench.writedata();
    
    }
  //  quench.~to_quench();
    return 0;
}