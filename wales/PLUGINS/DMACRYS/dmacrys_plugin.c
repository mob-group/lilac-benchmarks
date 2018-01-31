#include <potential_plugin.h>

// initialize the plugin
void wales_initialize()
{
    static int init=0;

    if(init==0)
        dmacrys_setup_();
    init = 1;
}

// get the number of atoms
int wales_get_dof()
{
    return dmacrys_get_dof_();
}

// calculate potential
double wales_potential(double *x, double *grad, double *hess, int do_grad, int do_hess)
{
   double energy;
   dmacrys_potential_(x, grad, &energy, &do_grad);
   return energy;;
}

// cleanup and close
void wales_cleanup()
{
}

