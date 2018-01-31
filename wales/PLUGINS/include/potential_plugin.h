#ifndef OPTIM_PLUGIN_H
#define OPTIM_PLUGIN_H

#ifdef __cplusplus
extern "C" {
#endif

// initialize the plugin
void wales_initialize_plugin();

// get the number of atoms
int wales_get_dof();

// calculate potential
double wales_potential(double *x, double *grad, double *hess, int do_grad, int do_hess);

// cleanup and close
void wales_cleanup();


#ifdef __cplusplus
}
#endif

#endif
