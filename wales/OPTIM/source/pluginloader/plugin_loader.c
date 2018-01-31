#include <stdio.h>
#include <stdlib.h>
#include <dlfcn.h>

static void *lib_handle = NULL;
static void (*lib_initialize)() = NULL;;
static int (*lib_get_dof)() = NULL;;
static double (*lib_potential)(double *x, double *grad, double *hess, int do_grad, int do_hess) = NULL;
static void (*lib_cleanup)()=NULL;

#define LOAD_FUNCTION(fn, name) \
    fn = dlsym(lib_handle, name); \
    if ((error = dlerror()) != NULL) { \
        fprintf(stderr, "%s\n", error); \
        exit(1); \
    }

void load_plugin_() 
{
   double (*fn)(int *);
   int x;
   char *error;
   int bytes_read;
   int nbytes=1023;
   char plugin_file[1024] = "./plugin.so";
   
   //bytes_read = getline (plugin_file, &nbytes, stdin);
   
   printf("loading plugin %s\n", plugin_file);
   
   lib_handle = dlopen(plugin_file, RTLD_LAZY | RTLD_LOCAL);

   if (!lib_handle) {
      fprintf(stderr, "%s\n", dlerror());
      exit(1);
   }
   
   LOAD_FUNCTION(lib_initialize, "wales_initialize");
   LOAD_FUNCTION(lib_get_dof, "wales_get_dof");
   LOAD_FUNCTION(lib_potential, "wales_potential");
   LOAD_FUNCTION(lib_cleanup, "wales_cleanup");

   (*lib_initialize)();
}

void unload_plugin_()
{
    (*lib_cleanup)();
    dlclose(lib_handle);
}

void plugin_get_natoms_(int *natoms) 
{
    natoms = (*lib_get_dof)()/3;
}

void plugin_potential_(double *x, double *energy, double *grad, double *hess,
    int *do_grad, int *do_hess) 
{
    *energy = (*lib_potential)(x, grad, hess, *do_grad, *do_hess);
}

