#ifndef RSCFT_H
#define RSCFT_H

#include <fftw3.h>
#include <iostream>
using namespace std;

class RScft
{
private:

    int nx, ny, nn, ns, mid;  // ns, nx, and ny are odd
                              // (purely for speed in the cosine transform)
              // mid is f*(ns-1).  mid must be less than ns.
              // mid must be even.  If the copolymer is all A, mid is ns-1.
              // If the copolmer is all B, mid is 0.
    double len, Q, *wA, *wB, *phiA, *phiB, *ksi, **q, **qdag, chiN;
    double ds, dx, freeE;

    // To be used with fftw
    fftw_plan plan;
    double *intermed, *factor;

    // Now the grand canonical variables.
    int alphaA, alphaB;  // These relate to ns in the same way as mid.  So
                         // if the copolymers and homopolymers have the
                         // same length, then alphaA = alphaB = ns-1.
                         // alphaA and alphaB must be even, just like mid.
    double QhA, QhB, *phiAH, *phiBH, zA, zB, **qAH, **qBH;

    bool gc, surf;  // A flags to indicate which ensemble is being used.

    // Additional variables.
    double *dens, *pot, eps, delta, dy;

    void set_up();

    void dealloc_mem();

    // Dynamically allocates memory for an ns by nx array.
    double** create2();

    // Releases dynamically allocated memory for an ns by nx array.
    void destroy2(double** a);

    double min_field();
    void adjust_field();

    void calc_phis_can();
    void calc_phis_gc();

    double free_energy_md();
    double ave_val(double val[]);

public:

    // For debugging purposes.
    void see_trans(ostream &stream);
    void write_q(ostream &stream);
    void write_q2(ostream &stream, int s);
    void write_qdag(ostream &stream);
    void write_ksi(ostream &stream);
    void write_field(ostream &stream);
    void write_dens(ostream &stream);

    void write_plot(ostream &stream);
    void read_data(istream &in);

    void set_field_rand(int seed);

    void set_field(double fA[], double fB[]);

    void set_nsx(int neu_ns, int neu_nx);
    void set_ns(int neu_ns);
    void set_nx(int neu_nx);
    void set_len(double nl);

    void calc_qs();
    //void calc_qs_fd();
    void calc_Q();
    void calc_phis();
    double free_energy();
    void calc_field(ostream &stream);
    void calc_dens();

    double ave_wA();
    double ave_wB();
    double ave_phiA();
    double ave_phiB();
    double ave_ksi();

    double max_phiA();
    double max_phiB();

    // Now the grand canonical specific functions.
    void set_alphaA(int naA);
    void set_alphaB(int naB);
    void set_zA(double nzA);
    void set_zB(double nzB);

    double ave_phiAH();
    double ave_phiBH();
    double ave_phiAC();
    double ave_phiBC();

    double hpA_ex(double bulk);
    double get_vol();
    double get_eff_vol();
    void profile_phiA_x(ostream &stream);
    void profile_phiA_y(ostream &stream);
    void profile_hpA_x(ostream &stream);
    void profile_hpA_y(ostream &stream);
    void profile_hp_x(ostream &stream);
    void profile_hp_y(ostream &stream);
    void fields2fs(ostream &stream);

    void set_pot(double gg);
    void set_pot3(double gg);

    // Canonical constructor.
    RScft(double segregation, double period_length, int mm, int xx, int ss);

    // Canonical constructor with surfaces.
    RScft(double segregation, double period_length, int mm, int xx, int ss,
                                                  bool dummy, double silon);

    // Grand canonical constructor.
    RScft(double segregation, double period_length, int mm, int xx, int ss,
                                                                  double z);

    // Grand canonical constructor with surfaces.
    RScft(double segregation, double period_length, int mm, int xx, int ss,
                                                    double z, double silon);

    // Two dimensions.
    RScft(double segregation, double period_length, double thickness,
                int mm, int xx, int yy, int ss, double z, double silon);

    ~RScft();
};

#endif  // RSCFT_H
