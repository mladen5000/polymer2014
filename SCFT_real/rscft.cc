#include "rscft.h"
#include <cmath>  // not needed on cae or che computers
#include <fftw3.h>
#include <iostream>
#include <cstdlib>
using namespace std;

// allocates memory and initializes variables
// nx, ny, nn, ns, len, and delta MUST already be set.  In the gc, also
// alphaA and alphaB.  With surfaces, eps and gg.  And of course, the gc
// and surf flags must be set.
void RScft::set_up()
{
    // Save time by doing certain calculations only once.
    ds = 1.0/(ns-1);
    if (ny > 1)
    {
        dx = len/(nx-1);
        dy = delta/(ny+1);
    }
    else if (surf)
        dx = len/(nx+1);
    else
        dx = len/(nx-1);

    // Memory allocations
    wA = new double[nn];
    wB = new double[nn];
    phiA = new double[nn];
    phiB = new double[nn];
    ksi = new double[nn];

    q = create2();
    qdag = create2();

    // Initial conditions of q and qdag.  They are the same for every
    // call to calc_qs.
    for (int i=0; i<nn; i++)
    {
        q[0][i] = 1.0;  // sin(dx*i);
        qdag[ns-1][i] = 1.0;  // sin(dx*i);
    }
    //for (int i=0; i<100; i++)
    //    q[0][i] = 0.0;
    //for (int i=201; i<301; i++)
    //    q[0][i] = 0.0;
    //if (surf)
    //    q[0][0] = q[0][1] = q[0][2] = q[0][3] =
    //    q[0][nx-4] = q[0][nx-3] = q[0][nx-2] = q[0][nx-1] = 0.0;

    // Get the forward and reverse FFT routines set up.
    intermed = (double*) fftw_malloc(sizeof(double) * nn);
    factor = new double[nn];
    //cout << "Initializing fftw plans.\n";
    if (ny > 1)
    {
        plan = fftw_plan_r2r_2d(nx, ny, intermed, intermed,
                  FFTW_REDFT00, FFTW_RODFT00, FFTW_EXHAUSTIVE);

        // To be used between forward and reverse FFT's.
        double tempc = -ds/6.0*M_PI*M_PI/len/len;  // These factors are used
        double temps = -ds/6.0*M_PI*M_PI/delta/delta;  // over and over.
        double temp, div = 4.0*(nx-1)*(ny+1);
        for (int i=0; i<nn; /* update i elsewhere */)
        {
            temp = tempc*(i/ny)*(i/ny);
            for (int j=1; j<=ny; j++, i++)
                factor[i] = exp(temp+temps*j*j)/div;
        }
    }
    else if (surf)
    {
        plan = fftw_plan_r2r_1d(nx, intermed, intermed, FFTW_RODFT00,
                                                FFTW_EXHAUSTIVE);

        // To be used between forward and reverse FFT's.
        double temp = -ds/6.0*M_PI*M_PI/len/len;  // This factor is used
                                                  // over and over.
        for (int i=1; i<=nx; i++)
            factor[i-1] = exp(temp*i*i)/(2.0*(nx+1));
    }
    else
    {
        plan = fftw_plan_r2r_1d(nx, intermed, intermed, FFTW_REDFT00,
                                                FFTW_EXHAUSTIVE);

        // To be used between forward and reverse FFT's.
        double temp = -ds/6.0*M_PI*M_PI/len/len;  // This factor is used
                                                  // over and over.
        for (int i=0; i<nx; i++)
            factor[i] = exp(temp*i*i)/(2.0*(nx-1));
    }
    //cout << "Done!\n";

    // More memory is needed for the grand canonical ensemble.
    if (gc)
    {
        phiAH = new double[nn];
        phiBH = new double[nn];

        qAH = new double*[alphaA+1];
        // The loop only runs if (alphaA > mid)
        for (int i=mid+1; i<=alphaA; i++)
            qAH[i] = new double[nn];

        qBH = new double*[alphaB+1];
        // The loop only runs if (alphaB > ns-1-mid)
        for (int i=ns-mid; i<=alphaB; i++)
            qBH[i] = new double[nn];
    }

    dens = new double[nn];
    pot = new double[nn];
    if (ny > 1)
    {
        for (int i=0; i<ny; i++)
        {
            if (dy*(i+1) < eps)
                dens[i] = (1.0-cos(M_PI*dy*(i+1)/eps))/2.0;
            else if (dy*(i+1) > delta-eps)
                dens[i] = (1.0-cos(M_PI*(delta-dy*(i+1))/eps))/2.0;
            else
                dens[i] = 1.0;
            pot[i] = 0.0;
            for (int j=1; j<nx; j++)
            {
                dens[ny*j+i] = dens[i];
                pot[ny*j+i] = pot[i];
            }
        }
    }
    else if (surf)
        for (int i=0; i<nn; i++)
        {
            if (dx*(i+1) < eps)
                dens[i] = (1.0-cos(M_PI*dx*(i+1)/eps))/2.0;
            else if (dx*(i+1) > len-eps)
                dens[i] = (1.0-cos(M_PI*(len-dx*(i+1))/eps))/2.0;
            else
                dens[i] = 1.0;
            pot[i] = 0.0;
        }
    else
        for (int i=0; i<nn; i++)
        {
            dens[i] = 1.0;
            pot[i] = 0.0;
        }
}

void RScft::dealloc_mem()
{
    if (gc)
    {
        delete[] phiAH;
        delete[] phiBH;

        // The loop only runs if (alphaA > mid)
        for (int i=mid+1; i<=alphaA; i++)
            delete[] qAH[i];
        delete[] qAH;

        // The loop only runs if (alphaB > ns-1-mid)
        for (int i=ns-mid; i<=alphaB; i++)
            delete[] qBH[i];
        delete[] qBH;
    }

    delete[] dens;
    delete[] pot;

    delete[] wA;
    delete[] wB;
    delete[] phiA;
    delete[] phiB;
    delete[] ksi;

    destroy2(q);
    destroy2(qdag);

    fftw_destroy_plan(plan);
    fftw_free(intermed);

    delete[] factor;
}

// Dynamically allocates memory for an ns by nn array.
double** RScft::create2()
{
    double** a;
    a = new double*[ns];
    for (int i=0; i<ns; i++)
        a[i] = new double[nn];
    return a;
}

// Releases dynamically allocated memory for an ns by nn array.
void RScft::destroy2(double** a)
{
    for (int i=0; i<ns; i++)
        delete[] a[i];
    delete[] a;
}

double RScft::min_field()
{
    double min = wA[0];
    if (wB[0] < min)
        min = wB[0];
    for (int i=1; i<nn; i++)
    {
        if (wA[i] < min)
            min = wB[i];
        if (wB[i] < min)
            min = wB[i];
    }
    return min;
}

void RScft::adjust_field()
{
    double min = min_field();
    for (int i=0; i<nn; i++)
    {
        wA[i] -= min;
        wB[i] -= min;
    }
}

void RScft::see_trans(ostream &stream)
{
    double exp_wA[nn];
    for (int i=0; i<nn; i++)
        exp_wA[i] = exp(-wA[i]*ds/2.0);

    for (int i=0; i<nn; i++)
        stream << q[ns-1][i] << ' ';
    stream << endl << endl;

    for (int i=0; i<nn; i++)
        intermed[i] = exp_wA[i]*q[ns-1][i];

    for (int i=0; i<nn; i++)
        stream << intermed[i] << ' ';
    stream << endl << endl;

    // FFT on intermed
    fftw_execute(plan);

    for (int i=0; i<nn; i++)
        stream << intermed[i] << ' ';
    stream << endl << endl;

    for (int i=0; i<nn; i++)
        intermed[i] *= factor[i];

    for (int i=0; i<nn; i++)
        stream << intermed[i] << ' ';
    stream << endl << endl;

    // Inverse FFT on intermed
    fftw_execute(plan);

    for (int i=0; i<nn; i++)
        stream << intermed[i] << ' ';
    stream << endl << endl;

    for (int i=0; i<nn; i++)
        exp_wA[i]*intermed[i];

    for (int i=0; i<nn; i++)
        stream << exp_wA[i]*intermed[i] << ' ';
    stream << endl << endl;
}

void RScft::write_q(ostream &stream)
{
    for (int s=0; s<ns; s++)
    {
        for (int i=0; i<nn; i++)
            stream << q[s][i] << " ";
        stream << endl;
    }
    stream << endl << endl;
}

void RScft::write_q2(ostream &stream, int s)
{
    if (ny > 1)
        for (int i=0; i<nn; i++)
            stream << dx*(i/ny) << ' ' << dy*(i%ny+1) << ' ' << q[s][i] <<endl;
    else if (surf)
        for (int i=0; i<nn; i++)
            stream << dx*(i+1) << ' ' << q[s][i] <<endl;
    else
        for (int i=0; i<nn; i++)
            stream << dx*i << ' ' << q[s][i] <<endl;
}

void RScft::write_qdag(ostream &stream)
{
    for (int s=0; s<ns; s++)
    {
        for (int i=0; i<nn; i++)
            stream << qdag[s][i] << " ";
        stream << endl;
    }
    stream << endl << endl;
}

void RScft::write_ksi(ostream &stream)
{
    for (int i=0; i<nn; i++)
        stream << ksi[i] << " ";
    stream << endl;
}

void RScft::write_field(ostream &stream)
{
    for (int i=0; i<nn; i++)
        stream << wA[i] << " ";
    stream << endl;
    for (int i=0; i<nn; i++)
        stream << wB[i] << " ";
    stream << endl;
}

void RScft::write_dens(ostream &stream)
{
    for (int i=0; i<nn; i++)
        stream << phiA[i] << " ";
    stream << endl;
    for (int i=0; i<nn; i++)
        stream << phiB[i] << " ";
    stream << endl;
}

void RScft::write_plot(ostream &stream)
{
    if (ny > 1)
        for (int i=0; i<nn; i++)
            stream << dx*(i/ny) << ' ' << dy*(i%ny+1) << ' ' << wA[i] << ' '
                   << wB[i] << ' ' << phiA[i] << ' ' << phiB[i] << ' '
                   << ksi[i] << ' ' << phiAH[i] << ' ' << phiBH[i] << ' '
                   << phiA[i]-phiAH[i] << ' ' << phiB[i]-phiBH[i] << ' '
                   << dens[i] << ' ' << phiA[i]+phiB[i] << ' '<< pot[i]
                   << endl;
    else if (gc)
        if (surf)
            for (int i=0; i<nn; i++)
                stream << dx*(i+1) << ' ' << wA[i] << ' ' << wB[i] << ' '
                       << phiA[i] << ' ' << phiB[i] << ' ' << ksi[i] << ' '
                       << phiAH[i] << ' ' << phiBH[i] << ' '
                       << phiA[i]-phiAH[i] << ' ' << phiB[i]-phiBH[i] << ' '
                       << dens[i] << ' ' << phiA[i]+phiB[i] << ' '<< pot[i]
                       << endl;
        else
            for (int i=0; i<nn; i++)
                stream << dx*i << ' ' << wA[i] << ' ' << wB[i] << ' '
                       << phiA[i] << ' ' << phiB[i] << ' ' << ksi[i] << ' '
                       << phiAH[i] << ' ' << phiBH[i] << ' '
                       << phiA[i]-phiAH[i] << ' ' << phiB[i]-phiBH[i] << ' '
                       << dens[i] << ' ' << phiA[i]+phiB[i] << ' '<< pot[i]
                       << endl;
    else
        if (surf)
            for (int i=0; i<nx; i++)
                stream << dx*(i+1) << ' ' << wA[i] << ' ' << wB[i] << ' '
                       << phiA[i] << ' ' << phiB[i] << ' ' << ksi[i] << ' '
                       << dens[i] << ' ' << phiA[i]+phiB[i] << ' '<< pot[i]
                       << endl;
        else
            for (int i=0; i<nx; i++)
                stream << dx*i << ' ' << wA[i] << ' ' << wB[i] << ' '
                       << phiA[i] << ' ' << phiB[i] << ' ' << ksi[i] << ' '
                       << dens[i] << ' ' << phiA[i]+phiB[i] << ' '<< pot[i]
                       << endl;
}

// Inputs what was previously output.  Assumes the number of data points match.
// Most importantly, the fields are set.  Also set: phiA, phiB, ksi, dens, pot,
// phiAH (gc), and phiBH (gc).
// May want to remove setting dens and pot.
void RScft::read_data(istream &in)
{
    // Move past the comments.
    char ch;
    in.get(ch);
    while (ch == '#')
    {
        while (ch != '\n')
            in.get(ch);
        in.get(ch);
    }
    in.putback(ch);

    // Now extract the data.
    double dd;
    if (ny > 1)
        for (int i=0; i<nn; i++)
            in >> dd >> dd >> wA[i] >> wB[i] >> phiA[i] >> phiB[i]
               >> ksi[i] >> phiAH[i] >> phiBH[i] >> dd >> dd
               >> dens[i] >> dd >> pot[i];
    else if (gc)
        for (int i=0; i<nn; i++)
            in >> dd >> wA[i] >> wB[i] >> phiA[i] >> phiB[i]
               >> ksi[i] >> phiAH[i] >> phiBH[i] >> dd >> dd
               >> dens[i] >> dd >> pot[i];
    else
        for (int i=0; i<nn; i++)
            in >> dd >> wA[i] >> wB[i] >> phiA[i] >> phiB[i]
               >> ksi[i] >> dens[i] >> dd >> pot[i];
}

// The rand numbers are between 0 and 1.
void RScft::set_field_rand(int seed)
{
    srand(seed);
    for (int i=0; i<nn; i++)
    {
        wA[i] = double(rand())/double(RAND_MAX);
        wB[i] = double(rand())/double(RAND_MAX);
    }
}

// An inefficient function, but it works.
// 1d only.
void RScft::set_nsx(int neu_ns, int neu_nx)
{
    dealloc_mem();
    ns = neu_ns;
    nx = neu_nx;
    set_up();
}

// An inefficient function, but it works.
// 1d only.
void RScft::set_ns(int neu_ns)
{
    destroy2(q);
    destroy2(qdag);
    if (gc)
        // The loop only runs if (alphaB > ns-1-mid)
        for (int i=ns-mid; i<=alphaB; i++)
            delete[] qBH[i];

    ns = neu_ns;

    ds = 1.0/(ns-1);
    if (surf)
    {
        // To be used between forward and reverse FFT's.
        double temp = -ds/6.0*M_PI*M_PI/len/len;  // This factor is used
                                                  // over and over.
        for (int i=1; i<=nx; i++)
            factor[i-1] = exp(temp*i*i)/(2.0*(nx+1));
    }
    else
    {
        double temp = -ds*2.0/3.0*M_PI*M_PI/len/len;  // This factor is used
                                                      // over and over.
        for (int i=0; i<=nx/2; i++)
            factor[i] = exp(temp*i*i)/nx;
        for (int i=nx/2+1; i<nx; i++)
            factor[i] = exp(temp*(nx-i)*(nx-i))/nx;
    }

    q = create2();
    qdag = create2();

    // Initial conditions of q and qdag.  They are the same for every
    // call to calc_qs.
    for (int i=0; i<nx; i++)
    {
        q[0][i] = 1.0;  // sin(dx*i);
        qdag[ns-1][i] = 1.0;  // sin(dx*i);
    }

    if (gc)
        // The loop only runs if (alphaB > ns-1-mid)
        for (int i=ns-mid; i<=alphaB; i++)
            qBH[i] = new double[nx];
}

// An inefficient function, but it works.
// 1d only.
void RScft::set_nx(int neu_nx)
{
    dealloc_mem();
    nx = neu_nx;
    set_up();  // The only thing done in set_up() that doesn't need to be
               // done is setting ds.
}

// 1d only.
void RScft::set_len(double nl)
{
    len = nl;

    dx = len/nx;

    if (surf)
    {
        // To be used between forward and reverse FFT's.
        double temp = -ds/6.0*M_PI*M_PI/len/len;  // This factor is used
                                                  // over and over.
        for (int i=1; i<=nx; i++)
            factor[i-1] = exp(temp*i*i)/(2.0*(nx+1));
    }
    else
    {
        double temp = -ds*2.0/3.0*M_PI*M_PI/len/len;  // This factor is used
                                                      // over and over.
        for (int i=0; i<=nx/2; i++)
            factor[i] = exp(temp*i*i)/nx;
        for (int i=nx/2+1; i<nx; i++)
            factor[i] = exp(temp*(nx-i)*(nx-i))/nx;
    }
}

void RScft::set_field(double fA[], double fB[])
{
    for (int i=0; i<nn; i++)
    {
        wA[i] = fA[i];
        wB[i] = fB[i];
    }
}

// Pseudo-spectral method.  Uses R_e as the frame of reference.
void RScft::calc_qs()
{
    double exp_wA[nn];
    double exp_wB[nn];
    for (int i=0; i<nn; i++)
    {
        exp_wA[i] = exp(-wA[i]*ds/2.0);
        exp_wB[i] = exp(-wB[i]*ds/2.0);
    }

    // Figure out q.
    for (int s=0; s<mid; s++)
    {
        for (int i=0; i<nn; i++)
            intermed[i] = exp_wA[i]*q[s][i];

        // FFT on intermed
        fftw_execute(plan);

        for (int i=0; i<nn; i++)
            intermed[i] *= factor[i];

        // Inverse FFT on intermed
        fftw_execute(plan);

        for (int i=0; i<nn; i++)
            q[s+1][i] = exp_wA[i]*intermed[i];

        //for (int i=0; i<100; i++)
        //{
        //    q[s+1][199-i] -= q[s+1][i];
        //    q[s+1][i] = 0.0;
        //}
        //for (int i=201; i<301; i++)
        //{
        //    q[s+1][401-i] -= q[s+1][i];
        //    q[s+1][i] = 0.0;
        //}
        //if (surf)
        //    q[s+1][0] = q[s+1][nx-1] = 0.0;
        //if (surf)
        //    q[s+1][0] = q[s+1][1] = q[s+1][2] = q[s+1][3] =
        //    q[s+1][nx-4] = q[s+1][nx-3] = q[s+1][nx-2] = q[s+1][nx-1] = 0.0;
    }

    for (int s=mid; s<ns-1; s++)
    {
        for (int i=0; i<nn; i++)
            intermed[i] = exp_wB[i]*q[s][i];

        // FFT on intermed
        fftw_execute(plan);

        for (int i=0; i<nn; i++)
            intermed[i] *= factor[i];

        // Inverse FFT on intermed
        fftw_execute(plan);

        for (int i=0; i<nn; i++)
            q[s+1][i] = exp_wB[i]*intermed[i];

        //for (int i=0; i<100; i++)
        //{
        //    q[s+1][199-i] -= q[s+1][i];
        //    q[s+1][i] = 0.0;
        //}
        //for (int i=201; i<301; i++)
        //{
        //    q[s+1][401-i] -= q[s+1][i];
        //    q[s+1][i] = 0.0;
        //}
        //if (surf)
        //    q[s+1][0] = q[s+1][nx-1] = 0.0;
        //if (surf)
        //    q[s+1][0] = q[s+1][1] = q[s+1][2] = q[s+1][3] =
        //    q[s+1][nx-4] = q[s+1][nx-3] = q[s+1][nx-2] = q[s+1][nx-1] = 0.0;
    }

/*
    double err = 0.0;
    for (int i=0; i<nx; i++)
        err += fabs(q[ns-1][i]-exp(-1.0)*sin(dx*i));
    err /= nx;
    cout << 1.0/ds << '\t' << 1.0/dx << '\t' << err << endl;
*/

    // Now figure out qdag.
    for (int s=ns-1; s>mid; s--)
    {
        for (int i=0; i<nn; i++)
            intermed[i] = exp_wB[i]*qdag[s][i];

        // FFT on intermed
        fftw_execute(plan);

        for (int i=0; i<nn; i++)
            intermed[i] *= factor[i];

        // Inverse FFT on intermed
        fftw_execute(plan);

        for (int i=0; i<nn; i++)
            qdag[s-1][i] = exp_wB[i]*intermed[i];

        //if (surf)
        //    qdag[s-1][0] = qdag[s-1][nx-1] = 0.0;
    }
    for (int s=mid; s>0; s--)
    {
        for (int i=0; i<nn; i++)
            intermed[i] = exp_wA[i]*qdag[s][i];

        // FFT on intermed
        fftw_execute(plan);

        for (int i=0; i<nn; i++)
            intermed[i] *= factor[i];

        // Inverse FFT on intermed
        fftw_execute(plan);

        for (int i=0; i<nn; i++)
            qdag[s-1][i] = exp_wA[i]*intermed[i];

        //if (surf)
        //    qdag[s-1][0] = qdag[s-1][nx-1] = 0.0;
    }

    // Done if in the canonical ensemble.
    if (!gc)
        return;

    // Don't create new entries for duplicate values.  Just point the
    // homopolymer q's to the copolymer values.  Let's start with qAH
    // Actually, this only needs to be done once in set_up, set_ns,
    // set_alphaA, and set_alphaB.
    if (alphaA < mid)
        for (int s=0; s<=alphaA; s++)
            qAH[s] = q[s];
    else
    {
        for (int s=0; s<=mid; s++)
            qAH[s] = q[s];
        for (int s=mid; s<alphaA; s++)
        {
            for (int i=0; i<nn; i++)
                intermed[i] = exp_wA[i]*qAH[s][i];

            // FFT on intermed
            fftw_execute(plan);

            for (int i=0; i<nn; i++)
                intermed[i] *= factor[i];

            // Inverse FFT on intermed
            fftw_execute(plan);

            for (int i=0; i<nn; i++)
                qAH[s+1][i] = exp_wA[i]*intermed[i];

            //if (surf)
            //    qAH[s+1][0] = qAH[s+1][nx-1] = 0.0;
        }
    }

    // Now qBH.  This time we copy from qdag.
    if (alphaB < ns-1-mid)
        for (int s=ns-1; s>=ns-1-alphaB; s--)
            qBH[ns-1-s] = qdag[s];
    else
    {
        for (int s=ns-1; s>=ns-1-mid; s--)
            qBH[ns-1-s] = qdag[s];

        for (int s=ns-1-mid; s<alphaB; s++)
        {
            for (int i=0; i<nn; i++)
                intermed[i] = exp_wB[i]*qBH[s][i];

            // FFT on intermed
            fftw_execute(plan);

            for (int i=0; i<nn; i++)
                intermed[i] *= factor[i];

            // Inverse FFT on intermed
            fftw_execute(plan);

            for (int i=0; i<nn; i++)
                qBH[s+1][i] = exp_wB[i]*intermed[i];

            //if (surf)
            //    qBH[s+1][0] = qBH[s+1][nx-1] = 0.0;
        }
    }
}


// Finite difference method.  Uses R_e as the frame of reference.
/* Replaced by the pseudo-spectral method.  Never extended into gc.
 * Numrec no longer necessary.
void RScft::calc_qs_fd()
{
    double coef1 = 12.0*dx*dx-2.0*ds;     // R_g => 2.0*dx*dx-2.0*ds
    double coef2 = 6.0*ds*dx*dx;          // R_g => ds*dx*dx
    double coef3 = 12.0*dx*dx+2.0*ds;     // R_g => 2.0*dx*dx+2.0*ds

    int i;
    double a[nx], b_wA[nx], b_wB[nx], c[nx], r[nx], val_wA[nx], val_wB[nx];

    // First q is calculated.
    for (i=0; i<nx; i++)
    {
        a[i] = c[i] = -ds;
        b_wA[i] = coef3+coef2*wA[i];
        b_wB[i] = coef3+coef2*wB[i];
        val_wA[i] = coef1-coef2*wA[i];
        val_wB[i] = coef1-coef2*wB[i];
    }

    for (int s=0; s<mid; s++)
    {
        r[0] = ds*(q[s][1]+q[s][nx-1])+val_wA[0]*q[s][0];
        for (i=1; i<nx-1; i++)
            r[i] = ds*(q[s][i+1]+q[s][i-1])+val_wA[i]*q[s][i];
        // i == nx-1
        r[i] = ds*(q[s][0]+q[s][i-1])+val_wA[i]*q[s][i];

        Numrec::cyclic(a, b_wA, c, r, q[s+1], nx);
    }

    // Switch from wA to wB in the middle.
    for (int s=mid; s<ns-1; s++)
    {
        r[0] = ds*(q[s][1]+q[s][nx-1])+val_wB[0]*q[s][0];
        for (i=1; i<nx-1; i++)
            r[i] = ds*(q[s][i+1]+q[s][i-1])+val_wB[i]*q[s][i];
        // i == nx-1
        r[i] = ds*(q[s][0]+q[s][i-1])+val_wB[i]*q[s][i];

        Numrec::cyclic(a, b_wB, c, r, q[s+1], nx);
    }

    // For finding the error relation between s and x.
    //double err = 0.0;
    //for (i=0; i<nx; i++)
    //    err += fabs(q[ns-1][i]-exp(-1.0)*sin(dx*i));
    //err /= nx;
    //cout << 1.0/ds << '\t' << 1.0/dx << '\t' << err << endl;

    // Now qdag is calculated.
    for (int s=ns-1; s>mid; s--)
    {
        r[0] = ds*(qdag[s][1]+qdag[s][nx-1])+val_wB[0]*qdag[s][0];
        for (i=1; i<nx-1; i++)
            r[i] = ds*(qdag[s][i+1]+qdag[s][i-1])+val_wB[i]*qdag[s][i];
        // i == nx-1
        r[i] = ds*(qdag[s][0]+qdag[s][i-1])+val_wB[i]*qdag[s][i];

        Numrec::cyclic(a, b_wB, c, r, qdag[s-1], nx);
    }

    // Switch from wB to wA in the middle.
    for (int s=mid; s>0; s--)
    {
        r[0] = ds*(qdag[s][1]+qdag[s][nx-1])+val_wA[0]*qdag[s][0];
        for (i=1; i<nx-1; i++)
            r[i] = ds*(qdag[s][i+1]+qdag[s][i-1])+val_wA[i]*qdag[s][i];
        // i == nx-1
        r[i] = ds*(qdag[s][0]+qdag[s][i-1])+val_wA[i]*qdag[s][i];

        Numrec::cyclic(a, b_wA, c, r, qdag[s-1], nx);
    }
}
*/

void RScft::calc_Q()
{
    double vol = get_vol();  // either len or len*delta

    Q = vol*ave_val(q[ns-1]);

    // Done if in the canonical ensemble.
    if (!gc)
        return;

    QhA = vol*ave_val(qAH[alphaA]);
    QhB = vol*ave_val(qBH[alphaB]);
}

// 1d only (len is used for volume)
void RScft::calc_phis_can()
{
    double sum1, sum2;
    for (int i=0; i<nx; i++)
    {
        sum1 = 0.0;
        sum2 = 0.0;
        for (int s=1; s<mid-1; /* s updated in loop */)
        {
            sum1 += q[s][i]*qdag[s][i];
            s++;
            sum2 += q[s][i]*qdag[s][i];
            s++;
        }
        sum1 += q[mid-1][i]*qdag[mid-1][i];
        phiA[i] = ((q[0][i]*qdag[0][i]+q[mid][i]*qdag[mid][i])
                           + (2.0*sum1+sum2)*2.0)/3.0*ds*len/Q;

        sum1 = 0.0;
        sum2 = 0.0;
        for (int s=mid+1; s<ns-2; /* s updated in loop */)
        {
            sum1 += q[s][i]*qdag[s][i];
            s++;
            sum2 += q[s][i]*qdag[s][i];
            s++;
        }
        sum1 += q[ns-1][i]*qdag[ns-1][i];
        phiB[i] = ((q[mid][i]*qdag[mid][i]+q[ns-1][i]*qdag[ns-1][i])
                             + (2.0*sum1+sum2)*2.0)/3.0*ds*len/Q;
    }
}

// 1d and 2d.
void RScft::calc_phis_gc()
{
    double sum1, sum2;
    int s;
    for (int i=0; i<nn; i++)
    {
        // The copolymer part of A.  This is done just as in the canonical
        // case, except the final integral isn't multiplied by V/Q.
        sum1 = 0.0;
        sum2 = 0.0;
        for (s=1; s<mid-1; /* s updated in loop */)
        {
            sum1 += q[s][i]*qdag[s][i];
            s++;
            sum2 += q[s][i]*qdag[s][i];
            s++;
        }
        sum1 += q[mid-1][i]*qdag[mid-1][i];
        phiA[i] = ((q[0][i]*qdag[0][i]+q[mid][i]*qdag[mid][i])
                           + (2.0*sum1+sum2)*2.0)/3.0*ds;

        // The homopolymer part of A.
        sum1 = 0.0;
        sum2 = 0.0;
        for (s=1; s<alphaA/2-1; /* s updated in loop */)
        {
            sum1 += qAH[s][i]*qAH[alphaA-s][i];
            s++;
            sum2 += qAH[s][i]*qAH[alphaA-s][i];
            s++;
        }
        if ((alphaA/2)%2)
        {
            // Currently s==alphaA/2
            sum2 += qAH[s][i]*qAH[s][i];
            phiAH[i] = (2.0*qAH[0][i]*qAH[alphaA][i]
                              + (2.0*sum1+sum2)*4.0)/3.0*ds*zA;
        }
        else
        {
            // Currently s==alphaA/2-1
            sum1 += qAH[s][i]*qAH[s+2][i];
            s++;  // Now s==alphaA/2
            phiAH[i] = (2.0*(qAH[0][i]*qAH[alphaA][i]+qAH[s][i]*qAH[s][i])
                              + (2.0*sum1+sum2)*4.0)/3.0*ds*zA;
        }
        phiA[i] += phiAH[i];

        // The copolymer part of B.  This is done just as in the canonical
        // case, except the final integral isn't multiplied by V/Q.
        sum1 = 0.0;
        sum2 = 0.0;
        for (s=mid+1; s<ns-2; /* s updated in loop */)
        {
            sum1 += q[s][i]*qdag[s][i];
            s++;
            sum2 += q[s][i]*qdag[s][i];
            s++;
        }
        sum1 += q[ns-1][i]*qdag[ns-1][i];
        phiB[i] = ((q[mid][i]*qdag[mid][i]+q[ns-1][i]*qdag[ns-1][i])
                             + (2.0*sum1+sum2)*2.0)/3.0*ds;

        // The homopolymer part of B.
        sum1 = 0.0;
        sum2 = 0.0;
        for (s=1; s<alphaB/2-1; /* s updated in loop */)
        {
            sum1 += qBH[s][i]*qBH[alphaB-s][i];
            s++;
            sum2 += qBH[s][i]*qBH[alphaB-s][i];
            s++;
        }
        if ((alphaB/2)%2)
        {
            // Currently s==alphaB/2
            sum2 += qBH[s][i]*qBH[s][i];
            phiBH[i] = (2.0*qBH[0][i]*qBH[alphaB][i]
                              + (2.0*sum1+sum2)*4.0)/3.0*ds*zB;
        }
        else
        {
            // Currently s==alphaB/2-1
            sum1 += qBH[s][i]*qBH[s+2][i];
            s++;  // Now s==alphaB/2
            phiBH[i] = (2.0*(qBH[0][i]*qBH[alphaA][i]+qBH[s][i]*qBH[s][i])
                              + (2.0*sum1+sum2)*4.0)/3.0*ds*zB;
        }
        phiB[i] += phiBH[i];
    }
}

// This is where it is critical that mid be even.  The integration is
// tuned around this assumption.  In the gc, alphaA and alphaB must also
// be even.
// 1d only.
void RScft::calc_phis()
{
    // The density calculation is slightly different depending on the
    // ensemble.  Pick the appropriate version.
    if (gc)
        calc_phis_gc();
    else
        calc_phis_can();
}

double RScft::free_energy_md()
{
    double sum1, sum2, sum3, sum4;
    double sumfe = 0.0, sumvol = 0.0;
    // Trapezoid rule in outer loop, Simpson's rule in inner loop.
    for (int i=0; i<nn; /* i updated elsewhere */)
    {
        sum1 = 0.0;
        sum2 = 0.0;
        sum3 = 0.0;
        sum4 = 0.0;
        // Simpson's rule with periodic boundaries.
        for (int j=0; j<ny-1; j+=2 /* i updated in loop */)
        {
            sum1 += chiN*phiA[i]*phiB[i]-pot[i]*(phiA[i]-phiB[i])-wA[i]*phiA[i]
                             -wB[i]*phiB[i]-ksi[i]*(dens[i]-phiA[i]-phiB[i]);
            sum3 += phiA[i]+phiB[i];
            i++;
            sum2 += chiN*phiA[i]*phiB[i]-pot[i]*(phiA[i]-phiB[i])-wA[i]*phiA[i]
                             -wB[i]*phiB[i]-ksi[i]*(dens[i]-phiA[i]-phiB[i]);
            sum4 += phiA[i]+phiB[i];
            i++;
        }
        sum1 += chiN*phiA[i]*phiB[i]-pot[i]*(phiA[i]-phiB[i])-wA[i]*phiA[i]
                             -wB[i]*phiB[i]-ksi[i]*(dens[i]-phiA[i]-phiB[i]);
        sum3 += phiA[i]+phiB[i];
        i++;
        if (i==ny || i==nn)
        {
            sumfe += (2.0*sum1+sum2)/3.0*dy;
            sumvol += (2.0*sum3+sum4)/3.0*dy;
        }
        else
        {
            sumfe += (2.0*sum1+sum2)*2.0/3.0*dy;
            sumvol += (2.0*sum3+sum4)*2.0/3.0*dy;
        }
    }

    double len_eff = sumvol*dx;

    if (gc)
        return freeE = (sumfe*dx -Q-zA*QhA-zB*QhB)/len_eff;
    else
        return freeE = sumfe*dx/len_eff - log(Q/len_eff);
}

// In the canonical case, F/(n*k*T) is returned.  Due to homopolymers with
// varible lengths, the grand canonical free energy is N*F/(k*T*rho0*V).
// (V=len in 1D).  When there are no homopolymers, the expression reduces
// to the canonical expression, except the F's are different.
double RScft::free_energy()
{
    if (ny > 1)
        return free_energy_md();

    if (surf)
    {
        double sum1 = 0.0;
        double sum2 = 0.0;
        double sum3 = 0.0;
        double sum4 = 0.0;
        int i;
        // Simpson's rule with the first and last points implicit and
        // equal to zero.
        for (i=0; i<nx-1; /* i updated in loop */)
        {
            sum1 += chiN*phiA[i]*phiB[i]-pot[i]*(phiA[i]-phiB[i])-wA[i]*phiA[i]
                             -wB[i]*phiB[i]-ksi[i]*(dens[i]-phiA[i]-phiB[i]);
            sum3 += phiA[i]+phiB[i];
            i++;
            sum2 += chiN*phiA[i]*phiB[i]-pot[i]*(phiA[i]-phiB[i])-wA[i]*phiA[i]
                              -wB[i]*phiB[i]-ksi[i]*(dens[i]-phiA[i]-phiB[i]);
            sum4 += phiA[i]+phiB[i];
            i++;
        }
        // i is nx-1
        sum1 += chiN*phiA[i]*phiB[i]-pot[i]*(phiA[i]-phiB[i])-wA[i]*phiA[i]
                             -wB[i]*phiB[i]-ksi[i]*(dens[i]-phiA[i]-phiB[i]);
        sum3 += phiA[i]+phiB[i];

        double len_eff = (2.0*sum3+sum4)*2.0/3.0*dx;

        cout << Q << ' ' << QhA << ' ' << QhB <<  ' '
             << (2.0*sum1+sum2)*2.0/3.0*dx << ' ' << len_eff << endl;

        if (gc)
            return freeE = ((2.0*sum1+sum2)*2.0/3.0*dx -Q-zA*QhA-zB*QhB)
                                                                /len_eff;
        else
            return freeE = (2.0*sum1+sum2)*2.0/3.0*dx/len_eff - log(Q/len_eff);
    }
    else
    {
        double sum1, sum3;
        // Trapezoid rule.
        sum1 = chiN*phiA[0]*phiB[0]-pot[0]*(phiA[0]-phiB[0])-wA[0]*phiA[0]
                             -wB[0]*phiB[0]-ksi[0]*(dens[0]-phiA[0]-phiB[0]);
        sum1 += chiN*phiA[nx-1]*phiB[nx-1]-pot[nx-1]*(phiA[nx-1]-phiB[nx-1])
                             -wA[nx-1]*phiA[nx-1]-wB[nx-1]*phiB[nx-1]
                             -ksi[nx-1]*(dens[nx-1]-phiA[nx-1]-phiB[nx-1]);
        sum1 /= 2.0;
        sum3 = phiA[0]+phiB[0];
        sum3 += phiA[nx-1]+phiB[nx-1];
        sum3 /= 2.0;
        for (int i=1; i<nx-1; i++)
        {
            sum1 += chiN*phiA[i]*phiB[i]-pot[i]*(phiA[i]-phiB[i])-wA[i]*phiA[i]
                             -wB[i]*phiB[i]-ksi[i]*(dens[i]-phiA[i]-phiB[i]);
            sum3 += phiA[i]+phiB[i];
        }

        double len_eff = sum3*dx;

        if (gc)
            return freeE = (sum1*dx -Q-zA*QhA-zB*QhB)/len_eff;
        else
            return freeE = sum1*dx/len_eff - log(Q/len_eff);
    }
}

// An initial value for field must be set before this function is called.
// The initial value may be a random number, or a good guess.
void RScft::calc_field(ostream &stream)
{
    // A trick from the Rasmussen and Kaloskas paper.  It didn't help me.
    //double temp1, temp2;
    //double sum[nx];
    //for (int i=0; i<nx; i++)
    //    sum[i] = dens[i];

    // For the calculation of lambda.
    double diffwA_raw[nn], diffwB_raw[nn];
    double temp, num = 0.0, denom, lambda = 0.1;
    bool stable = false;
    for (int i=0; i<nn; i++)
        diffwA_raw[i] = diffwB_raw[i] = 0.0;

    double resid1, resid2;
    for (int iter=0; iter<10001; iter++)
    {
        for (int i=0; i<nn; i++)
            ksi[i] = (wA[i]+wB[i]-chiN*dens[i])/2.0;
        calc_qs();
        calc_Q();
        calc_phis();

        resid1 = 0.0, resid2 = 0.0;
        for (int i=0; i<nn; i++)
        {
            resid1 += fabs(dens[i]-phiA[i]-phiB[i]);
            resid2 += fabs(wA[i]-wB[i] + 2.0*pot[i] + chiN*(phiA[i]-phiB[i]));
        }
        if (iter%100 == 0)
            stream << iter << ' ' << resid1/nn << ' ' << resid2/nn << ' '
                   << free_energy() << endl;
        if (resid1+resid2 < 1.0e-9)
            stable = true;

        // Calculate denom, then lambda.  num was calculated during the
        // last iteration.
        denom = 0.0;
        for (int i=0; i<nn; i++)
        {
            //temp1 = phiA[i]+phiB[i];
            //temp2 = temp1-sum[i];
            //sum[i] = temp1;

            temp = chiN*phiB[i]-pot[i]+ksi[i]-wA[i] - diffwA_raw[i] /*+temp2*/;
            diffwA_raw[i] += temp;
            temp *= temp;
            denom += temp;

            temp = chiN*phiA[i]+pot[i]+ksi[i]-wB[i] - diffwB_raw[i] /*+temp2*/;
            diffwB_raw[i] += temp;
            temp *= temp;
            denom += temp;
        }
        if (denom==0.0)
            lambda = 0.1;
        else
        {
            lambda = sqrt(num/denom);
            if (lambda < 0.1)
                lambda = 0.1;
            else if (lambda > 0.3) // 0.5
                lambda = 0.3;
        }
        //cout << num << ' ' << denom << ' ' << iter << ' ' << lambda << endl;
        //if (stable)
            lambda = 0.1;

        // While calculating the new wA's and wB's, find num, which will
        // be used to find lambda in the next iteration.
        num = 0.0;
        for (int i=0; i<nn; i++)
        {
            temp = wA[i];
            wA[i] = wA[i] + lambda*(diffwA_raw[i] /*+temp2*/ );
            temp -= wA[i];
            temp *= temp;
            num += temp;

            temp = wB[i];
            wB[i] = wB[i] + lambda*(diffwB_raw[i] /*+temp2*/ );
            temp -= wB[i];
            temp *= temp;
            num += temp;
        }
        if (num == 0.0)
        {
            cout << "Converged to machine precision.\n";
            return;
        }
        if (!gc)
            adjust_field();
        // No field adjustment in the grand canonical case.
    }
}

void RScft::calc_dens()
{
    for (int i=0; i<nn; i++)
        ksi[i] = (wA[i]+wB[i]-chiN*dens[i])/2.0;
    calc_qs();
    calc_Q();
    calc_phis();
}

double RScft::ave_wA()
{
    return ave_val(wA);
}

double RScft::ave_wB()
{
    return ave_val(wB);
}

double RScft::ave_phiA()
{
    return ave_val(phiA);
}

double RScft::ave_phiB()
{
    return ave_val(phiB);
}

double RScft::ave_ksi()
{
    return ave_val(ksi);
}

double RScft::max_phiA()
{
    double max = phiA[0];
    for (int i=1; i<nn; i++)
        if (phiA[i] > max)
            max = phiA[i];
    return max;
}

double RScft::max_phiB()
{
    double max = phiB[0];
    for (int i=1; i<nn; i++)
        if (phiB[i] > max)
            max = phiB[i];
    return max;
}

// An inefficient function, but it works.
void RScft::set_alphaA(int naA)
{
    if (!gc)
    {
        cerr << "No alphas in canonical ensemble!\n";
        return;
    }

    // The loop only runs if (alphaA > mid)
    for (int i=mid+1; i<=alphaA; i++)
        delete[] qAH[i];
    delete[] qAH;

    alphaA = naA;
    qAH = new double*[alphaA+1];
    // The loop only runs if (alphaA > mid)
    for (int i=mid+1; i<=alphaA; i++)
        qAH[i] = new double[nn];
}

// An inefficient function, but it works.
void RScft::set_alphaB(int naB)
{
    if (!gc)
    {
        cerr << "No alphas in canonical ensemble!\n";
        return;
    }

    // The loop only runs if (alphaB > ns-1-mid)
    for (int i=ns-mid; i<=alphaB; i++)
        delete[] qBH[i];
    delete[] qBH;

    alphaB = naB;
    qBH = new double*[alphaB+1];
    // The loop only run if (alphaB > ns-1-mid)
    for (int i=ns-mid; i<=alphaB; i++)
        qBH[i] = new double[nn];
}

void RScft::set_zA(double nzA)
{
    if (!gc)
    {
        cerr << "No z's in canonical ensemble!\n";
        return;
    }

    zA = nzA;
}

void RScft::set_zB(double nzB)
{
    if (!gc)
    {
        cerr << "No z's in canonical ensemble!\n";
        return;
    }

    zB = nzB;
}

double RScft::ave_phiAH()
{
    if (!gc)
    {
        cerr << "No homopolymers in canonical ensemble!\n";
        return 0.0;
    }

    return ave_val(phiAH);
}

double RScft::ave_phiBH()
{
    if (!gc)
    {
        cerr << "No homopolymers in canonical ensemble!\n";
        return 0.0;
    }

    return ave_val(phiBH);
}

double RScft::ave_phiAC()
{
    if (!gc)
    {
        cerr << "No homopolymers in canonical ensemble!\n";
        return ave_phiA();
    }

    return ave_val(phiA) - ave_val(phiAH);
}

double RScft::ave_phiBC()
{
    if (!gc)
    {
        cerr << "No homopolymers in canonical ensemble!\n";
        return ave_phiB();
    }

    return ave_val(phiB) - ave_val(phiBH);
}

// Old version still needs to be deleted.
double RScft::hpA_ex(double bulk)
{
    return get_vol()*ave_val(phiAH) - bulk*get_eff_vol();

// The below is limited to 1d situations with surfaces.
//    double sum1 = 0.0;
//    double sum2 = 0.0;
//    // Simpson's rule with periodic boundaries.
//    for (int i=0; i<nx-1; /* i updated in loop */)
//    {
//        sum1 += phiAH[i]-bulk*dens[i];
//        i++;
//        sum2 += phiAH[i]-bulk*dens[i];
//        i++;
//    }
//    if (surf)
//        sum1 += phiAH[nx-1]-bulk*dens[nx-1];
//    return (2.0*sum1+sum2)*2.0/3.0*dx;
}

double RScft::get_vol()
{
    if (ny > 1)
        return len*delta;
    else
        return len;
}

double RScft::get_eff_vol()
{
    return get_vol()*(ave_val(phiA)+ave_val(phiB));
}

// Computes the appropriate weighted average of the parameter sent.
// To get the integral over all space, the returned average merely has to
// be multipled by the actual volume (not the effective volume).
double RScft::ave_val(double val[])
{
    if (ny > 1)
    {
        double sum1, sum2;
        double sumt = 0.0;
        // Trapezoid rule in outer loop, Simpson's rule in inner loop.
        for (int i=0; i<nn; /* i updated elsewhere */)
        {
            sum1 = 0.0;
            sum2 = 0.0;
            // Simpson's rule with the first and last points implicit and
            // equal to zero.
            for (int j=0; j<ny-1; j+=2 /* i updated in loop */)
            {
                sum1 += val[i++];
                sum2 += val[i++];
            }
            sum1 += val[i++];
            if (i==ny || i==nn)
                sumt += (2.0*sum1+sum2)/3.0/(ny+1);
            else
                sumt += (2.0*sum1+sum2)*2.0/3.0/(ny+1);
        }
        return sumt/(nx-1);
    }
    else if (surf)
    {
        double sum1 = 0.0;
        double sum2 = 0.0;
        // Simpson's rule with the first and last points implicit and
        // equal to zero.
        for (int i=0; i<nx-1; /* i updated in loop */)
        {
            sum1 += val[i++];
            sum2 += val[i++];
        }
        sum1 += val[nx-1];
        return (2.0*sum1+sum2)*2.0/3.0/(nx+1);
    }
    else
    {
        //double sum = 0.0;
        //for (int i=0; i<nx; i++)
        //    sum += val[i];
        //return sum/nx;

        double sum = 0.0;
        for (int i=1; i<nx-1; i++)
            sum += val[i];
        return ((val[0]+val[nx-1])/2.0 + sum)/(nx-1);
    }
}

void RScft::profile_phiA_x(ostream &stream)
{
    if (ny <= 1)
    {
        cout << "Profile functions only valid in 2d.\n";
        return;
    }

    for (int i=0; i<nn; /* i updated elsewhere */)
    {
        stream << dx*(i/ny);
        for (int j=0; j<ny; j++, i++)
            stream << ' ' << phiA[i]/(phiA[i]+phiB[i]);
        stream << endl;
    }
}

void RScft::profile_phiA_y(ostream &stream)
{
    if (ny <= 1)
    {
        cout << "Profile functions only valid in 2d.\n";
        return;
    }

    for (int i=0; i<ny; i++)
    {
        stream << dy*(i+1);
        for (int j=0; j<nx; j++)
        {
            int indx = i+j*ny;
            stream << ' ' << phiA[indx]/(phiA[indx]+phiB[indx]);
        }
        stream << endl;
    }
}

void RScft::profile_hpA_x(ostream &stream)
{
    if (ny <= 1)
    {
        cout << "Profile functions only valid in 2d.\n";
        return;
    }

    for (int i=0; i<nn; /* i updated elsewhere */)
    {
        stream << dx*(i/ny);
        for (int j=0; j<ny; j++, i++)
            stream << ' ' << phiAH[i]/(phiA[i]+phiB[i]);
        stream << endl;
    }
}

void RScft::profile_hpA_y(ostream &stream)
{
    if (ny <= 1)
    {
        cout << "Profile functions only valid in 2d.\n";
        return;
    }

    for (int i=0; i<ny; i++)
    {
        stream << dy*(i+1);
        for (int j=0; j<nx; j++)
        {
            int indx = i+j*ny;
            stream << ' ' << phiAH[indx]/(phiA[indx]+phiB[indx]);
        }
        stream << endl;
    }
}

void RScft::profile_hp_x(ostream &stream)
{
    if (ny <= 1)
    {
        cout << "Profile functions only valid in 2d.\n";
        return;
    }

    for (int i=0; i<nn; /* i updated elsewhere */)
    {
        stream << dx*(i/ny);
        for (int j=0; j<ny; j++, i++)
            stream << ' ' << (phiAH[i]+phiBH[i])/(phiA[i]+phiB[i]);
        stream << endl;
    }
}

void RScft::profile_hp_y(ostream &stream)
{
    if (ny <= 1)
    {
        cout << "Profile functions only valid in 2d.\n";
        return;
    }

    for (int i=0; i<ny; i++)
    {
        stream << dy*(i+1);
        for (int j=0; j<nx; j++)
        {
            int indx = i+j*ny;
            stream << ' ' << (phiAH[indx]+phiBH[indx])/(phiA[indx]+phiB[indx]);
        }
        stream << endl;
    }
}

// Takes the real space data and transforms it into Fourier space, in a
// style that matches the spectral program's input/output.
void RScft::fields2fs(ostream &stream)
{
    double div[nn];

    if (ny > 1)
    {
        stream << chiN << ' ' << 2.0*len << ' ' << delta << ' ' << nx << ' '
               << ny << ' ' << zA << ' ' << (double)alphaA/(ns-1) << endl;
        div[0] = (2.0*(nx-1))*(sqrt(2.0)*(ny+1));
        for (int i=1; i<ny; i++)
            div[i] = div[0];
        div[ny] = 2.0*(nx-1)*(ny+1);
        for (int i=ny+1; i<nn-ny; i++)
            div[i] = div[ny];
        div[nn-ny] = 4.0*(nx-1)*(ny+1);
        for (int i=nn-ny+1; i<nn; i++)
            div[i] = div[nn-ny];
    }
    else if (surf)
    {
        stream << chiN << ' ' << len << ' ' << nx << ' ' << zA << ' '
               << (double)alphaA/(ns-1) << endl;
        div[0] = sqrt(2.0)*(nx+1);
        for (int i=1; i<nn; i++)
            div[i] = div[0];
    }
    else
    {
        stream << chiN << ' ' << 2.0*len << ' ' << nx << ' ' << zA << ' '
               << (double)alphaA/(ns-1) << endl;
        div[0] = 2.0*(nx-1);
        div[1] = sqrt(2.0)*(nx-1);
        for (int i=2; i<nn-1; i++)
            div[i] = div[1];
        div[nn-1] = 2.0*sqrt(2.0)*(nx-1);
    }

    for (int i=0; i<nn; i++)
        intermed[i] = wA[i];
    // FFT on intermed
    fftw_execute(plan);
    for (int i=0; i<nn; i++)
        stream << intermed[i]/div[i] << ' ';
    stream << endl;

    for (int i=0; i<nn; i++)
        intermed[i] = wB[i];
    // FFT on intermed
    fftw_execute(plan);
    for (int i=0; i<nn; i++)
        stream << intermed[i]/div[i] << ' ';
    stream << endl;

    for (int i=0; i<nn; i++)
        intermed[i] = phiA[i];
    // FFT on intermed
    fftw_execute(plan);
    for (int i=0; i<nn; i++)
        stream << intermed[i]/div[i] << ' ';
    stream << endl;

    for (int i=0; i<nn; i++)
        intermed[i] = phiB[i];
    // FFT on intermed
    fftw_execute(plan);
    for (int i=0; i<nn; i++)
        stream << intermed[i]/div[i] << ' ';
    stream << endl;

    if (gc)
    {
        for (int i=0; i<nn; i++)
            intermed[i] = phiAH[i];
        // FFT on intermed
        fftw_execute(plan);
        for (int i=0; i<nn; i++)
            stream << intermed[i]/div[i] << ' ';
        stream << endl;

        for (int i=0; i<nn; i++)
            intermed[i] = phiBH[i];
        // FFT on intermed
        fftw_execute(plan);
        for (int i=0; i<nn; i++)
            stream << intermed[i]/div[i] << ' ';
        stream << endl;
    }
    stream << free_energy() << endl;
}

void RScft::set_pot(double gg)
{
    for (int i=0; i<ny; i++)
        pot[i] = -gg/eps*exp(-0.5*(i+1)*(i+1)*dy*dy/eps/eps);
    for (int i=1; i<nx/2; i++)
        for (int j=0; j<ny; j++)
            pot[ny*i+j] = pot[j];
    for (int j=0; j<ny; j++)
        pot[ny*(nx/2)+j] = 0.0;
    for (int i=nx/2+1; i<nx; i++)
        for (int j=0; j<ny; j++)
            pot[ny*i+j] = -pot[j];

    //for (int i=0; i<nn; i++)
    //{
    //    wA[i] -= pot[i];
    //    wB[i] += pot[i];
    //}
}

void RScft::set_pot3(double gg)
{
    for (int i=0; i<ny; i++)
        pot[i] = gg/eps*exp(-0.5*(i+1)*(i+1)*dy*dy/eps/eps);
    for (int i=1; i<=nx/6; i++)
        for (int j=0; j<ny; j++)
            pot[ny*i+j] = pot[j];
    for (int i=nx/6+1; i<nx/2; i++)
        for (int j=0; j<ny; j++)
            pot[ny*i+j] = -pot[j];
    for (int j=0; j<ny; j++)
        pot[ny*(nx/2)+j] = 0.0;
    for (int i=nx/2+1; i<nx*5/6; i++)
        for (int j=0; j<ny; j++)
            pot[ny*i+j] = pot[j];
    for (int i=nx*5/6; i<nx; i++)
        for (int j=0; j<ny; j++)
            pot[ny*i+j] = -pot[j];
}

RScft::RScft(double segregation, double period_length, int mm, int xx, int ss)
{
    chiN = segregation;
    len = period_length/2.0;

    mid = mm;
    nx = nn = xx;
    ny = 1;
    ns = ss;

    gc = false;
    surf = false;

    set_up();
}

// The sole purpose of dummy is to distinguish the constructor from the gc
// constructor.  The variable dummy is never read.
RScft::RScft(double segregation, double period_length, int mm, int xx, int ss,
                                                  bool dummy, double silon)
{
    chiN = segregation;
    len = period_length;

    mid = mm;
    nx = nn = xx;
    ny = 1;
    ns = ss;

    eps = silon;

    gc = false;
    surf = true;

    set_up();
}

RScft::RScft(double segregation, double period_length, int mm, int xx, int ss,
                                                                     double z)
{
    chiN = segregation;
    len = period_length/2.0;

    mid = mm;
    nx = nn = xx;
    ny = 1;
    ns = ss;

    alphaA = mid;
    alphaB = ns-1-mid;
    zA = zB = z;

    gc = true;
    surf = false;

    set_up();
}

RScft::RScft(double segregation, double period_length, int mm, int xx, int ss,
                                                       double z, double silon)
{
    chiN = segregation;
    len = period_length;

    mid = mm;
    nx = nn = xx;
    ny = 1;
    ns = ss;

    alphaA = mid;
    alphaB = ns-1-mid;
    zA = zB = z;

    eps = silon;

    gc = true;
    surf = true;

    set_up();
}

RScft::RScft(double segregation, double period_length, double thickness,
                int mm, int xx, int yy, int ss, double z, double silon)
{
    chiN = segregation;
    len = period_length/2.0;  // Only half a period is used.
    delta = thickness;

    mid = mm;
    nx = xx;
    ny = yy;
    nn = nx*ny;
    ns = ss;

    alphaA = mid;
    alphaB = ns-1-mid;
    zA = zB = z;

    eps = silon;

    gc = true;
    surf = true;  // Not used with ny > 1.

    set_up();
}

RScft::~RScft()
{
    dealloc_mem();
}
