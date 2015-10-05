#include "rscft.h"
//#include <ctime>
#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdlib>
using namespace std;

int main(int argc, char **argv)
{
/*
    ofstream outfile1(argv[1]);
    ofstream outfile2(argv[2]);
    ofstream outfile3(argv[3]);
    cout.precision(17);
    outfile1.precision(17);
    outfile2.precision(17);
    outfile3.precision(17);

    RScft obj(37.6, 4.0, 500, 9, 1001, 1.45);
    obj.set_alphaA(200);
    obj.set_alphaB(200);
    double wA[129], wB[129];
    for (int i=0; i<64; i++)
    {
        wA[i] = 35.0;
        wB[i] = -5.0;
    }
    wA[64] = wB[64] = 15.0;
    for (int i=65; i<129; i++)
    {
        wA[i] = -5.0;
        wB[i] = 35.0;
    }
    obj.set_field(wA, wB);
    //obj.set_field_rand(1);
    obj.calc_field(outfile1);
    obj.write_plot(outfile2);
    obj.fields2fs(outfile3);

    outfile1.close();
    outfile2.close();
    outfile3.close();

    return 0;
*
    ifstream infile(argv[1]);
    ofstream outfile(argv[2]);
    outfile.precision(17);
    RScft obj(37.6, 5.0, 1.43, 50, 65, 63, 101, 0.271, 0.15);
    obj.read_data(infile);

    outfile << "# chiN=" << 37.6 << ", period=" << 5.0 << ", z=" << 0.271
            << endl;
    outfile << "# thickness=" << 1.43 << ", eps=" << 0.15 << endl;
    outfile << "# phiA profile in x direction.\n";
    obj.profile_phiA_x(outfile);
    return 0;
*/
    if (argc < 10)
    {
        cout << "Usage: rscft <chiN> <len> <z> <alpha> <pot> <outfile1>"
             << " <outfile2> <infile> <thickness>\n";
        return 0;
    }

    //time_t seconds;

    //int s = atoi(argv[1]);
    double chiN = atof(argv[1]);
    double len = atof(argv[2]);
    double z = atof(argv[3]);
    int alpha = atoi(argv[4]);
    double pot = atof(argv[5]);
    ofstream outfile1(argv[6]);
    ofstream outfile2(argv[7]);
    ifstream infile(argv[8]);
    double thick = atof(argv[9]);

    //time(&seconds);
    //cout << seconds << endl;
    //RScft obj(chiN, len, 60, 401, 121);
    //RScft obj(chiN, len, 50, 129, 101, z);
    //RScft obj(chiN, len, 60, 399, 121, z, eps);

    RScft obj(chiN, len, thick, 250, 65, 63, 501, z, 0.15);
    //RScft obj(chiN, 1.43, 250, 63, 501, z, 0.15);

    //time(&seconds);
    //cout << seconds << endl;

    //RScft obj(chiN, len, 500, 129, 1001, z);
    obj.set_alphaA(alpha);
    obj.set_alphaB(alpha);
    //double wA[129], wB[129];
    //for (int i=0; i<64; i++)
    //{
    //    wA[i] = 35.0;
    //    wB[i] = -5.0;
    //}
    //wA[64] = wB[64] = 15.0;
    //for (int i=65; i<129; i++)
    //{
    //    wA[i] = -5.0;
    //    wB[i] = 35.0;
    //}
    //obj.set_field(wA, wB);

    cout.precision(17);
    outfile1.precision(17);
    outfile2.precision(17);
    //outfile3.precision(17);

    //obj.read_data(infile);
    infile.close();
    //obj.calc_dens();
    //obj.free_energy();
    //return 0;
    //obj.fields2fs(outfile1);
    obj.set_pot3(pot);
    //obj.set_pot(pot);
    //obj.set_pot(0.09375);

    double wA[65*63], wB[65*63];
    for (int i=0; i<32*63; i++) // 150, 250 -> asym
    {
        wA[i] = 33.3;  // -2.006;
        wB[i] = -3.3;  // 6.4706;
    }
    for (int i=32*63; i<33*63; i++)
    {
        wA[i] = 15.0; // log(z); // 6.4706;
        wB[i] = 15.0; // log(z)+chiN; // -2.006;
    }
    for (int i=33*63; i<65*63; i++)
    {
        wA[i] = -3.3;  // -2.006;
        wB[i] = 33.3;  // 6.4706;
    }

    //for (int i=0; i<65*63; i++)
    //{
    //    wA[i] = log(z);
    //    wB[i] = log(z)+chiN;
    //}

    //for (int i=0; i<1001; i++)
    //    wA[i] = wB[i] = 0.0;
    obj.set_field(wA, wB);
    //obj.set_field_rand(1);

    outfile1 << "# chiN=" << chiN << ", period=" << len << ", z=" << z << endl;
    outfile1 << "# thickness=" << thick << ", eps=" << 0.15 << endl;
    cout << chiN << ' ' << len << ' ' << z << endl;
    obj.calc_field(outfile1);
    outfile2 << "# chiN=" << chiN << ", period=" << len << ", z=" << z << endl;
    outfile2 << "# thickness=" << thick << ", eps=" << 0.15 << endl;
    outfile2 << "# FE: " << obj.free_energy() << endl;
    outfile2 << "# Ave wA: " << obj.ave_wA() << endl;
    outfile2 << "# Ave wB: " << obj.ave_wB() << endl;
    outfile2 << "# Ave ksi: " << obj.ave_ksi() << endl;
    outfile2 << "# Ave phiA: " << obj.ave_phiA() << ", max: "
             << obj.max_phiA() << endl;
    outfile2 << "# Ave phiAH: " << obj.ave_phiAH() << endl;
    outfile2 << "# Ave phiAC: " << obj.ave_phiAC() << endl;
    outfile2 << "# Ave phiB: " << obj.ave_phiB() << ", max: "
             << obj.max_phiB() << endl;
    outfile2 << "# Ave phiBH: " << obj.ave_phiBH() << endl;
    outfile2 << "# Ave phiBC: " << obj.ave_phiBC() << endl;
    //outfile2 << "# hpA excess at chiN=2, z=0.2, and epsilon=" << eps << ": "
    //         << obj.hpA_ex(0.090498756211223261) << endl;
    double vol = obj.get_vol();
    double eff_vol = obj.get_eff_vol();
    double ratio = vol/eff_vol;
    outfile2 << "# Actual volume: " << vol << endl;
    outfile2 << "# Effective volume: " << eff_vol << endl;
    outfile2 << "# Ave phiH: " << (obj.ave_phiAH()+obj.ave_phiBH())*ratio
             << endl;
    outfile2 << "# Ave phiC: " << (obj.ave_phiAC()+obj.ave_phiBC())*ratio
             << endl;
    outfile2 << "# The following data is formatted as in the write_plot "
             << "function.\n";
    obj.write_plot(outfile2);

    //obj.write_q2(outfile1, 4800);

    //obj.fields2fs(outfile3);

    outfile1.close();
    outfile2.close();
    //outfile3.close();

    //obj.calc_qs();

    //time(&seconds);
    //cout << seconds << endl << endl;

    //obj.see_trans(cout);

    return 0;
}
