#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <vector>
#include <array>
#include <cmath>
#include <limits>
#include <cstdlib>
#include <algorithm>
#include <string>

#include "dodecahedron.h"

using namespace std;

bool file_exists (std::string& name)
{
    ifstream f(name.c_str());
    return f.good();
}



string getFileName(int i)
{
    std::stringstream ss;
    string str;

    ss << "md000" << i << ".xtc";
    ss >> str;

    return str;
}

string getFileName_res(int i)
{
    std::stringstream ss;
    string str;

    if( i<10 )
    {
        ss << "res.000" << i << ".xtc";
        ss >> str;
    }
    else
    {
        ss << "res.00" << i << ".xtc";
        ss >> str;
    }

    return str;
}

vector<string> get_file_vector(int i, int j)
{
    vector<string> vec;
    string str;
    for(int ii=i; ii<=j; ++ii)
    {
        str = getFileName(i);
        //str = getFileName_res(i);
        if( file_exists(str) )
        {
            vec.push_back(str);
        }
    }
}

void analyzeAll(int start, int stop, double force_aver)
{
    Dodecahedron dode;
    XTCAnalysis xtc_analysis;
    string first, second;
    vector<string> xtc_infile;
    //xtc_infile = get_file_vector(1, 50);
    xtc_infile.push_back("file.xtc");

    /*for(double i = 16.5652; i< 19.7162+5.0; i+=0.01)
    {
       cout << i-19.7162 << " " << slowAnalys.potential(i, 7.2, 5.0, 19.7162 ) << endl;
    }*/

    xtc_analysis.fluctuationAnalyse(xtc_infile[0], start, stop, force_aver);

    return;

    //
    // Locate the xtc file where release happened
    //
    for(string str : xtc_infile)
    {
        second = first;
        if( dode.slowRele == numeric_limits<int>::max() )
        {
            first = str;
            dode.slowRele = xtc_analysis.slowReleAnalyse(str,1); // look at first frame, is genome detected outside capsid
        }
        if(dode.slowRele != numeric_limits<int>::max()) { break; }
    }

    //
    //
    //
    if(dode.slowRele != numeric_limits<int>::max() || first.compare( getFileName(1) ) == 0 || first.compare( second ) == 0 ) // analyze xtc where release occured
    {
        if( !second.empty() )
        {
            dode.slowRele = xtc_analysis.slowReleAnalyse(second);
        }
        else
        {
            dode.slowRele = xtc_analysis.slowReleAnalyse(first);
        }
    }

    dode.isfold2 = xtc_analysis.isfold2;
    dode.isfold3 = xtc_analysis.isfold3;
    dode.ismulti_slow = xtc_analysis.ismultiple_slow;
    dode.rim_time = xtc_analysis.rim_time;

    if( dode.loadAnalyze("thermo", start, stop)) // load and at the same time analyze -> several GB worth of data, must be done at the same time
    {
        double stop = xtc_analysis.rim_time;
        stop = ( (stop > dode.getRuptureTime()) ? dode.getRuptureTime() : stop );
        stop = ( (stop < 11000) ? 1000000 : stop );
        xtc_analysis.fluctuationAnalyse(getFileName(1), stop-10000);

        double average = 0.0, variance = 0.0;
        for(int i=0; i<240; ++i)
        {
            //cout << slowAnalys.aver[i].getAverage() << " " << slowAnalys.aver[i].getVariance() << endl;
            dode.average += xtc_analysis.aver[i].getAverage();
            dode.variance += xtc_analysis.aver[i].getVariance();
        }
        dode.average /= 240.0;
        dode.variance /= 240.0;

        dode.info();
    }


}

double analyzeThermo(int start, int stop)
{
    Dodecahedron dode;
    if( dode.loadAnalyze("thermo", start, stop) ) // load and at the same time analyze -> several GB worth of data, must be done at the same time
    {
        ofstream average_force;
        average_force.open("average_penta_penta_force");
        average_force << dode.penta_penta_force.getAverage() << endl;
        average_force.close();
        dode.info();
        return dode.penta_penta_force.getAverage();
    }
    return 0.0;
}

/**
 *
 * Analyze capsid genome release
 *
 * Rupture - one pentamer-pentamer bond is 0, getRuptureTime()
 * slowRelease - genome outside of capsid and for 100ns - (100 000 steps) no bonds at 0
 *
 */
int main()
{    
    int start=0;//10*1000;
    int stop=550*1000*1000;
    double force_aver = analyzeThermo(start, stop);
    //analyzeAll(start, stop,force_aver);
    /*Dodecahedron dode;
    dode.mollweide.genDiagram("fragment.svg", "fragment");
    dode.mollweide.genDiagram("rupture.svg", "rupture");
    dode.mollweide.genDiagram("peel.svg", "peel");

    dode.mollweide.genDiagram("fragment_hamming.svg", "fragment_hamming");
    dode.mollweide.genDiagram("rupture_hamming.svg", "rupture_hamming");
    dode.mollweide.genDiagram("peel_hamming.svg", "peel_hamming");*/

    return 0;
}

