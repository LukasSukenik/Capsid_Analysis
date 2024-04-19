#ifndef HISTOGRAM_H
#define HISTOGRAM_H


#include <string>
#include <limits>
#include <cmath>
#include <iostream>
#include <vector>
#include "welford.h"



using namespace std;

class Histogram
{
    Welford_Algo val_average;
    ofstream histo;
    ofstream force_value;
    vector<int> bins;
    vector<int> bins_neg;
    double factor;
    int count=0;
    double sigma;
    double bead_radius;
    double negative_limit;
public:
    Histogram(){}

    void init(double sigm, int size, double rm, int extra_size)
    {
        histo.open ("histogram");
        force_value.open ("force_value");
        bins.resize(size);
        for(int i=0; i<bins.size(); ++i)
        {
            bins[i] = 0;
        }
        bins_neg.resize(extra_size);
        for(int i=0; i<bins_neg.size(); ++i)
        {
            bins_neg[i] = 0;
        }
        factor = (double)size / sigm;
        sigma = sigm;
        bead_radius = rm;
        negative_limit = extra_size / factor;
    }

    void add(double dist)
    {
        count++;
        if(dist < 0.0 && dist > -1.0*negative_limit)
        {
            bins_neg[(int)(-1.0*dist*factor)]++;
        }
        if(dist > 0.0 && dist < sigma)
        {
          bins[(int)(dist*factor)]++;
        }
        val_average.push(dist);
    }

    void print()
    {
        for(int i=0; i<bins_neg.size(); ++i)
        {
            histo << -i/factor << " " << ( (double)bins_neg[i] ) /( ( (double)count ) / 240.0 ) << "\n";
        }
        for(int i=0; i<bins.size(); ++i)
        {
            histo << i/factor << " " << ( (double)bins[i] ) /( ( (double)count ) / 240.0 ) << "\n";
        }
        force_value << val_average.getAverage() << "\n"; // still a little off
    }

    ~Histogram()
    {
        histo.close();
        force_value.close();
    }
};

#endif // HISTOGRAM_H
