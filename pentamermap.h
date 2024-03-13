#ifndef PENTAMERMAP_H
#define PENTAMERMAP_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <array>
#include <cmath>
#include <limits>
#include <algorithm>

using namespace std;

class PentamerMap {
public:
    array< array<int, 12>, 12 > interactions;
    int time=-1;

    void setAllValues(int value)
    {
        for(int i=0;i<12;i++)
            for(int j=0;j<12;j++)
                this->interactions[i][j] = value;
    }

    void show()
    {
        for(int i=0;i<12;i++)
            for(int j=0;j<12;j++)
                cout << this->interactions[i][j] << " ";
    }

    int getNumberOfInteractions()
    {
        int sum=0;
        for(int i=0;i<12;i++)
            for(int j=0;j<12;j++)
                if( this->interactions[i][j] == 1 )
                    ++sum;
        return sum;
    }
};

#endif // PENTAMERMAP_H
