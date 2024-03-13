#ifndef DIAGRAM_H
#define DIAGRAM_H

//#include <iostream>
//#include <fstream>
//#include <sstream>
//#include <vector>
#include <array>
//#include <cmath>

#include "pentamermap.h"

using namespace std;

/**
 * @brief The Edge class
 * - Edge in mollweide graph
 * - represents interface between two pentamers
 */
class Edge {
public:
    int N;
    int pentamerA;
    int pentamerB;
};


class Diagram {
public:
    Diagram(){}

    const int height=500;
    const int width=700;
    const int distortion_x=-200;

    //
    // Map of nodes from mollweide projection to real Lammps data
    //
    int pentamerMap[12]; // nodeMap[mollweide penta] = lammps penta

    //
    // Map of connection of lammps data before any breaking
    //
    PentamerMap fullMap;

    vector< array<int,30> > save_state;

    const string connected = string("rgb(0,255,0)");
    const string disconnected = string("rgb(255,0,0)");

    int getRupturePentas(array<int, 2>& penta, vector< PentamerMap >& map, int ruptureTime)
    {
        int connect[12] = {0};
        int rupture_id = -1;
        int count=0;
        int total=0;

        for(int id=0; id<map.size(); ++id)
        {
            if(map[id].time >= ruptureTime) // get first rupture -> identify 2 pentamers between the broken bond is
            {
                rupture_id = id;
                for(int i=0; i<12; ++i)
                {
                    for(int j=i+1; j<12; ++j)
                    {
                        connect[i] += map[id].interactions[i][j];
                        connect[j] += map[id].interactions[i][j];
                        total +=  map[id].interactions[i][j] +  map[id].interactions[i][j];
                    }
                }
                break;
            }
        }

        if(total == 29) {
            for(int i=0; i<12; ++i)
            {
                if( connect[i] == 4 ) {

                    penta[count] = i;
                    ++count;
                }
            }
        }
        else
        {
            penta[0] = pentamerMap[0];
            penta[1] = pentamerMap[1];
        }

        return rupture_id;
    }

    void getHalfCapsidPentamers(array<int, 6>& penta, PentamerMap& map, int fP)
    {
        // identify surrounding pentamers of firstRupturePenta
        int countt = 1;
        penta[0] = fP;
        for(int j=0; j<12; ++j)
        {
            if( map.interactions[fP][j] == 1 || map.interactions[j][fP] == 1 )
            {
                penta[countt] = j;
                ++countt;
            }
        }
    }

    array<int, 6> getOtherCapsid(array<int, 6>& hcapsid)
    {
        array<int, 6> other;
        bool o = true;
        int countt = 0;

        for(int i=0; i<12; ++i)
        {
            o=true;
            for(int j=0; j<6; ++j)
            {
                if( i == hcapsid[j] )
                {
                    o = false;
                }
            }
            if(o)
            {
                other[countt]=i;
                ++countt;
            }
        }

        return other;
    }

    /**
     * @brief getConnectionDistrib
     * @param inside of half-capsid
     * @param border of half-capsid
     * @param outside of half-capsid
     */
    void getConnectionDistrib(int& inside, int& border, int& outside, array<int, 6>& hcapsid, PentamerMap& map)
    {
        array<int, 6> other = getOtherCapsid(hcapsid);

        for(int i=0; i<6; ++i)
        {
            for(int j=0; j<6; ++j)
            {
                inside += map.interactions[ hcapsid[i] ][ hcapsid[j] ] + map.interactions[ hcapsid[j] ][ hcapsid[i] ];
                border += map.interactions[ hcapsid[i] ][ other[j] ] + map.interactions[ other[j] ][ hcapsid[i] ];
                outside += map.interactions[ other[i] ][ other[j] ] + map.interactions[ other[j] ][ other[i] ];
            }
        }

        inside /= 2;
        outside /= 2;
    }

    bool isPeelingRelease(vector< PentamerMap >& map, int ruptureTime, int seq)
    {
        int count=0;
        int rupture_id=0;
        int inside=0, border=0, outside=0;
        int inside2=0, border2=0, outside2=0;

        array<int, 2> pentas = {-1,-1};
        array<int, 6> fir = {-1,-1,-1,-1,-1,-1};
        array<int, 6> sec = {-1,-1,-1,-1,-1,-1};

        rupture_id = getRupturePentas(pentas, map, ruptureTime);
        getHalfCapsidPentamers(fir, map[0], pentas[0]);
        getHalfCapsidPentamers(sec, map[0], pentas[1]);

        getConnectionDistrib(inside, border, outside, fir, map[seq]);
        getConnectionDistrib(inside2, border2, outside2, sec, map[seq]);

        if(outside == 10 || outside2 == 10)
            return true;

        return false;
    }


    /**
     * @brief genDiagram
     * @param out
     * @param map - map from lammps data
     */
    void genDiagram(string out, PentamerMap& map)
    {
        std::fstream fs (out, std::fstream::out);
        std::fstream fs2 ("state", std::fstream::out);

        vector<double> point_x;
        vector<double> point_y;

        beginHtml(fs);
        string diagramVertices = generateVertices(point_x, point_y);

        for(int i=0; i<20; ++i) // draw connections between vertices -> 3-fold axis symetry points
        {
          for(int j=i+1; j<20; ++j)
          {
            if(isDraw(i,j))
            {
              drawEdge(fs, i, j, map, point_x, point_y);
              fs2 << isConnected(i, j, map) << " ";
            }
          }
        }

        fs << diagramVertices; // draw Vertices

        endHtml(fs);

        fs.close();
        fs2.close();
    }

    string getColor(int max, int color)
    {
        //string("diagram") + std::to_string(states[seq].time) + string(".svg")
        //string("rgb(0,255,0)"); // connected
        //string("rgb(255,0,0)"); // disconected
        if(color > max/2)
        {
            return string("rgb(") + std::to_string(2*(255.0-255.0*color/max)) + string(",255,0)");
        }
        else
        {
            return string("rgb(255,") + std::to_string(2*(255.0*color/max)) + string(",0)");
        }
    }

    void genDiagram(string out, string state)
    {
        std::fstream fs (out, std::fstream::out);
        std::ifstream fs2 (state);

        vector<double> point_x;
        vector<double> point_y;
        string rgb;
        int num[30];
        int ii=0;
        int max=0;

        beginHtml(fs);
        string diagramVertices = generateVertices(point_x, point_y);

        for(int i=0; i<30; ++i)
            fs2 >> num[i];

        for(int i=0; i<30; ++i)
        {
            if(max < num[i])
                max = num[i];
        }

        for(int i=0; i<20; ++i) // draw connections between vertices -> 3-fold axis symetry points
        {
          for(int j=i+1; j<20; ++j)
          {
            if(isDraw(i,j))
            {
              drawEdge(fs, i, j, getColor(max, num[ii]), point_x, point_y);
              ++ii;
            }
          }
        }

        fs << diagramVertices; // draw Vertices

        endHtml(fs);

        fs.close();
        fs2.close();
    }

    int getConnectionNumber(PentamerMap map)
    {
        const int connectedd = 1;
        int count = 0;

        for(int i=0; i<20; ++i)
        {
            for(int j=i+1; j<20; ++j)
            {
                if( isDraw(i,j) )
                {
                if( getColor(i, j, map) == connected )
                    ++count;
                }
            }
        }

        return count;
    }

    void buildMollweideConnectivity(array<int, 12>& nei)
    {
        // build mollweide connectivity
        nei[0] = 0;
        nei[1] = getNeighbor(nei[0], 0);
        nei[2] = getNeighbor(nei[0], nei[1], true); // 0,1,2 - around 3 fold axis

        // around pentamer nei[0]:
        nei[3] = getNeighbor(nei[0], nei[1], nei[2]); // 0,1,3 - around 3 fold axis
        nei[4] = getNeighbor(nei[0], nei[3], nei[1]); // 0,3,4
        nei[5] = getNeighbor(nei[0], nei[4], nei[3]); // 0,4,5 // or getNeighbor(nei[0], nei[2], nei[1]);

        // opposite pentamer nei[0]
        nei[6] = getNeighbor(nei[1], nei[2], nei[0]); // 1,2,6
        nei[7] = getNeighbor(nei[1], nei[3], nei[0]); // 1,3,7
        nei[8] = getNeighbor(nei[2], nei[5], nei[0]); // 2,5,8
        nei[9] = getNeighbor(nei[3], nei[4], nei[0]); // 3,4,9
        nei[10] = getNeighbor(nei[4], nei[5], nei[0]);// 4,5,10

        nei[11] = getNeighbor(nei[6], nei[7], nei[1]); // opposite pentamer nei[1] // 6,7,11
    }


    /**
     * @brief setNodeMap
     * @param edge
     * @param i - pentamer id
     * @param j
     * @param k
     */
    void setNodeMap(vector<Edge> edgeFull, int i=0, int j=0, int k=0)
    {
        array<int, 12> nei;
        array<int, 12> neiLammps;
        buildMollweideConnectivity(nei);

        neiLammps[0] = i;
        neiLammps[1] = getNeighbor(edgeFull, neiLammps[0], j);
        neiLammps[2] = getNeighbor(edgeFull, neiLammps[0], neiLammps[1], (k == 0));

        // around pentamer nei[0]:
        neiLammps[3] = getNeighbor(edgeFull, neiLammps[0], neiLammps[1], neiLammps[2]); // 0,1,3 - around 3 fold axis
        neiLammps[4] = getNeighbor(edgeFull, neiLammps[0], neiLammps[3], neiLammps[1]); // 0,3,4
        neiLammps[5] = getNeighbor(edgeFull, neiLammps[0], neiLammps[4], neiLammps[3]); // 0,4,5 // or getNeighbor(nei[0], nei[2], nei[1]);

        // opposite pentamer nei[0]
        neiLammps[6] = getNeighbor(edgeFull, neiLammps[1], neiLammps[2], neiLammps[0]); // 1,2,6
        neiLammps[7] = getNeighbor(edgeFull, neiLammps[1], neiLammps[3], neiLammps[0]); // 1,3,7
        neiLammps[8] = getNeighbor(edgeFull, neiLammps[2], neiLammps[5], neiLammps[0]); // 2,5,8
        neiLammps[9] = getNeighbor(edgeFull, neiLammps[3], neiLammps[4], neiLammps[0]); // 3,4,9
        neiLammps[10] = getNeighbor(edgeFull, neiLammps[4], neiLammps[5], neiLammps[0]);// 4,5,10

        neiLammps[11] = getNeighbor(edgeFull, neiLammps[6], neiLammps[7], neiLammps[1]); // opposite pentamer nei[1] // 6,7,11

        for(int i=0; i<12; ++i)
        {
            pentamerMap[nei[i]] = neiLammps[i];
        }

        /*pentamerMap[0] = 0;
        pentamerMap[1] = 8;
        pentamerMap[2] = 2;
        pentamerMap[3] = 4;
        pentamerMap[4] = 5;

        pentamerMap[5] = 9;
        pentamerMap[6] = 7;
        pentamerMap[7] = 6;
        pentamerMap[8] = 10;
        pentamerMap[9] = 3;

        pentamerMap[10] = 11;
        pentamerMap[11] = 1;

        for(int i=0; i<12; ++i)
        {
            cout << i << " " << nei[i] << " == " << neiLammps[i] << endl;
        }*/
    }

    int getNeighbor(vector<Edge> edgeFull, int pentamer, int other)
    {
        int ii=0;
        for(int i=0; i<30; ++i)
        {
            if(edgeFull[i].pentamerA == pentamer || edgeFull[i].pentamerB == pentamer)
            {
                 if(ii==other)
                 {
                     if( edgeFull[i].pentamerA == pentamer )
                     {
                         return edgeFull[i].pentamerB;
                     }
                     else
                     {
                         return edgeFull[i].pentamerA;
                     }
                 }
                 else
                 {
                     ++ii;
                 }
            }
        }
        cerr << "diagram.h::getNeighbor(vector<Edge> edgeFull, int pentamer, int other), too many neighbors" << endl;
        return -1;
    }

    int getNeighbor(vector<Edge> edgeFull, int pentamerA, int pentamerB, bool other)
    {
        int ii=0;
        int neighbor[2] = {-2,-3};
        for(int i=0; i<5; ++i)
        {
            for(int j=0; j<5; ++j)
            {
                if( getNeighbor(edgeFull, pentamerA, i) == getNeighbor(edgeFull, pentamerB, j))
                {
                    neighbor[ii] = getNeighbor(edgeFull, pentamerA, i);
                    ++ii;
                }
            }
        }
        if( ii > 2 )
        {
            cerr << "ERROR: diagram.h::getNeighbor(vector<Edge> edgeFull, int pentamerA, int pentamerB, int other), 3 pentamer neighbors found" << endl;
            return -1;
        }
        if(other)
            return neighbor[0];
        else
            return neighbor[1];
    }

    int getNeighbor(vector<Edge> edgeFull, int pentamerA, int pentamerB, int pentamerC)
    {
        if( getNeighbor(edgeFull, pentamerA, pentamerB, true) == pentamerC)
        {
            return getNeighbor(edgeFull, pentamerA, pentamerB, false);
        }
        else
        {
            return getNeighbor(edgeFull, pentamerA, pentamerB, true);
        }
    }

    int getNeighbor(int mollweidePentamer, int other)
    {
        return mollweideNeighbor[mollweidePentamer][other];
    }

    int getNeighbor(int mollweidePentamerA, int mollweidePentamerB, int mollweidePentamerC)
    {
        if( getNeighbor(mollweidePentamerA, mollweidePentamerB, true) == mollweidePentamerC)
        {
            return getNeighbor(mollweidePentamerA, mollweidePentamerB, false);
        }
        else
        {
            return getNeighbor(mollweidePentamerA, mollweidePentamerB, true);
        }
    }

    int getNeighbor(int mollweidePentamerA, int mollweidePentamerB, bool other)
    {
        int ii=0;
        int neighbor[2] = {-2,-3};
        for(int i=0; i<5; ++i)
        {
            for(int j=0; j<5; ++j)
            {
                if( getNeighbor(mollweidePentamerA,i) == getNeighbor(mollweidePentamerB,j) )
                {
                    neighbor[ii] = getNeighbor(mollweidePentamerA,i);
                    ++ii;
                }
            }
        }
        if( ii > 2 )
        {
            cerr << "ERROR: diagram.h::getNeighbor(int mollweidePentamerA, int mollweidePentamerB, int other), 3 pentamer neighbors found" << endl;
            return -1;
        }
        if(other)
            return neighbor[0];
        else
            return neighbor[1];
    }

    int isConnected(int verticeA, int verticeB, PentamerMap map)
    {
        string con = getColor(verticeA, verticeB, map);
        if( con == connected)
            return 1;
        else
            return 0;
    }

//private:
    void beginHtml (std::fstream& fs)
    {
        fs << "<!DOCTYPE html>" << endl;
        fs << "<html>" << endl;
        fs << "<body>" << endl;
        fs << endl;
        fs << "<svg height=\"" << height <<"\" width=\"" << width << "\">" << endl;
    }

    void endHtml (std::fstream& fs)
    {
        fs << "</svg>" << endl;
        fs << endl;
        fs << "</body>" << endl;
        fs << "</html>" << endl;
    }

    void drawEdge (std::fstream& fs, int i, int j, PentamerMap& map, vector<double>& point_x, vector<double>& point_y)
    {
        string rgb;

        rgb=getColor(i, j, map); //`analyze $i $j $thermo`
        if( i == 0 && j == 1 ) {
          fs << "  <path d=\"M " << point_x[i] << " " << point_y[i] << " C " << distortion_x << " -30, " << width - distortion_x
             << " -30, " << point_x[j] << " " << point_y[j] << "\" style=\"stroke:rgb(0,0,0);stroke-width:10\" fill=\"transparent\" fill-opacity=\"0.0\" />" << endl;
          fs << "  <path d=\"M " << point_x[i] << " " << point_y[i] << " C " << distortion_x << " -30, " << width - distortion_x
             << " -30, " << point_x[j] << " " << point_y[j] << "\" style=\"stroke:" <<  rgb << ";stroke-width:6\" fill=\"transparent\" fill-opacity=\"0.0\" />" << endl;
        } else {
          fs << "  <line x1=\"" << point_x[i] << "\" y1=\"" << point_y[i] << "\" x2=\"" << point_x[j] << "\" y2=\"" << point_y[j]
                << "\" style=\"stroke:rgb(0,0,0);stroke-width:10\" />" << endl;
          fs << "  <line x1=\"" << point_x[i] << "\" y1=\"" << point_y[i] << "\" x2=\"" << point_x[j] << "\" y2=\"" << point_y[j]
                << "\" style=\"stroke:" << rgb << ";stroke-width:6\" />" << endl;
        }
    }

    void drawEdge (std::fstream& fs, int i, int j, string rgb, vector<double>& point_x, vector<double>& point_y)
    {
        if( i == 0 && j == 1 ) {
          fs << "  <path d=\"M " << point_x[i] << " " << point_y[i] << " C " << distortion_x << " -30, " << width - distortion_x
             << " -30, " << point_x[j] << " " << point_y[j] << "\" style=\"stroke:rgb(0,0,0);stroke-width:10\" fill=\"transparent\" fill-opacity=\"0.0\" />" << endl;
          fs << "  <path d=\"M " << point_x[i] << " " << point_y[i] << " C " << distortion_x << " -30, " << width - distortion_x
             << " -30, " << point_x[j] << " " << point_y[j] << "\" style=\"stroke:" <<  rgb << ";stroke-width:6\" fill=\"transparent\" fill-opacity=\"0.0\" />" << endl;
        } else {
          fs << "  <line x1=\"" << point_x[i] << "\" y1=\"" << point_y[i] << "\" x2=\"" << point_x[j] << "\" y2=\"" << point_y[j]
                << "\" style=\"stroke:rgb(0,0,0);stroke-width:10\" />" << endl;
          fs << "  <line x1=\"" << point_x[i] << "\" y1=\"" << point_y[i] << "\" x2=\"" << point_x[j] << "\" y2=\"" << point_y[j]
                << "\" style=\"stroke:" << rgb << ";stroke-width:6\" />" << endl;
        }
    }

    /**
     * @brief getColor - return collor of Edge in diagram based on PentamerMap
     * @param verticeA - vertices in mollweide == 3-fold axis symetry points
     * @param verticeB
     * @param map
     * @return
     */
    string getColor(int verticeA, int verticeB, PentamerMap map)
    {
        array<int, 2> pentamer = getPentamer(verticeA, verticeB);
        if( map.interactions[ pentamer[0] ][ pentamer[1] ] == 1) {
            return connected;
        }
        return disconnected;
    }

    string generateVertices(vector<double>& point_x, vector<double>& point_y)
    {
        std::stringstream ss;
        double x,y,z,hack_y, xx, yy;
        double offset_x = width/2;
        double offset_y = height/2;
        double dot_size=10.0;

        for(int i=0;i<20;++i) {

            hack_y=0.0;

            x=vertices[i][0];
            y=vertices[i][1];
            z=vertices[i][2];

            if( i == 10 ) {
                x=vertices[7][0];
                y=vertices[7][1];
                z=vertices[7][2];
                hack_y=40.0;
            }

            if( i == 2 ) {
                x=vertices[13][0];
                y=vertices[13][1];
                z=vertices[13][2];
                hack_y=-40.0;
            }

            mollweideProject(x, y, z, xx, yy);

            xx += offset_x;
            yy += offset_y + hack_y;

            point_x.push_back(xx);
            point_y.push_back(yy);

            ss << "  <circle cx=\"" << xx << "\" cy=\"" << yy << "\" r=\"" << dot_size << "\" stroke=\"black\" stroke-width=\"2\" fill=\"grey\" />" << endl;

            /*xx += -3.0;
            yy += 5.0;
            ss << "  <text x=\"" << xx << "\" y=\"" << yy << "\" fill=\"black\">" << i << "</text>" << endl;*/
        }

        return ss.str();
    }


    /**
     * @brief getPentamer
     * @param nodeA
     * @param nodeB
     * @return
     */
    array<int, 2> getPentamer(int verticeA, int verticeB)
    {
        //
        // verticeA + verticeB -> edge between 2 pentamers
        //
        array<int, 2> out;
        int ii=0;
        bool verA=false;
        bool verB=false;
        for(int i=0; i<12; ++i) // 12 pentamers
        {
            verA=false;
            verB=false;
            for(int j=0; j<5; ++j) // each 5 vertices
            {
                if( verticeA == pentamers_mollweide[i][j] )
                {
                    verA=true;
                }
                if( verticeB == pentamers_mollweide[i][j] )
                {
                    verB=true;
                }
            }
            if( verA && verB )
            {
                out[ii] = i;
                ++ii;
            }
        }

        out[0] = pentamerMap[out[0]];
        out[1] = pentamerMap[out[1]];

        if(out[1] < out[0]) {
            int temp = out[1];
            out[1] = out[0];
            out[0] = temp;
        }

        return out;
    }

    void mollweideProject(double x, double y, double z, double& xx, double& yy)
    {
        double R=100.0;
        double lambda_0 = 0.0;
        double sigma = asin(z / R);
        double lambda = atan2(y,x);

        xx= R * 2 * sqrt(2) / 3.14159265359 * (lambda-lambda_0) * cos(sigma);
        yy= R * sqrt(2) *sin(sigma);
    }


    /**
     * @brief isDraw
     * @param i - node from mollweide projection
     * @param j - node from mollweide projection
     * @return
     */
    bool isDraw(int i, int j)
    {
        if( i == 0  && ( j == 1  || j == 3 || j == 11 ) ) { return true; }
        if( i == 1  && ( j == 4  || j == 12 ) ) { return true; }
        if( i == 2  && ( j == 3  || j == 13 || j == 4 ) ) { return true; }
        if( i == 3  && ( j == 16 ) ) { return true; }
        if( i == 4  && ( j == 17 ) ) { return true; }
        if( i == 5  && ( j == 6  || j == 8  || j == 14 ) ) { return true; }
        if( i == 6  && ( j == 9  || j == 15 ) ) { return true; }
        if( i == 7  && ( j == 8  || j == 9  || j == 10 ) ) { return true; }
        if( i == 8  && ( j == 18 ) ) { return true; }
        if( i == 9  && ( j == 19 ) ) { return true; }
        if( i == 10 && ( j == 11 || j == 12 ) ) { return true; }
        if( i == 11 && ( j == 18 ) ) { return true; }
        if( i == 12 && ( j == 19 ) ) { return true; }
        if( i == 13 && ( j == 14 || j == 15 ) ) { return true; }
        if( i == 14 && ( j == 16 ) ) { return true; }
        if( i == 15 && ( j == 17 ) ) { return true; }
        if( i == 16 && ( j == 18 ) ) { return true; }
        if( i == 17 && ( j == 19 ) ) { return true; }

        return false;
    }

protected:
    //
    // Order set partially by generation program and manual rotation in vmd
    //
    const double vertices[20][3] = {
        {-87.274673, -33.041664,  -3.031746},
        {-87.283112,  33.624931,  -3.073639},
        {-33.165180,   0.244000, -90.209511},
        {-53.827793, -53.671818, -56.884789},
        {-53.841656,  54.196720, -56.952446},
        { 87.260803, -33.019405,  -2.673706},
        { 87.252487,  33.646938,  -2.715416},
        { 33.143063,   0.361232,  84.462440},
        { 53.819481, -53.591534,  51.204964},
        { 53.805992,  54.277565,  51.137688},
        {-33.523731,   0.352736,  84.325867},
        {-54.049580, -53.604828,  50.984009},
        {-54.063107,  54.263775,  50.916424},
        { 33.501408,   0.252555, -90.072762},
        { 54.041134, -53.658440, -56.663689},
        { 54.027405,  54.210564, -56.731102},
        {  0.068443, -86.985161, -36.152351},
        {  0.046722,  87.550453, -36.261459},
        { -0.068176, -86.943314,  30.514156},
        { -0.090363,  87.592201,  30.405018}
    };

    //
    // The order is set by vertices. Vertices order essentially random. [mollweide_pentamer][vertice_count] = vertice_index
    //
    const int pentamers_mollweide[12][5] = {
        {5,6,13,14,15},
        {5,6,7,8,9},
        {6,9,15,17,19},
        {5,8,14,16,18},
        {2,3,13,14,16},

        {2,4,13,15,17},
        {1,4,12,17,19},
        {7,9,10,12,19},
        {7,8,10,11,18},
        {0,3,11,16,18},

        {0,1,2,3,4},
        {0,1,10,11,12}
    };

    // mollweide pentamers
    // 0 - 1,2,3,4,5     6  - 2,5,7,10,11
    // 1 - 0,2,3,7,8     7  - 1,2,6,8,11
    // 2 - 0,1,5,6,7     8  - 1,3,7,9,11
    // 3 - 0,1,4,8,9     9  - 3,4,8,10,11
    // 4 - 0,3,5,9,10    10 - 4,5,6,9,11
    // 5 - 0,2,4,6,10    11 - 6,7,8,9,10
    int mollweideNeighbor[12][5] = {{1,2,3,4,5},   // 0
                                    {0,2,3,7,8},   // 1
                                    {0,1,5,6,7},   // 2
                                    {0,1,4,8,9},   // 3
                                    {0,3,5,9,10},  // 4
                                    {0,2,4,6,10},  // 5
                                    {2,5,7,10,11}, // 6
                                    {1,2,6,8,11},  // 7
                                    {1,3,7,9,11},  // 8
                                    {3,4,8,10,11}, // 9
                                    {4,5,6,9,11},  // 10
                                    {6,7,8,9,10},};// 11
};

#endif // DIAGRAM_H
