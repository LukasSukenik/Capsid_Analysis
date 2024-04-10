#ifndef DODECAHEDRON_H
#define DODECAHEDRON_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <array>
#include <cmath>
#include <limits>
#include <algorithm>

#include "diagram.h"
#include "pentamermap.h"
#include "lammpscapsidparam.h"
#include "xtcanalysis.h"

class MapKey {
public:
    MapKey(int i, int j, int k) : i(i), j(j), k(k) {}
    int i;
    int j;
    int k;
};

class Dodecahedron {
private:
    const double interaction_factor = 0.0; // threshold of pentamer interface interaction below which the interaction is consider broken
    const double pressure_limit = 1e-5;
    const int later = 5434782; // 500*1000*1000/92 = 500 ns in steps, 1 step 92 fs

public:
    int slowRele = numeric_limits<int>::max();
    double fold2_dist = 0.0;
    double fold3_dist = 0.0;

    bool isfold2 = false;
    bool isfold3 = false;
    bool ismulti_slow = false;

    int rim_time = -1;

    double average=0.0;
    double variance=0.0;

    Welford_Algo energy[84];
    Welford_Algo penta_penta_force;

    ofstream plot_data;
    ofstream plot_halfs;

    Dodecahedron() {
        for(int i=0;i<12;i++) {
            for(int j=0;j<12;j++) {
                cluster.interactions[i][j] = -1;
            }
        }
        plot_data.open("plot_data");
        plot_halfs.open("plot_halfs");
    }

    ~Dodecahedron()
    {
        plot_data.close();
        plot_halfs.close();
    }

    Diagram mollweide;

    LammpsCapsidParam param_zero_press;
    LammpsCapsidParam param_half_press;
    LammpsCapsidParam param_500ns;
    LammpsCapsidParam param_later;
    LammpsCapsidParam param_full_press;

    double interaction_min;
    int map_correction=0;
    bool cluster_gen=true;
    PentamerMap cluster;
    PentamerMap connectionMap; // per state-change PentamerMap
    vector< PentamerMap > states;
    array< array<vector<int>, 12>, 12 > timed_connectionMap_break; // time of bond breaking
    array< array<vector<int>, 12>, 12 > timed_connectionMap_form; // time of bond breaking
    vector<Edge> edgesFull;

    vector<double> perStep_Aver_penta_penta_E;
    vector<double> perStep_Aver_chain_penta_E;

    /**
     * @brief load - Load interaction data from log.lammps
     * @return
     */
    bool loadAnalyze(string inName, int start, int stop)
    {
        ifstream input (inName, std::ifstream::in);
        string line;
        stringstream ss;
        int timestep;
        LammpsCapsidParam param_current;
        bool first = true;

        while( !input.eof() ) // Lines in input
        {
            getline(input, line);
            ss.str(line);
            ss >> timestep;

            //
            // Store thermo capsid parameters
            //
            for(int i=2; i<param_current.param.size(); ++i)
            {
                ss >> param_current.param[i]; // store data from log
            }

            //
            // Set Interaction thresholds
            //
            if( first )
            {
                setInteractionThreshold(param_current);
                start += timestep;
                stop += timestep;
                first = false;
            }

            //
            // Analyze each step
            //
            if( stepAnalyze(timestep, param_current, start, stop) ) // analyze the step
            {
                break;
            }
            ss.flush();
            ss.clear();
        }
        input.close();

        return true;
    }

    /**
     * @brief info
     * step_break - set by low interaction between genome and capsid
     *
     * For analysis
     * X - before and after time in steps
     * : - before and after pressure in kT
     * % - before and after energy in kT
     *
     * STATE TH time_half_press TZ time_zero_press TR time_rupture AE aver_energy PF press_full PH press_half PZ press_zero
     *
     */
    void info()
    {
        if(cluster_gen) // analyze state of the capsid - either 500ns after zero pressure or on last analyzed step
        {
            if( param_500ns.set )
            {
                cluster = param_500ns.getClusters(interaction_min);
            }
            else
            {
                cluster = param_later.getClusters(interaction_min);
            }
        }

        string release_state = setReleaseState();

        cout << release_state << " TH-X" << param_half_press.step << "X TZ-X" << param_zero_press.step << "X TR-X" << getRuptureTime() << "X AE-:" << getAverageE() << ": PF-%"
             << getAverageP() << "% PH-%" << param_half_press.get_pentamer_chain_energy() << "% PZ-%" << param_zero_press.get_pentamer_chain_energy() << "% " << "TS-X" << slowRele << "X " << "TD-Q" << rim_time << "Q "
             << "AF-Q" << average << "Q AV-Q" << variance << "Q ";

        for(int i=0;i<12;i++) {
            for(int j=11;j>=0;j--) {
                if(cluster.interactions[i][j] != -1) {
                    cout << j+1 << "-mer ";
                    break;
                }
            }
        }

        //
        // Generate diagram
        //
        bool peel = true;
        if( true || release_state == string("Rupture") )
        {
            map_correction=0;
            int min=999;
            bool once = true;
            for(unsigned int seq=0; seq<states.size(); ++seq)
            {
                if( min > states[seq].getNumberOfInteractions() )
                    min = states[seq].getNumberOfInteractions();
            }
            for(unsigned int seq=0; seq<states.size(); ++seq)
            {
                if( states[seq].time >= getRuptureTime())
                {
                    if( map_correction == 0 )
                    {
                        init(seq);
                        map_correction = map(seq);
                        mollweide.genDiagram(string("map.svg"), states[seq]);
                    }

                    if( !mollweide.isPeelingRelease(states, getRuptureTime(), seq) )
                        peel = false;

                    if( min == states[seq].getNumberOfInteractions() && once ) // minimal number of interactions
                    {
                        mollweide.genDiagram(string("diagram") + std::to_string(states[seq].time) + string(".svg"), states[seq]);
                        once = false;
                    }
                }
            }
            if(getRuptureTime() != numeric_limits<int>::max())
            {
                if(peel)
                {
                    cout << " peel " << endl;
                }
                else
                {
                    cout << " meh " << endl;
                }
            }
        }

        cout << endl;
    }

    /**
     * @brief getRuptureTime - time when the a bond broken for more then 20000 steps broke
     * @return time of rupture
     */
    int getRuptureTime()
    {
        int ruptureTime=numeric_limits<int>::max();

        for(int i=0; i<12; ++i)
        {
            for(int j=i+1; j<12; ++j) // loop over pentamer-pentamer bonds
            {
                for( int k=0; k < timed_connectionMap_break[i][j].size(); ++k ) // loop over all stored break times
                {
                    if( getDuration(i,j,k) > 20000 && ruptureTime > timed_connectionMap_break[i][j][k] )
                    {
                        ruptureTime = timed_connectionMap_break[i][j][k];
                    }
                }
            }
        }

        return ruptureTime;
    }

    int thermo_step=10;
    int block_size=1000;
    double pen_chain_E=0.0; double pen_chain_F=0.0;
    double pen_pen_E=0.0; double pen_pen_F=0.0;
    double chain_chain_E=0.0; double chain_chain_F=0.0;

    double apen_pen_E[12]={0.0};
    double apen_pen_F[12]={0.0};

    double halfs[6] = {0.0};
    double half_chain[6] = {0.0};

    double chain_wall_F=0.0;

private:

    /**
     * @brief analyze - per step analysis
     * @param step
     * @return
     */
    bool stepAnalyze(const int timestep, LammpsCapsidParam& param, int start, int stop)
    {
        /*for(int i=18; i<=83; ++i)
        {
            energy[i].push(param.param[i]);
        }*/

        if(timestep > start && timestep < stop)
        {
            pen_chain_E   += param.get_pentamer_chain_energy();    pen_chain_F   += param.get_pentamer_chain_force();
            pen_pen_E     += param.get_pentamer_pentamer_energy(); pen_pen_F     += param.get_pentamer_pentamer_force();
            chain_chain_E += param.get_chain_chain_energy();       chain_chain_F += param.get_chain_chain_force();

            chain_wall_F = param.get_chain_wall_F();

            for(int i=0; i<12; ++i)
            {
                apen_pen_E[i] += param.get_pentamer_pentamer_energy(i);
                apen_pen_F[i] += param.get_pentamer_pentamer_force(i);
            }
            for(int i=0; i<6; ++i)
            {
                halfs[i] += param.get_half_half_F(i);
                half_chain[i] += param.get_half_genome_F(i);
            }

            if(timestep%(block_size*thermo_step) == 0)
            {

                plot_data << timestep << " " << pen_chain_E/block_size << " " << pen_chain_F/block_size // penta-chain
                    << " " << pen_pen_E/block_size << " " << pen_pen_F/block_size                   // penta-penta
                    << " " << chain_chain_E/block_size << " " << chain_chain_F/block_size;          // chain-chain
                for(int i=0; i<12; ++i)
                {
                    plot_data << " " << apen_pen_E[i]/block_size << " " << apen_pen_F[i]/block_size;
                    apen_pen_E[i]=0.0;
                    apen_pen_F[i]=0.0;
                }
                plot_data << "\n";

                plot_halfs << timestep << " ";
                for(int i=0; i<6; ++i)
                {
                    plot_halfs << " " << halfs[i]/block_size << " " << half_chain[i]/block_size;
                    halfs[i]=0.0;
                    half_chain[i]=0.0;
                }
                plot_halfs << " " << chain_wall_F << "\n";

                pen_chain_E=0.0; pen_chain_F=0.0;
                pen_pen_E=0.0; pen_pen_F=0.0;
                chain_chain_E=0.0; chain_chain_F=0.0;
                chain_wall_F=0.0;
            }
            penta_penta_force.push(param.get_pentamer_pentamer_force());
        }

        //
        // Set average penta-penta and chain-penta energy
        //
        if(timestep < 20000)
        {
            perStep_Aver_penta_penta_E.push_back( param.get_pentamer_pentamer_energy()/30 );
            perStep_Aver_chain_penta_E.push_back( param.get_pentamer_chain_energy() );
        }

        //
        // identify genome release -> zero pressure exerted by genome on capsid
        //
        if( timestep > 20000 && !param_zero_press.set && param.get_pentamer_chain_energy() < pressure_limit)  // pressure_limit basically 0, set the time release finished
        {
            param_zero_press.set = true;
            param_zero_press.step = timestep;
            param_zero_press.param = param.param;
            param_zero_press.pressure = param.get_pentamer_chain_energy();
        }

        //
        // identify ongoing genome release -> minimal pressure exerted by genome on capsid
        //
        if( timestep > 20000 && !param_half_press.set && param.get_pentamer_chain_energy() < 0.5 * getAverageP() )
        {
            param_half_press.set = true;
            param_half_press.step = timestep;
            param_half_press.param = param.param;
            param_half_press.pressure = param.get_pentamer_chain_energy();
        }

        //
        // save state after 500 ns after detection of release process start
        //
        if(param_zero_press.set && timestep > later + param_zero_press.step)
        {
            param_500ns.param = param.param;
            param_500ns.set = true;
        }

        //
        // save any changes in pentamer conectivity for mollweide diagram
        //
        if( stateChange(timestep, param) )
        {
            connectionMap.interactions = param.getConnectionMap( interaction_min ); // refresh connection map, next time function call we are looking at difference vs the last change
            connectionMap.time=timestep;
            states.push_back(connectionMap);
        }

        param_later.param = param.param;

        return false;
    }


    string setReleaseState()
    {
        string release_state = "confined";

        if( slowRele != -1 ) {
            if( getRuptureTime() != numeric_limits<int>::max() ) // rupture occurs before slow release. bias by 109000 steps -> 10ns
            {
                release_state = "Rupture";
            }
            else
            {
                if( isfold2 )
                {
                    release_state = "slow2fold";
                }
                if( isfold3 )
                {
                    release_state = "slow3fold";
                }
                if( ismulti_slow )
                {
                    release_state = "slowmulti";
                }
            }
        }

        return release_state;
    }







    /**
     * @brief Dodecahedron::stateChange - determine if new pentamer-pentamer bonds were created or broken
     * - record the change in timed_connectionMap_break or timed_connectionMap_form
     * @param step - time of state change
     * @return 0 = no change, 1 - change, 2 no change after breaking for 100000 steps
     */
    bool stateChange(int step, LammpsCapsidParam &param)
    {
        PentamerMap aktual;
        aktual.interactions = param.getConnectionMap(interaction_min);
        bool stateChange = false;

        for(int i=0; i<12; ++i) {
            for(int j=i+1; j<12; ++j) {
                if(connectionMap.interactions[i][j] != aktual.interactions[i][j]) { // diference vs last step
                    if(aktual.interactions[i][j] == 0) { // Bond disconnected
                        timed_connectionMap_break[i][j].push_back(step);
                    }
                    if(aktual.interactions[i][j] == 1) { // Bond formed
                        timed_connectionMap_form[i][j].push_back(step);
                    }
                    stateChange = true;
                }
            }
        }

        return stateChange;
    }


    /**
     * @brief init - Find the Maximal amount of connected pentamer-pentamer interfaces for a given state is sequence [seq]
     * @param seq - index of state in sequence
     * @return
     */
    int init(int seq)
    {
        int num = states[seq].getNumberOfInteractions();
        for(int i=0; i<12; ++i)
        {
            for(int j=i+1; j<5; ++j)
            {
                for(int k=0; k<2; ++k)
                {
                    mollweide.setNodeMap(edgesFull, i, j, k);

                    if( mollweide.getConnectionNumber(states[seq]) == num )
                    {
                        return num;
                    }
                }
            }
        }
        return -1;
    }

    vector<MapKey> key;

    int map(int seq)
    {
        const int disconnected = 0;
        int connections = mollweide.getConnectionNumber(states[seq]);
        //cout << "Number of connected: " << connections << endl;

        if( connections == 29 && map_correction == 0) // one missing connections
        {
            for(int i=0; i<12; ++i) {
                for(int j=0; j<5; ++j) {
                    for(int k=0; k<2; ++k)
                    {
                        mollweide.setNodeMap(edgesFull, i, j, k);
                        if( mollweide.getConnectionNumber(states[seq]) == 29 && mollweide.isConnected(5,6,states[seq]) == disconnected)
                        {
                            key.push_back(MapKey(i,j,k));
                        }
                    }
                }
            }
            map_correction=1;
            mollweide.setNodeMap(edgesFull, key[0].i, key[0].j, key[0].k);
            if( mollweide.getConnectionNumber(states[seq]) == 29 && mollweide.isConnected(5,6,states[seq]) == disconnected)
            {
                return map(seq+1); // map 2+ bonds
            }
        }

        // two missing connections - in sequence
        if( connections == 28 && map_correction==1) {
            for(int i=0; i<key.size(); ++i)
            {
                mollweide.setNodeMap(edgesFull, key[i].i, key[i].j, key[i].k);
                if( mollweide.getConnectionNumber(states[seq]) == 28 &&
                        mollweide.isConnected(5,6,states[seq]) == disconnected &&
                        mollweide.isConnected(5,14,states[seq]) == disconnected )
                {
                    map_correction=2;
                    return 2;
                }

            }
        }


        // two missing and more connections - at once
        if( connections <= 28 && map_correction==0)
        {
            for(int i=0; i<12; ++i) {
                for(int j=0; j<5; ++j) {
                    for(int k=0; k<2; ++k)
                    {
                        if(map_correction == 0)
                        {
                            mollweide.setNodeMap(edgesFull, i, j, k);
                            if( mollweide.getConnectionNumber(states[seq]) <= 28 &&
                                    mollweide.isConnected(5,6,states[seq]) == disconnected
                                && mollweide.isConnected(5,14,states[seq]) == disconnected )
                            {
                                map_correction=2;
                                return 2;
                            }
                        }
                    }
                }
            }
        }
        return -1;
    }

    int getDuration (int i, int j, int k)
    {
        if( k < timed_connectionMap_form[i][j].size() )
        {
            return timed_connectionMap_form[i][j][k] - timed_connectionMap_break[i][j][k];
        }
        else
        {
            return numeric_limits<int>::max();
        }
        return -1;
    }


    void setInteractionThreshold(LammpsCapsidParam& param)
    {
        double interaction_sum=0.0; // interaction sum in kT
        int count=0; // count of pentamer-pentamer interactions == number of pentamer-pentamer interfaces

        for(int i=18; i<=83; ++i)
        {
            interaction_sum += param.param[i];
            if(param.param[i] != 0.0)
                count++;
        }

        interaction_min = interaction_factor*interaction_sum/count;

        mollweide.fullMap.interactions = param.getConnectionMap(interaction_min);

        param_full_press.set = true;
        param_full_press.pressure = param.get_pentamer_chain_energy();
        param_full_press.param = param.param;

        edgesFull = param.getEdges(interaction_min);
        connectionMap.interactions = param.getConnectionMap(interaction_min);
        connectionMap.time=0;
        states.push_back(connectionMap);
    }

    double getAverageE()
    {
        if( !perStep_Aver_penta_penta_E.empty() )
        {
            double sum=0.0;
            for(unsigned int i=0; i<perStep_Aver_penta_penta_E.size(); ++i)
            {
                sum += perStep_Aver_penta_penta_E[i];
            }
            return sum/perStep_Aver_penta_penta_E.size();
        }
        return 0.0;
    }

    double getAverageP()
    {
        if( !perStep_Aver_chain_penta_E.empty() )
        {
            double sum=0.0;
            for(unsigned int i=0; i<perStep_Aver_chain_penta_E.size(); ++i)
            {
                sum += perStep_Aver_chain_penta_E[i];
            }
            return sum/perStep_Aver_chain_penta_E.size();
        }
        return 0.0;
    }
};

#endif // DODECAHEDRON_H
