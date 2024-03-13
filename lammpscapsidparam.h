#ifndef LAMMPSCAPSIDPARAM_H
#define LAMMPSCAPSIDPARAM_H

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
#include "welford.h"

/**
 * @brief The LammpsCapsidParam class
 * \$1 step
 * \$2 temp
 * \$3 pe
 * \$4 ke
 *
 * Temperature
 * \$5 temp_penta_0  \$6 temp_penta_1  \$7 t_penta_2 \$8 t_penta_3 \$9 t_penta_4 \$10 t_penta_5 \$11 t_penta_6 \$12 t_penta_7 \$13 t_penta_8 \$14 t_penta_9 \$15 t_penta_10 \$16 t_penta_11
 * \$17 temp_chain
 *
 * Interaction energy columns: Pentamer-Pentamer
 *      #            1     2     3     4     5     6     7     8     9    10    11
 *      # 0       \$18  \$19  \$20  \$21  \$22  \$23  \$24  \$25  \$26  \$27  \$28
 *      # 1             \$29  \$30  \$31  \$32  \$33  \$34  \$35  \$36  \$37  \$38
 *      # 2                   \$39  \$40  \$41  \$42  \$43  \$44  \$45  \$46  \$47
 *      # 3                         \$48  \$49  \$50  \$51  \$52  \$53  \$54  \$55
 *      # 4                               \$56  \$57  \$58  \$59  \$60  \$61  \$62
 *      # 5                                     \$63  \$64  \$65  \$66  \$67  \$68
 *      # 6                                           \$69  \$70  \$71  \$72  \$73
 *      # 7                                                 \$74  \$75  \$76  \$77
 *      # 8                                                       \$78  \$79  \$80
 *      # 9                                                             \$81  \$82
 *      # 10                                                                  \$83
 *      # 11
 *
 * Interaction energy columns: Chain-Pentamer
 * \$84 \$85 \$86 \$87 \$88 \$89 \$90 \$91 \$92 \$93 \$94 \$95
 *
 * Interaction energy columns: Chain-Chain
 * \$96
 *
 * Force columns (x,y,z): Pentamer-Pentamer
 *      #            1                 2                  3                  4                  5                  6                  7                  8                  9                 10                 11
 *      # 0   \$97+\$97+\$98  \$100+\$100+\$101  \$103+\$103+\$104  \$106+\$106+\$107  \$109+\$109+\$110  \$112+\$112+\$113  \$115+\$115+\$116  \$118+\$118+\$119  \$121+\$121+\$122  \$124+\$124+\$125  \$127+\$127+\$128
 *      # 1                   \$130+\$130+\$131  \$133+\$133+\$134  \$136+\$136+\$137  \$139+\$139+\$140  \$142+\$142+\$143  \$145+\$145+\$146  \$148+\$148+\$149  \$151+\$151+\$152  \$154+\$154+\$155  \$157+\$157+\$158
 *      # 2                                      \$160+\$160+\$161  \$162+\$163+\$164  \$165+\$166+\$167  \$168+\$169+\$170  \$171+\$172+\$173  \$174+\$175+\$176  \$177+\$178+\$179  \$180+\$181+\$182  \$184+\$184+\$185
 *      # 3                                                         \$186+\$187+\$188  \$189+\$190+\$191  \$192+\$193+\$194  \$195+\$196+\$197  \$198+\$199+\$200  \$201+\$202+\$203  \$204+\$205+\$206  \$208+\$208+\$209
 *      # 4                                                                            \$210+\$211+\$212  \$213+\$214+\$215  \$216+\$217+\$218  \$219+\$220+\$221  \$222+\$223+\$224  \$225+\$226+\$227  \$229+\$229+\$230
 *      # 5                                                                                               \$231+\$232+\$233  \$234+\$235+\$236  \$237+\$238+\$239  \$240+\$241+\$242  \$243+\$244+\$245  \$247+\$247+\$248
 *      # 6                                                                                                                  \$249+\$250+\$251  \$252+\$253+\$254  \$255+\$256+\$257  \$258+\$259+\$260  \$262+\$262+\$263
 *      # 7                                                                                                                                     \$264+\$265+\$266  \$267+\$268+\$269  \$270+\$271+\$272  \$274+\$274+\$275
 *      # 8                                                                                                                                                        \$276+\$277+\$278  \$279+\$280+\$281  \$283+\$283+\$284
 *      # 9                                                                                                                                                                           \$285+\$286+\$287  \$289+\$289+\$290
 *      # 10                                                                                                                                                                                             \$292+\$293+\$294
 *      # 11
 *
 * Force columns (x,y,z): Chain-Pentamer
 * Chain - Penta_0: \$295 \$296 \$297
 * Chain - Penta_1: \$298 \$299 \$300
 * Chain - Penta_2: \$301 \$302 \$303
 * Chain - Penta_3: \$304 \$305 \$306
 * Chain - Penta_4: \$307 \$308 \$309
 * Chain - Penta_5: \$310 \$311 \$312
 * Chain - Penta_6: \$313 \$314 \$315
 * Chain - Penta_7: \$316 \$317 \$318
 * Chain - Penta_8: \$319 \$320 \$321
 * Chain - Penta_9: \$322 \$323 \$324
 * Chain - Penta_10: \$325 \$326 \$327
 * Chain - Penta_11: \$328 \$329 \$330
 *
 * Force columns (x,y,z): Chain-Chain
 * \$331 \$332 \$333
 *
 *  half_force column:
 *  half-half(x,y,z)      half_A-chain(x,y,z)    half_B-chain(x,y,z)
 *  \$334+\$335+\$336     \$337+\$338+\$339      \$340+\$341+\$342
 *  \$343+\$344+\$345     \$346+\$347+\$348      \$349+\$350+\$351
 *  \$352+\$353+\$354     \$355+\$356+\$357      \$358+\$359+\$360
 *  \$361+\$362+\$363     \$364+\$365+\$366      \$367+\$368+\$369
 *  \$370+\$371+\$372     \$373+\$374+\$375      \$376+\$377+\$378
 *  \$379+\$380+\$381     \$382+\$383+\$384      \$385+\$386+\$387
 */
class LammpsCapsidParam {
public:
    LammpsCapsidParam(){}
    bool set = false;
    array<double, 400> param;
    double pressure = 1e5;
    int step;

    //
    // 12 pentamers, Each 11 potential partners, order set by lammps output setting
    //
    const int interaction[12][12] = {
        {0, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28},
        {18, 0, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38},
        {19, 29, 0, 39, 40, 41, 42, 43, 44, 45, 46, 47},
        {20, 30, 39, 0, 48, 49, 50, 51, 52, 53, 54, 55},
        {21, 31, 40, 48, 0, 56, 57, 58, 59, 60, 61, 62},

        {22, 32, 41, 49, 56, 0, 63, 64, 65, 66, 67, 68},
        {23, 33, 42, 50, 57, 63, 0, 69, 70, 71, 72, 73},
        {24, 34, 43, 51, 58, 64, 69, 0, 74, 75, 76, 77},
        {25, 35, 44, 52, 59, 65, 70, 74, 0, 78, 79, 80},
        {26, 36, 45, 53, 60, 66, 71, 75, 78, 0, 81, 82},

        {27, 37, 46, 54, 61, 67, 72, 76, 79, 81, 0, 83},
        {28, 38, 47, 55, 62, 68, 73, 77, 80, 82, 83, 0},
    };

    const int force[12][12] = {
        {  0,  97, 100, 103, 106, 109, 112, 115, 118, 121, 124, 127},
        { 97,   0, 130, 133, 136, 139, 142, 145, 148, 151, 154, 157},
        {100, 130,   0, 160, 163, 166, 169, 172, 175, 178, 181, 184},
        {103, 133, 160,   0, 187, 190, 193, 196, 199, 202, 205, 208},
        {106, 136, 163, 187,   0, 211, 214, 217, 220, 223, 226, 229},

        {109, 139, 166, 190, 211,   0, 232, 235, 238, 241, 244, 247},
        {112, 142, 169, 193, 214, 232,   0, 250, 253, 256, 259, 262},
        {115, 145, 172, 196, 217, 235, 250,   0, 265, 268, 271, 274},
        {118, 148, 175, 199, 220, 238, 253, 265,   0, 277, 280, 283},
        {121, 151, 178, 202, 223, 241, 256, 268, 277,   0, 286, 289},

        {124, 154, 181, 205, 226, 244, 259, 271, 280, 286,   0, 292},
        {127, 157, 184, 208, 229, 247, 262, 274, 283, 289, 292,   0},
    };



    double get_pentamer_pentamer_energy()
    {
        double sum=0.0;
        for(int i=18; i<=83; ++i)
        {
            sum += param[i];
        }
        return sum;
    }

    double get_pentamer_pentamer_force()
    {
        double sum=0.0;
        for(int i=97; i<=294; i+=3)
        {
            sum += sqrt(param[i]*param[i] + param[i+1]*param[i+1] + param[i+2]*param[i+2]);
        }
        return sum;
    }

    /**
     * @brief get_pentamer_pentamer_energy
     * @param i - from 0 to 11
     * @return
     */
    double get_pentamer_pentamer_energy(int i)
    {
        double sum=0.0;
        for(int j=0; j<12; ++j)
        {
            sum += param[ interaction[i][j] ];
        }
        return sum;
    }

    double get_pentamer_pentamer_force(int i)
    {
        double sum=0.0, x,y,z;
        for(int j=0; j<12; ++j)
        {
            x=param[ force[i][j] ];
            y=param[ force[i][j] +1 ];
            z=param[ force[i][j] +2 ];
            sum += sqrt(x*x + y*y + z*z);
        }
        return sum;
    }



    double get_pentamer_chain_energy()
    {
        double energy_sum=0.0; // energy sum in kT
        for(int i = 84; i<=95; ++i) // energy of genome exerted on 12 pentamers
            energy_sum += param[i];
        return energy_sum;
    }

    double get_pentamer_chain_force()
    {
        double force_sum=0.0; // energy sum in kT
        for(int i = 295; i<=328; i+=3) // energy of genome exerted on 12 pentamers
            force_sum += sqrt(param[i]*param[i] + param[i+1]*param[i+1] + param[i+2]*param[i+2]);
        return force_sum;
    }

    /**
     * @brief get_pentamer_chain_energy
     * @param i - from 0 to 11
     * @return
     */
    double get_pentamer_chain_energy(int i)
    {
        return param[84+i];
    }

    double get_pentamer_chain_force(int i)
    {
        return sqrt(param[295+(3*i)]*param[295+(3*i)] + param[295+(3*i)+1]*param[295+(3*i)+1] + param[295+(3*i)+2]*param[295+(3*i)+2]);;
    }



    double get_chain_wall_F()
    {
        return ( sqrt( param[388]*param[388] + param[389]*param[389] + param[390]*param[390]) +
               sqrt( param[391]*param[391] + param[392]*param[392] + param[393]*param[393] ) ) / 2.0;
    }

    double get_chain_chain_energy()
    {
        return param[96];
    }

    double get_chain_chain_force()
    {
        return sqrt(param[331]*param[331] + param[332]*param[332] + param[333]*param[333]);
    }

    /**
     * @brief get_halfs
     * @return
     *  \$334+\$335+\$336     \$337+\$338+\$339      \$340+\$341+\$342
     *  \$343+\$344+\$345     \$346+\$347+\$348      \$349+\$350+\$351
     *  \$352+\$353+\$354     \$355+\$356+\$357      \$358+\$359+\$360
     *  \$361+\$362+\$363     \$364+\$365+\$366      \$367+\$368+\$369
     *  \$370+\$371+\$372     \$373+\$374+\$375      \$376+\$377+\$378
     *  \$379+\$380+\$381     \$382+\$383+\$384      \$385+\$386+\$387
     */
    string get_halfs()
    {
        stringstream ss;
        for(int i=334; i<388; i+=9)
        {
            ss << sqrt( param[i]*param[i] + param[i+1]*param[i+1] + param[i+2]*param[i+2] ) << " ";
            ss << 0.5*sqrt( param[i+3]*param[i+3] + param[i+4]*param[i+4] + param[i+5]*param[i+5] ) + 0.5*sqrt( param[i+6]*param[i+6] + param[i+7]*param[i+7] + param[i+8]*param[i+8] ) << " ";
        }
        return ss.str();
    }

    /**
     * @brief get_half_genome_F
     * @param i : 0, 1, 2, 3, 4, 5
     * @return
     */
    double get_half_half_F(int i)
    {
        i*=9;
        i+=334;
        return sqrt( param[i]*param[i] + param[i+1]*param[i+1] + param[i+2]*param[i+2] );
    }

    double get_half_genome_F(int i)
    {
        i*=9;
        i+=334;
        return 0.5*sqrt( param[i+3]*param[i+3] + param[i+4]*param[i+4] + param[i+5]*param[i+5] ) + 0.5*sqrt( param[i+6]*param[i+6] + param[i+7]*param[i+7] + param[i+8]*param[i+8] );
    }


    void get_interaction_partners(int number, int& a, int& b)
    {
        for(int i=0; i<12; ++i) {
            for(int j=0; j<12; ++j) {
                if(interaction[i][j] == number) {
                    a=i;
                    b=j;
                    return;
                }
            }
        }
    }

    array< array<int, 12>, 12 > getConnectionMap(double interaction_min)
    {
        array< array<int, 12>, 12 > map = {{0}};
        vector<Edge> edge = getEdges(interaction_min);

        for(unsigned int i=0; i<edge.size(); ++i ) {
            map[edge[i].pentamerA][edge[i].pentamerB] = 1;
        }

        return map;
    }

    vector<Edge> getEdges(double interaction_min)
    {
        vector<Edge> edge;
        Edge temp;
        int a=-1,b=-1;

        int count=0;
        for(int i=18;i<=83;i++)
        {
            if(param[i] < interaction_min) { // energy in negative
                get_interaction_partners(i,a,b);
                temp.N=count;
                temp.pentamerA=a;
                temp.pentamerB=b;
                edge.push_back(temp);
                ++count;
            }
        }

        return edge;
    }

    PentamerMap getClusters(double interaction_min)
    {
        PentamerMap cluster;
        vector<Edge> edge = getEdges(interaction_min);

        cluster.setAllValues(-1);

        vector<int> pentas;
        for(int i=0; i<12; i++)
            pentas.push_back(i);

        int num;
        for(int i=0;i<12;i++)
        { // for possible cluster
            // initialize first cluster
            if(!edge.empty()) {
                cluster.interactions[i][0]=edge[0].pentamerA;
                cluster.interactions[i][1]=edge[0].pentamerB;
                edge.erase(edge.begin());
                num=2;
                erase(cluster.interactions[i][0], pentas);
                erase(cluster.interactions[i][1], pentas);
            } else {
                if(!pentas.empty()) {
                    cluster.interactions[i][0]= pentas.front();
                    pentas.erase(pentas.begin());
                }
            }

            for(int j=0;j<12;j++) // for possible pentamers in cluster i
            {
                auto it = edge.begin();
                while( it != edge.end() ) // for all edges,  // edge, [0] and [1] are nodes
                {
                    if( match(it->pentamerA, cluster.interactions[i]) && match(it->pentamerB, cluster.interactions[i]) ) {
                        erase(it->pentamerA, pentas);
                        erase(it->pentamerB, pentas);
                        edge.erase(it);
                        it = edge.begin();
                        continue;
                    }
                    if(match(it->pentamerA, cluster.interactions[i]) && !match(it->pentamerB, cluster.interactions[i]))
                    {
                        cluster.interactions[i][num]=it->pentamerB;
                        erase(it->pentamerB, pentas);
                        num++;
                        edge.erase(it);
                        it = edge.begin();
                        continue;
                    }
                    if(match(it->pentamerB, cluster.interactions[i]) && !match(it->pentamerA, cluster.interactions[i]))
                    {
                        cluster.interactions[i][num]=it->pentamerA;
                        erase(it->pentamerA, pentas);
                        num++;
                        edge.erase(it);
                        it = edge.begin();
                        continue;
                    }
                    it++;
                }
            }
        }
        return cluster;
    }

    bool match(int graph, array<int, 12>& cluster)
    {
        for(int i=0; i<12; ++i) {
            if(cluster[i] == graph)
                return true;
        }
        return false;
    }

    void erase(int node, vector<int>& pentas)
    {
        for(unsigned int k=0; k<pentas.size(); ++k)
        {
            if(pentas[k] == node)
            {
                pentas.erase( pentas.begin() + k );
                break;
            }
        }

    }
};

#endif // LAMMPSCAPSIDPARAM_H
