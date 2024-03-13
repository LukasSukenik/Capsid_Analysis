#ifndef XTCANALYSIS_H
#define XTCANALYSIS_H

#include <string>
#include <limits>
#include <cmath>
#include <iostream>
#include <vector>

#include "welford.h"
#include "histogram.h"

#include "xdrfile-1.1.4/include/xdrfile.h"
#include "xdrfile-1.1.4/include/xdrfile_xtc.h"

using namespace std;




class XTCAnalysis
{
private:
/*
 3-fold axis:
 name 4
 serial 170 811 2706 177 804 3030 205 1445 1755 212 1438 2713 234 1762 3023 487 1128 3340 494 1121 3664 522 2079 2389 529 2072 3347 2474 3114 3742 551 2396 3657 846
 2107 2431 839 2114 2735 874 2424 3058 1163 1473 1797 1156 1480 3369 1191 1790 3692 1523 2769 3397 1846 3080 3720 2163 2791 3431

 2-fold axis:
 name 3
 serial 173 807 208 1441 237 1758 258 2709 279 3026 490 1124 525 2075 554 2392 575 3343 596 3660 842 2110 870 2427 892 2738 913 3054 1159 1476 1187
 1793 1209 3372 1230 3688 1498 1815 1526 2765 1547 3400 1842 3083 1864 3716 2132 2449 2159 2794 2181 3427 2477 3110 2498 3745 2815 3449 3132 3766

 name 2
 serial 174 808 209 1442 238 1759 259 2710 280 3027 491 1125 526 2076 555 2393 576 3344 597 3661 843 2111 871 2428 893 2739 914 3055 1160 1477
 1188 1794 1210 3373 1231 3689 1499 1816 1527 2766 1548 3401 1843 3084 1865 3717 2133 2450 2160 2795 2182 3428 2478 3111 2499 3746 2816 3450 3133 3767
*/

    int fold3[60] = {170, 811, 2706, 177, 804, 3030, 205, 1445, 1755, 212, 1438, 2713, 234, 1762, 3023, 487, 1128, 3340, 494, 1121, 3664, 522, 2079, 2389, 529, 2072, 3347, 2474, 3114, 3742, 551, 2396, 3657, 846, 2107, 2431, 839, 2114, 2735, 874, 2424, 3058, 1163, 1473, 1797, 1156, 1480, 3369, 1191, 1790, 3692, 1523, 2813, 3397, 1846, 3080, 3720, 2163, 2791, 3431};
    int fold2_1[60] = {173, 807, 208, 1441, 237, 1758, 258, 2709, 279, 3026, 490, 1124, 525, 2075, 554, 2392, 575, 3343, 596, 3660, 842, 2110, 870, 2427, 892, 2738, 913, 3054, 1159, 1476, 1187, 1793, 1209, 3372, 1230, 3688, 1498, 1815, 1526, 2765, 1547, 3400, 1842, 3083, 1864, 3716, 2132, 2449, 2159, 2794, 2181, 3427, 2477, 3110, 2498, 3745, 2815, 3449, 3132, 3766};
    int fold2_2[60] = {174, 808, 209, 1442, 238, 1759, 259, 2710, 280, 3027, 491, 1125, 526, 2076, 555, 2393, 576, 3344, 597, 3661, 843, 2111, 871, 2428, 893, 2739, 914, 3055, 1160, 1477, 1188, 1794, 1210, 3373, 1231, 3689, 1499, 1816, 1527, 2766, 1548, 3401, 1843, 3084, 1865, 3717, 2133, 2450, 2160, 2795, 2182, 3428, 2478, 3111, 2499, 3746, 2816, 3450, 3133, 3767};

    int rim[420] = {811, 890, 2706, 894, 893, 2738, 892, 2739, 2740, 891, 2741, 2707, 2708, 2709, 2114, 2157, 2113, 2797, 2796, 841, 842, 839, 840, 2737, 2736, 2735, 895, 2111, 2112, 2159, 2158, 2795, 2794, 845, 844, 2110, 843, 870, 871, 872, 2764, 2765, 2160, 2161, 2162, 2792, 2791, 2793, 2818, 2817, 2431, 2107, 2108, 2109, 2816, 2815, 2814, 869, 2430, 2429, 868, 846, 2768, 2769, 2813, 2428, 2766, 2767, 2426, 2427, 1527, 1526, 2180, 2179, 2163, 3431, 3447, 2134, 3448, 3449, 3450, 2447, 2448, 2135, 3451, 3452, 3397, 1525, 1523, 1524, 2181, 3429, 3430, 2133, 2132, 2450, 2449, 3399, 3398, 1550, 1549, 170, 808, 809, 810, 171, 261, 260, 2710, 807, 806, 172, 173, 259, 258, 805, 174, 175, 176, 2711, 2712, 2713, 2763, 914, 915, 257, 256, 212, 916, 3053, 804, 177, 3052, 3030, 3029, 277, 278, 279, 911, 874, 1529, 1528, 1438, 912, 3056, 913, 211, 1439, 1440, 3055, 3054, 210, 209, 3028, 3027, 280, 281, 2424, 2425, 873, 3108, 3058, 3057, 1442, 1441, 208, 207, 1443, 3026, 238, 239, 240, 205, 206, 3025, 3024, 282, 234, 235, 236, 237, 2480, 2479, 3110, 3109, 1444, 1756, 1755, 1496, 1445, 3023, 3085, 3086, 1840, 1762, 1761, 1760, 1759, 1758, 1757, 2478, 3111, 3112, 1497, 1498, 1818, 1817, 3084, 1842, 1841, 2183, 2182, 3427, 3428, 2131, 2452, 2130, 2451, 3400, 3401, 1548, 2077, 2076, 2075, 2074, 2184, 2072, 2073, 3426, 3425, 3347, 2389, 2079, 2078, 2390, 2391, 3403, 3402, 1547, 1546, 1545, 526, 527, 528, 529, 573, 3346, 3345, 522, 523, 524, 525, 2392, 557, 556, 2393, 2394, 3370, 3369, 3371, 2498, 1156, 1479, 1480, 574, 575, 576, 3344, 3343, 555, 3342, 554, 3374, 3373, 3372, 2395, 2396, 553, 552, 1210, 1211, 2496, 2497, 3748, 1157, 1212, 577, 578, 487, 3341, 3340, 488, 3375, 1207, 1208, 551, 599, 598, 597, 1209, 3657, 3658, 3659, 3660, 489, 490, 1128, 1127, 492, 491, 1125, 1126, 596, 595, 494, 594, 1122, 493, 1123, 1124, 2476, 2477, 3113, 1499, 1816, 3082, 3083, 1843, 1844, 2500, 2475, 2474, 2501, 3114, 3130, 3742, 3131, 3132, 1500, 1501, 3133, 3134, 3135, 1814, 1815, 3080, 3081, 1846, 1862, 1845, 2499, 1478, 3745, 3744, 1477, 1476, 3769, 3743, 1475, 1474, 3768, 3767, 1795, 1813, 1797, 1796, 1163, 1473, 3766, 3764, 3765, 3719, 3720, 1864, 1863, 3747, 3746, 1158, 1159, 1160, 1161, 1162, 1187, 1185, 1186, 3718, 3717, 3716, 1866, 1865, 1867, 1792, 1793, 1794, 3691, 3692, 3714, 1228, 1191, 1190, 1189, 1188, 3715, 1791, 1790, 3661, 3662, 3663, 3664, 3686, 1233, 3687, 3688, 1232, 1231, 3689, 3690, 1230, 1229, 1121};

    int bead_i_general[240] = {287, 288, 288, 289, 286, 286, 285, 285, 297, 296, 296, 295, 294, 294, 293, 293, 305, 304, 304, 303, 302, 302, 301, 301, 309, 310, 310, 311, 308, 308, 307, 307, 317, 315, 316, 316, 313, 313, 314, 314, 604, 605, 605, 606, 602, 602, 603, 603, 614, 613, 613, 612, 611, 611, 610, 610, 622, 620, 621, 621, 619, 619, 618, 618, 626, 627, 627, 628, 625, 625, 624, 624, 632, 633, 633, 634, 630, 630, 631, 631, 929, 930, 930, 931, 927, 927, 928, 928, 937, 937, 938, 936, 934, 934, 935, 935, 943, 944, 944, 945, 941, 941, 942, 942, 951, 950, 950, 949, 948, 948, 947, 947, 1246, 1248, 1247, 1247, 1245, 1245, 1244, 1244, 1253, 1254, 1254, 1255, 1251, 1251, 1252, 1252, 1262, 1261, 1261, 1260, 1259, 1259, 1258, 1258, 1266, 1267, 1267, 1268, 1264, 1264, 1265, 1265, 1571, 1571, 1570, 1572, 1568, 1568, 1569, 1569, 1579, 1577, 1578, 1578, 1576, 1576, 1575, 1575, 1585, 1584, 1584, 1583, 1581, 1581, 1582, 1582, 1893, 1894, 1894, 1895, 1891, 1891, 1892, 1892, 1901, 1901, 1900, 1902, 1899, 1899, 1898, 1898, 2206, 2205, 2205, 2204, 2203, 2203, 2202, 2202, 2212, 2210, 2211, 2211, 2208, 2208, 2209, 2209, 2217, 2219, 2218, 2218, 2215, 2215, 2216, 2216, 2530, 2529, 2529, 2528, 2527, 2527, 2526, 2526, 2536, 2535, 2535, 2534, 2532, 2532, 2533, 2533, 2851, 2852, 2852, 2853, 2849, 2849, 2850, 2850, 3170, 3168, 3169, 3169, 3167, 3167, 3166, 3166};
    int bead_j_general[240] = {920, 920, 919, 919, 922, 921, 923, 922, 1553, 1553, 1554, 1554, 1555, 1556, 1556, 1557, 1870, 1870, 1871, 1871, 1872, 1873, 1873, 1874, 2822, 2822, 2821, 2821, 2823, 2824, 2825, 2824, 3138, 3139, 3139, 3138, 3142, 3141, 3140, 3141, 1237, 1237, 1236, 1236, 1240, 1239, 1239, 1238, 2187, 2187, 2188, 2188, 2189, 2190, 2190, 2191, 2504, 2505, 2504, 2505, 2506, 2507, 2507, 2508, 3456, 3455, 3456, 3455, 3457, 3458, 3459, 3458, 3773, 3773, 3772, 3772, 3775, 3776, 3775, 3774, 2196, 2196, 2195, 2195, 2199, 2198, 2198, 2197, 2513, 2512, 2512, 2513, 2516, 2515, 2515, 2514, 2830, 2829, 2830, 2829, 2833, 2832, 2831, 2832, 3145, 3145, 3146, 3146, 3147, 3148, 3148, 3149, 1562, 1561, 1562, 1561, 1563, 1564, 1565, 1564, 1879, 1879, 1878, 1878, 1882, 1881, 1881, 1880, 3463, 3463, 3464, 3464, 3465, 3466, 3466, 3467, 3780, 3780, 3779, 3779, 3783, 3782, 3781, 3782, 1886, 1885, 1886, 1885, 1889, 1888, 1887, 1888, 2835, 2836, 2835, 2836, 2838, 2837, 2838, 2839, 3470, 3471, 3470, 3471, 3473, 3474, 3472, 3473, 3154, 3154, 3153, 3153, 3157, 3156, 3156, 3155, 3786, 3787, 3787, 3786, 3789, 3788, 3789, 3790, 2519, 2520, 2519, 2520, 2521, 2522, 2522, 2523, 2843, 2844, 2843, 2844, 2847, 2846, 2845, 2846, 3477, 3476, 3477, 3476, 3480, 3479, 3479, 3478, 3159, 3159, 3160, 3160, 3161, 3162, 3163, 3162, 3794, 3794, 3795, 3795, 3797, 3798, 3796, 3797, 3484, 3484, 3483, 3483, 3487, 3486, 3485, 3486, 3800, 3801, 3800, 3801, 3802, 3803, 3803, 3804};
    const double center=750.0;

    char fileName[64];
    int natoms;

    int step;
    int release_step = numeric_limits<int>::max();
    float time;
    matrix box;
    float prec;

    double dist_3fold = 999;
    double dist_2fold = 999;

public:
    double fold2_dist = 0.0;
    double fold3_dist = 0.0;

    bool isfold2 = false;
    bool isfold3 =false;
    bool ismultiple_slow = false;

    int rim_time;

    Welford_Algo aver[240];

    XTCAnalysis()
    {
        // minus 1 for all fold indexes
        for(int i=0; i<60; ++i)
        {
            --fold2_1[i];
            --fold2_2[i];
            --fold3[i];
        }
        for(int i=0; i<240; ++i)
        {
            --bead_i_general[i];
            --bead_j_general[i];
        }
        for(int i=0; i<420; ++i)
        {
            --rim[i];
        }
    }

    double get_3fold_center_dist(rvec k[], int natoms)
    {
        double dist = 0.0;

        for(int j=0; j<60; ++j)
        {
            dist += (k[fold3[j]][0] - center)*(k[fold3[j]][0] - center)
                  + (k[fold3[j]][1] - center)*(k[fold3[j]][1] - center)
                  + (k[fold3[j]][2] - center)*(k[fold3[j]][2] - center);
        }

        return sqrt( dist/60.0 );
    }

    double calc_3fold_dist(rvec k[], int natoms)
    {
        double dist_3fold1, dist_3fold2, dist_3fold3;
        double min = 999999999.0;

        for(int i=3804; i<7509; ++i)
        {
            for(int j=0; j<60; j+=3)
            {
                dist_3fold1 = (k[i][0]-k[fold3[j]][0])*(k[i][0]-k[fold3[j]][0])
                            + (k[i][1]-k[fold3[j]][1])*(k[i][1]-k[fold3[j]][1])
                            + (k[i][2]-k[fold3[j]][2])*(k[i][2]-k[fold3[j]][2]);
                dist_3fold2 = (k[i][0]-k[fold3[j+1]][0])*(k[i][0]-k[fold3[j+1]][0])
                            + (k[i][1]-k[fold3[j+1]][1])*(k[i][1]-k[fold3[j+1]][1])
                            + (k[i][2]-k[fold3[j+1]][2])*(k[i][2]-k[fold3[j+1]][2]);
                dist_3fold3 = (k[i][0]-k[fold3[j+2]][0])*(k[i][0]-k[fold3[j+2]][0])
                            + (k[i][1]-k[fold3[j+2]][1])*(k[i][1]-k[fold3[j+2]][1])
                            + (k[i][2]-k[fold3[j+2]][2])*(k[i][2]-k[fold3[j+2]][2]);

                if(min > dist_3fold1) { min = dist_3fold1; }
                if(min > dist_3fold2) { min = dist_3fold2; }
                if(min > dist_3fold3) { min = dist_3fold3; }
            }
        }
        return sqrt(min);
    }

    double calc_2fold_dist(rvec k[], int natoms)
    {
        double dist_2fold1, dist_2fold2, dist_2fold3, dist_2fold4;
        double min = 999999999.0;

        for(int i=3804; i<7509; ++i)
        {
            for(int j=0; j<60; j+=2)
            {
                dist_2fold1 = (k[i][0]-k[fold2_1[j]][0])*(k[i][0]-k[fold2_1[j]][0])
                            + (k[i][1]-k[fold2_1[j]][1])*(k[i][1]-k[fold2_1[j]][1])
                            + (k[i][2]-k[fold2_1[j]][2])*(k[i][2]-k[fold2_1[j]][2]);
                dist_2fold2 = (k[i][0]-k[fold2_1[j+1]][0])*(k[i][0]-k[fold2_1[j+1]][0])
                            + (k[i][1]-k[fold2_1[j+1]][1])*(k[i][1]-k[fold2_1[j+1]][1])
                            + (k[i][2]-k[fold2_1[j+1]][2])*(k[i][2]-k[fold2_1[j+1]][2]);

                dist_2fold3 = (k[i][0]-k[fold2_2[j]][0])*(k[i][0]-k[fold2_2[j]][0])
                            + (k[i][1]-k[fold2_2[j]][1])*(k[i][1]-k[fold2_2[j]][1])
                            + (k[i][2]-k[fold2_2[j]][2])*(k[i][2]-k[fold2_2[j]][2]);
                dist_2fold4 = (k[i][0]-k[fold2_2[j+1]][0])*(k[i][0]-k[fold2_2[j+1]][0])
                            + (k[i][1]-k[fold2_2[j+1]][1])*(k[i][1]-k[fold2_2[j+1]][1])
                            + (k[i][2]-k[fold2_2[j+1]][2])*(k[i][2]-k[fold2_2[j+1]][2]);

                if(min > dist_2fold1) { min = dist_2fold1; }
                if(min > dist_2fold2) { min = dist_2fold2; }
                if(min > dist_2fold3) { min = dist_2fold3; }
                if(min > dist_2fold4) { min = dist_2fold4; }
            }
        }

        return sqrt(min);
    }

    double get_2fold_center_dist(rvec k[], int natoms)
    {
        double dist = 0.0;

        for(int j=0; j<60; ++j)
        {
            dist += (k[fold2_1[j]][0] - center)*(k[fold2_1[j]][0] - center)
                  + (k[fold2_1[j]][1] - center)*(k[fold2_1[j]][1] - center)
                  + (k[fold2_1[j]][2] - center)*(k[fold2_1[j]][2] - center);
            dist += (k[fold2_2[j]][0] - center)*(k[fold2_2[j]][0] - center)
                  + (k[fold2_2[j]][1] - center)*(k[fold2_2[j]][1] - center)
                  + (k[fold2_2[j]][2] - center)*(k[fold2_2[j]][2] - center);
        }

        return sqrt( dist/120.0 );
    }

    /**
     * @brief isReleased
     * @param k
     * @param natoms
     * @return
     */
    bool isReleased(double max, rvec k[], int natoms)
    {
        return ( (max > get_2fold_center_dist(k, natoms)) || (max > get_3fold_center_dist(k, natoms)) );
    }

    //
    // pair_coeff 2 3 lj/cut 1.0 18.7373 21.0319
    // pair_coeff 6 7 lj/cut 1.0 16.3931 18.4006
    // pair_coeff 9 10 lj/cut 1.0 17.5652 19.7162
    //
    // pair_coeff 2 3 cosatt ${coeff_out} 21.0319 ${intr_dist_out}
    // pair_coeff 6 7 cosatt ${coeff_in} 18.4006 ${intr_dist_in}
    // pair_coeff 9 10 cosatt ${coeff_mid} 19.7162 ${intr_dist_mid}
    //
    double potential(double dist, double epsilon, double sigma, double rm)
    {
        if(dist > rm+sigma)
            return 0.0;
        if(dist > rm)
        {
            const double pih = std::atan(1.0)*2;
            return -epsilon * cos( pih * (dist - rm) / sigma ) * cos( pih * (dist - rm) / sigma );
        }
        else
        {
            // https://en.wikipedia.org/wiki/Lennard-Jones_potential
            // -epsilon from cos^2 potential, eps for LJ is 1.0
            double ssigma = rm / pow(2.0, 1.0/6.0);
            return 4*epsilon*( pow(ssigma/dist, 12) - pow(ssigma/dist,6) );
        }
    }

    double force(double dist, double epsilon, double w_c, double r_c)
    {
        double pi=3.14159265359;
        double arg= pi/2.0 * (dist - r_c) / w_c;
        double sinn=sin(arg);
        double coss=cos(arg);
        double f = epsilon*pi*coss*sinn/(w_c);
        double e = -epsilon*coss*coss;

        if( dist < r_c )
        {
            return 0.0;
        }
        if( dist > r_c + w_c )
        {
            return 0.0;
        }

        return f;
    }

    void fluctuationAnalyse(string inName, int start, int stop=numeric_limits<int>::max(), double force_av=0.0 )
    {
        Welford_Algo force_aver;
        double rm = 19.7162;
        double sigma = 0.0;
        double eps = 0.0;

        ifstream input;
        input.open("ff_values");
        input >> sigma >> eps;
        input.close();

        double sumE=0.0;
        double sumF=0.0;
        Histogram histo;
        histo.init(sigma, 100, rm, 20);

        int status=exdrOK;
        double dist;

        //
        // Add end symbol to xtc filename
        //
        int i;
        for (i = 0; inName[i] != '\0'; ++i) {
          fileName[i] = inName[i];
        }
        fileName[i] = '\0';

        //
        // read xtc file
        //
        status = read_xtc_natoms(fileName, &natoms);
        if (status == exdrOK)
        {
            XDRFILE* xfp = xdrfile_open(fileName, "r");
            if (xfp != NULL)
            {
                rvec k[natoms];
                status = read_xtc(xfp, natoms, &step, &time, box, k, &prec);

                //
                // loop over frames
                //
                while(status == exdrOK && step < stop)
                {
                    status = read_xtc(xfp, natoms, &step, &time, box, k, &prec);

                    if(step > start)
                    {
                        //
                        // analyse fluctuation of bead-bead distance of interaction pairs
                        //
                        sumE = 0.0;
                        sumF = 0.0;
                        for(int i=0; i<240; ++i)
                        {
                            dist = pow( ( k[bead_i_general[i]][0] - k[bead_j_general[i]][0] ), 2)
                                    + pow( ( k[bead_i_general[i]][1] - k[bead_j_general[i]][1] ), 2)
                                    + pow( ( k[bead_i_general[i]][2] - k[bead_j_general[i]][2] ), 2);
                            sumE += potential(sqrt(dist), eps, sigma, rm);
                            sumF += force(sqrt(dist), eps, sigma, rm);

                            //
                            // Add value to Histogram
                            //
                            histo.add( sqrt(dist) - rm );
                        }
                        //cout << step << " " << sumF << endl;
                        force_aver.push( sumF );
                        //cout << step << " " << sumF/240.0 << "  " << endl;
                        //histo << step << " " << average.getAverage() << " " << sumE << "\n";
                    }
                }

                //
                // Print out histogram
                //
                histo.print();
                xdrfile_close(xfp);
            }
            else
            {
                cout << "File not opened:" << fileName << endl;
            }  
        }
        else
        {
            fprintf(stderr, "read_xtc_natoms failure; return code %d", status);
        }
    }

    int slowReleAnalyse(string inName, int stop=-1)
    {
        int status=exdrOK;
        double dist,max=0;

        int first_step;

        Welford_Algo fold2aver;
        Welford_Algo fold3aver;

        int i;
        for (i = 0; inName[i] != '\0'; ++i) {
          fileName[i] = inName[i];
        }
        fileName[i] = '\0';

        status = read_xtc_natoms(fileName, &natoms);
        if (status == exdrOK)
        {
            XDRFILE* xfp = xdrfile_open(fileName, "r");
            if (xfp != NULL)
            {
                rvec k[natoms];
                status = read_xtc(xfp, natoms, &first_step, &time, box, k, &prec);

                while(status == exdrOK && ( stop == -1 || step < first_step+stop) ) // Analyze frame
                {
                    status = read_xtc(xfp, natoms, &step, &time, box, k, &prec);

                    for(int i=3804; i<7509; ++i)
                    {
                        dist = (k[i][0]-center)*(k[i][0]-center) + (k[i][1]-center)*(k[i][1]-center) + (k[i][2]-center)+(k[i][2]-center);
                        if(max < dist)
                        {
                            max=dist;
                        }
                    }
                    max = sqrt(max);
                    if( max > 150.0)
                    {
                        if(rim_time > step) {
                            rim_time = step;
                        }
                    }

                    if( isReleased( max, k, natoms ) )
                    {
                        fold2_dist = calc_2fold_dist(k, natoms);
                        fold3_dist = calc_3fold_dist(k, natoms);

                        fold2aver.push( fold2_dist );
                        fold3aver.push( fold3_dist );

                        if(release_step > step) {
                            release_step = step;
                        }
                    }
                }
                isfold2 = (fold2aver.getAverage() < 16.5);
                isfold3 = (fold3aver.getAverage() < 16.5);
                ismultiple_slow = (isfold2 && isfold3);

                xdrfile_close(xfp);
                return release_step;
            }
            else
            {
                cout << "File not opened:" << fileName << endl;
            }
        }
        else
        {
            fprintf(stderr, "read_xtc_natoms failure; return code %d", status);
        }
        return -1;
    }
};

#endif // XTCANALYSIS_H
