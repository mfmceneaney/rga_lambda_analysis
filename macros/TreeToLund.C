/**
* Description: Convert a ROOT TTree to Lund format.
* Author: Matthew F. McEneaney
* Date: 16 September 2025
*/

float getCharge(int pid) {
    switch (pid) {

        // Lund particles
        case 91:  // Lund string
            return 0;
        case 92:  // Lund cluster
            return 0;

        /*
        
        DIQUARKS
        (dd)1 1103
        (ud)0 2101
        (ud)1 2103
        (uu)1 2203
        (sd)0 3101
        (sd)1 3103
        (su)0 3201
        (su)1 3203
        (ss)1 3303
        
        */

        // Diquarks // Just set to 0 charge for now
        case 1103: // dd_1 diquark -2/3
            return 0;
        case 2101: // ud_0 diquark 1/3
            return 0;
        case 2103: // ud_1 diquark 1/3
            return 0;
        case 2203: // uu_1 diquark 4/3
            return 0;
        case 3101: // sd_0 diquark -2/3
            return 0;
        case 3103: // sd_1 diquark -2/3
            return 0;
        case 3201: // su_0 diquark 1/3
            return 0;
        case 3203: // su_1 diquark 1/3
            return 0;
        case 3303: // ss_1 diquark -2/3
            return 0;

        // Leptons
        case 11:   // electron
            return -1;
        case -11:  // positron
            return 1;
        case 13:   // muon-
            return -1;
        case -13:  // muon+
            return 1;

        // Bosons
        case 21:   // gluon
            return 0;
        case 22:   // photon
            return 0;

        // Light I=1 Mesons
        case 111:  // pi0
            return 0;
        case -111: // anti-pi0
            return 0;
        case 211:  // pi+
            return 1;
        case -211: // pi-
            return -1;
        case 213:  // rho+
            return 1;
        case -213: // rho-
            return -1;
        case 113:  // rho0
            return 0;
        case -113: // anti-rho0
            return 0;

        // Light I=0 Mesons
        case 221:  // eta
            return 0;
        case 331:  // eta'
            return 0;
        case 223:  // omega
            return 0;
        case 333:  // phi
            return 0;
        
        // Strange Mesons
        case 130:  // K0_L
            return 0;
        case 310:  // K0_S
            return 0;
        case 311:  // K0
            return 0;
        case -311: // anti-K0
            return 0;
        case 321:  // K+
            return 1;
        case -321: // K-
            return -1;
        case 313:  // K*0
            return 0;
        case -313: // anti-K*0
            return 0;
        case 323:  // K*+
            return 1;
        case -323: // K*-
            return -1;

        // Light Baryons
        case 2212: // proton
            return 1;
        case -2212:// anti-proton
            return -1;
        case 2112: // neutron
            return 0;
        case -2112:// anti-neutron
            return 0;
        case 2208: // Delta++
            return 2;
        case -2208:// anti-Delta++
            return -2;
        case 2214: // Delta+
            return 1;
        case -2214:// anti-Delta+
            return -1;
        case 2224: // Delta0
            return 0;
        case -2224:// anti-Delta0
            return 0;
        case 2234: // Delta-
            return -1;
        case -2234:// anti-Delta-
            return 1;

        // Strange Baryons
        case 3122: // Lambda
            return 0;
        case -3122:// anti-Lambda
            return 0;
        case 3222: // Sigma+
            return 1;
        case -3222:// anti-Sigma+
            return -1;
        case 3212: // Sigma0
            return 0;
        case -3212:// anti-Sigma0
            return 0;
        case 3112: // Sigma-
            return -1;
        case -3112:// anti-Sigma-
            return 1;
        case 3322: // Xi0
            return 0;
        case -3322:// anti-Xi0
            return 0;
        case 3312: // Xi-
            return -1;
        case -3312:// anti-Xi-
            return 1;
        case 3334: // Omega-
            return -1;
        case -3334:// anti-Omega-
            return 1;

        // Default
        default:
            return 0; // neutral or unknown
    }
}

/**
CURRENT -> -12, 21, -12, -21, -21, 1, -23, -62, -63, -83, 1, -84, 1, 1, 1, 1
1	-1.00000	-12	11	0	4	-0.00000	-0.00000	11.98460	11.98460	-0.00000	0.00000	0.00000	0.00000
2	0.00000	21	22	1	0	-0.00000	0.00000	-0.03200	-0.03200	-0.00000	0.00000	0.00000	0.00000
3	1.00000	-12	2212	0	5	-0.00000	0.00000	-0.00120	0.93827	0.93827	0.00000	0.00000	0.00000
4	-1.00000	-21	11	1	6	0.00000	-0.00000	12.01660	12.01660	-0.00000	0.00000	0.00000	0.00000
5	0.00000	-21	2	2	6	0.00000	-0.00000	-0.03209	0.03209	-0.00000	0.00000	0.00000	0.00000
6	-1.00000	1	11	4	0	-0.09960	-0.52539	3.34761	3.39006	0.00755	0.00000	0.00000	0.00000
7	0.00000	-23	2	4	8	0.09960	0.52539	9.07344	9.09517	0.32992	0.00000	0.00000	0.00000
8	0.00000	-62	2	7	10	0.09622	0.50758	8.76569	8.78710	0.33001	0.00000	0.00000	0.00000
9	0.00000	-63	2103	2	10	0.00338	0.01781	-0.09789	0.77773	0.77133	0.00000	0.00000	0.00000
10	0.00000	-83	113	8	13	-0.08312	-0.03175	4.54481	4.59112	0.64434	0.00000	0.00000	0.00000
11	1.00000	1	211	8	0	0.24974	0.42440	3.44304	3.48088	0.13969	0.00000	0.00000	0.00000
12	0.00000	-84	2114	8	15	-0.06702	0.13273	0.67995	1.49284	1.32066	0.00000	0.00000	0.00000
13	1.00000	1	211	10	0	0.20396	0.13439	2.38503	2.40156	0.13950	0.00000	0.00000	0.00000
14	-1.00000	1	-211	10	0	-0.28708	-0.16613	2.15978	2.18956	0.13967	0.00000	0.00000	0.00000
15	0.00000	1	2112	12	0	-0.30382	0.26474	0.47601	1.12773	0.93957	0.00000	0.00000	0.00000
16	0.00000	1	111	12	0	0.23680	-0.13200	0.20393	0.36511	0.13498	0.00000	0.00000	0.00000

16	1	1	0.00000	1.00000	11	10.60000	2212	1	0.00000
1	-1.00000	12	11	0	4	-0.00000	0.00000	11.92820	11.92820	-0.00000	0.00000	0.00000	0.00000
2	0.00000	21	22	1	0	0.00000	0.00000	-0.14960	-0.14960	-0.00000	0.00000	0.00000	0.00000
3	1.00000	12	2212	0	5	-0.00000	0.00000	-0.00563	0.93829	0.93827	0.00000	0.00000	0.00000
4	-1.00000	21	11	1	6	-0.00000	0.00000	12.07780	12.07780	-0.00000	0.00000	0.00000	0.00000
5	0.00000	21	1	2	6	-0.00000	0.00000	-0.15141	0.15141	-0.00000	0.00000	0.00000	0.00000
6	-1.00000	1	11	4	0	0.54202	-1.05325	8.99844	9.07607	0.00290	0.00000	0.00000	0.00000
7	0.00000	23	1	4	8	-0.54202	1.05325	3.24471	3.46989	0.32997	0.00000	0.00000	0.00000
8	0.00000	62	1	7	10	-0.49168	0.95544	2.94255	3.14994	0.33001	0.00000	0.00000	0.00000
9	0.00000	63	2203	2	10	-0.05034	0.09781	0.13121	0.79010	0.77133	0.00000	0.00000	0.00000
10	-1.00000	83	-213	8	12	-0.40470	0.63771	2.14567	2.39994	0.76509	0.00000	0.00000	0.00000
11	0.00000	84	2224	8	14	-0.13731	0.41554	0.92809	1.54011	1.14850	0.00000	0.00000	0.00000
12	-1.00000	1	-211	10	0	-0.56961	0.42202	1.14111	1.35062	0.13958	0.00000	0.00000	0.00000
13	0.00000	1	111	10	0	0.16491	0.21569	1.00456	1.04932	0.13496	0.00000	0.00000	0.00000
14	1.00000	1	2212	11	0	0.01385	0.38961	0.82346	1.30784	0.93828	0.00000	0.00000	0.00000
15	1.00000	1	211	11	0	-0.15116	0.02593	0.10463	0.23227	0.13957	0.00000	0.00000	0.00000


13	1	1	0.00000	1.00000	11	10.60000	2212	1	0.00000
1	-1.00000	12	11	0	4	-0.00000	-0.00000	11.96700	11.96700	-0.00000	0.00000	0.00000	0.00000
2	0.00000	21	22	1	0	0.00000	-0.00000	-0.06870	-0.06870	-0.00000	0.00000	0.00000	0.00000
3	1.00000	12	2212	0	5	-0.00000	0.00000	-0.00259	0.93827	0.93827	0.00000	0.00000	0.00000
4	-1.00000	21	11	1	6	-0.00000	-0.00000	12.03570	12.03570	-0.00000	0.00000	0.00000	0.00000
5	0.00000	21	2	2	6	-0.00000	0.00000	-0.06910	0.06910	-0.00000	0.00000	0.00000	0.00000
6	-1.00000	1	11	4	0	-0.09912	0.84263	7.96210	8.00718	0.00625	0.00000	0.00000	0.00000
7	0.00000	23	2	4	8	0.09912	-0.84263	4.40363	4.49674	0.32994	0.00000	0.00000	0.00000
8	0.00000	62	2	7	10	0.09599	-0.81607	4.26447	4.35543	0.32998	0.00000	0.00000	0.00000
9	0.00000	63	2101	2	10	0.00312	-0.02657	-0.19346	0.61136	0.57933	0.00000	0.00000	0.00000
10	0.00000	1	111	8	0	0.30877	-1.02875	3.39005	3.55870	0.13506	0.00000	0.00000	0.00000
11	1.00000	1	211	8	0	-0.17657	0.12305	0.29211	0.38875	0.13957	0.00000	0.00000	0.00000
12	0.00000	1	2112	8	0	-0.03309	0.06307	0.38885	1.01935	0.93957	0.00000	0.00000	0.00000
*/

/**
NEW CURRENT:
16	1	1	0.00000	1.00000	11	10.60000	2212	1	0.00000
1	-1.00000	21	11	0	4	-0.00000	-0.00000	11.98460	11.98460	-0.00000	0.00000	0.00000	0.00000
2	1.00000	21	2212	0	5	-0.00000	0.00000	-0.00120	0.93827	0.93827	0.00000	0.00000	0.00000
3	0.00000	21	22	1	0	0.09960	0.52539	8.63699	8.59454	1.00868	0.00000	0.00000	0.00000
4	-1.00000	1	11	1	0	-0.09960	-0.52539	3.34761	3.39006	0.00755	0.00000	0.00000	0.00000
5	0.00000	21	2	2	4	0.00000	-0.00000	-0.03209	0.03209	-0.00000	0.00000	0.00000	0.00000
6	0.00000	13	2	4	8	0.09622	0.50758	8.76569	8.78710	0.33001	0.00000	0.00000	0.00000
7	0.00000	13	2103	2	8	0.00338	0.01781	-0.09789	0.77773	0.77133	0.00000	0.00000	0.00000
8	0.00000	11	113	6	11	-0.08312	-0.03175	4.54481	4.59112	0.64434	0.00000	0.00000	0.00000
9	1.00000	1	211	6	0	0.24974	0.42440	3.44304	3.48088	0.13969	0.00000	0.00000	0.00000
10	0.00000	11	2114	6	13	-0.06702	0.13273	0.67995	1.49284	1.32066	0.00000	0.00000	0.00000
11	1.00000	1	211	8	0	0.20396	0.13439	2.38503	2.40156	0.13950	0.00000	0.00000	0.00000
12	-1.00000	1	-211	8	0	-0.28708	-0.16613	2.15978	2.18956	0.13967	0.00000	0.00000	0.00000
13	0.00000	1	2112	10	0	-0.30382	0.26474	0.47601	1.12773	0.93957	0.00000	0.00000	0.00000
14	0.00000	12	0	0	0	0.00000	0.00000	0.00000	0.00000	-0.00000	0.00000	0.00000	0.00000
15	0.00000	12	0	0	0	0.00000	0.00000	0.00000	0.00000	-0.00000	0.00000	0.00000	0.00000
*/

/**
TARGET -> 21, 1, 13, 12, 11
1 -1.   21     11  0  0    0.0000    0.0000   10.6000   10.6000    0.0005     -0.0068   -0.0031    0.0000
2  1.   21   2212  0  0    0.0000    0.0000    0.0000    0.9383    0.9383     -0.0068   -0.0031    0.0000
3  0.   21     22  1  0   -0.8889    1.3687    4.6345    4.4152   -2.1558     -0.0068   -0.0031    0.0000
4 -1.    1     11  1  0    0.8889   -1.3687    5.9655    6.1848    0.0005     -0.0068   -0.0031    0.0000
5  1.   13      2  0  6   -1.1985    1.4051    4.1313    4.5253    0.0056     -0.0068   -0.0031    0.0000
6  0.   13   2101  2  0    0.3096   -0.0364    0.5032    0.8282    0.5793     -0.0068   -0.0031    0.0000
7  1.   12      2  5  9   -1.1985    1.4051    4.1313    4.5253    0.0056     -0.0068   -0.0031    0.0000
8  0.   11   2101  6  9    0.3096   -0.0364    0.5032    0.8282    0.5793     -0.0068   -0.0031    0.0000
9  0.   11     92  7 10   -0.8889    1.3687    4.6345    5.3535    2.1257     -0.0068   -0.0031    0.0000
10  0.   11    113  9 12   -1.1537    1.0718    3.0055    3.4320    0.5157     -0.0068   -0.0031    0.0000
11  1.    1   2212  9  0    0.2648    0.2969    1.6290    1.9215    0.9383     -0.0068   -0.0031    0.0000
12  1.    1    211 10  0   -0.4846    0.2310    0.7846    0.9608    0.1396     -0.0068   -0.0031    0.0000
13 -1.    1   -211 10  0   -0.6691    0.8408    2.2209    2.4712    0.1396     -0.0068   -0.0031    0.0000
*/

bool isQuark(int pid) {
    return abs(pid)<10 && pid!=0;
}

bool isDiQuark(int pid) {
    int abs_pid = abs(pid);
    return (abs_pid == 1103 || 
        abs_pid == 2101 || abs_pid == 2103 || 
        abs_pid == 2203 || 
        abs_pid == 3101 || abs_pid == 3103 || 
        abs_pid == 3201 || abs_pid == 3203 || 
        abs_pid == 3303
    );
}

bool isFinal(int pid, int daughter) {
    int abs_pid = abs(pid);
    return (abs_pid==11 || abs_pid==211 || abs_pid==321 || abs_pid==2212 || abs_pid==2112 && daughter==0);
}

int setStatus(int original_status, int lund_idx, int daughter, int pid) {

    // Incoming electron, proton, and virtual photon
    if (lund_idx<4) return 21;

    // Outgoing electron
    if (lund_idx==4) return 1;

    // Incoming quark
    if (lund_idx==5) return 21;

    // Quarks
    if (isQuark(pid)) {
        if (lund_idx<=6) return 13;
        else return 12;
    }

    // Quarks
    if (isDiQuark(pid)) {
        if (lund_idx<=7) return 13;
        else return 11;
    }

    // Final state particles
    if (isFinal(pid,daughter)) return 1;

    // Default for intermediate state particles
    return 11;
}

std::tuple<int, int, int, int, double, double, double, double> popParticleMapEntry(
        std::map<int, std::tuple<int, int, int, int, double, double, double, double>>& particleMap, // key: particle index[1:nParticles], value: (status, pid, mother, daughter, px, py, pz, E)
        int idx,
        int new_mother = 0,
        int new_daughter = 0
    ) {

    // Grab number of particles
    int nparticles = particleMap.size();

    // Replace incoming beam with duplicate incoming beam
    std::tuple<int, int, int, int, double, double, double, double> tuple_at_idx = particleMap[idx];

    // Reset all entries non-destructively
    for (int i=1; i<=nparticles; i++) {
        
        // Grab particle data
        auto& tup = particleMap[i];
        int status = std::get<0>(tup);
        int pid = std::get<1>(tup);
        int mother = std::get<2>(tup);
        int daughter = std::get<3>(tup);
        double px = std::get<4>(tup);
        double py = std::get<5>(tup);
        double pz = std::get<6>(tup);
        double E = std::get<7>(tup);

        // Reset mother
        if (mother>0 && mother==idx) mother = new_mother;
        else if (mother>0 && mother>idx) mother--;

        // Reset daughter
        if (daughter>0 && daughter==idx) daughter = new_daughter;
        else if (daughter>0 && daughter>idx) daughter--;

        // Reset entries
        if ( i==idx ) {
            tuple_at_idx = std::make_tuple(status, pid, mother, daughter, px, py, pz, E);
        } else {
            particleMap[i>idx ? i-1 : i] = std::make_tuple(status, pid, mother, daughter, px, py, pz, E);
        }
    }

    // Remove final entry
    particleMap.erase(nparticles);

    return tuple_at_idx;
}

void insertParticleMapEntry(
        std::map<int, std::tuple<int, int, int, int, double, double, double, double>>& particleMap, // key: particle index[1:nParticles], value: (status, pid, mother, daughter, px, py, pz, E)
        int idx,
        std::tuple<int, int, int, int, double, double, double, double> tuple_at_idx
    ) {

    // Grab number of particles
    int nparticles = particleMap.size();

    // Reset all entries non-destructively
    for (int i=nparticles+1; i>=1; i--) {
        
        // Grab particle data
        auto& tup = i==idx? tuple_at_idx : particleMap[i>idx ? i-1 : i];
        int status = std::get<0>(tup);
        int pid = std::get<1>(tup);
        int mother = std::get<2>(tup);
        int daughter = std::get<3>(tup);
        double px = std::get<4>(tup);
        double py = std::get<5>(tup);
        double pz = std::get<6>(tup);
        double E = std::get<7>(tup);

        // Reset mother and daughter
        if (mother>0 && mother>=idx) mother++;
        if (daughter>0 && daughter>=idx) daughter++;

        // Reset entries
        particleMap[i] = std::make_tuple(status, pid, mother, daughter, px, py, pz, E);
    }
}

void squeezeParticleMapEntries(
        std::map<int, std::tuple<int, int, int, int, double, double, double, double>>& particleMap, // key: particle index[1:nParticles], value: (status, pid, mother, daughter, px, py, pz, E)
        int idx1,
        int idx2
    ) {

    // Pull out earlier particle
    auto& tup1 = particleMap[idx1];
    int pid1 = std::get<1>(tup1);
    int mother1 = std::get<2>(tup1);
    int daughter1 = std::get<3>(tup1);

    // Pull out later particle
    auto& tup2 = particleMap[idx2];
    int status = std::get<0>(tup2);
    int pid = std::get<1>(tup2);
    int mother2 = std::get<2>(tup2);
    int daughter2 = std::get<3>(tup2);
    double px = std::get<4>(tup2);
    double py = std::get<5>(tup2);
    double pz = std::get<6>(tup2);
    double E = std::get<7>(tup2);

    if (daughter1!=idx2 && mother2!=idx1) {
        std::cerr << "ERROR: Cannot squeeze particles that are not in a mother-daughter relationship." << std::endl;
        std::cerr << "       idx1: " << idx1 << " (pid: " << pid1 << ", mother: " << mother1 << ", daughter: " << daughter1 << ")" << std::endl;
        std::cerr << "       idx2: " << idx2 << " (pid: " << pid << ", mother: " << mother2 << ", daughter: " << daughter2 << ")" << std::endl;
        return;
    }

    // Create squeezed particle
    auto tup_squeezed = std::make_tuple(status, pid, mother1, daughter2, px, py, pz, E);

    // Insert squeezed particle
    particleMap[idx1] = tup_squeezed;

    // pop later particle
    auto late = popParticleMapEntry(particleMap, idx2, idx1, idx1);
}

void createVirtualPhoton(
    std::map<int, std::tuple<int, int, int, int, double, double, double, double>>& particleMap,
    int incoming_beam_idx = 1,
    int outgoing_beam_idx = 3,
    int virtual_photo_idx = 3
) {

    //----- Manually insert virtual photon from incoming and scattered electrons -----//

    // Grab incoming electron
    auto &tup_e_0 = particleMap[incoming_beam_idx];
    double px_e_0 = std::get<4>(tup_e_0);
    double py_e_0 = std::get<5>(tup_e_0);
    double pz_e_0 = std::get<6>(tup_e_0);
    double E_e_0  = std::get<7>(tup_e_0);

    // Grab outgoing electron
    auto &tup_e_3 = particleMap[outgoing_beam_idx];
    double px_e_3 = std::get<4>(tup_e_3);
    double py_e_3 = std::get<5>(tup_e_3);
    double pz_e_3 = std::get<6>(tup_e_3);
    double E_e_3  = std::get<7>(tup_e_3);

    // Set virtual photon with momentum and energy conservation
    int status_ph    = 21;
    int pid_ph       = 22;
    int mother_ph    = incoming_beam_idx;
    int daughter_ph  = 0;
    double px_ph     = px_e_0-px_e_3;
    double py_ph     = py_e_0-py_e_3;
    double pz_ph     = pz_e_0-pz_e_3;
    double E_ph      = E_e_0-E_e_3;

    insertParticleMapEntry(
        particleMap,
        virtual_photo_idx,
        std::make_tuple(status_ph, pid_ph, mother_ph, daughter_ph, px_ph, py_ph, pz_ph, E_ph)
    );
}

void writeLundEvent(std::ofstream &out,
        std::map<int, std::tuple<int, int, int, int, double, double, double, double>>& particleMap, // key: particle index[1:nParticles], value: (status, pid, mother, daughter, px, py, pz, E)
        int nTargetNucleons = 1, // assuming hydrogen target
        int nTargetProtons = 1,  // assuming hydrogen target
        float targetPol = 0,       // unpolarized target
        float beamPol = 1,         // assuming positive helicity
        int beamPid = 11,        // beam particle id
        float beamEnergy = 10.6, // GeV
        int targetPid = 2212,    // proton
        int processId = 1,       // assuming DIS
        float integratedXS = 0.0, // not provided
        std::string delimiter = " ",
        int precision = 5
    ) {

    if (particleMap.size()<7) {
        std::cerr << "ERROR: particleMap.size() = " << particleMap.size() << " but need at least 7 particles to write LUND event." << std::endl;
        return;
    }

    // Modify particle map to comply with gemc lund format
    squeezeParticleMapEntries(particleMap, 3, 5); //NOTE: Squeeze scattered electron into duplicate incoming electron
    createVirtualPhoton(particleMap, 1, 3, 3); //NOTE: Create virtual photon from incoming and scattered electrons
    squeezeParticleMapEntries(particleMap, 6, 7); //NOTE: Squeeze scattered quarks AFTER squeezing scattered electron AND inserting virtual photon

    // Get number of particles
    int nparticles = particleMap.size();

    // Write event header
    // nparticles, nTargetNucleons(PARL1), nTargetProtons(PARL2),
    // targetPol, beamPol
    // beamPid, beamEnergy, targetPid, processId, integratedXS(PARL23)
    out << std::fixed << std::setprecision(precision) 
        << nparticles << delimiter.c_str()
        << nTargetNucleons << delimiter.c_str() << nTargetProtons << delimiter.c_str() 
        << targetPol << delimiter.c_str() << beamPol << delimiter.c_str()
        << beamPid << delimiter.c_str() << beamEnergy << delimiter.c_str()
        << targetPid << delimiter.c_str() << processId << delimiter.c_str()
        << integratedXS << std::endl;

    // //NOTE: ASSUME FIRST ENTRIES ARE:
    // - LUND STRING, (SKIPPED MANUALLY IN EVENT LOOP)
    // - INITIAL ELECTRON, 
    // - TARGET PROTON,
    // - VIRTUAL PHOTON, 
    // - SCATTERED ELECTRON,
    // - INITIAL QUARK
    // - INITIAL DIQUARK
    // - SCATTERED QUARK
    // - SCATTERED DIQUARK
    // - OTHERS...

    // Loop over particles
    for (int j = 1; j <=nparticles; j++) { //NOTE: LUND INDEX BEGINS AT 1

        // Set particle data
        auto& tup = particleMap[j];
        int status = std::get<0>(tup);
        int pid = std::get<1>(tup);
        int mother = std::get<2>(tup);
        int daughter = std::get<3>(tup);
        double px = std::get<4>(tup);
        double py = std::get<5>(tup);
        double pz = std::get<6>(tup);
        double E = std::get<7>(tup);

        // // Sanity checks
        // if (j==1 && pid!=beamPid)        cout << "ERROR: j=1 but pid="<<pid<<endl;
        // else if (j==2 && pid!=targetPid) cout << "ERROR: j=2 but pid="<<pid<<endl;
        // else if (j==3 && pid!=22)        cout << "ERROR: j=3 but pid="<<pid<<endl;
        // else if (j==4 && pid!=beamPid)   cout << "ERROR: j=4 but pid="<<pid<<endl;

        // Set particle status to comply with gemc
        status = setStatus(status, j, daughter, pid);

        // Compute mass and set vertices and get charge
        float charge = getCharge(pid);
        float m2 = E*E - (px*px + py*py + pz*pz);
        float m = (m2 > 0) ? sqrt(m2) : sqrt(-m2);
        float vx = 0.0, vy = 0.0, vz = 0.0; //NOTE: DO NOT RASTER FOR NOW

        // Write particle data
        out << std::fixed << std::setprecision(precision)
            << j << delimiter.c_str() //IMPORTANT: Lund index starts at 1!
            << charge << delimiter.c_str()
            << status << delimiter.c_str() << pid << delimiter.c_str()
            << mother << delimiter.c_str()
            << daughter << delimiter.c_str()
            << px << delimiter.c_str() << py << delimiter.c_str() << pz << delimiter.c_str()
            << E << delimiter.c_str() << m << delimiter.c_str()
            << vx << delimiter.c_str() << vy << delimiter.c_str() << vz
            << std::endl;
    }

    return;
}

void TreeToLund(const char* filename = "/Users/mfm45/Downloads/Tree-ep-at-JLAB12-posHelicity.root",
                const char* treename = "tree",
                const char* outname = "output.lund",
                int batch_idx = 1, // batch index in [1,nbatches]
                int batch_size = 10000, // batch size
                int nTargetNucleons = 1, // assuming hydrogen target
                int nTargetProtons = 1,  // assuming hydrogen target
                float targetPol = 0,       // unpolarized target
                float beamPol = 1,         // assuming positive helicity
                int beamPid = 11,        // beam particle id
                float beamEnergy = 10.6, // GeV
                int targetPid = 2212,    // proton
                int processId = 1        // assuming DIS
            ) {

    // Open input file and get TTree
    TFile* file = TFile::Open(filename);
    if (!file || file->IsZombie()) {
        std::cerr << "Cannot open input file!" << std::endl;
        return;
    }

    // Open TTree
    TTree* tree = (TTree*)file->Get(treename);
    if (!tree) {
        std::cerr << "Tree not found!" << std::endl;
        return;
    }

    // Set branches
    Int_t iEvent, iPos, status, pid, mother, daughter;
    Double_t px, py, pz, E;

    // Set branch addresses
    tree->SetBranchAddress("iEvent", &iEvent);
    tree->SetBranchAddress("iPos", &iPos);
    tree->SetBranchAddress("status", &status);
    tree->SetBranchAddress("id", &pid);
    tree->SetBranchAddress("iMother1", &mother);
    tree->SetBranchAddress("iDaughter1", &daughter);
    tree->SetBranchAddress("Px", &px);
    tree->SetBranchAddress("Py", &py);
    tree->SetBranchAddress("Pz", &pz);
    tree->SetBranchAddress("E", &E);

    // Open output file
    std::ofstream out(outname);
    if (!out.is_open()) {
        std::cerr << "Failed to open output file!" << std::endl;
        return;
    }

    // Initialize looping variables
    Long64_t nentries = tree->GetEntries();
    Int_t currentEvent = -1;
    std::map<int, std::tuple<int, int, int, int, double, double, double, double>> particleMap;

    // Loop TTree
    cout << "Processing entries min , max: " << (batch_idx-1)*batch_size << " , " << batch_idx*batch_size << endl;
    int eventCounter = 0;
    int pcounter = 0;
    for (Long64_t i = 0; i < nentries; ++i) {
        tree->GetEntry(i);

        // Assume events are grouped by iEvent in the TTree
        if (iEvent != currentEvent && currentEvent != -1) {

            eventCounter++;

            // Check if minimum event number has been reached
            if (currentEvent > (batch_idx-1)*batch_size) {

                // Write previous event
                writeLundEvent(
                    out,
                    particleMap,
                    nTargetNucleons, // nTargetNucleons
                    nTargetProtons, // nTargetProtons
                    targetPol, // targetPol (unpolarized)
                    beamPol, // beamPol (positive helicity)
                    beamPid, // beamPid (electron)
                    beamEnergy, // beamEnergy (GeV)
                    targetPid, // targetPid (proton)
                    processId, // processId (DIS)
                    0.0, // integratedXS
                    "\t", // delimiter
                    5 // precision
                );

                // Clear for next event
                particleMap.clear();
                pcounter = 0;
            }
        }

        // Set event number
        currentEvent = iEvent;

        // Break if maximum number of events for this batch_idx has been reached
        if (currentEvent >= batch_idx*batch_size) {
            break;
        }

        // Check if minimum index has been reached
        if (currentEvent >= (batch_idx-1)*batch_size) {

            // Skip initial LUND string entry
            if (pcounter==0 && pid!=beamPid) {
                pcounter++;
                continue; 
            }

            // Store particle data by particle counter [0:nParticles-1]
            particleMap[pcounter++] = std::make_tuple(status, pid, mother, daughter, px, py, pz, E);
        }
    }

    cout << "Total events processed: " << eventCounter << endl;

    // Write last event
    if (!particleMap.empty()) {
        writeLundEvent(out, particleMap);
    }

    // Close files
    out.close();
    file->Close();
    std::cout << "Lund file written to: " << outname << std::endl;
}
