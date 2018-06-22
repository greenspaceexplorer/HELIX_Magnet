#ifndef HELIX_MAGNET_H
#define HELIX_MAGNET_H

#include "coil.h"

class helix{
    public:
        // constructor
        //  the stationary coil is the last one in the helix_coils vector
        helix(vector<coil> helix_coils, bool contract_al = true);
        // destructor
        ~helix();
        // Returns B field at given position
        Vector3d B(Vector3d &position);
        // Print stored magnet information
        void print_magnet_info();
        // vector containing helix coil objects
        vector<coil> coil_vec;
    private:
        // magnitude of thermal contraction for copper
        double cont_cu = 0.99674;
        // magnitude of thermal contraction for aluminum
        double cont_al = 0.99585;
        bool apply_cont_al;

};
#endif // HELIX_MAGNET_H
