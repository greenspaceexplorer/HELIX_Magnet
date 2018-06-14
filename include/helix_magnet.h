#ifndef HELIX_MAGNET_H
#define HELIX_MAGNET_H

#include "coil.h"

class helix{
    public:
        // constructor
        helix(vector<coil> helix_coils);
        // destructor
        ~helix();
        // Returns B field at given position
        Vector3d B(Vector3d &position);
        // Print stored magnet information
        void print_magnet_info();
    private:
        // vector containing helix coil objects
        vector<coil> coil_vec;
        // magnitude of thermal contraction for copper
        double cont_cu = 0.99674;
        // magnitude of thermal contraction for aluminum
        double cont_al = 0.99585;

};
#endif // HELIX_MAGNET_H
