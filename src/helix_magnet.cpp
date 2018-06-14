#include "helix_magnet.h"

helix::helix(vector<coil> helix_coils){
    coil_vec = helix_coils;
    // apply thermal contraction
    for(int i = 0; i < coil_vec.size(); i++){
        // apply copper contraction to coil width
        double coil_w = coil_vec[i].get_width();
        coil_vec[i].set_width(coil_w*cont_cu);

        // apply aluminum contraction to coil separation
        Vector3d coil_origin = coil_vec[i].get_origin();
        coil_origin(2)*=cont_al;
        coil_vec[i].set_origin(coil_origin);

        // apply copper contraction to subcoil radii
        vector<double> coil_ir = coil_vec[i].get_inner_radius();
        vector<double> coil_or = coil_vec[i].get_outer_radius();
        for(int j = 0; j < coil_ir.size(); j++){
            coil_ir[j]*=cont_cu;
            coil_or[j]*=cont_cu;
        } 
        coil_vec[i].set_inner_radius(coil_ir);
        coil_vec[i].set_outer_radius(coil_or); 
    }  
}

helix::~helix(){};


Vector3d helix::B(Vector3d &position){
    Vector3d out(0.,0.,0.);
    for(int i=0; i < coil_vec.size(); i++){
        out += coil_vec[i].B(position);
    }
    return out;
}

void helix::print_magnet_info(){
    cout << "Note: thermal contraction applied to coil sizes and displacements." 
        << endl;
    cout << "Thermal contraction of copper = " << cont_cu << endl;
    cout << "Thermal contraction of aluminum = " << cont_al << endl;
    for(int i = 0; i < coil_vec.size(); i++){
        cout << "------Coil " << i << "------" << endl;
        coil_vec[i].print_coil_info();
    }
}
