#include "helix_magnet.h"
#include "coil.h"
#include <iostream>

using namespace std;

int main(){
    // magnet characteristics

    double radius = 0.255524; // for single coil test
    double current,width,rotation_magnitude1,rotation_magnitude2;

    current            = 93.0;
    width              = 0.075818;
    rotation_magnitude1 = 0.0;
    rotation_magnitude2 = 0.0;

    double inner_radius1_a[]{0.203517,0.215011,0.234632},
           outer_radius1_a[]{0.215011,0.234632,0.255524},
           inner_radius2_a[]{0.203492,0.214732,0.234709},
           outer_radius2_a[]{0.214731,0.234709,0.256222},
           subturns1_a[]{1995.9,5150.0,5489.5},
           subturns2_a[]{1976.0,5110.0,5486.7};

    vector<double> 
        inner_radius1(inner_radius1_a,inner_radius1_a+
                sizeof(inner_radius1_a)/sizeof(double)),
        outer_radius1(outer_radius1_a,outer_radius1_a+
                sizeof(outer_radius1_a)/sizeof(double)),
        inner_radius2(inner_radius2_a,inner_radius2_a+
                sizeof(inner_radius2_a)/sizeof(double)),
        outer_radius2(outer_radius2_a,outer_radius2_a+
                sizeof(outer_radius2_a)/sizeof(double)),
        subturns1(subturns1_a,subturns1_a+
                sizeof(subturns1_a)/sizeof(double)),
        subturns2(subturns2_a,subturns2_a+
                sizeof(subturns2_a)/sizeof(double));

    unsigned int elements_rho_a[]{8,8,8},
                 elements_z_a[]{52,32,32};

    vector<unsigned int>
        elements_rho(elements_rho_a,elements_rho_a+
                sizeof(elements_rho_a)/sizeof(unsigned int)),
        elements_z(elements_z_a,elements_z_a+
                sizeof(elements_z_a)/sizeof(unsigned int));

    Vector3d rotation_axis1,rotation_axis2,coil_origin1,coil_origin2;

    rotation_axis1 = Vector3d(0.0,0.0,1.0);
    rotation_axis2 = Vector3d(0.0,0.0,1.0);
    coil_origin1 = Vector3d(0.0,0.0,-0.359747);
    coil_origin2 = Vector3d(0.0,0.0,0.359747); 

    //  Single coil test
    cout << "Single coil test:" << endl;
    Vector3d b_field;
    cout.precision(4);
    cout << scientific;
    // test on axis
    Vector3d position(0.0,0.0,0.1);
    b_field = b_loop(current,radius,position);
    cout << 
        "\n\t( X  Y  Z) = ( " <<
        position(0) << "  " << 
        position(1) << "  " <<
        position(2) << " )" << endl;
    cout <<
        "\t(BX BY BZ) = ( " <<
        b_field(0) << "  " <<
        b_field(1) << "  " <<
        b_field(2) << " )" << endl;

    // test off axis
    position = Vector3d(0.05,0.07,0.1);
    b_field = b_loop(current,radius,position);
    cout << 
        "\n\t( X  Y  Z) = ( " <<
        position(0) << "  " << 
        position(1) << "  " <<
        position(2) << " )" << endl;
    cout <<
        "\t(BX BY BZ) = ( " <<
        b_field(0) << "  " <<
        b_field(1) << "  " <<
        b_field(2) << " )" << endl;

// Fixed magnet test
    cout << "\nFull magnet test, no change in coil position:" << endl;
    coil coil1(current,width,rotation_axis1,rotation_magnitude1,coil_origin1,
            inner_radius1,outer_radius1,subturns1,elements_rho,elements_z),
         coil2(current,width,rotation_axis2,rotation_magnitude2,coil_origin2,
                 inner_radius2,outer_radius2,subturns2,elements_rho,elements_z);
    vector<coil> helix_coils;
    helix_coils.push_back(coil1);
    helix_coils.push_back(coil2);
    helix magneto(helix_coils);

    // test on axis
    position = Vector3d(0.0,0.0,0.1);
    b_field = magneto.B(position);

    cout << 
        "\n\t( X  Y  Z) = ( " <<
        position(0) << "  " << 
        position(1) << "  " <<
        position(2) << " )" << endl;
    cout <<
        "\t(BX BY BZ) = ( " <<
        b_field(0) << "  " <<
        b_field(1) << "  " <<
        b_field(2) << " )" << endl;

    // test off axis
    position = Vector3d(0.05,0.07,0.1);
    b_field = magneto.B(position);
    cout << 
        "\n\t( X  Y  Z) = ( " <<
        position(0) << "  " << 
        position(1) << "  " <<
        position(2) << " )" << endl;
    cout <<
        "\t(BX BY BZ) = ( " <<
        b_field(0) << "  " <<
        b_field(1) << "  " <<
        b_field(2) << " )\n" << endl;

// Magnet position change testing
    cout << "Positioning test:" << endl;
    cout << "\n\t--Rotations about the z-axis should not change the field--\n" << endl;
    magneto.coil_vec[0].set_rotation(0.0,0.0,3.141593/2.);
    magneto.coil_vec[1].set_rotation(0.0,0.0,3.141593/2.);
    magneto.print_magnet_info();
    position = Vector3d(0.05,0.07,0.1);
    b_field = magneto.B(position);
    cout << 
        "\n\t( X  Y  Z) = ( " <<
        position(0) << "  " << 
        position(1) << "  " <<
        position(2) << " )" << endl;
    cout <<
        "\t(BX BY BZ) = ( " <<
        b_field(0) << "  " <<
        b_field(1) << "  " <<
        b_field(2) << " )" << endl;

    cout << "\n\t--Flipping both coils should reverse the field from its original values--\n" << endl;
    magneto.coil_vec[0].set_rotation(3.141593/2.0,0.0,3.141593);
    magneto.coil_vec[1].set_rotation(3.141593/2.0,0.0,3.141593);
    magneto.print_magnet_info();
    position = Vector3d(0.05,0.07,0.1);
    b_field = magneto.B(position);
    cout << 
        "\n\t( X  Y  Z) = ( " <<
        position(0) << "  " << 
        position(1) << "  " <<
        position(2) << " )" << endl;
    cout <<
        "\t(BX BY BZ) = ( " <<
        b_field(0) << "  " <<
        b_field(1) << "  " <<
        b_field(2) << " )" << endl;

    cout << "\n\t--Flipping only one coil should mostly negate the field near the origin--\n" << endl;
    magneto.coil_vec[0].set_rotation(0.0,0.0,0.0);
    magneto.coil_vec[1].set_rotation(3.141593/2.0,0.0,3.141593);
    magneto.print_magnet_info();
    double midpoint = (magneto.coil_vec[0].get_origin()(2)+magneto.coil_vec[1].get_origin()(2))/2.;
    position = Vector3d(0.0,0.0,midpoint);
    b_field = magneto.B(position);
    cout << 
        "\n\t( X  Y  Z) = ( " <<
        position(0) << "  " << 
        position(1) << "  " <<
        position(2) << " )" << endl;
    cout <<
        "\t(BX BY BZ) = ( " <<
        b_field(0) << "  " <<
        b_field(1) << "  " <<
        b_field(2) << " )" << endl;

    return 0;
}
