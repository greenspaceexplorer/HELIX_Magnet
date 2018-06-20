#include "helix_magnet.h"
#include "coil.h"
#include <iostream>

using namespace std;

int main(){
    
//  Single coil test
    double current = 93.0;
    double radius = 0.255224;
    Vector3d b_field;
    // first test on axis
    Vector3d on_axis(0.0,0.0,0.1);
    b_field = b_loop(current,radius,on_axis);
    cout << 
        "(X Y Z) = (" <<
        on_axis(0) << " " << 
        on_axis(1) << " " <<
        on_axis(2) << ")" << endl;
    cout <<
        "(BX BY BZ) = (" <<
        b_field(0) << " " <<
        b_field(1) << " " <<
        b_field(2) << ")" << endl;
    

    return 0;
}
