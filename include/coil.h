//------------------------------------------------------------------------------
// Class for FEA calculation of the magnetic field of a compound coil with
//  a rotation and displacement.
//
//  Author: Noah Green
//------------------------------------------------------------------------------
#ifndef HELIX_COIL_H
#define HELIX_COIL_H

#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <vector>
#include <Eigen/Geometry>
#include <cmath>
#include <gsl/gsl_sf_ellint.h>

using namespace Eigen;
using namespace std;

class coil{
    public:
        // constructor
        coil(
                double curr,
                double wid,
                double theta,
                double phi,
                double gamma,
                double x_origin,
                double y_origin,
                double z_origin,
                vector<double> inner_radius,
                vector<double> outer_radius,
                vector<double> sub_turns,
                vector<unsigned int> elements_rho,
                vector<unsigned int> elements_z);

        // alternate constructor using Eigen objects
        coil(
                double curr,
                double wid,
                Vector3d rotation_axis,
                double rotation_magnitude,
                Vector3d coil_origin,
                vector<double> inner_radius,
                vector<double> outer_radius,
                vector<double> sub_turns,
                vector<unsigned int> elements_rho,
                vector<unsigned int> elements_z);
        coil();
        // destructor
        ~coil();

        // set and get coil current
        double get_current();
        void set_current(double curr);

        // set and get coil width
        double get_width();
        void set_width(double wid);

        // set the coil's rotation
        void set_rotation(double theta, double phi, double gamma);
        void set_rotation(double gamma, Vector3d axis);
        // compose another rotation on top of the coil's current rotation
        void compose_rotation(double theta, double phi, double gamma);
        void compose_rotation(double gamma, Vector3d axis);
        // get rotation matrix of coil's rotation
        Matrix3d get_rotation();
        Matrix3d get_rotation_inv();

        // set and get coil origin
        void set_origin(Vector3d coil_origin);
        void set_origin(double x, double y, double z);
        Vector3d get_origin();

        // set and get subcoil properties
        void set_inner_radius(vector<double> &i_r);
        void set_outer_radius(vector<double> &o_r);
        void set_subcoil_turns(vector<double> &s_t);
        void set_rho_elements(vector<unsigned int> &e_rho);
        void set_z_elements(vector<unsigned int> &e_z);

        vector<double> get_inner_radius(){ return r_inner; }
        vector<double> get_outer_radius(){ return r_outer; }
        vector<double> get_subcoil_turns(){ return turns; }
        vector<unsigned int> get_rho_elements(){ return n_rho; }
        vector<unsigned int> get_z_elements(){ return n_z; }

        // Returns magnetic field of coil at the input position
        Vector3d B(Vector3d &position);

        // prints out stored coil information
        void print_coil_info();

        // Returns magnetic field of coil at position if it were positioned
        //  at the origin.
        Vector3d coil_at_origin(Vector3d &position); 
    private:				
        // coil current in amperes
        double current; 
        // coil width in meters
        double width; 
        // rotaton matrix for coil
        Matrix3d coil_rotation; 
        // inverse rotation matrix for coil
        Matrix3d coil_rotation_inv; 
        // coil origin
        Vector3d origin; 
        // number of subcoils
        unsigned int number_subcoils;
        // inner radius of each subcoil
        vector<double> r_inner; 
        // outer radius of each subcoil
        vector<double> r_outer; 
        // number of element divisions in the rho direction of each subcoil
        vector<unsigned int> n_rho; 
        // number of element divisions in the z direction of each subcoil
        vector<unsigned int> n_z; 		
        // number of turns in each subcoil
        vector<double> turns; 
};

// Returns the B field of an ideal current loop positioned at the origin.
Vector3d b_loop(double current, double radius, Vector3d &position);

#endif // HELIX_COIL_H

