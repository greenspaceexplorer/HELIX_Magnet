#include "coil.h"
//------------------------------------------------------------------------------
coil::coil(
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
        vector<unsigned int> elements_z){

    current = curr;
    width = wid;
    set_rotation(theta,phi,gamma);
    set_origin(x_origin,y_origin,z_origin);
    set_inner_radius(inner_radius);
    set_outer_radius(outer_radius);
    set_subcoil_turns(sub_turns);
    set_rho_elements(elements_rho);
    set_z_elements(elements_z);
}
//------------------------------------------------------------------------------
coil::coil(
        double curr,
        double wid,
        Vector3d rotation_axis,
        double rotation_magnitude,
        Vector3d coil_origin,
        vector<double> inner_radius,
        vector<double> outer_radius,
        vector<double> sub_turns,
        vector<unsigned int> elements_rho,
        vector<unsigned int> elements_z){
    set_current(curr);
    set_width(wid);
    set_rotation(rotation_magnitude,rotation_axis);
    set_origin(coil_origin);
    set_inner_radius(inner_radius);
    set_outer_radius(outer_radius);
    set_subcoil_turns(sub_turns);
    set_rho_elements(elements_rho);
    set_z_elements(elements_z);

}
//------------------------------------------------------------------------------
coil::coil(){

    set_current(93.);
    set_width(0.075818);
    set_rotation(0.,0.,0.);
    set_origin(0.,0.,0.);
    number_subcoils = 0;

    double inner_radius_a[] = {0.203517,0.215011,0.234632};
    double outer_radius_a[] = {0.215011,0.234632,0.255224};
    double sub_turns_a[] = {1995.9,5150.0,5489.5};
    unsigned int elements_rho_a[] = {8,8,8};
    unsigned int elements_z_a[] = {52,32,32};

    vector<double> inner_radius(inner_radius_a,inner_radius_a +\
            sizeof(inner_radius_a)/sizeof(double));
    vector<double> outer_radius(outer_radius_a,outer_radius_a +\
            sizeof(outer_radius_a)/sizeof(double));
    vector<double> sub_turns(sub_turns_a,sub_turns_a +\
            sizeof(sub_turns_a)/sizeof(double));
    vector<unsigned int> elements_rho(elements_rho_a,elements_rho_a +\
            sizeof(elements_rho_a)/sizeof(unsigned int));
    vector<unsigned int> elements_z(elements_z_a,elements_z_a +\
            sizeof(elements_z_a)/sizeof(unsigned int));

    set_inner_radius(inner_radius);
    set_outer_radius(outer_radius);
    set_subcoil_turns(sub_turns);
    set_rho_elements(elements_rho);
    set_z_elements(elements_z);
    // z-coordinate of coil origin is different though
}
//------------------------------------------------------------------------------
coil::~coil(){}
//------------------------------------------------------------------------------
Vector3d coil::B(Vector3d &position){
    Vector3d position_to_coil = position - origin;
    // rotate position to magnet coordinates
    position_to_coil = coil_rotation_inv*position_to_coil; 
    // calculate field
    Vector3d out = coil_at_origin(position_to_coil);
    // rotate field back to original orientation
    out = coil_rotation*out;
    // return field
    return out;
}
//------------------------------------------------------------------------------
void coil::set_current(double curr){ current = curr; }
//------------------------------------------------------------------------------
double coil::get_current(){ return current; }
//------------------------------------------------------------------------------
void coil::set_width(double wid){ width = wid; }
//------------------------------------------------------------------------------
double coil::get_width(){ return width; }
//------------------------------------------------------------------------------
void coil::set_rotation(double theta, double phi, double gamma){
    // convert rotation information to rotation matrix
    double ct,st,cp,sp;
    ct = cos(theta);
    st = sin(theta);
    cp = cos(phi);
    sp = sin(phi);
    Vector3d rotation_axis(st*cp,st*sp,ct);
    AngleAxisd coil_rotation_aa(gamma,rotation_axis);
    AngleAxisd coil_rotation_inv_aa = coil_rotation_aa.inverse();
    coil_rotation = coil_rotation_aa.toRotationMatrix();
    coil_rotation_inv = coil_rotation_inv_aa.toRotationMatrix();
}
//------------------------------------------------------------------------------
void coil::set_rotation(double gamma, Vector3d axis){
    AngleAxisd coil_rotation_aa(gamma,axis);
    AngleAxisd coil_rotation_inv_aa = coil_rotation_aa.inverse();
    coil_rotation = coil_rotation_aa.toRotationMatrix();
    coil_rotation_inv = coil_rotation_inv_aa.toRotationMatrix();
}
//------------------------------------------------------------------------------
void coil::compose_rotation(double theta, double phi, double gamma){
    double ct,st,cp,sp;
    ct = cos(theta);
    st = sin(theta);
    cp = cos(phi);
    sp = sin(phi);
    Vector3d rotation_axis(st*cp,st*sp,ct);
    AngleAxisd coil_rotation_aa(gamma,rotation_axis);
    AngleAxisd coil_rotation_inv_aa = coil_rotation_aa.inverse();
    Matrix3d next_rotation = coil_rotation_aa.toRotationMatrix();
    Matrix3d next_rotation_inv = coil_rotation_inv_aa.toRotationMatrix();
    coil_rotation = next_rotation*coil_rotation;
    coil_rotation_inv = coil_rotation_inv*next_rotation_inv;
}
//------------------------------------------------------------------------------
void coil::compose_rotation(double gamma, Vector3d axis){
    AngleAxisd coil_rotation_aa(gamma,axis);
    AngleAxisd coil_rotation_inv_aa = coil_rotation_aa.inverse();
    Matrix3d next_rotation = coil_rotation_aa.toRotationMatrix();
    Matrix3d next_rotation_inv = coil_rotation_inv_aa.toRotationMatrix();
    coil_rotation = next_rotation*coil_rotation;
    coil_rotation_inv = coil_rotation_inv*next_rotation_inv;
}
//------------------------------------------------------------------------------
Matrix3d coil::get_rotation(){ return coil_rotation; }
//------------------------------------------------------------------------------
Matrix3d coil::get_rotation_inv(){ return coil_rotation_inv; }
//------------------------------------------------------------------------------
void coil::set_origin(Vector3d coil_origin){ origin = coil_origin; }
//------------------------------------------------------------------------------
void coil::set_origin(double x, double y, double z){ 
    origin = Vector3d(x,y,z);
}
//------------------------------------------------------------------------------
Vector3d coil::get_origin(){ return origin; }
//------------------------------------------------------------------------------
void coil::set_inner_radius(vector<double> &i_r){
    if(number_subcoils){
        if(i_r.size() == number_subcoils){ 
            r_inner = vector<double>(i_r);
        }
        else{
            throw invalid_argument("Error: vector length must equal number of subcoils.");
        }
    }
    else{
        number_subcoils = i_r.size();
        r_inner = i_r;
    }
}
//------------------------------------------------------------------------------
void coil::set_outer_radius(vector<double> &o_r){
    if(number_subcoils){
        if(o_r.size() == number_subcoils){ r_outer = o_r; }
        else{
            throw invalid_argument("Error: vector length must equal number of subcoils.");
        }
    }
    else{
        number_subcoils = o_r.size();
        r_outer = o_r;
    }
}
//------------------------------------------------------------------------------
void coil::set_subcoil_turns(vector<double> &s_t){
    if(number_subcoils){
        if(s_t.size() == number_subcoils){ turns = s_t; }
        else{
            throw invalid_argument("Error: vector length must equal number of subcoils.");
        }
    }
    else{
        number_subcoils = s_t.size();
        turns = s_t;
    }
}
//------------------------------------------------------------------------------
void coil::set_rho_elements(vector<unsigned int> &e_rho){
    if(number_subcoils){
        if(e_rho.size() == number_subcoils){ n_rho = e_rho; }
        else{
            throw invalid_argument("Error: vector length must equal number of subcoils.");
        }
    }
    else{
        number_subcoils = e_rho.size();
        n_rho = e_rho;
    }
}
//------------------------------------------------------------------------------
void coil::set_z_elements(vector<unsigned int> &e_z){
    if(number_subcoils){
        if(e_z.size() == number_subcoils){ n_z = e_z; }
        else{
            throw invalid_argument("Error: vector length must equal number of subcoils.");
        }
    }
    else{
        number_subcoils = e_z.size();
        n_z = e_z;
    }
}
//------------------------------------------------------------------------------
void coil::print_coil_info(){
    AngleAxisd rotation_aa(coil_rotation);
    Vector3d rotation_axis = rotation_aa.axis();
    cout << "Current = " << current << " A" << endl;
    cout << "Width = " << width << " m" << endl;
    cout << "Rotation Angle = " << rotation_aa.angle() << " radians" << endl;
    cout << "Rotation Axis = (" 
        << rotation_axis(0)
        << ", " << rotation_axis(1)
        << ", " << rotation_axis(2) 
        << ") " << endl;
    cout << "Origin = (" 
        << origin(0)
        << ", " << origin(1)
        << ", " << origin(2) 
        << ") " << endl;

    const char fill = ' ';
    const int lblwidth = 15;
    const int valwidth = 15;

    cout << left << setw(lblwidth) << setfill(fill) << "Subcoil:";
    for(int i = 0; i < number_subcoils; i++){
        cout << left << setw(valwidth) << setfill(fill) << i;    
    }
    cout << endl;

    cout << left << setw(lblwidth) << setfill(fill) << "Turns:";
    for(int i = 0; i < number_subcoils; i++){
        cout << left << setw(valwidth) << setfill(fill) << turns[i];    
    }
    cout << endl;

    cout << left << setw(lblwidth) << setfill(fill) << "Inner Radius:";
    for(int i = 0; i < number_subcoils; i++){
        cout << left << setw(valwidth) << setfill(fill) << r_inner[i];    
    }
    cout << endl;
    cout << left << setw(lblwidth) << setfill(fill) << "Outer Radius:";
    for(int i = 0; i < number_subcoils; i++){
        cout << left << setw(valwidth) << setfill(fill) << r_outer[i];    
    }
    cout << endl;

    cout << left << setw(lblwidth) << setfill(fill) << "Rho Divisions:";
    for(int i = 0; i < number_subcoils; i++){
        cout << left << setw(valwidth) << setfill(fill) << n_rho[i];    
    }
    cout << endl;

    cout << left << setw(lblwidth) << setfill(fill) << "Z Divisions:";
    for(int i = 0; i < number_subcoils; i++){
        cout << left << setw(valwidth) << setfill(fill) << n_z[i];    
    }
    cout << endl;
}
//------------------------------------------------------------------------------
Vector3d coil::coil_at_origin(Vector3d &position){

    Vector3d B(0.,0.,0.);
    for(int i = 0; i < number_subcoils; i++){
        // current in ideal coils inside subcoils
        double current_ic = current*turns[i]/(n_rho[i]*n_z[i]);
        // width of subcoil in the rho direction
        double wrho = r_outer[i]-r_inner[i];
        // rho step size
        double drho = wrho/double(n_rho[i]);
        // z step size
        double dz = width/double(n_z[i]);
        // initialize subcoil field
        Vector3d bsub(0.,0.,0.);
        for(int j = 0; j < n_rho[i]; j++){
            // ideal coil radius
            double ric = r_outer[i]-drho*(double(j)+0.5);
            for(int k = 0; k < n_z[i]; k++){
                // ideal coil z-displacement
                double zic = (double(k)+0.5)*dz-width/2.;
                Vector3d pos_ic(position(0),position(1),position(2)-zic);
                bsub += b_loop(current_ic,ric,pos_ic);
            }
        }
        B += bsub;
    }	
    return B;
}
//------------------------------------------------------------------------------
Vector3d b_loop(double current, double radius, Vector3d &position){
    Vector3d B;
    double rhosq = position(0)*position(0)+position(1)*position(1);
    double rho = sqrt(rhosq);
    double d1 = (radius+rho)*(radius+rho)+position(2)*position(2);
    double bb = 2.e-7*current/sqrt(d1);

    if(rhosq == 0.){
        B(0) = 0.;
        B(1) = 0.;
        B(2) = bb*M_PI*radius*radius/d1;
    }
    else{

        double d2 = (radius-rho)*(radius-rho)+position(2)*position(2);
        double distsq = rhosq+position(2)*position(2);
        double n1 = radius*radius+distsq;
        double n2 = radius*radius - distsq;
        double kvar = sqrt(4.0*radius*rho/d1);

        double elk = gsl_sf_ellint_Kcomp(kvar,GSL_PREC_DOUBLE);
        double ele = gsl_sf_ellint_Ecomp(kvar,GSL_PREC_DOUBLE);

        double sint = position(1)/rho;
        double cost = position(0)/rho;

        double brho = (bb*position(2))/rho*(n1*ele/d2-elk);

        B(0) = brho*cost;
        B(1) = brho*sint;
        B(2) = bb*(elk+n2*ele/d2); 
    }
    return B;
}
