#include "helix_magnet.h"
#include "token_parser.h"
#include <iostream>
#include <gsl/gsl_multimin.h>

using namespace std;

// struct for passing parameters to fit function
struct fit_params{
    helix* magnet;
    vector< vector<Vector3d> >* data;
    unsigned int dilution;
};

double fit_func(const gsl_vector* v, void* params);
vector<coil> parse_magnet_csv(const string csv_file);
void field_to_csv(vector< vector<Vector3d> >& data_vec, string filename);
vector< vector<Vector3d> > parse_data_csv(const vector<string>& data_csv_vec);
//double chisq(helix* my_helix,vector<vector<Vector3d> >& data_vec);
//------------------------------------------------------------------------------
int main(int argc, char *argv[]){
    if(argc < 2){
        throw invalid_argument("HELIX_MAGNET takes 2 or more arguments: magnet_config.csv data1.csv (data2.csv ...)");
    }
// MAKE HELIX MAGNET OBJECT
    // first input is the magnet configuration
    string magnet_csv = argv[1];
    vector<coil> helix_coils = parse_magnet_csv(magnet_csv);  
    // instantiate the HELIX magnet
    helix my_helix(helix_coils);
    // print data to check it is correctly loaded
    cout << "Loaded magnet data from " << magnet_csv << " ..." << endl;
    my_helix.print_magnet_info();

// LOAD MAGNET SCAN DATA
    // second and subsequent inputs are magnet scan files
    vector<string> data_csv_vec;
    for(int i = 2; i < argc; i++){
        data_csv_vec.push_back(argv[i]);
    }
    // make a vector of coordinates and field measurements
    vector< vector<Vector3d> > data_vec = parse_data_csv(data_csv_vec);
    // vector for the origin of the data in measuring apparatus steps
    Vector3d scan_origin(9.1e3,1.2e4,-1.755e5);
    // vector to convert steps to inches
    Vector3d step_conversion(9.7960e-4,7.7690e-4,6.0000e-5);
    step_conversion*=2.54e-2; // convert inches to meters
    // Create rotation matrix to put magnet coils along the z-axis in scan data.
    Vector3d rotation_vector1(0.,0.,1.);
    double rotation_angle1 = -M_PI/2.;
    AngleAxisd rotation_aa1(rotation_angle1,rotation_vector1);
    Matrix3d rotation_matrix1 = rotation_aa1.toRotationMatrix();
    Vector3d rotation_vector2(0.,1.,0.);
    double rotation_angle2 = M_PI/2.;
    AngleAxisd rotation_aa2(rotation_angle2,rotation_vector2);
    Matrix3d rotation_matrix2 = rotation_aa2.toRotationMatrix();
    Matrix3d rotation_matrix = rotation_matrix2*rotation_matrix1;
    // change coordinate system and unit of data for consistency with model
    for(int i = 0; i < data_vec[0].size(); i++){
        // convert b field to teslas
        data_vec[1][i]*=1.e-4;
        // shift origin of coordinates
        data_vec[0][i]-=scan_origin;
        // convert units to meters
        for(int j = 0; j < 3; j++){data_vec[0][i](j)*=step_conversion(j);}
        // apply rotation to put coils along z-axis
        data_vec[0][i] = rotation_matrix*data_vec[0][i];
        data_vec[1][i] = rotation_matrix*data_vec[1][i];
    }


    fit_params helix_fit_params;
    helix_fit_params.magnet = &my_helix;
    helix_fit_params.data = &data_vec;
    helix_fit_params.dilution = 100; // only evaluate every nth point
    gsl_vector *v;
    v = gsl_vector_alloc(4);
    gsl_vector_set(v,0,0.);
    gsl_vector_set(v,1,0.);
    gsl_vector_set(v,2,0.);
    gsl_vector_set(v,3,0.);

    double test_sum = fit_func(v,&helix_fit_params);
    cout << "test_sum = " << test_sum << endl;

// DETERMINE CALCULATED MAGNETIC FIELD AT MAGNET SCAN POINTS
/*
    cout << "CALCULATING FIELD FOR " << data_vec[0].size() << " POINTS" << endl;
    vector<Vector3d> b_calc;
    for(int i = 0; i < data_vec[0].size(); i++){
        if(i%500 == 0){
            cout << data_vec[0].size()-i << " points remaining..." << endl;
        }
        Vector3d calc = my_helix.B(data_vec[0][i]);
        b_calc.push_back(calc);
    }
    vector< vector<Vector3d> > b_calc_and_coords;
    b_calc_and_coords.push_back(data_vec[0]);
    b_calc_and_coords.push_back(b_calc);

// SAVE MEASURED AND CALCULATED FIELD DATA TO FILE
    string calc_name = "calculated_field.csv";
    string meas_name = "measured_field.csv";
    field_to_csv(b_calc_and_coords,calc_name);
    field_to_csv(data_vec,meas_name);
    cout << "Measured field output file: " << meas_name << endl;
    cout << "Calculated field output file: " << calc_name << endl;
*/
    gsl_vector_free(v);
    return 0;
}

//------------------------------------------------------------------------------
// Minimizing this function will give the best fit b/w calculated and measured
//  magnetic fields.
double fit_func(const gsl_vector* v, void* params){
    // get parameters from void pointer
    fit_params* par;
    par = (fit_params*)params;
    helix* magnet = par->magnet;
    vector< vector<Vector3d> >* data = par->data;
    
    // position magnet coils from values in v
    magnet->coil_vec[0].set_rotation(0.,gsl_vector_get(v,0),gsl_vector_get(v,1));
    magnet->coil_vec[1].set_rotation(0.,gsl_vector_get(v,2),gsl_vector_get(v,3));

    // Sum over the modulus squared of the difference 
    //  in measured and calculated fields.
    double sum = 0.0;
    for(int i = 0; i < (*data)[0].size(); i++){
        if(i%par->dilution == 0){
            Vector3d diff = magnet->B((*data)[0][i]) - (*data)[1][i];
            sum += diff.dot(diff);
        }
    }
    return sum;

}
//------------------------------------------------------------------------------
// Saves magnetic field data to a csv file
void field_to_csv(vector< vector<Vector3d> >& data_vec, string filename){
    ofstream ofs(filename,ofstream::out);
    for(int i = 0; i < data_vec[0].size(); i++){
        for(int j = 0; j < 3; j++){ ofs << data_vec[0][i](j) << ",";}
        for(int k = 0; k < 2; k++){ ofs << data_vec[1][i](k) << ",";}
        ofs << data_vec[1][i](2) << endl;
    }
    ofs.close();
}
//------------------------------------------------------------------------------
// Returns a vector of the form 
//  (vector<Vector3d> coordinates, vector<Vector3d> B_field)
//  from a vector of magnet scan csv files
vector< vector<Vector3d> > parse_data_csv(const vector<string>& data_csv_vec){
    vector<Vector3d> coordinates,bfield;
    vector<vector<Vector3d> > out;
    out.push_back(coordinates);
    out.push_back(bfield);
    int count;
    for(int i = 0; i < data_csv_vec.size(); i++){
        ifstream csv_stream(data_csv_vec[i]);
        token_parser csv_row(',');
        count = 0;
        while(csv_stream >> csv_row){
            Vector3d coord;
            Vector3d field;
            coord(0) = stod(csv_row[1]);
            coord(1) = stod(csv_row[2]);
            coord(2) = stod(csv_row[3]);
            field(0) = stod(csv_row[4]);
            field(1) = stod(csv_row[5]);
            field(2) = stod(csv_row[6]);
            out[0].push_back(coord);
            out[1].push_back(field);
            count++;
        }
        cout << "Loaded " << count << " points from " 
            << data_csv_vec[i] << " ..." << endl;
    }
    return out;
}

//------------------------------------------------------------------------------
// Returns a vector of coil objects that can instantiate a helix object given
//  a properly formatted magnet configuration csv.
vector<coil> parse_magnet_csv(const string csv_file){
    ifstream csv_stream(csv_file);
    vector<coil> out;
    token_parser csv_row(',');

    int count = -1;
    while(csv_stream >> csv_row){
        if(csv_row[0] == "new_coil"){
            count++;
            coil new_coil = coil();
            out.push_back(new_coil);
        }
        else if(csv_row[0] == "current"){
            double current = stod(csv_row[1]);
            out[count].set_current(current);
        }
        else if(csv_row[0] == "width"){
            double width = stod(csv_row[1]);
            out[count].set_width(width);
        }
        else if(csv_row[0] == "rotation_angleaxis"){
            double rotation_angle = stod(csv_row[1]);
            Vector3d rotation_axis;
            for(int i = 2; i < csv_row.size(); i++){
                rotation_axis(i-2) = stod(csv_row[i]);
            }
            out[count].set_rotation(rotation_angle,rotation_axis);
        }
        else if (csv_row[0] == "origin"){
            Vector3d origin;
            for(int i = 1; i < csv_row.size(); i++){
                origin(i-1) = stod(csv_row[i]);
            }
            out[count].set_origin(origin);
        }
        else if(csv_row[0] == "inner_radius"){
            vector<double> inner_radius;
            for(int i = 1; i < csv_row.size(); i++){
                double ir_d = stod(csv_row[i]);
                inner_radius.push_back(ir_d);
            }
            out[count].set_inner_radius(inner_radius);
        }
        else if(csv_row[0] == "outer_radius"){
            vector<double> outer_radius;
            for(int i = 1; i < csv_row.size(); i++){
                outer_radius.push_back(stod(csv_row[i]));
            }
            out[count].set_outer_radius(outer_radius);
        }
        else if(csv_row[0] == "subcoil_turns"){
            vector<double> subcoil_turns;
            for(int i = 1; i < csv_row.size(); i++){
                subcoil_turns.push_back(stod(csv_row[i]));
            }
            out[count].set_subcoil_turns(subcoil_turns);
        }
        else if(csv_row[0] == "subcoil_rho_div"){
            vector<unsigned int> subcoil_rho_div;
            for(int i = 1; i < csv_row.size(); i++){
                subcoil_rho_div.push_back(stoi(csv_row[i]));
            }
            out[count].set_rho_elements(subcoil_rho_div);
        }
        else if(csv_row[0] == "subcoil_z_div"){
            vector<unsigned int> subcoil_z_div;
            for(int i = 1; i < csv_row.size(); i++){
                subcoil_z_div.push_back(stoi(csv_row[i]));
            }
            out[count].set_z_elements(subcoil_z_div);
        }
        else{
            throw "Error: invalid row encountered while parsing magnet data.";
        }
    }
    return out; 
}
