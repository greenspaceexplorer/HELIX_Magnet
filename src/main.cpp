#include "helix_magnet.h"
#include "token_parser.h"
#include <iostream>

using namespace std;

vector<coil> parse_magnet_csv(const string csv_file);
void field_to_csv(vector< vector<Vector3d> >& data_vec, string filename);
vector< vector<Vector3d> > parse_data_csv(const vector<string>& data_csv_vec);
//double chisq(helix* my_helix,vector<vector<Vector3d> >& data_vec);
//------------------------------------------------------------------------------
int main(int argc, char *argv[]){
    if(argc < 2){
        throw invalid_argument("HELIX_MAGNET takes 2 or more arguments: magnet_config.csv data1.csv (data2.csv ...)");
    }
    // first input is the magnet configuration
    string magnet_csv = argv[1];
    vector<coil> helix_coils = parse_magnet_csv(magnet_csv);  

    // second and subsequent inputs are magnet scan files
    vector<string> data_csv_vec;
    for(int i = 2; i < argc; i++){
        data_csv_vec.push_back(argv[i]);
    }
    // make a vector of coordinates and field measurements
    vector< vector<Vector3d> > data_vec = parse_data_csv(data_csv_vec);
    // vector for the origin of the data in steps
    Vector3d scan_origin(9.1e3,1.2e4,1.83e5);
    // stage step distance in inches
    Vector3d step_conversion(9.7960e-4,7.7690e-4,6.0000e-5);
    // convert step distance to meters
    step_conversion*=2.54e-2;
    // rotation matrix for rotating coordinate system
    Vector3d rotation_vector(0.,1.,0.);
    double rotation_angle = M_PI/2.;
    AngleAxisd rotation_aa(rotation_angle,rotation_vector);
    Matrix3d rotation_matrix = rotation_aa.toRotationMatrix();
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

    // instantiate the HELIX magnet
    helix my_helix(helix_coils);

    // print data to check it is correctly loaded
    cout << "Loaded magnet data from " << magnet_csv << " ..." << endl;
    my_helix.print_magnet_info();

    // Calculate field and output to file
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
    string calc_name = "calculated_field.csv";
    string meas_name = "measured_field.csv";
    field_to_csv(b_calc_and_coords,calc_name);
    field_to_csv(data_vec,meas_name);
    cout << "Measured field output file: " << meas_name << endl;
    cout << "Calculated field output file: " << calc_name << endl;

    return 0;
}

//------------------------------------------------------------------------------
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
