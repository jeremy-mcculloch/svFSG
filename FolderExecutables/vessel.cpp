// vessel.h
#include <iostream>
#include <fstream>
#include <cstring>
#include <cmath>
#include <string>
#include <vector>
#include <iomanip>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_multiroots.h>
#include "json.hpp"
#include "vessel.h"
#include "functions.h"
#include "layer.h"
#include "constituent.h"

using std::string;
using std::vector;
using std::cout;

using json = nlohmann::json;

template <typename T>
// C++ template to print vessel contents
std::ostream& operator<<(std::ostream& os, const vector<T>& v) {
    for (int i = 0; i < v.size(); ++i) { 
        os << v[i]; 
        if (i != v.size() - 1) 
            os << " "; 
    }
    return os; 
}

template <typename T>
// C++ template to read vector container elements
std::ifstream& operator>>(std::ifstream& in, vector<T>& v) {
    T x;
    char next;
    int i=0;
    while(in.get(next))
    {
        if (next == '\n')  
        {    break; }
    }
    while ((in.peek()!='\n') && (in>>x))
    {
        v[i]=x;
        i++;
    }
    return in; 
}

vessel::vessel() { 
    vessel_name = "default";
    file_name = "Vs_out";
    gnr_name = "GnR_out";
    equil_gnr_name = "Equil_GnR_out";
    exp_name = "Exp_out";
    hns_name = "HnS_out";
}
void vessel::initializeJSON(string json_name, double n_days_inp, double dt_inp) {
    
    //Set native vessel time parameters
    double n_days = n_days_inp;
    dt = dt_inp;
    nts = int(n_days / dt);
    sn = 0;
    s = 0.0;

    double mu = 0; //apparent viscosity

    std::ifstream f(json_name);
    json json_in = json::parse(f);
    vessel_name = json_in["vessel_name"];
    lambda_z_h = json_in["in_vivo_z_stretch"];
    P_h = json_in["pressure"]; 
    Q_h = json_in["flow_rate"];
   
    // Computed values
    f_h = 0;
    P = P_h;
    Q = Q_h;
    P_prev = P_h;     
//$    f = 0.0;
    mb_equil = 0;
    mb_equil_e = 0;
    f_z_e = 0;

    pol_only_flag = 0;
    wss_calc_flag = 0;
    app_visc_flag = 0;
    num_exp_flag = 0;

    // Load constituents
    json layers_in = json_in["layers"];
    layers = {};
    int n_layers = layers_in.size();
    for (int i = 0; i < n_layers; i++) {
        layer next_layer;
        layers.push_back(next_layer);
        layers[i].initializeJSON(layers_in[i], this);
    } 

    // Set up appropriate parent relationships
    for (int layer = 0; layer < layers.size(); layer ++) {
        for (int alpha = 0; alpha < layers[layer].constituents.size(); alpha++) {
            layers[layer].constituents[alpha].parent_layer = &layers[layer];
        }
        layers[layer].parent_vessel = this;
    }
    
    double a_h = layers[0].a_h; 
    mu = get_app_visc(this, sn);
    bar_tauw_h = 4*mu*Q_h/(3.14159265*pow(a_h*100, 3));
    bar_tauw = bar_tauw_h;    

    // int equil_check = 0;
    // int num_act_incr = 10;
    // num_exp_flag = 1;
    // for (int i = 1; i < num_act_incr; i++){
    //     T_act = T_act_h * (1 - static_cast<double>(i) / num_act_incr);
    //     equil_check = find_iv_geom(this);
    // }
    // num_exp_flag = 0;


    //Set reference values
    for (int layer = 0; layer < layers.size(); layer ++) {

        //Update stress to in vivo stress for homeostatic stress calc
        update_sigma(&layers[layer]);
        layers[layer].sigma_h = layers[layer].sigma;
        layers[layer].sigma_inv_h = layers[layer].sigma_inv;


    }

}

