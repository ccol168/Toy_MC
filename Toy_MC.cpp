#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <random>


using namespace std;

struct Photon {
    double theta;
    double phi;
    double gen_time;
    double arr_time;
    bool type;
};

struct Direction {
    double theta;
    double phi;
};

//parameters (some should be taken from a configuration file)

double Juno_radius = 40.1 ; //m
int Effective_Rate = 1600; //scintillation photons @ 1 MeV
double Cher_Fraction = 0.01; // Cherenkov photons produced as a fraction of scintilation photons
double n = 1.55 ; //refraction index
double c = 299792458 ; // m/s
double tau = 4*pow(10,-9); // s
double m_e = 0.51099895; //MeV    electron mass
double Be7_energy = 0.862; //MeV    enegy of a 7Be neutrino

//random number generator
double seed = 56;
mt19937_64 rnd (seed) ; //stdlib mersenne twister
uniform_real_distribution <> flat(0,1);
exponential_distribution <> expo (1./tau) ;


//to keep phi and theta in their intended range
double Pbc_phi (double in) {
    if (in<0) return in + M_PI;
    else if (in>M_PI) return in - M_PI;
    else return in;
}

double Pbc_theta (double in) {
    if (in<0) return in + 2*M_PI;
    else if (in>2*M_PI) return in - 2*M_PI;
    else return in;
}

void PrintPhoton (ofstream& fileout , Photon in ) {

    fileout << in.theta << "  " << in.phi<< "  " << in.gen_time << "  " << in.arr_time << "  "<< in.type << endl;

    return;
}

Direction Generate_Cone (double theta_0, double phi_0, double angle) {

    Direction dir_out;
    double theta, phi, cos_th1, cos_ph1, sin_th1, sin_ph_sin_th, sin_ph_cos_th;

    theta = 2*M_PI*flat(rnd);
    phi = angle; 

    if (phi_0 == 0) {

        dir_out.phi = phi;
        dir_out.theta = theta;
        return dir_out;

    } else if (phi_0 == M_PI) {

        dir_out.phi = M_PI - phi;
        dir_out.theta = theta;
        return dir_out;

    } else {

        cos_ph1 = -sin(phi_0)*cos(theta)*sin(phi) + cos(phi_0)*cos(phi);
        sin_ph_cos_th = cos(phi_0)*cos(theta_0)*cos(theta)*sin(phi) - sin(theta_0)*sin(theta)*sin(phi) + cos(theta_0)*sin(phi_0)*cos(phi);
        sin_ph_sin_th = sin(theta_0)*cos(phi_0)*cos(theta)*sin(phi) + cos(theta_0)*sin(theta)*sin(phi) + sin(theta_0)*sin(phi_0)*cos(phi);        

        dir_out.phi = acos(cos_ph1);
        cos_th1 = sin_ph_cos_th/sin(dir_out.phi);
        sin_th1 = sin_ph_sin_th/sin(dir_out.phi);

        if (sin_th1 >= 0) {
            dir_out.theta = acos(cos_th1) ;
        } else if (sin_th1 < 0) {
            dir_out.theta = Pbc_theta(-acos(cos_th1));
        }


    //cout<<dir_out.phi<<"  "<<cos_th1<<"   "<<sin_th1<<"   "<<cos_th1*cos_th1 + sin_th1*sin_th1<<endl;


        return dir_out;

    }

    
}

void Generate_Central_Photons (ofstream& fileout , double Event_Energy) {

    int N_Scint_Photons = Effective_Rate * Event_Energy;
    int N_Cher_Photons = Cher_Fraction * N_Scint_Photons;
    Photon Scint_Photon;

    //scintillation photons
    for (int i = 0; i < N_Scint_Photons; i++)  {
        Scint_Photon.theta = flat(rnd) * 2 * M_PI ;
        Scint_Photon.phi = flat(rnd) * M_PI ;
        Scint_Photon.gen_time = expo(rnd);
        Scint_Photon.arr_time = Scint_Photon.gen_time + Juno_radius*n/c ;
        Scint_Photon.type = 0;

        PrintPhoton (fileout,Scint_Photon);
    }

    double theta_e = acos((1+m_e/Be7_energy)*pow(Event_Energy/(Event_Energy+2*m_e),0.5)); //angle between the solar-nu and the electron scattered (assuming 7Be-nu)
    double beta = pow(1-(pow(m_e/(Event_Energy+m_e),2)),0.5) ; //beta of the electron generated
    double theta_Cher = acos(1/(beta*n)); //Cherenkov angle

    //cout<<theta_e<<"   "<<beta<<"   "<<theta_Cher<<endl;
    Direction Cher_electron_dir, Cher_photon_dir;
    Photon Cher_photon;

    for (int i = 0; i < N_Cher_Photons; i++) {
        Cher_electron_dir = Generate_Cone(M_PI-1,M_PI_2,theta_e);
        Cher_photon_dir = Generate_Cone(Cher_electron_dir.theta,Cher_electron_dir.phi,theta_Cher);
        Cher_photon.phi = Cher_photon_dir.phi;
        Cher_photon.theta = Cher_photon_dir.theta;
        Cher_photon.gen_time = 0.;
        Cher_photon.arr_time = Juno_radius*n/c;
        Cher_photon.type = 1;

        PrintPhoton (fileout,Cher_photon);

    }
    
    return;
}


int main() {

    ofstream fileout ("output.txt");
    double N_events = 1000;
    double Event_Energy = 0.5;

    for (int i=0;i<N_events;i++) {
        Generate_Central_Photons (fileout,Event_Energy);
    }

    fileout.close();

    return 0;
}


