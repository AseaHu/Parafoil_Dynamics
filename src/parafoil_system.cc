#include <iostream>
#include <vector>
#include "Parafoil_Aerodynamics.h"

cv::RNG rng((unsigned)time(NULL));
double sig_wind = 0.3;

int main(){
    
    Parafoil_Aerodynamics PAD;

    PAD.SystemInitialization();

    //PAD.Sim_With_DesignedInputs();
    //PAD.Save_State();
    
    PVA pva;
    Vector3d v_wind;
    double time = PAD.Get_Time();
    for(int i = 0; i < 300; i ++){
        v_wind = Vector3d(rng.gaussian(sig_wind), rng.gaussian(sig_wind), rng.gaussian(sig_wind));
        pva = PAD.Sim_With_SingleStep(10, 10, v_wind);
    }
    for(int i = 0; i < 500; i ++){
        v_wind = Vector3d(rng.gaussian(sig_wind), rng.gaussian(sig_wind), rng.gaussian(sig_wind));
        pva = PAD.Sim_With_SingleStep(10, 20, v_wind);
    }
    for(int i = 0; i < 900; i ++){
        v_wind = Vector3d(rng.gaussian(sig_wind), rng.gaussian(sig_wind), rng.gaussian(sig_wind));
        pva = PAD.Sim_With_SingleStep(15, 15, v_wind);
    }
    for(int i = 0; i < 2000; i ++){
        v_wind = Vector3d(rng.gaussian(sig_wind), rng.gaussian(sig_wind), rng.gaussian(sig_wind));
        pva = PAD.Sim_With_SingleStep(15, 10, v_wind);
    }

    PAD.Save_State();

    return 0;
}
