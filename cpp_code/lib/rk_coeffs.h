#ifndef _COEFFS_H_
#define _COEFFS_H_

/*
 * This Section describes the classical RK4 Butcher Tableau
 *
const int s = 4;

double bE [] = {1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0};
double c  [] = {0.0, 0.5, 0.5, 1.0};

double AE [][s] = {
   {0.0, 0.0, 0.0, 0.0},
   {0.5, 0.0, 0.0, 0.0}, 
   {0.0, 0.5, 0.0, 0.0},
   {0.0, 0.0, 1.0, 0.0}
};

double AI [][s] = {
   {0.0, 0.0, 0.0 ,0.0},
   {0.0, 0.0, 0.0 ,0.0},
   {0.0, 0.0, 0.0 ,0.0},
   {0.0, 0.0, 0.0 ,0.0}
};
*/


/* IMEX Coefficients */

    const int s        = 5;
    const double g     = 8.717330430169640e-01;

    const double bE [] = {
        4.128980428124743e-01,
        0.0,
        1.973398903788692e-01,
        -4.610445469982568e-02,
        4.358665215084820e-01
    };

    const double bI [] = {
        4.128980428124743e-01,
        0.000000000000000e+00,
        1.973398903788692e-01,
        -4.610445469982568e-02,
        4.358665215084820e-01
    };


    // Intermediate time values
    const double c [] = { 
        0.0, 
        8.717330430169640e-01, 
        8.717330430168252e-01, 
        2.340212575108680e+00, 
        1.0, 
    };

    // Explicit Matrix //
    const double AE [][s] = 
    {
        {0.0, 0.0, 0.0, 0.0, 0.0},
        {8.717330430169640e-01, 0.0, 0., 0., 0.},
        {0.435866521508482, 0.435866521508482, 0., 0., 0.},
        {
            -8.009984530656288e-01,
            0.,
            3.141211028174309e+00,
            0.,
            0.
        },
        {
            3.567532077796395e-01,
            -1.973398903788692e-01,
            8.819488413937908e-01,
            -4.136215879462423e-02,
            0.0
        }
    };

    // implicit matrix //
    const double AI [][s] = 
    {
        {0., 0., 0., 0., 0.},
        {g, g,   0., 0., 0.},
        {g, -6.963508421466583e-14, g, 0., 0., },
        {-6.675868701958376e-02, -6.129659659417483e-13, 1.971104740620395e+00, g, 0.},
        {   4.128980428124743e-01,
            0.000000000000000e+00,
            1.973398903788692e-01,
            -4.610445469982568e-02,
            4.358665215084820e-01}
    };

#endif
