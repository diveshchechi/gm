/*************************************************************************************/
/*                                                                                    */
/*  Visualization Library                                                             */
/*  http://visualizationlibrary.org                                                   */
/*                                                                                    */
/*  Copyright (c) 2005-2020, Michele Bosi                                             */
/*  All rights reserved.                                                              */
/*                                                                                    */
/*  Redistribution and use in source and binary forms, with or without modification,  */
/*  are permitted provided that the following conditions are met:                     */
/*                                                                                    */
/*  - Redistributions of source code must retain the above copyright notice, this     */
/*  list of conditions and the following disclaimer.                                  */
/*                                                                                    */
/*  - Redistributions in binary form must reproduce the above copyright notice, this  */
/*  list of conditions and the following disclaimer in the documentation and/or       */
/*  other materials provided with the distribution.                                   */
/*                                                                                    */
/*  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND   */
/*  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED     */
/*  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE            */
/*  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR  */
/*  ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES    */
/*  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;      */
/*  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON    */
/*  ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT           */
/*  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS     */
/*  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                      */
/*                                                                                    */
/*************************************************************************************/

#include "BaseDemo.hpp"
#include <vlVolume/RaycastVolume.hpp>
#include <vlVolume/VolumeUtils.hpp>
#include <vlGraphics/Light.hpp>
#include <vlGraphics/Text.hpp>
#include <vlGraphics/FontManager.hpp>
#include <vlGraphics/GLSL.hpp>
#include <vlGraphics/GeometryPrimitives.hpp>
#include <vlCore/Object.hpp>
#include <stdlib.h>
#include <time.h> 
#include <math.h>
#include <string.h>
#include <string>
#define DRND(x) ((double)(x)/RAND_MAX*rand())




#define NDX 100 //64
#define NDY 100 //128
#define NDZ 100 //64



#define N 16
#define GNP 16
#define INX 400						% Drawing window 1 Pixel size of side x
#define INY 800						% Drawing window 1 Pixel size of side y
#define INZ 400						% Drawing window 1 Pixel size of side z

int ndx = NDX, ndxm = NDX - 1;
int ndy = NDY, ndym = NDY - 1;
int ndz = NDZ, ndzm = NDZ - 1;
int nm = N - 1, nmm = N - 2, GN = GNP - 1;
double PI = 3.141592, time1;
double RR = 8.3145;
int totalPts = 0;
double ph[N][NDX][NDY][NDZ];		//% Ph [1] [i] [j] to ph [n] [i] [j]: pf in position [i] [j]
int qh[N][NDX][NDY][NDZ];	 		//% Qh [1] [i] [j] to qh [n] [i] [j]: grain number in position [i] [j]
int n00[NDX][NDY][NDZ];				//% The number of cases where ph is not 0 in position [i] [j]
int n00p[NDX][NDY][NDZ];            //% The number of cases when ph in position [i] [j] and ph around it are not 0
int is[NDX][NDY][NDZ];
double ph2[N][NDX][NDY][NDZ];
int qh2[N][NDX][NDY][NDZ];
//int inside[NDZ][NDX][NDY];
int inside[NDX][NDY][NDZ];

int intermediate[NDX / 2][NDY / 2][NDZ / 2];
//double boxCoordinates[NDX][NDY][NDZ][3];
int seeds[NDX][NDY][NDZ];
int centroid[N][3];
int seedVolume[N];
int lloyd = 1;

double aij[GNP][GNP], wij[GNP][GNP], tij[GNP][GNP], eij[GNP][GNP];


int vol[NDX][NDY][NDZ];

double time1max;
double delt, L, b1, t, s;
double dG, M1, W1, K1, E1;
double mu_chem, mu_grad;
double temp;
double sum1, sum2, sum3, sxs;
double pddtt;

double gamma, delta, amobi;
double aaa, vm0;
using namespace vl;

/*
  This example demonstrates how to implement raycasting using the vl::RaycastVolume class and a few ad-hoc fragment shaders.

  The user has several controls as well:

  Left/Right arrow keys changes the raycasting mode:
  - Isosurface_Mode
  - Isosurface_Transp_Mode,
  - MIP_Mode
  - RaycastBrightnessControl_Mode
  - RaycastDensityControl_Mode
  - RaycastColorControl_Mode

  Mouse wheel:
  - In Isosurface_Mode controls the iso-value of the isosurface
  - In Isosurface_Transp_Mode controls the iso-value of the isosurface
  - In MIP_Mode all the volume values less than this are discarded
  - In RaycastBrightnessControl_Mode controls the brightness of the voxels
  - In RaycastDensityControl_Mode controls the density of the voxels
  - In RaycastColorControl_Mode controls the color-bias of the voxels

  The Up/Down arrow keys are used to higher/lower the ray-advancement precision.

  The 'L' key toggles the dynamic and colored lights.
*/



class App_VolumeRaycast : public BaseDemo
{
    /* ----- raycast volume rendering options ----- */

    /* The sample step used to render the volume, the smaller the number
       the  better ( and slower ) the rendering will be. */
    float SAMPLE_STEP;

    /* volume visualization mode */
    enum RaycastMode {
        Isosurface_Mode,
        Isosurface_Transp_Mode,
        MIP_Mode,
        RaycastBrightnessControl_Mode,
        RaycastDensityControl_Mode,
        RaycastColorControl_Mode
    } MODE;

    /* If enabled, renders the volume using 3 animated lights. */
    bool DYNAMIC_LIGHTS;

    /* If enabled 3 colored lights are used to render the volume. */
    bool COLORED_LIGHTS;

    /* Use a separate 3d texture with a precomputed gradient to speedup the fragment shader.
       Requires more memory ( for the gradient texture ) but can speedup the rendering. */
    bool PRECOMPUTE_GRADIENT;

public:
    virtual String appletInfo()
    {
        return BaseDemo::appletInfo() +
            "- Left/Right Arrow: change raycast technique.\n" +
            "- Up/Down Arrow: changes SAMPLE_STEP.\n" +
            "- L: toggles lights (useful only for isosurface).\n" +
            "- Mouse Wheel: change the bias used to render the volume.\n" +
            "\n" +
            "- Drop inside the window a set of 2D files or a DDS or DAT volume to display it.\n" +
            "\n";
    }

    App_VolumeRaycast()
    {
        SAMPLE_STEP = 512.0f;
        MODE = Isosurface_Mode;
        DYNAMIC_LIGHTS = false;
        COLORED_LIGHTS = false;
        PRECOMPUTE_GRADIENT = false;
    }




    /* initialize the applet with a default volume */
    virtual void initEvent()
    {
        // loading meshin
        Log::debug(Say("Init Event called \n"));

        FILE* datin0;
        int n0p, n1, n2, n3;
        int i, j, k, l, it;
        int ii, jj, kk, kkk, iGN;
        int ip, im, jp, jm, kp, km;
        double sum1, t;
        printf("India");
        double r;
        int x1, y1, z1, r0;

        delt = 2.0;

        temp = 1000.0;
        L = 500.0;
        b1 = L / (double)NDX * 1.0e-9;

        vm0 = 7.0e-6;
        gamma = 0.5 * vm0 / RR / temp / b1;
        delta = 3.0;

        aaa = 2.0 / PI * sqrt(2.0 * delta * gamma);
        K1 = aaa * aaa;
        W1 = 4.0 * gamma / delta;
        amobi = 1.;
        M1 = amobi * PI * PI / (8.0 * delta);
        E1 = 500.0 / RR / temp;
        //E1=0.0;

        time1 = 0.;  time1max = 1.0e+20 + 1.0;

        for (i = 1; i <= GN; i++) {
            for (j = 1; j <= GN; j++) {
                wij[i][j] = W1;
                aij[i][j] = K1;
                tij[i][j] = M1;
                eij[i][j] = 0.0;
                if ((i == GN) || (j == GN)) { eij[i][j] = E1; }
                if (i > j) { eij[i][j] = -eij[i][j]; }
                if (i == j) { wij[i][j] = 0.0;  aij[i][j] = 0.0;  tij[i][j] = 0.0; eij[i][j] = 0.0; }
            }
        }


        //datin0 = fopen("C:\\Users\\dk18\\Work\\hiwi\\files\\6July\\eightparam\\eighty.txt", "r");
        //datin0 = fopen("C:\\Users\\dk18\\Work\\hiwi\\files\\dancer\\dancers32mod.txt", "r");
        //datin0 = fopen("C:\\Users\\dk18\\Work\\hiwi\\files\\pumpkin32\\pumpkin32.txt", "r");
        //datin0 = fopen("C:\\Users\\dk18\\Work\\hiwi\\files\\mushroom100\\mushroom100.txt", "r");
        datin0 = fopen("C:\\Users\\dk18\\Work\\hiwi\\files\\tweety\\tweety100half.txt", "r");
        int x, y, z;
        x = NDX;
        y = NDY;
        z = NDZ;
        // For eightparam 
        //z = 64;
        int factor = 1;
        for (i = 0; i < NDX / factor; i++) {
            for (j = 0; j < NDY / factor; j++) {
                for (k = 0; k < NDZ / factor; k++) {
                    fscanf(datin0, "%d  ", &inside[i][j][k]);
                    //fscanf(datin0, "%d  ", &intermediate[i][j][k]);
                    //inside[2 * i + 1][2 * j + 1][2 * k + 1] = inside[2 * i][2 * j][2 * k];
                }
            }
        }
        /*
        for (i = 0; i < NDX; i++) {
            for (j = 0; j < NDY; j++) {
                for (k = 0; k < NDZ; k++) {

                    inside[i][j][k] = intermediate[ i/factor][ j/factor][ k/factor];
                }
            }
        }
        */

        // TO get real coordinates uncomment this
        /*
        for (i = 0; i < NDX; i++) {
            for (j = 0; j < NDY; j++) {
                for (k = 0; k < NDZ; k++) {
                    fscanf(datin0, "%lf %lf %lf", &boxCoordinates[i][j][k][0], &boxCoordinates[i][j][k][1], &boxCoordinates[i][j][k][2]);
                }
            }
        }
        */
        fclose(datin0);


        // mesh loaded


        // ------------  added code 



            // Creating tape like torus
            /*
            int startX = 0;
            int endX = NDX;
            int startY = 0;
            int endY = NDY;
            int cX, cY, cZ;
            cX = NDX / 2;
            cY = NDY / 2;
            cZ = NDZ / 2;
            int radiusInner = 20;
            int radiusOuter = 30;
            for (k = 15; k >= 0; k--) {
                for (i = startX; i < endX; i++) {
                    for (j = startY; j < endY; j++) {
                        r = sqrt((double(i - cX)) * (double(i - cX))
                            + (double(j - cY)) * (double(j - cY)));
                        if (r < radiusOuter && r> radiusInner) {
                            inside[NDZ/2 - k][i][j] = 1;
                            inside[NDZ/2 + k][i][j] = 1;
                        }
                    }
                }
                //radiusInner--;
                //radiusOuter++;
            }
            */

            // Creating Torus
            /*
            int startX = 0;
            int endX = NDX;
            int startY = 0;
            int endY = NDY;
            int cX, cY, cZ;
            cX = NDX / 2;
            cY = NDY / 2;
            cZ = NDZ / 2;
            int radiusInner = 20;
            int radiusOuter = 35;
            for (k = 25; k >= 0; k--) {
                for (i = startX; i < endX; i++) {
                    for (j = startY; j < endY; j++) {
                        r = sqrt((double(i - cX)) * (double(i - cX))
                            + (double(j - cY)) * (double(j - cY)));
                        if (r < radiusOuter && r> radiusInner) {
                            inside[NDZ/2 - k][i][j] = 1;
                            inside[NDZ/2 + k][i][j] = 1;
                        }
                    }
                }
                radiusInner--;
                radiusOuter++;
            }
            */

            // Creating Sphere
            /*
            int startX = 0;
            int endX = NDX;
            int startY = 0;
            int endY = NDY;
            int cX, cY, cZ;
            cX = NDX / 2;
            cY = NDY / 2;
            cZ = NDZ / 2;
            int radius = NDX / 2.2;
            for (k = 0; k < NDZ; k++) {
                for (i = startX; i < endX; i++) {
                    for (j = startY; j < endY; j++) {
                        r = sqrt((double(i - cX)) * (double(i - cX))
                            + (double(j - cY)) * (double(j - cY))
                            + (double(k - cZ)) * (double(k - cZ)));
                        if (r < radius) {
                            inside[k][i][j] = 1;
                            totalPts++;
                        }
                    }
                }
            }
            */

            //Creating Cylinder
            /*
            int startX = 0;
            int endX = NDX;
            int startY = 0;
            int endY = NDY;
            int cX, cY, cZ;
            cX = NDX / 2;
            cY = NDY / 2;
            cZ = NDZ / 2;
            int radius = NDX / 3;
            for (k = 0; k < NDZ; k++) {
                for (i = startX; i < endX; i++) {
                    for (j = startY; j < endY; j++) {
                        r = sqrt((double(i - cX)) * (double(i - cX))
                            + (double(j - cY)) * (double(j - cY)));
                        if (r < radius) {
                            inside[k][i][j] = 1;
                            totalPts++;
                        }
                    }
                }
            }
            */

            //Creating a cube
            /*
            int startX = 1;
            int endX = NDX-1;
            int startY = 1;
            int endY = NDY-1;
            for (k = 0; k < NDZ; k++) {
                for (i = startX; i < endX; i++) {
                    for (j = startY; j < endY; j++) {
                        inside[k][i][j] = 1;
                    }
                }
            }
            */

            // Creating Diamond
            /*

            int startX = 0;
            int endX = NDX;
            int startY = 0;
            int endY = NDY;
            for (k = NDZ/2; k >= 0; k--) {
                if (startX > endX) {
                    break;
                }
                for (i = startX; i < endX; i++) {
                    for (j = startY; j < endY; j++) {
                        inside[k][i][j] = 1;
                    }
                }
                startX++;
                startY++;
                endX--;
                endY--;
            }

            startX = 0;
            endX = NDX;
            startY = 0;
            endY = NDY;
            for (k = NDZ/2; k <NDZ; k++) {
                if (startX > endX) {
                    break;
                }
                for (i = startX; i < endX; i++) {
                    for (j = startY; j < endY; j++) {
                        inside[k][i][j] = 1;
                    }
                }
                startX++;
                startY++;
                endX--;
                endY--;
            }

            */

            // Creating up down pyramid
            /*
            for (k = 0; k < NDZ; k++) {
                if (startX > endX) {
                    break;
                }
                for (i = startX; i < endX; i++) {
                    for (j = startY; j < endY; j++) {
                        inside[k][i][j] = 1;
                    }
                }
                startX++;
                startY++;
                endX--;
                endY--;
            }

            startX = 0;
            endX = NDX;
            startY = 0;
            endY = NDY;
            for (k = NDZ - 1; k >= 0; k--) {
                if (startX > endX) {
                    break;
                }
                for (i = startX; i < endX; i++) {
                    for (j = startY; j < endY; j++) {
                        inside[k][i][j] = 1;
                    }
                }
                startX++;
                startY++;
                endX--;
                endY--;
            }
            */

        for (i = 0; i <= ndxm; i++) {
            for (j = 0; j <= ndym; j++) {
                for (k = 0; k <= ndzm; k++) {
                    if (inside[i][j][k] == 1) {
                        ph[1][i][j][k] = 1.0;  //qh[1][i][j][k]=GN;
                        n00[i][j][k] = 1;      n00p[i][j][k] = 1;
                        is[i][j][k] = GN;
                    }
                }
            }
        }
        r0 = 2;
        //%Write a sphere with a radius 5 with different grain number (grain number written in circle)
        for (iGN = 1; iGN <= GN - 1; iGN++) {
            x1 = ndx * DRND(1); y1 = ndy * DRND(1); z1 = ndz * DRND(1);
            while (inside[x1][y1][z1] != 1) { x1 = ndx * DRND(1); y1 = ndy * DRND(1); z1 = ndz * DRND(1); }
            for (i = -r0; i <= (ndxm + r0); i++) {
                ii = i;
                //if (i < 0) { ii = ndx + i; }  if (i > ndxm) { ii = i - ndx; }
                if (i < 0) { ii = 0; } if (i > ndxm) { ii = ndxm; }
                for (j = -r0; j <= (ndym + r0); j++) {
                    jj = j;
                    //if (j < 0) { jj = ndy + j; }  if (j > ndym) { jj = j - ndy; }
                    if (j < 0) { jj = 0; } if (j > ndym) { jj = ndym; }
                    for (k = -r0; k <= (ndzm + r0); k++) {
                        kk = k;
                        //if (k < 0) { kk = ndz + k; }  if (k > ndzm) { kk = k - ndz; }
                        if (k < 0) { kk = 0; } if (k > ndzm) { kk = ndzm; }
                        r = sqrt((double(i - x1)) * (double(i - x1))
                            + (double(j - y1)) * (double(j - y1))
                            + (double(k - z1)) * (double(k - z1)));
                        if (r <= r0 && inside[ii][jj][kk] == 1) {
                            if (!(i<0 || j<0 || k <0 || i >ndxm || j> ndym || k>ndzm)) {
                                is[ii][jj][kk] = iGN;
                            }
                        }
                    }
                }
            }
        }

        /*
        r0 = 5.0; x1 = ndx / 2; y1 = ndy / 2; z1 = ndz / 2;
        for (i = 0; i <= ndxm; i++) {
            for (j = 0; j <= ndym; j++) {
                for (k = 0; k <= ndzm; k++) {
                    r = sqrt((double(i - x1)) * (double(i - x1))
                        + (double(j - y1)) * (double(j - y1))
                        + (double(k - z1)) * (double(k - z1)));
                    //if(r<=r0){ is[i][j][k]=1; }
                }
            }
        }
        */

        for (i = 0; i <= ndxm; i++) {
            for (j = 0; j <= ndym; j++) {
                for (k = 0; k <= ndzm; k++) {
                    if (inside[i][j][k] == 1) {

                        qh[1][i][j][k] = is[i][j][k];
                    }
                }
            }
        }

        //--- ph=0 Special case  (if there are other orientations around) added -------------------------
        for (i = 0; i <= ndxm; i++) {
            for (j = 0; j <= ndym; j++) {
                for (k = 0; k <= ndzm; k++) {

                    ip = i + 1; im = i - 1;  jp = j + 1;  jm = j - 1;  kp = k + 1;  km = k - 1;
                    /*
                    if (i == ndxm) { ip = 0; }  if (i == 0) { im = ndxm; }
                    if (j == ndym) { jp = 0; }  if (j == 0) { jm = ndym; }
                    if (k == ndzm) { kp = 0; }  if (k == 0) { km = ndzm; }
                    */

                    if (i == ndxm) { ip = ndxm; }  if (i == 0) { im = 0; }
                    if (j == ndym) { jp = ndym; }  if (j == 0) { jm = 0; }
                    if (k == ndzm) { kp = ndzm; }  if (k == 0) { km = 0; }

                    for (kk = 1; kk <= n00[ip][j][k]; kk++) {
                        kkk = qh[kk][ip][j][k]; //The kk-th  in position [ip] [j] [k]
                        ii = 1; for (l = 1; l <= n00p[i][j][k]; l++) { if (qh[l][i][j][k] == kkk) { ii = 0; } }
                        //% If direction kkk lies in the direction of position [i] [j] [k], ii = 0, if none, ii = 1
                        if (ii == 1) {
                            n00p[i][j][k] += 1;
                            ph[n00p[i][j][k]][i][j][k] = 0.0; qh[n00p[i][j][k]][i][j][k] = kkk;
                        }
                    }

                    for (kk = 1; kk <= n00[im][j][k]; kk++) {
                        kkk = qh[kk][im][j][k];
                        ii = 1; for (l = 1; l <= n00p[i][j][k]; l++) { if (qh[l][i][j][k] == kkk) { ii = 0; } }
                        if (ii == 1) {
                            n00p[i][j][k] += 1;
                            ph[n00p[i][j][k]][i][j][k] = 0.0; qh[n00p[i][j][k]][i][j][k] = kkk;
                        }
                    }

                    for (kk = 1; kk <= n00[i][jp][k]; kk++) {
                        kkk = qh[kk][i][jp][k];
                        ii = 1; for (l = 1; l <= n00p[i][j][k]; l++) { if (qh[l][i][j][k] == kkk) { ii = 0; } }
                        if (ii == 1) {
                            n00p[i][j][k] += 1;
                            ph[n00p[i][j][k]][i][j][k] = 0.0; qh[n00p[i][j][k]][i][j][k] = kkk;
                        }
                    }

                    for (kk = 1; kk <= n00[i][jm][k]; kk++) {
                        kkk = qh[kk][i][jm][k];
                        ii = 1; for (l = 1; l <= n00p[i][j][k]; l++) { if (qh[l][i][j][k] == kkk) { ii = 0; } }
                        if (ii == 1) {
                            n00p[i][j][k] += 1;
                            ph[n00p[i][j][k]][i][j][k] = 0.0; qh[n00p[i][j][k]][i][j][k] = kkk;
                        }
                    }

                    for (kk = 1; kk <= n00[i][j][kp]; kk++) {
                        kkk = qh[kk][i][j][kp];
                        ii = 1; for (l = 1; l <= n00p[i][j][k]; l++) { if (qh[l][i][j][k] == kkk) { ii = 0; } }
                        if (ii == 1) {
                            n00p[i][j][k] += 1;
                            ph[n00p[i][j][k]][i][j][k] = 0.0; qh[n00p[i][j][k]][i][j][k] = kkk;
                        }
                    }

                    for (kk = 1; kk <= n00[i][j][km]; kk++) {
                        kkk = qh[kk][i][j][km];
                        ii = 1; for (l = 1; l <= n00p[i][j][k]; l++) { if (qh[l][i][j][k] == kkk) { ii = 0; } }
                        if (ii == 1) {
                            n00p[i][j][k] += 1;
                            ph[n00p[i][j][k]][i][j][k] = 0.0; qh[n00p[i][j][k]][i][j][k] = kkk;
                        }
                    }

                }
            }
        }

        // -----------------

        vl::Log::notify(appletInfo());

        if (!Has_GLSL)
        {
            vl::Log::error("OpenGL Shading Language not supported.\n");
            vl::Time::sleep(5000);
            exit(1);
        }

        mLight0 = new Light;
        mLight1 = new Light;
        mLight2 = new Light;

        mLight0Tr = new Transform;
        mLight1Tr = new Transform;
        mLight2Tr = new Transform;
        rendering()->as<Rendering>()->transform()->addChild(mLight0Tr.get());
        rendering()->as<Rendering>()->transform()->addChild(mLight1Tr.get());
        rendering()->as<Rendering>()->transform()->addChild(mLight2Tr.get());

        // volume transform
        mVolumeTr = new Transform;

        // val_threshold: manipulated via mouse wheel
        // - In Isosurface_Mode controls the iso-value of the isosurface
        // - In Isosurface_Transp_Mode controls the iso-value of the isosurface
        // - In MIP_Mode all the volume values less than this are discarded
        // - In RaycastBrightnessControl_Mode controls the brightness of the voxels
        // - In RaycastDensityControl_Mode controls the density of the voxels
        // - In RaycastColorControl_Mode controls the color-bias of the voxels
        mValThreshold = new Uniform("val_threshold");
        mValThreshold->setUniformF(0.5f);

        // default volume image
        //FILE* datin0;
        //mVolumeImage = loadImage( "/volume/VLTest.dat" );
        Log::debug(Say("Opening file realCordinates100.txt: \n"));
        int w, h, d;
        w = NDX;
        h = w;
        printf("data loaded");
        d = w;
        mVolumeImage = new Image(w, h, d, 1, IF_RGBA, IT_UNSIGNED_BYTE);

        datin0 = fopen("C:\\Users\\dk18\\Work\\hiwi\\files\\mushroom100\\txtFilesSeed\\realCordinates1125.txt", "r");
        Log::debug(Say("-------------Opened file realCordinates1000.txt: \n"));



        for (int i = 0; i < NDX; i++) {
            for (int j = 0; j < NDY; j++) {
                for (int k = 0; k < NDZ; k++) {
                    fscanf(datin0, "%d ", &vol[i][j][k]);
                }
            }
        }

        fclose(datin0);

        Log::debug(Say("Closing file realCordinates100.txt: \n"));

        Log::debug(Say("Volume: %n %n %n\n") << mVolumeImage->width() << mVolumeImage->height() << mVolumeImage->depth());

        ubvec4* rgba_px = (ubvec4*)mVolumeImage->pixels();


        for (int z = 0; z < d; ++z) {
            for (int y = 0; y < h; ++y) {
                for (int x = 0; x < w; ++x, ++rgba_px) {
                    rgba_px->r() = (unsigned char)vol[x][y][z] * 100.0f;
                    rgba_px->g() = (unsigned char)vol[x][y][z] * 100.0f;
                    rgba_px->b() = (unsigned char)vol[x][y][z] * 100.0f;
                    rgba_px->a() = (unsigned char)(200 * 1.0f);
                }
            }
        }


        Log::debug(Say("Updated rgba_px values  \n"));

        setupScene();

        vl::Time::sleep(2000);

        // trying to do computation here and calling setupScene everytime 



    }

    void setupScene()
    {

        Log::debug(Say("Inside setup scene  \n"));
        mRaycastVolume = new vl::RaycastVolume;
        mVolumeAct = new Actor;

        // scrap previous scene
        sceneManager()->tree()->eraseAllChildren();
        sceneManager()->tree()->actors()->clear();
        mLight0->bindTransform(NULL);
        mLight1->bindTransform(NULL);
        mLight2->bindTransform(NULL);

        ref<Effect> volume_fx = new Effect;
        // we don't necessarily need this:
        // volume_fx->shader()->enable( EN_BLEND );
        volume_fx->shader()->enable(EN_CULL_FACE);
        volume_fx->shader()->enable(EN_DEPTH_TEST);

        // NOTE
        // in these cases we render the back faces and raycast in back to front direction
        // in the other cases we render the front faces and raycast in front to back direction
        if (MODE == RaycastBrightnessControl_Mode || MODE == RaycastDensityControl_Mode || MODE == RaycastColorControl_Mode)
        {
            volume_fx->shader()->enable(vl::EN_CULL_FACE);
            volume_fx->shader()->gocCullFace()->set(vl::PF_FRONT);
        }

        mRaycastVolume->lights().push_back(mLight0);

        // light bulbs
        if (DYNAMIC_LIGHTS)
        {
            // you can color the lights!
            if (COLORED_LIGHTS)
            {
                mLight0->setAmbient(fvec4(0.1f, 0.1f, 0.1f, 1.0f));
                mLight1->setAmbient(fvec4(0.1f, 0.1f, 0.1f, 1.0f));
                mLight2->setAmbient(fvec4(0.1f, 0.1f, 0.1f, 1.0f));
                mLight0->setDiffuse(vl::gold);
                mLight1->setDiffuse(vl::green);
                mLight2->setDiffuse(vl::royalblue);
            }

            // add the other two lights
            mRaycastVolume->lights().push_back(mLight1);
            mRaycastVolume->lights().push_back(mLight2);

            // animate the three lights
            mLight0->bindTransform(mLight0Tr.get());
            mLight1->bindTransform(mLight1Tr.get());
            mLight2->bindTransform(mLight2Tr.get());

            // add also a light bulb actor
            ref<Effect> fx_bulb = new Effect;
            fx_bulb->shader()->enable(EN_DEPTH_TEST);
            ref<Geometry> light_bulb = vl::makeIcosphere(vec3(0, 0, 0), 1, 1);
            sceneManager()->tree()->addActor(light_bulb.get(), fx_bulb.get(), mLight0Tr.get());
            sceneManager()->tree()->addActor(light_bulb.get(), fx_bulb.get(), mLight1Tr.get());
            sceneManager()->tree()->addActor(light_bulb.get(), fx_bulb.get(), mLight2Tr.get());
        }

        // the GLSL program that performs the actual raycasting
        mGLSL = volume_fx->shader()->gocGLSLProgram();
        mGLSL->gocUniform("sample_step")->setUniformF(1.0f / SAMPLE_STEP);

        // attach vertex shader (common to all the raycasting techniques)
        mGLSL->attachShader(new GLSLVertexShader("/glsl/volume_luminance_light.vs"));

        // attach fragment shader implementing the specific raycasting tecnique
        if (MODE == Isosurface_Mode)
            mGLSL->attachShader(new GLSLFragmentShader("/glsl/volume_raycast_isosurface.fs"));
        else
            if (MODE == Isosurface_Transp_Mode)
                mGLSL->attachShader(new GLSLFragmentShader("/glsl/volume_raycast_isosurface_transp.fs"));
            else
                if (MODE == MIP_Mode)
                    mGLSL->attachShader(new GLSLFragmentShader("/glsl/volume_raycast_mip.fs"));
                else
                    if (MODE == RaycastBrightnessControl_Mode)
                        mGLSL->attachShader(new GLSLFragmentShader("/glsl/volume_raycast01.fs"));
                    else
                        if (MODE == RaycastDensityControl_Mode)
                            mGLSL->attachShader(new GLSLFragmentShader("/glsl/volume_raycast02.fs"));
                        else
                            if (MODE == RaycastColorControl_Mode)
                                mGLSL->attachShader(new GLSLFragmentShader("/glsl/volume_raycast03.fs"));

        // manipulate volume transform with the trackball
        trackball()->setTransform(mVolumeTr.get());

        // volume actor
        mVolumeAct->setEffect(volume_fx.get());
        mVolumeAct->setTransform(mVolumeTr.get());
        sceneManager()->tree()->addActor(mVolumeAct.get());
        // bind val_threshold uniform to the volume actor
        mVolumeAct->setUniform(mValThreshold.get());

        // RaycastVolume will generate the actual actor's geometry upon setBox() invocation.
        // The geometry generated is actually a simple box with 3D texture coordinates.
        mRaycastVolume->bindActor(mVolumeAct.get());
        AABB volume_box(vec3(-10, -10, -10), vec3(+10, +10, +10));
        mRaycastVolume->setBox(volume_box);

        // val_threshold text
        mValThresholdText = new Text;
        mValThresholdText->setFont(defFontManager()->acquireFont("/font/bitstream-vera/VeraMono.ttf", 12));
        mValThresholdText->setTextAlignment(TextAlignCenter);
        mValThresholdText->setAlignment(AlignHCenter | AlignBottom);
        mValThresholdText->setViewportAlignment(AlignHCenter | AlignBottom);
        mValThresholdText->translate(0, 5, 0);
        mValThresholdText->setBackgroundEnabled(true);
        mValThresholdText->setBackgroundColor(fvec4(0, 0, 0, 0.75));
        mValThresholdText->setColor(vl::white);
        ref<Effect> effect = new Effect;
        effect->shader()->enable(EN_BLEND);
        updateText();
        // Don't add text in GL core profile - not supported yet.
        if (Has_Fixed_Function_Pipeline) {
            sceneManager()->tree()->addActor(mValThresholdText.get(), effect.get());
        }

        // let's visualize the volume!

        Log::debug(Say("Closing setup scene  \n"));
        setupVolume();
    }

    /* visualize the given volume */
    void setupVolume()
    {


        //Log::debug(Say("Inside setup volume  \n"));
        Effect* volume_fx = mVolumeAct->effect();

        volume_fx->shader()->enable(EN_BLEND);

        //Log::debug(Say("shader enabled  \n"));
        // for semplicity we don't distinguish between different image formats, i.e. IF_LUMINANCE, IF_RGBA etc.

        ref<Image> gradient;
        if (PRECOMPUTE_GRADIENT)
        {

            Log::debug(Say("gradient starting  \n"));
            // note that this can take a while...
            gradient = vl::genGradientNormals(mVolumeImage.get());

            Log::debug(Say("gradient computed  \n"));
        }

        Log::debug(Say("Printing volume  \n"));
        //Log::debug( Say("Volume: %n %n %n\n") << mVolumeImage->width() << mVolumeImage->height() << mVolumeImage->depth() );

        // volume image textue must be on sampler #0
        //mVolumeImage->dimension() = ID_3D;
        vl::ref< vl::Texture > vol_tex = new vl::Texture(mVolumeImage.get(), TF_RED, false, false);
        volume_fx->shader()->gocTextureImageUnit(0)->setTexture(vol_tex.get());
        vol_tex->getTexParameter()->setMagFilter(vl::TPF_LINEAR);
        vol_tex->getTexParameter()->setMinFilter(vl::TPF_LINEAR);
        vol_tex->getTexParameter()->setWrap(vl::TPW_CLAMP_TO_EDGE);
        volume_fx->shader()->gocUniform("volume_texunit")->setUniformI(0);
        mRaycastVolume->generateTextureCoordinates(ivec3(mVolumeImage->width(), mVolumeImage->height(), mVolumeImage->depth()));

        Log::debug(Say("texture cordinates generated  \n"));
        // generate a simple colored transfer function
        ref<Image> trfunc;
        if (COLORED_LIGHTS && DYNAMIC_LIGHTS) {
            trfunc = vl::makeColorSpectrum(128, vl::white, vl::white); // let the lights color the volume
        }
        else {
            trfunc = vl::makeColorSpectrum(128, vl::blue, vl::royalblue, vl::green, vl::yellow, vl::crimson);
        }

        // installs the transfer function as texture #1
        vl::ref< vl::Texture > trf_tex = new Texture(trfunc.get(), vl::TF_RGBA, false, false);
        trf_tex->getTexParameter()->setMagFilter(vl::TPF_LINEAR);
        trf_tex->getTexParameter()->setMinFilter(vl::TPF_LINEAR);
        trf_tex->getTexParameter()->setWrap(vl::TPW_CLAMP_TO_EDGE);
        volume_fx->shader()->gocTextureImageUnit(1)->setTexture(trf_tex.get());
        volume_fx->shader()->gocUniform("trfunc_texunit")->setUniformI(1);

        // gradient computation, only use for isosurface methods
        if (MODE == Isosurface_Mode || MODE == Isosurface_Transp_Mode)
        {
            volume_fx->shader()->gocUniform("precomputed_gradient")->setUniformI(PRECOMPUTE_GRADIENT ? 1 : 0);
            if (PRECOMPUTE_GRADIENT)
            {
                volume_fx->shader()->gocUniform("gradient_texunit")->setUniformI(2);
                vl::ref< vl::Texture > tex = new Texture(gradient.get(), TF_RGBA, false, false);
                tex->getTexParameter()->setMagFilter(vl::TPF_LINEAR);
                tex->getTexParameter()->setMinFilter(vl::TPF_LINEAR);
                tex->getTexParameter()->setWrap(vl::TPW_CLAMP_TO_EDGE);
                volume_fx->shader()->gocTextureImageUnit(2)->setTexture(tex.get());

                Log::debug(Say("Inside isosurface mode  \n"));
            }
        }

        //Log::debug(Say("closing setup volume  \n"));
        // update text
        updateText();

        // refresh window
        openglContext()->update();

    }

    /* load files drag&dropped in the window */
    void fileDroppedEvent(const std::vector<String>& files)
    {
        mVolumeImage = NULL;

        if (files.size() == 1) // if there is one file load it directly
        {
            if (files[0].endsWith(".dat") || files[0].endsWith(".dds") || files[0].endsWith(".mhd"))
            {
                mVolumeImage = loadImage(files[0]);

                // CT volumes are often in Hounsfield units saved as signed 16 bits ints ranging from -1000 to +3000.
                // You can keep the values as they are but keep in mind that when OpenGL converts that image to a texture
                // the values are mapped to a 0 ... 1 range following the OpenGL rules defined in the manual, so to render
                // them properly you'll need to do an appropriate scaling operation in the GLSL shader.
                // In the example below we prefer to scale/contrast the image values directly so we can reuse the usual
                // shaders.

                // If image format is SHORT we assume it contains Hounsfield units so we rescale its values for
                // optimal visibility. Note that you can do this in the shader as well.
                if (mVolumeImage->type() == vl::IT_SHORT) {
                    mVolumeImage->setType(vl::IT_UNSIGNED_SHORT);
                    float h_start = -1000;
                    float h_end = +1500;
                    vl::i16* ival = (vl::i16*)mVolumeImage->pixels();
                    vl::u16* uval = (vl::u16*)mVolumeImage->pixels();
                    for (int i = 0; i < mVolumeImage->width() * mVolumeImage->height() * mVolumeImage->depth(); ++i, ++ival, ++uval) {
                        float t = ((*ival) - h_start) / (h_end - h_start);
                        t = vl::clamp(t, 0.0f, 1.0f);
                        *uval = (u16)(t * 65535);
                    }
                }

                if (mVolumeImage) {
                    setupVolume();
                }
            }
        }
        else // if there is more than one file load and sort them and assemble a 3D image
        {
            // sort files by their name
            std::vector<String> files_sorted = files;
            std::sort(files_sorted.begin(), files_sorted.end());
            // load the files
            std::vector< ref<Image> > images;
            for (unsigned int i = 0; i < files_sorted.size(); ++i)
            {
                images.push_back(loadImage(files_sorted[i]));
                if (files_sorted[i].endsWith(".dcm"))
                    images.back()->contrastHounsfieldAuto();
            }
            // assemble the volume
            mVolumeImage = assemble3DImage(images);
            // set the volume
            if (mVolumeImage)
                setupVolume();
        }

        if (!mVolumeImage)
            Log::error("Error loading volume data!\n");
    }

    void updateText()
    {
        // update the val_threshold value and the raycast technique name
        String technique_name;
        switch (MODE)
        {
        case Isosurface_Mode: technique_name = "raycast isosurface >"; break;
        case Isosurface_Transp_Mode: technique_name = "< raycast transparent isosurface >"; break;
        case MIP_Mode: technique_name = "< raycast maximum intensity projection >"; break;
        case RaycastBrightnessControl_Mode: technique_name = "< raycast brightness control >"; break;
        case RaycastDensityControl_Mode: technique_name = "< raycast density control >"; break;
        case RaycastColorControl_Mode: technique_name = "< raycast color control"; break;
        };

        float val_threshold = 0;
        mValThreshold->getUniform(&val_threshold);
        mValThresholdText->setText(Say("val_threshold = %n\n" "sample_step = 1.0 / %.0n\n" "%s") << val_threshold << SAMPLE_STEP << technique_name);
    }

    void updateValThreshold(int val)
    {
        float val_threshold = 0.0f;
        mValThreshold->getUniform(&val_threshold);
        val_threshold += val * 0.01f;
        val_threshold = clamp(val_threshold, 0.0f, 1.0f);
        mValThreshold->setUniformF(val_threshold);

        updateText();
        openglContext()->update();
    }

    void mouseWheelEvent(int val)
    {
        updateValThreshold(val);
    }

    /* animate the lights */
    virtual void updateScene()
    {
        if (DYNAMIC_LIGHTS)
        {
            mat4 mat;
            // light 0 transform.
            mat = mat4::getRotation(Time::currentTime() * 43, 0, 1, 0) * mat4::getTranslation(20, 20, 20);
            mLight0Tr->setLocalMatrix(mat);
            // light 1 transform.
            mat = mat4::getRotation(Time::currentTime() * 47, 0, 1, 0) * mat4::getTranslation(-20, 0, 0);
            mLight1Tr->setLocalMatrix(mat);
            // light 2 transform.
            mat = mat4::getRotation(Time::currentTime() * 47, 0, 1, 0) * mat4::getTranslation(+20, 0, 0);
            mLight2Tr->setLocalMatrix(mat);
        }
    }

    virtual void keyPressEvent(unsigned short, EKey key)
    {
        int i, j, k, l, n1, n2, n3, im, ip, it, jp, jm, kp, km, ii, jj, kk, kkk, n0p;
        int ctr = 0;

        double sum1, t;
        //printf("India");
        double r;
        int x1, y1, z1, r0;

        delt = 2.0;

        temp = 1000.0;
        L = 500.0;
        b1 = L / (double)NDX * 1.0e-9;

        vm0 = 7.0e-6;
        gamma = 0.5 * vm0 / RR / temp / b1;
        delta = 3.0;

        aaa = 2.0 / PI * sqrt(2.0 * delta * gamma);
        K1 = aaa * aaa;
        W1 = 4.0 * gamma / delta;
        amobi = 1.;
        M1 = amobi * PI * PI / (8.0 * delta);
        E1 = 500.0 / RR / temp;
        //E1=0.0;

        time1 = 0.;  time1max = 1.0e+20 + 1.0;

        for (i = 1; i <= GN; i++) {
            for (j = 1; j <= GN; j++) {
                wij[i][j] = W1;
                aij[i][j] = K1;
                tij[i][j] = M1;
                eij[i][j] = 0.0;
                if ((i == GN) || (j == GN)) { eij[i][j] = E1; }
                if (i > j) { eij[i][j] = -eij[i][j]; }
                if (i == j) { wij[i][j] = 0.0;  aij[i][j] = 0.0;  tij[i][j] = 0.0; eij[i][j] = 0.0; }
            }
        }
        ubvec4* rgba_px = (ubvec4*)mVolumeImage->pixels();
        RaycastMode modes[] = { Isosurface_Mode, Isosurface_Transp_Mode, MIP_Mode, RaycastBrightnessControl_Mode, RaycastDensityControl_Mode, RaycastColorControl_Mode };
        int mode = MODE;
        if (key == vl::Key_K) {
            Log::debug(Say("Right key pressed \n"));

            while (ctr < 5) {
                ctr++;

                //if ((((int)(time1) % 100) == 0)) { datsave(); }
                //datin();
                //datsave();
                //setupScene();

                for (i = 0; i <= ndxm; i++) {
                    for (j = 0; j <= ndym; j++) {
                        for (k = 0; k <= ndzm; k++) {
                            ip = i + 1; im = i - 1;  jp = j + 1;  jm = j - 1;  kp = k + 1;  km = k - 1;
                            // Picking neighbours and cyclic picking
                            /*
                            if (i == ndxm) { ip = 0; }  if (i == 0) { im = ndxm; }
                            if (j == ndym) { jp = 0; }  if (j == 0) { jm = ndym; }
                            if (k == ndzm) { kp = 0; }  if (k == 0) { km = ndzm; }
                            */

                            if (i == ndxm) { ip = ndxm; }  if (i == 0) { im = 0; }
                            if (j == ndym) { jp = ndym; }  if (j == 0) { jm = 0; }
                            if (k == ndzm) { kp = ndzm; }  if (k == 0) { km = 0; }

                            n0p = n00p[i][j][k];
                            if (n0p >= 2) {
                                for (n1 = 1; n1 <= n0p; n1++) {
                                    ii = qh[n1][i][j][k];  pddtt = 0.0;
                                    for (n2 = 1; n2 <= n0p; n2++) {
                                        jj = qh[n2][i][j][k];  sum1 = 0.0;
                                        for (n3 = 1; n3 <= n0p; n3++) {
                                            kk = qh[n3][i][j][k]; sum2 = 0.0;
                                            // Trying to calculate Laplacian here from neighbours
                                            for (l = 1; l <= n00[ip][j][k]; l++) { if (qh[l][ip][j][k] == kk) { sum2 += ph[l][ip][j][k]; } }
                                            for (l = 1; l <= n00[im][j][k]; l++) { if (qh[l][im][j][k] == kk) { sum2 += ph[l][im][j][k]; } }
                                            for (l = 1; l <= n00[i][jp][k]; l++) { if (qh[l][i][jp][k] == kk) { sum2 += ph[l][i][jp][k]; } }
                                            for (l = 1; l <= n00[i][jm][k]; l++) { if (qh[l][i][jm][k] == kk) { sum2 += ph[l][i][jm][k]; } }
                                            for (l = 1; l <= n00[i][j][kp]; l++) { if (qh[l][i][j][kp] == kk) { sum2 += ph[l][i][j][kp]; } }
                                            for (l = 1; l <= n00[i][j][km]; l++) { if (qh[l][i][j][km] == kk) { sum2 += ph[l][i][j][km]; } }
                                            // Formula on line 13 in paper in algorithm 1
                                            sum1 += 0.5 * (aij[ii][kk] - aij[jj][kk]) * (sum2 - 6.0 * ph[n3][i][j][k])
                                                + (wij[ii][kk] - wij[jj][kk]) * ph[n3][i][j][k];
                                        }
                                        // Formula on line 14 in paper in algorithm 1
                                        pddtt += -2.0 * tij[ii][jj] / double(n00p[i][j][k])
                                            * (sum1 - 8.0 / PI * eij[ii][jj] * sqrt(ph[n1][i][j][k] * ph[n2][i][j][k]));
                                    }
                                    if (inside[i][j][k] == 1) {
                                        // Formula on line 15 in paper in algorithm 1
                                        ph2[n1][i][j][k] = ph[n1][i][j][k] + pddtt * delt;
                                        // commeting qh2
                                        qh2[n1][i][j][k] = qh[n1][i][j][k];
                                        if (ph2[n1][i][j][k] >= 1.0) { ph2[n1][i][j][k] = 1.0; }
                                        if (ph2[n1][i][j][k] <= 0.0) { ph2[n1][i][j][k] = 0.0; }
                                        //if(ph2[n1][i][j][k]<=1.0e-4){ph2[n1][i][j][k]=0.0;}
                                    }
                                }
                            }
                            else {
                                if (inside[i][j][k] == 1) {
                                    ph2[1][i][j][k] = ph[1][i][j][k];
                                    // commeting qh2
                                    qh2[1][i][j][k] = qh[1][i][j][k];
                                }
                            }//if
                        }//k
                    }//j
                }//i

                //%--- ph2?ph, qh2?qh and ph2=0,qh2=0 -------------------------
                for (i = 0; i <= ndxm; i++) {
                    for (j = 0; j <= ndym; j++) {
                        for (k = 0; k <= ndzm; k++) {
                            if (inside[i][j][k] == 1) {
                                for (kk = 1; kk <= nm; kk++) {
                                    ph[kk][i][j][k] = ph2[kk][i][j][k];
                                    // commeting qh2
                                    qh[kk][i][j][k] = qh2[kk][i][j][k];
                                    ph2[kk][i][j][k] = 0.0;
                                }
                            }
                        }
                    }
                }

                //%--- Sort (descending) -------------------------
                for (i = 0; i <= ndxm; i++) {
                    for (j = 0; j <= ndym; j++) {
                        for (k = 0; k <= ndzm; k++) {
                            if (inside[i][j][k] == 1) {
                                for (kk = 1; kk <= nm - 2; kk++) { //Sort from kk = 1 to kk = nm in descending order (large to small) (kk = 0 is ignored)
                                    for (l = nm - 1; l > kk; l--) {
                                        if (ph[l][i][j][k] > ph[l - 1][i][j][k]) {
                                            t = ph[l][i][j][k];  ph[l][i][j][k] = ph[l - 1][i][j][k]; ph[l - 1][i][j][k] = t;
                                            it = qh[l][i][j][k]; qh[l][i][j][k] = qh[l - 1][i][j][k]; qh[l - 1][i][j][k] = it;
                                        }
                                    }
                                }

                                //%--- Standardization-------------------------
                                // Normalization line 20 in paper in algorithm 1
                                sum1 = 0.0; ii = 0;
                                for (kk = 1; kk <= nm; kk++) { if (ph[kk][i][j][k] > 0.0) { ii++; sum1 += ph[kk][i][j][k]; } }
                                n00[i][j][k] = ii; n00p[i][j][k] = ii;
                                for (kk = 1; kk <= n00[i][j][k]; kk++) { ph[kk][i][j][k] = ph[kk][i][j][k] / sum1; }
                                //for(kk=1;kk<=n00[i][j][k];kk++){ if(sum1>=1.0){ph[kk][i][j][k]=ph[kk][i][j][k]/sum1;} }
                            }
                        }
                    }
                }


                //%--- Added special case of ph = 0 (if there are other orientations around)-------------------------
                for (i = 0; i <= ndxm; i++) {
                    for (j = 0; j <= ndym; j++) {
                        for (k = 0; k <= ndzm; k++) {

                            ip = i + 1; im = i - 1;  jp = j + 1;  jm = j - 1;  kp = k + 1;  km = k - 1;
                            /*
                            if (i == ndxm) { ip = 0; }  if (i == 0) { im = ndxm; }
                            if (j == ndym) { jp = 0; }  if (j == 0) { jm = ndym; }
                            if (k == ndzm) { kp = 0; }  if (k == 0) { km = ndzm; }
                            */

                            if (i == ndxm) { ip = ndxm; }  if (i == 0) { im = 0; }
                            if (j == ndym) { jp = ndym; }  if (j == 0) { jm = 0; }
                            if (k == ndzm) { kp = ndzm; }  if (k == 0) { km = 0; }

                            for (kk = 1; kk <= n00[ip][j][k]; kk++) {
                                kkk = qh[kk][ip][j][k]; //%The kth grain number in position [ip] [j] [k]
                                ii = 1; for (l = 1; l <= n00p[i][j][k]; l++) { if (qh[l][i][j][k] == kkk) { ii = 0; } }
                                //%If direction kkk lies in the direction of position [i] [j] [k], ii = 0, if none, ii = 1
                                if (ii == 1) {
                                    n00p[i][j][k] += 1;
                                    ph[n00p[i][j][k]][i][j][k] = 0.0; qh[n00p[i][j][k]][i][j][k] = kkk;
                                }
                            }

                            for (kk = 1; kk <= n00[im][j][k]; kk++) {
                                kkk = qh[kk][im][j][k];
                                ii = 1; for (l = 1; l <= n00p[i][j][k]; l++) { if (qh[l][i][j][k] == kkk) { ii = 0; } }
                                if (ii == 1) {
                                    n00p[i][j][k] += 1;
                                    ph[n00p[i][j][k]][i][j][k] = 0.0; qh[n00p[i][j][k]][i][j][k] = kkk;
                                }
                            }

                            for (kk = 1; kk <= n00[i][jp][k]; kk++) {
                                kkk = qh[kk][i][jp][k];
                                ii = 1; for (l = 1; l <= n00p[i][j][k]; l++) { if (qh[l][i][j][k] == kkk) { ii = 0; } }
                                if (ii == 1) {
                                    n00p[i][j][k] += 1;
                                    ph[n00p[i][j][k]][i][j][k] = 0.0; qh[n00p[i][j][k]][i][j][k] = kkk;
                                }
                            }

                            for (kk = 1; kk <= n00[i][jm][k]; kk++) {
                                kkk = qh[kk][i][jm][k];
                                ii = 1; for (l = 1; l <= n00p[i][j][k]; l++) { if (qh[l][i][j][k] == kkk) { ii = 0; } }
                                if (ii == 1) {
                                    n00p[i][j][k] += 1;
                                    ph[n00p[i][j][k]][i][j][k] = 0.0; qh[n00p[i][j][k]][i][j][k] = kkk;
                                }
                            }

                            for (kk = 1; kk <= n00[i][j][kp]; kk++) {
                                kkk = qh[kk][i][j][kp];
                                ii = 1; for (l = 1; l <= n00p[i][j][k]; l++) { if (qh[l][i][j][k] == kkk) { ii = 0; } }
                                if (ii == 1) {
                                    n00p[i][j][k] += 1;
                                    ph[n00p[i][j][k]][i][j][k] = 0.0; qh[n00p[i][j][k]][i][j][k] = kkk;
                                }
                            }

                            for (kk = 1; kk <= n00[i][j][km]; kk++) {
                                kkk = qh[kk][i][j][km];
                                ii = 1; for (l = 1; l <= n00p[i][j][k]; l++) { if (qh[l][i][j][k] == kkk) { ii = 0; } }
                                if (ii == 1) {
                                    n00p[i][j][k] += 1;
                                    ph[n00p[i][j][k]][i][j][k] = 0.0; qh[n00p[i][j][k]][i][j][k] = kkk;
                                }
                            }

                        }
                    }
                }



                //setupVolume();
            }
            // visualizing data
            for (i = 0; i <= ndxm; i++) {
                for (j = 0; j <= ndym; j++) {
                    for (k = 0; k <= ndzm; k++, ++rgba_px) {
                        /*
                        fprintf(stream, "%d  \n", n00[i][j][k]);
                        for (kk = 1; kk <= n00[i][j][k]; kk++) {
                            fprintf(stream, "%d , %e  ", qh[kk][i][j][k], ph[kk][i][j][k]);
                        }*/
                        //fprintf(stream, "\n");
                        //float v = 0.23 * k;
                        double col = 0.; for (kk = 1; kk <= n00[i][j][k]; kk++) { col += ph[kk][i][j][k] * ph[kk][i][j][k]; }


                        double max = -1;
                        int seed = 0;
                        for (kk = 1; kk <= n00[i][j][k]; kk++) {
                            if (max < ph[kk][i][j][k]) {
                                max = ph[kk][i][j][k];
                                seed = qh[kk][i][j][k];
                                seeds[i][j][k] = seed;
                            }
                            //fprintf(stream, "%d  %e\n ", qh[kk][i][j][k], ph[kk][i][j][k]);
                        }
                        seedVolume[seed] ++;


                        // commeting out seed visualization
                        rgba_px->r() = (unsigned char)(col) * 1000.0f;
                        rgba_px->g() = (unsigned char)(col) * 1000.0f;
                        rgba_px->b() = (unsigned char)(col) * 1000.0f;
                        rgba_px->a() = (unsigned char)(200 * 1.0f);
                        /*
                        if (inside[i][j][k] && seed != GNP - 1 && seed != 0) {
                            // seed
                            if (seed % 2 == 0 ) {
                                rgba_px->r() = (unsigned char)(seed*seed) * 100.0f;
                                rgba_px->g() = (unsigned char)(seed*seed) * 100.0f;
                                rgba_px->b() = (unsigned char)(seed*seed) * 100.0f;
                                rgba_px->a() = (unsigned char)(200 * 1.0f);
                            }
                            else {
                                rgba_px->r() = (unsigned char)(seed) * 100.0f;
                                rgba_px->g() = (unsigned char)(seed) * 100.0f;
                                rgba_px->b() = (unsigned char)(seed) * 100.0f;
                                rgba_px->a() = (unsigned char)(200 * 1.0f);
                            }
                        }
                        else if (inside[i][j][k] && seed == GNP - 1) {
                            rgba_px->r() = (unsigned char)2900.0f;
                            rgba_px->g() = (unsigned char)2900.0f;
                            rgba_px->b() = (unsigned char)2900.0f;
                            rgba_px->a() = (unsigned char)(200 * 1.0f);
                        }

                        else {
                            rgba_px->r() = (unsigned char)0.0f;
                            rgba_px->g() = (unsigned char)0.0f;
                            rgba_px->b() = (unsigned char)0.0f;
                            rgba_px->a() = (unsigned char)(200 * 1.0f);
                        }
                        */
                    }
                }
            }

            Log::debug(Say("Visualizing after computation %n iter  \n") << ctr);
            for (int i = 0; i < N; i++) {
                Log::debug(Say("Seed %n Volume is %n  \n") << i << seedVolume[i]);
                seedVolume[i] = 0;
            }
            setupVolume();
            //setupScene();
        }


        if (key == vl::Key_R) {
            r0 = 2;
            //%Write a sphere with a radius 5 with different grain number (grain number written in circle)
            for (int iGN = 1; iGN <= GN - 1; iGN++) {
                x1 = ndx * DRND(1); y1 = ndy * DRND(1); z1 = ndz * DRND(1);
                while (inside[x1][y1][z1] != 1) { x1 = ndx * DRND(1); y1 = ndy * DRND(1); z1 = ndz * DRND(1); }
                for (i = -r0; i <= (ndxm + r0); i++) {
                    ii = i;
                    //if (i < 0) { ii = ndx + i; }  if (i > ndxm) { ii = i - ndx; }
                    if (i < 0) { ii = 0; } if (i > ndxm) { ii = ndxm; }
                    for (j = -r0; j <= (ndym + r0); j++) {
                        jj = j;
                        //if (j < 0) { jj = ndy + j; }  if (j > ndym) { jj = j - ndy; }
                        if (j < 0) { jj = 0; } if (j > ndym) { jj = ndym; }
                        for (k = -r0; k <= (ndzm + r0); k++) {
                            kk = k;
                            //if (k < 0) { kk = ndz + k; }  if (k > ndzm) { kk = k - ndz; }
                            if (k < 0) { kk = 0; } if (k > ndzm) { kk = ndzm; }
                            r = sqrt((double(i - x1)) * (double(i - x1))
                                + (double(j - y1)) * (double(j - y1))
                                + (double(k - z1)) * (double(k - z1)));
                            if (r <= r0 && inside[ii][jj][kk] == 1) {
                                if (!(i<0 || j<0 || k <0 || i >ndxm || j> ndym || k>ndzm)) {
                                    is[ii][jj][kk] = iGN;
                                }
                            }
                        }
                    }
                }
            }

            /*
            r0 = 5.0; x1 = ndx / 2; y1 = ndy / 2; z1 = ndz / 2;
            for (i = 0; i <= ndxm; i++) {
                for (j = 0; j <= ndym; j++) {
                    for (k = 0; k <= ndzm; k++) {
                        r = sqrt((double(i - x1)) * (double(i - x1))
                            + (double(j - y1)) * (double(j - y1))
                            + (double(k - z1)) * (double(k - z1)));
                        //if(r<=r0){ is[i][j][k]=1; }
                    }
                }
            }
            */

            for (i = 0; i <= ndxm; i++) {
                for (j = 0; j <= ndym; j++) {
                    for (k = 0; k <= ndzm; k++) {
                        if (inside[i][j][k] == 1) {

                            qh[1][i][j][k] = is[i][j][k];
                        }
                    }
                }
            }
        }
        if (key == vl::Key_Right)
            mode++;
        else
            if (key == vl::Key_Left)
                mode--;
        MODE = modes[vl::clamp(mode, 0, 5)];

        // up/down changes SAMPLE_STEP
        if (key == vl::Key_Up)
        {
            SAMPLE_STEP += 64; // more precision
        }
        else
            if (key == vl::Key_Down)
            {
                SAMPLE_STEP -= 64; // less precision
            }

        // L key toggles lights (useful only for isosurface)
        if (key == vl::Key_L)
        {
            if (!DYNAMIC_LIGHTS)
            {
                DYNAMIC_LIGHTS = true;
                COLORED_LIGHTS = false;
            }
            else
                if (DYNAMIC_LIGHTS && !COLORED_LIGHTS)
                {
                    DYNAMIC_LIGHTS = true;
                    COLORED_LIGHTS = true;
                }
                else
                {
                    DYNAMIC_LIGHTS = false;
                    COLORED_LIGHTS = false;
                }
        }

        setupScene();
    }

    virtual void destroyEvent() {
        BaseDemo::destroyEvent();
        // We need to release objects containing OpenGL objects before the OpenGL context is destroyed.
        mRaycastVolume = NULL;
        mVolumeAct = NULL;
        mGLSL = NULL;
        mLight0 = mLight1 = mLight2 = NULL;
        mValThreshold = NULL;
    }

private:
    int getColorR(int seed) {
        switch (seed) {
        case 1:
            return 255;

        case 2:
            return 255;

        case 3:
            return 0;

        case 4:
            return 170;

        case 5:
            return 255;

        case 6:
            return 191;

        case 7:
            return 0;

        case 8:
            return 255;

        case 9:
            return 255;

        case 10:
            return 106;

        case 11:
            return 0;

        case 12:
            return 237;

        case 13:
            return 185;

        case 14:
            return 231;

        case 15:
            return 220;

        case 16:
            return 185;

        case 17:
            return 143;

        case 18:
            return 35;

        case 19:
            return 143;

        case 20:
            return 107;

        case 21:
            return 79;

        case 22:
            return 10;

        case 23:
            return 115;

        case 24:
            return 204;

        case 25:
            return 128;

        case 26:
            return 128;

        case 27:
            return 0;

        case 28:
            return 128;

        case 29:
            return 0;

        case 30:
            return 0;

        case 31:
            return 139;

        case 32:
            return 244;

        case 33:
            return 188;


        default:
            return 255;
        }
    }
    int getColorG(int seed) {
        switch (seed) {
        case 1:
            return 0;

        case 2:
            return 255;

        case 3:
            return 234;

        case 4:
            return 0;

        case 5:
            return 127;

        case 6:
            return 255;

        case 7:
            return 149;

        case 8:
            return 0;

        case 9:
            return 212;

        case 10:
            return 255;

        case 11:
            return 64;

        case 12:
            return 185;

        case 13:
            return 215;

        case 14:
            return 233;

        case 15:
            return 185;

        case 16:
            return 237;

        case 17:
            return 35;

        case 18:
            return 98;

        case 19:
            return 106;

        case 20:
            return 35;

        case 21:
            return 143;

        case 22:
            return 10;

        case 23:
            return 115;

        case 24:
            return 204;

        case 25:
            return 0;

        case 26:
            return 128;

        case 27:
            return 128;

        case 28:
            return 0;

        case 29:
            return 128;

        case 30:
            return 0;

        case 31:
            return 69;

        case 32:
            return 164;

        case 33:
            return 143;


        default:
            return 255;
        }
    }
    int getColorB(int seed) {
        switch (seed) {
        case 1:
            return 0;

        case 2:
            return 0;

        case 3:
            return 255;

        case 4:
            return 255;

        case 5:
            return 0;

        case 6:
            return 0;

        case 7:
            return 255;

        case 8:
            return 170;

        case 9:
            return 0;

        case 10:
            return 0;

        case 11:
            return 255;

        case 12:
            return 185;

        case 13:
            return 237;

        case 14:
            return 185;

        case 15:
            return 237;

        case 16:
            return 224;

        case 17:
            return 35;

        case 18:
            return 143;

        case 19:
            return 35;

        case 20:
            return 143;

        case 21:
            return 35;

        case 22:
            return 0;

        case 23:
            return 115;

        case 24:
            return 204;

        case 25:
            return 0;

        case 26:
            return 0;

        case 27:
            return 0;

        case 28:
            return 128;

        case 29:
            return 128;

        case 30:
            return 128;

        case 31:
            return 19;

        case 32:
            return 96;

        case 33:
            return 143;

        default:
            return 255;
        }
    }
    ref<Transform> mVolumeTr;
    ref<Transform> mLight0Tr;
    ref<Transform> mLight1Tr;
    ref<Transform> mLight2Tr;
    ref<Uniform> mValThreshold;
    ref<Text> mValThresholdText;
    ref<Light> mLight0;
    ref<Light> mLight1;
    ref<Light> mLight2;
    ref<GLSLProgram> mGLSL;
    ref<Actor> mVolumeAct;
    ref<vl::RaycastVolume> mRaycastVolume;
    ref<Image> mVolumeImage;
};
// Have fun!

BaseDemo* Create_App_VolumeRaycast() { return new App_VolumeRaycast; }
