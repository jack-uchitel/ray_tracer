/* Jack Uchitel and Jade Kandel*/

/* On macOS, compile with...
    clang 610mainTexturing.c 000pixel.o -lglfw3 -framework OpenGL -framework Cocoa -framework IOKit
On Ubuntu, compile with...
    cc 600mainSpheres.c 000pixel.o -lglfw -lGL -lm -ldl
*/
#include <stdio.h>
#include <math.h>
#include <GLFW/glfw3.h>
#include "000pixel.h"

#include "610vector.c"
#include "040texture.c"
#include "140matrix.c"
#include "600isometry.c"
#include "600camera.c"
#include "620ray.c"

#define SCREENWIDTH 512
#define SCREENHEIGHT 512
#define TEXR 0
#define TEXG 1
#define TEXB 2



/*** Lights, camera ***/

camCamera camera;
double cameraTarget[3] = {0.0, 0.0, 0.0};
double cameraRho = 10.0, cameraPhi = M_PI / 3.0, cameraTheta = -M_PI / 3.0;

void initializeLightsCamera(void) {
    camSetProjectionType(&camera, camPERSPECTIVE);
    camSetFrustum(&camera, M_PI / 6.0, cameraRho, 10.0, SCREENWIDTH,
        SCREENHEIGHT);
    camLookAt(&camera, cameraTarget, cameraRho, cameraPhi, cameraTheta);
}

void destroyLightsCamera(void) {
    return;
}



/*** Scene ***/

isoIsometry isomA, isomB;
double radiusA = 1.0, radiusB = 1.5;
double colorA[3] = {1.0, 0.0, 1.0}, colorB[3] = {1.0, 1.0, 0.0};
texTexture tex;


void initializeScene(void) {
    double center[3] = {0.0, 0.0, 0.0}, axis[3] = {1.0, 0.0, 0.0}, r[3][3];
    mat33AngleAxisRotation(0.0, axis, r);
    isoSetTranslation(&isomA, center);
    isoSetRotation(&isomA, r);
    vec3Set(1.0, 0.0, 1.0, center);
    isoSetTranslation(&isomB, center);
    isoSetRotation(&isomB, r);
}

void destroyScene(void) {
    return;
}



/*** Rendering ***/

/* Given the isometry and radius of a spherical body. Given a ray query. Returns information about how the ray first intersects the sphere. */
rayResponse intersectSphere(
        const isoIsometry *iso, double radius, rayQuery query) {
    rayResponse result;
    double d[3] = {query.d[0],query.d[1] ,query.d[2]};
    double e[3] = {query.e[0],query.e[1] ,query.e[2]};
    double center[3] = {iso->translation[0],iso->translation[1],iso->translation[2] };
    double twoTimesD[3];
    double eMinusCenter[3];
    vecScale(3,2,d,twoTimesD);
    vecSubtract(3,e,center,eMinusCenter);
    double dInA;
    dInA = vecLength(3,d);
    //printf("dInA: %f\n", dInA);
    double a = pow(dInA,2) ;
    double b = vecDot(3,twoTimesD,eMinusCenter);
    double magEMinusCenter;
    magEMinusCenter= vecLength(3,eMinusCenter);
    double rSquared = pow(radius,2);
    double magEMinusCenterSquared = pow(magEMinusCenter,2);
    double c = magEMinusCenterSquared - rSquared;

   /* printf("a: %f ", a);
    printf("b: %f ", b);
    printf("c: %f\n", c);*/
    if ((pow(b,2) - (4*a*c))<=0){
        result.intersected = rayNONE;
    } else {
        double t1, t2;
        t2 = (-b + sqrt(pow(b,2)-4*a*c))/(2*a);
        t1 = (-b - sqrt(pow(b,2)-4*a*c))/(2*a);
        /*printf("a: %f ", a);
        printf("t1: %f ", t1);
        printf("t2: %f\n", t2);*/
        /*if ( query.tStart< t1 && t1 < query.tEnd ){
            if(query.tStart< t2 && t2 < query.tEnd  ){
                if (t1<t2){
                    query.tEnd= t1;
                    result.t=t1;
                    result.intersected = rayEXIT;

                } else {
                    query.tEnd = t2;
                    result.t=t2;
                    result.intersected = rayENTER;

                }
            }
            query.tEnd = t1;
            result.t=t1;
            result.intersected = rayEXIT;
        }
        
        if (query.tStart< t2 && t2 < query.tEnd){
            query.tEnd = t2;
            result.t =t2;
            result.intersected = rayENTER;
        } else {
            result.intersected = rayNONE;

        }*/
        if ( query.tStart< t1 && t1 < query.tEnd ){
            result.t=t1;
            result.intersected = rayENTER;
        } else if (query.tStart< t2 && t2 < query.tEnd) {
            result.t =t2;
            result.intersected = rayEXIT;
        } else {
            result.intersected = rayNONE;

        }
    }
    
    return result;
}

/* Given a sphere described by an isometry, radius, color, and texture. Given
the query and response that prompted this request for color. Outputs the RGB
color. */
void colorSphere(
        const isoIsometry *iso, double radius, const double color[3],
        const texTexture *tex, rayQuery query, rayResponse resp,
                 double rgb[3]) {
    double localCoords[3];
    double tD[3];
    vecScale(3, resp.t, query.d, tD);
    double worldCoords[3];
    vecAdd(3, query.e, tD, worldCoords);

    isoUntransformPoint(&(camera.isometry), worldCoords, localCoords);

    double rho;
    double phi;
    double theta;
    vec3Rectangular(localCoords, &rho, &phi, &theta);
    double s = phi/M_PI;
    double t = theta/(2*M_PI);

    double texColor[3];
    texSample(tex, s, t, texColor);
    rgb[0] = texColor[0] * color[0];
    rgb[1] = texColor[1] * color[1];
    rgb[2] = texColor[2] * color[2];
    
    
}

void render(void) {
    /* Get camera world position e and transformation from screen to world. */
    rayQuery query;
    double trans4mation[4][4];
    camWorldFromScreenHomogeneous(&camera, SCREENWIDTH, SCREENHEIGHT, trans4mation);
    vec3Set(camera.isometry.translation[0],
            camera.isometry.translation[1],
            camera.isometry.translation[2],
            query.e);
    /*printf("query.e: %f,%f,%f ",query.e[0], query.e[1], query.e[2]);
    printf("translation: %f,%f,%f ",camera.isometry.translation[0],
           camera.isometry.translation[1],
           camera.isometry.translation[2]);
*/

    /* Each screen point is arbitrarily chosen on the near plane. */
    double screen[4];
    screen[2] = 0.0;
    screen[3] = 1.0;
    for (int i = 0; i < SCREENWIDTH; i += 1) {
        screen[0] = i;
        for (int j = 0; j < SCREENHEIGHT; j += 1) {
            screen[1] = j;
            /* Compute the direction d from the camera to the pixel. */
            double pixelWorld[4];
            mat441Multiply(trans4mation, screen, pixelWorld);
            double worldCoord[3] = {pixelWorld[0]/pixelWorld[3],pixelWorld[1]/pixelWorld[3],pixelWorld[2]/pixelWorld[3]};
            vecSubtract(3,worldCoord, query.e ,query.d);
           // printf("worldCoord: %f, %f, %f, %f \n",pixelWorld[0],pixelWorld[1],pixelWorld[2],pixelWorld[3] );
            //!!Student code goes here.
            /* Prepare to loop over all bodies. */
            query.tStart = rayEPSILON;
            query.tEnd = rayINFINITY;
            /* Test the first sphere. */
            //call the funciton on isometry A get rayResponseA
            rayResponse responseA = intersectSphere(&isomA, radiusA, query);
            if (responseA.intersected != rayNONE){
                query.tEnd = responseA.t;
            }

            //change query tEnd here
            //!!Student code goes here.
            /* Test the second sphere. */
            //call the funciton on isometry B get rayResponseB
            rayResponse responseB = intersectSphere(&isomB, radiusB, query);
            //!!Student code goes here.
            /* Choose the winner. */
            

            double rgb[3];
        
            //!!Student code goes here.
            //need to check if rayNone. Then can compare t values.
            if (responseA.intersected == rayNONE){
                if (responseB.intersected == rayNONE) {
                    rgb[0] = 0;
                    rgb[1] = 0;
                    rgb[2] = 0;
                } else {
                    colorSphere(&(camera.isometry), radiusB, colorB,
                            &tex, query, responseB,
                                rgb);

                }
            } else {
                if (responseB.intersected == rayNONE) {
                    colorSphere(
                            &(camera.isometry), radiusA, colorA,
                            &tex, query, responseA,
                                rgb);

                }
                else {
                    if (responseA.t < responseB.t) {
                        colorSphere(
                                &(camera.isometry), radiusA, colorA,
                                &tex, query, responseA,
                                    rgb);
                    }
                    else {
                        colorSphere(
                                &(camera.isometry), radiusB, colorB,
                                &tex, query, responseB,
                                    rgb);
                    }
                }
            }
            
            //printf("rgb: %f,%f,%f", rgb[0],rgb[1],rgb[2] );
            pixSetRGB(i, j, rgb[0], rgb[1], rgb[2]);
        }
    }
}



/*** User interface ***/

void handleKey(
        int key, int shiftIsDown, int controlIsDown, int altOptionIsDown,
        int superCommandIsDown) {
    texSetFiltering(&tex,0);
    if (key == GLFW_KEY_W)
        cameraPhi -= 0.1;
    else if (key == GLFW_KEY_A)
        cameraTheta -= 0.1;
    else if (key == GLFW_KEY_S)
        cameraPhi += 0.1;
    else if (key == GLFW_KEY_D)
        cameraTheta += 0.1;
    else if (key == GLFW_KEY_E)
        cameraRho *= 0.9;
    else if (key == GLFW_KEY_C)
        cameraRho *= 1.1;
    camSetFrustum(&camera, M_PI / 6.0, cameraRho, 10.0, SCREENWIDTH,
        SCREENHEIGHT);
    camLookAt(&camera, cameraTarget, cameraRho, cameraPhi, cameraTheta);
}

void handleTimeStep(double oldTime, double newTime) {
    if (floor(newTime) - floor(oldTime) >= 1.0)
        printf("handleTimeStep: %f frames/s\n", 1.0 / (newTime - oldTime));
    render();
}

int main(void) {
    if (pixInitialize(SCREENWIDTH, SCREENHEIGHT, "600mainSpheres") != 0)
        return 1;
    else {
        texInitializeFile(&tex, "foliage.jpg");
        initializeLightsCamera();
        initializeScene();
        pixSetKeyDownHandler(handleKey);
        pixSetKeyRepeatHandler(handleKey);
        pixSetTimeStepHandler(handleTimeStep);
        pixRun();
        destroyScene();
        destroyLightsCamera();
        texDestroy(&tex);
        return 0;
    }
}

