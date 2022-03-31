/* Jack Uchitel and Jade Kandel*/

/* On macOS, compile with...
 clang 630mainShadows.c 000pixel.o -lglfw3 -framework OpenGL -framework Cocoa -framework IOKit
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
void debug(int integer, double real, const char *message) {
    fprintf(stderr, "%d %f %s", integer, real, message);
    fflush(stderr);
}

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
/* Intersects the ray with every object in the scene. If no object is hit in
[tStart, tEnd], then the response signals rayNONE. If any object is hit, then
the response is the response from the nearest object, and index is the index of
that object in the scene's list of bodies (0, 1, 2, ...). */
rayResponse intersectScene(rayQuery query, int *index){
    //calling for sphereA
    //debug(0,0,"beginning of intersectScene \n");
    rayResponse result;
    rayResponse responseA = intersectSphere(&isomA, radiusA, query);
    if (responseA.intersected != rayNONE){
        query.tEnd = responseA.t;
    }
    rayResponse responseB = intersectSphere(&isomB, radiusB, query);
    
    if (responseA.intersected == rayNONE){
        if (responseB.intersected == rayNONE) {
            result.intersected = rayNONE;
            
        } else {
            result = responseB;
            *index = 1;

        }
    } else {
        if (responseB.intersected == rayNONE) {
            result = responseA;
            *index = 0;

        }
        else {
            if (responseA.t < responseB.t) {
                result = responseA;
                *index = 0;
            }
            else {
                result = responseB;
                *index = 1;

            }
        }
    }
    //debug(0,0,"end of intersectScene \n");

    return result;
}

/* Given a sphere described by an isometry, radius, color, and texture. Given
the query and response that prompted this request for color. Outputs the RGB
color. */
void colorSphere(
        const isoIsometry *iso, double radius, const double color[3],
        const texTexture *tex, rayQuery query, rayResponse resp,
                 double rgb[3]) {
    //texSampling and creating rgb
    
    double newCoords[3];
    double localCoords[3];
    double tD[3];
    vecScale(3, resp.t, query.d, tD);
    //double worldCoords[3];
    vecAdd(3, query.e, tD, newCoords);



    
    isoUntransformPoint(&(camera.isometry), newCoords, localCoords);

    double rho;
    double phi;
    double theta;
    vec3Rectangular(localCoords, &rho, &phi, &theta);
    double s = phi/M_PI;
    double t = theta/(2*M_PI);
    
    double cDiff[3];
    
    double texColor[3];
    texSample(tex, s, t, texColor);
    cDiff[0] = texColor[0] * color[0];
    cDiff[1] = texColor[1] * color[1];
    cDiff[2] = texColor[2] * color[2];
    
    //dNormal
    double dNormal[3];
    double xMinusCenter[3];
    double center[3] = {iso->translation[0],iso->translation[1],iso->translation[2]};
    vecSubtract(3,newCoords,center, xMinusCenter);
    vecUnit(3,xMinusCenter, dNormal);
    
    //dCam
    double dCamUnscaled[3];
    double dCam[3];
    vecUnit(3,query.d,dCamUnscaled);
    vecScale(3,-1,dCamUnscaled,dCam);
    
    //dRefl
    double dRefl[3];
    rayReflect(3,dCam,dNormal,dRefl);
    
    
    double dLight[3] = {0,0,1};
    double white[3] = {1,1,1};
    double shininess = 5;
    
    rayQuery shadowQuery;
    //double shadowVector[3];
    //vecSubtract(3,newCoords,dLight,shadowVector);
    
    vec3Set(newCoords[0],
            newCoords[1],
            newCoords[2],
            shadowQuery.e);
    
    vec3Set(dLight[0],
            dLight[1],
            dLight[2],
            shadowQuery.d);
    
    shadowQuery.tStart = rayEPSILON;
    shadowQuery.tEnd = rayINFINITY;
    int index;
    
    rayResponse result = intersectScene(shadowQuery, &index);
    if (result.intersected == rayNONE) {
        rayDiffuseAndSpecular(
                              dNormal, dLight, dRefl,
                              cDiff, white, shininess,
                              white, rgb);
    } else {
        rgb[0] = 0;
        rgb[1] = 0;
        rgb[2] = 0;
        
    }
    
    
   
    //printf("rgb: %f,%f,%f", rgb[0],rgb[1],rgb[2]);
    
    rgb[0] = rgb[0] +(0.1 * cDiff[0]);
    rgb[1] = rgb[1] +(0.1 * cDiff[1]);
    rgb[2] = rgb[2] +(0.1 * cDiff[2]);
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
    int index;
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
            /*rayResponse responseA = intersectSphere(&isomA, radiusA, query);
            if (responseA.intersected != rayNONE){
                query.tEnd = responseA.t;
            }*/

            //change query tEnd here
            //!!Student code goes here.
            /* Test the second sphere. */
            //call the funciton on isometry B get rayResponseB
            /*rayResponse responseB = intersectSphere(&isomB, radiusB, query);*/
            //!!Student code goes here.
            /* Choose the winner. */
            

            double rgb[3];
            
            rayResponse result = intersectScene(query, &index);
            
            if (result.intersected == rayNONE) {
                rgb[0] = 0;
                rgb[1] = 0;
                rgb[2] = 0;
            }
            else if (index == 0) {
                colorSphere(
                            &isomA, radiusA, colorA,
                        &tex, query, result,
                            rgb);
            }
            else if (index == 1) {
                colorSphere(&isomB, radiusB, colorB,
                        &tex, query, result,
                            rgb);
            }
        
            //!!Student code goes here.
            //need to check if rayNone. Then can compare t values.
            /*if (responseA.intersected == rayNONE){
                if (responseB.intersected == rayNONE) {
                    rgb[0] = 0;
                    rgb[1] = 0;
                    rgb[2] = 0;
                } else {
                    colorSphere(&isomB, radiusB, colorB,
                            &tex, query, responseB,
                                rgb);

                }
            } else {
                if (responseB.intersected == rayNONE) {
                    colorSphere(
                                &isomA, radiusA, colorA,
                            &tex, query, responseA,
                                rgb);

                }
                else {
                    if (responseA.t < responseB.t) {
                        colorSphere(
                                     &isomA, radiusA, colorA,
                                &tex, query, responseA,
                                    rgb);
                    }
                    else {
                        colorSphere(
                                    &isomB, radiusB, colorB,
                                &tex, query, responseB,
                                    rgb);
                    }
                }
            }*/
            
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

