/* Jack Uchitel and Jade Kandel*/

/* On macOS, compile with...
 clang 670mainAbstracted.c 000pixel.o -lglfw3 -framework OpenGL -framework Cocoa -framework IOKit
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
#include "660ray.c"

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
    double center[3] = {-1, -0.5, -1}, axis[3] = {1.0, 0.0, 0.0}, r[3][3];
    mat33AngleAxisRotation(0.0, axis, r);
    isoSetTranslation(&isomA, center);
    isoSetRotation(&isomA, r);
    vec3Set(1.2, 0.0, 1.2, center);
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
    double a = pow(dInA,2) ;
    double b = vecDot(3,twoTimesD,eMinusCenter);
    double magEMinusCenter;
    magEMinusCenter= vecLength(3,eMinusCenter);
    double rSquared = pow(radius,2);
    double magEMinusCenterSquared = pow(magEMinusCenter,2);
    double c = magEMinusCenterSquared - rSquared;
    
    
   
    
    if ((pow(b,2) - (4*a*c))<=0){
        result.intersected = rayNONE;
    } else {
        double t1, t2;
        t2 = (-b + sqrt(pow(b,2)-4*a*c))/(2*a);
        t1 = (-b - sqrt(pow(b,2)-4*a*c))/(2*a);
        
        if ( query.tStart< t1 && t1 < query.tEnd ){
            
            result.t=t1;
            result.intersected = rayENTER;
            
            //get worldCoords
            double newCoords[3];
            double localCoords[3];
            double tD[3];
            vecScale(3, t1, query.d, tD);
            vecAdd(3, query.e, tD, newCoords);
            
            ////dNormal
            double dNormal[3];
            double xMinusCenter[3];
            vecSubtract(3,newCoords,center, xMinusCenter);
            vecUnit(3,xMinusCenter, dNormal);
            
            //get tex Coords
            isoUntransformPoint(&(camera.isometry), newCoords, localCoords);

            double rho;
            double phi;
            double theta;
            vec3Rectangular(localCoords, &rho, &phi, &theta);

            
            
            vec3Set(dNormal[0], dNormal[1], dNormal[2], result.normal);
            result.texCoords[0] = phi/M_PI;
            result.texCoords[1] = theta/(2*M_PI);
            
        } else if (query.tStart< t2 && t2 < query.tEnd) {
            result.t =t2;
            result.intersected = rayEXIT;
            
            //get worldCoords
            double newCoords[3];
            double localCoords[3];
            double tD[3];
            vecScale(3, t2, query.d, tD);
            vecAdd(3, query.e, tD, newCoords);
            
            ////dNormal
            double dNormal[3];
            double xMinusCenter[3];
            vecSubtract(3,newCoords,center, xMinusCenter);
            vecUnit(3,xMinusCenter, dNormal);
            
            //get tex Coords
            isoUntransformPoint(&(camera.isometry), newCoords, localCoords);

            double rho;
            double phi;
            double theta;
            vec3Rectangular(localCoords, &rho, &phi, &theta);
            double s = phi/M_PI;
            double t = theta/(2*M_PI);
            
            
            vec3Set(dNormal[0], dNormal[1], dNormal[2], result.normal);
            result.texCoords[0] = phi/M_PI;
            result.texCoords[1] = theta/(2*M_PI);
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

void colorScene(rayQuery query, double rgb[3]);


/* Given a sphere described by an isometry, radius, color, and texture. Given
the query and response that prompted this request for color. Outputs the RGB
color. */

void rayColor(
              rayQuery query, rayResponse resp, rayMaterial mat, double rgb[3]) {
    
    double newCoords[3];
    double localCoords[3];
    double tD[3];
    vecScale(3, resp.t, query.d, tD);
    //double worldCoords[3];
    vecAdd(3, query.e, tD, newCoords);
    
    //add ambient
    
   
    rgb[0] = (0.2 * mat.cDiff[0]);
    rgb[1] = (0.2 * mat.cDiff[1]);
    rgb[2] = (0.2 * mat.cDiff[2]);
    
    //dCam
    double dCamUnscaled[3];
    double dCam[3];
    vecUnit(3,query.d,dCamUnscaled);
    vecScale(3,-1,dCamUnscaled,dCam);
    
    //dRefl
    double dRefl[3];
    rayReflect(3,dCam,resp.normal,dRefl);
 
    
    if (mat.hasDiffAndSpec != 0) {
        double diffAndSpec[3];
        double dLight[3] = {0,0,1};
        
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
                                  resp.normal, dLight, dRefl,
                                  mat.cDiff, mat.cSpec, mat.shininess,
                                  mat.cSpec, diffAndSpec);
        rgb[0] = rgb[0] + diffAndSpec[0];
        rgb[1] = rgb[1] + diffAndSpec[1];
        rgb[2] = rgb[2] + diffAndSpec[2];
        }
    }
        
    if (mat.hasMirror != 0) {
        double mirror[3];
        rayQuery mirrorQuery;
        vec3Set(newCoords[0], newCoords[1], newCoords[2], mirrorQuery.e);
        vec3Set(dRefl[0], dRefl[1], dRefl[2], mirrorQuery.d);
        mirrorQuery.tStart = rayEPSILON;
        mirrorQuery.tEnd = rayINFINITY;
        colorScene(mirrorQuery, mirror);
        rgb[0] = rgb[0]+ mirror[0]* mat.cMirror[0];
        rgb[1] = rgb[1] + mirror[1]* mat.cMirror[1];
        rgb[2] = rgb[2] + mirror[2]* mat.cMirror[2];
    }
}

rayMaterial getMaterialA(
        const isoIsometry *iso, double radius, const double color[3],
        const texTexture *tex, rayQuery query, rayResponse resp,
                 double rgb[3]) {
    
    //diffuse and specular
    rayMaterial mat;
    
    
    mat.hasDiffAndSpec = 1;
    
    double cDiff[3];
    double texColor[3];
    texSample(tex, resp.texCoords[0], resp.texCoords[1], texColor);
    cDiff[0] = texColor[0] * color[0];
    cDiff[1] = texColor[1] * color[1];
    cDiff[2] = texColor[2] * color[2];
    
    vec3Set(cDiff[0],cDiff[1],cDiff[2], mat.cDiff);
    vec3Set(1,1,1, mat.cSpec);
    mat.shininess = 5;
    
    //mirror
    mat.hasMirror = 0;

    return mat;
    //rayColor(query,resp,mat,rgb);
}

rayMaterial getMaterialMirror(
        const isoIsometry *iso, double radius, const double color[3],
        const texTexture *tex, rayQuery query, rayResponse resp,
                 double rgb[3]) {
    
    //diffuse and specular
    //printf("%f,%f,%f", rgb[0],rgb[1],rgb[2]);
    rayMaterial mat;
    
    double cDiff[3];
    double texColor[3];
    texSample(tex, resp.texCoords[0], resp.texCoords[1], texColor);
    cDiff[0] = texColor[0] * color[0];
    cDiff[1] = texColor[1] * color[1];
    cDiff[2] = texColor[2] * color[2];
    
    vec3Set(cDiff[0],cDiff[1],cDiff[2], mat.cDiff);
    
    mat.hasDiffAndSpec = 0;

    mat.hasMirror = 1;
    vec3Set(1,1,1, mat.cMirror);
    //mirror
    
    return mat;
    //rayColor(query,resp,mat,rgb);
    //printf("%f,%f,%f", rgb[0],rgb[1],rgb[2]);

}



/*void getMaterialMirror(
        const isoIsometry *iso, double radius, const double color[3],
        const texTexture *tex, rayQuery query, rayResponse resp,
                 double rgb[3]) {


}
 */
/* Given a query, the resulting response, and the resulting material, computes
the fragment's final RGB color. As of the 660 version, this function is allowed
to access the scene (bodies and lights) through global variables. */

/* takes in material and slaps on ambient lighting. then, if it also has diffuse ans specular, add that on, then, if it implements mirroring, add that on. */


/* Given a ray query, tests the entire scene, setting rgb to be the color of
whatever object that query hits (or the background). */
void colorScene(rayQuery query, double rgb[3]){
    int index;
    rayResponse result = intersectScene(query, &index);
    
    if (result.intersected == rayNONE) {
        rgb[0] = 0;
        rgb[1] = 0;
        rgb[2] = 0;
    }
    else if (index == 0) {
        rayMaterial mat = getMaterialA(
                    &isomA, radiusA, colorA,
                &tex, query, result,
                    rgb);
        rayColor(query,result,mat,rgb);
        
    }
    else if (index == 1) {
        rayMaterial mat = getMaterialMirror(&isomB, radiusB, colorB,
                &tex, query, result,
                    rgb);
        rayColor(query,result,mat,rgb);
    }
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
           
            //!!Student code goes here.
            /* Prepare to loop over all bodies. */
            query.tStart = rayEPSILON;
            query.tEnd = rayINFINITY;

            double rgb[3] = {0,0,0};
            colorScene(query,rgb);
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

