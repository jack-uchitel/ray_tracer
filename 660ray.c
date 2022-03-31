


/* In principle, these constants might have to be tuned to your scene. In
practice, probably you can just leave them as they are. */
#define rayEPSILON 0.00000001
#define rayINFINITY 100000000

/* These constants help make rayResponse less mysterious. */
#define rayENTER -1
#define rayEXIT 1
#define rayNONE 0

typedef struct rayQuery rayQuery;
struct rayQuery {
    /* The ray is parametrized as x(t) = e + t d for t in [tStart, tEnd]. The
    direction d must be non-zero but not necessarily unit. */
    double e[3], d[3], tStart, tEnd;
};

typedef struct rayResponse rayResponse;
struct rayResponse {
    /* Records the nature of the first intersection of the ray with the body,
    considering only times t in [tStart, tEnd]: rayNONE for no intersection (or
    just tangency), rayENTER for entering the body, or rayEXIT for exiting the
    body. */
    int intersected;
    /* If the intersection code is non-zero, then contains the first time of
    intersection. */
    double t;
    
    /* If the intersection code is not rayNONE, then contains the unit
    outward-pointing normal vector at the intersection point, in world
    coordinates. */
    double normal[3];
    /* If the intersection code is not rayNONE, then contains the texture
    coordinates at the intersection point. */
    double texCoords[2];
};

typedef struct rayMaterial rayMaterial;
struct rayMaterial {
    /* In every case, the ray tracer should use cDiff to add some ambient light
    to the fragment color. If hasDiffAndSpec != 0, then the ray tracer should
    also calculate the rest of the Phong lighting model and add the result into
    the fragment color. If hasDiffAndSpec == 0, then cSpec is effectively black,
    but cDiff still matters for the ambient lighting. */
    int hasDiffAndSpec;
    double cDiff[3], cSpec[3], shininess;
    /* If hasMirror != 0, then the ray tracer should send out a mirror ray, get
    its color, modulate that color by cMirror, and add the result into the
    fragment color. If hasMirror == 0, then cMirror is effectively black. */
    int hasMirror;
    double cMirror[3];
};


/* Calculates the sum of the diffuse and specular parts of the Phong reflection
model, storing that sum in the rgb argument. Directions are unit vectors. They
can be either local or global, as long as they're consistent. dRefl is the
reflected camera direction. */
void rayDiffuseAndSpecular(
        const double dNormal[3], const double dLight[3], const double dRefl[3],
        const double cDiff[3], const double cSpec[3], double shininess,
                           const double cLight[3], double rgb[3]){
    //calculate iDiff
    double iDiff;
    //printf("dNormal: %f %f %f \n", dNormal[0], dNormal[1], dNormal[2]);
    //printf("dLight: %f %f %f \n", dLight[0], dLight[1], dLight[2]);
    iDiff= vecDot(3,dNormal, dLight);
    //printf("iDiff = %f\n", iDiff);
    if(iDiff<0){
        iDiff=0;
    }
    
    //calculate iSpec
    double iSpec;
    if (iDiff==0){
        iSpec=0;
    }else{
        iSpec = vecDot(3,dRefl,dLight);
        /*printf("iSpec: %f\n", iSpec);
        printf("dLight: %f, %f,%f\n", dLight[0],dLight[1],dLight[2]);
        printf("dRefl: %f, %f,%f\n", dRefl[0],dRefl[1],dRefl[2]);*/
        if (iSpec<0){
            iSpec=0;
        }
    }
    
    //diffuse and specular
    double diffuse[3] = {iDiff*cDiff[0]*cLight[0],
                         iDiff*cDiff[1]*cLight[1],
                         iDiff*cDiff[2]*cLight[2]};
    //printf("diffuse: %f,%f,%f \n", diffuse[0],diffuse[1],diffuse[2]);

    //specular = pow(iSpec,shininess) * cSpec * cLight
    double specular[3] = {pow(iSpec,shininess) * cSpec[0] * cLight[0],
                          pow(iSpec,shininess) * cSpec[1] * cLight[1],
                          pow(iSpec,shininess) * cSpec[2] * cLight[2]};
    //printf("specular: %f,%f,%f \n", specular[0],specular[1],specular[2]);

    
    rgb[0] = diffuse[0]+specular[0];
    
    rgb[1] = diffuse[1]+specular[1];
    rgb[2] = diffuse[2]+specular[2];
    
}

/* Given a (not necessarily unit) vector v and a unit vector n, reflects v
across n, storing the result in refl. The output can safely alias the input. */
void rayReflect(int dim, const double v[], const double n[], double refl[]){
    //2*dot(v,n)*n-v

    double scalar = 2*vecDot(dim, v, n);
    double bigboi[dim];
    vecScale(dim,scalar,n, bigboi);
    vecSubtract(dim,bigboi,v,refl);
    
}

/*
void getMaterialA(
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

    
    rayColor(query,resp,mat,rgb);
}

void getMaterialMirror(
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
    
    rayColor(query,resp,mat,rgb);
    //printf("%f,%f,%f", rgb[0],rgb[1],rgb[2]);

}
*/
