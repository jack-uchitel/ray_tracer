/*** General ***/

/* Feel free to read from this struct's members, but don't write to them, except through the accessors below. */
typedef struct bodyBody bodyBody;
struct bodyBody {
    rayResponse (*getIntersection)(const bodyBody *body, rayQuery query);
    rayMaterial (*getMaterial)(
        const bodyBody *body, rayQuery query, rayResponse response);
    isoIsometry isometry;
    int auxNum, texNum;
    double *auxiliaries;
    const texTexture **textures;
    const void *data;
};

/* Each auxiliary is one double. Returns 0 on success and non-zero on error.
Upon success, don't forgot to call bodyDestroy when you are done with the body.
*/
int bodyInitialize(
        bodyBody *body, int auxNum, int texNum,
        rayResponse (*getIntersection)(const bodyBody *body, rayQuery query),
        rayMaterial (*getMaterial)(
            const bodyBody *body, rayQuery query, rayResponse response)) {
    double rotation[3][3] = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
    double translation[3] = {0.0, 0.0, 0.0};
    isoSetRotation(&(body->isometry), rotation);
    isoSetTranslation(&(body->isometry), translation);
    body->auxiliaries = (double *)malloc(auxNum * sizeof(double) +
            texNum * sizeof(texTexture *));
    if (body->auxiliaries == NULL)
        return 1;
    body->textures = (const texTexture **)&(body->auxiliaries[auxNum]);
    body->auxNum = auxNum;
    body->texNum = texNum;
    body->getIntersection = getIntersection;
    body->getMaterial = getMaterial;
    return 0;
}

void bodyDestroy(bodyBody *body) {
    if (body->auxiliaries != NULL) {
        free(body->auxiliaries);
        body->auxiliaries = NULL;
    }
}

/* Sets some of the body's auxiliaries. Namely, starting at
body->auxiliaries[index], copies dim numbers from a into the auxiliaries. */
void bodySetAuxiliaries(bodyBody *body, int index, int dim, const double *a) {
    if (0 <= index && index <= body->auxNum - dim)
        vecCopy(dim, a, &(body->auxiliaries[index]));
    else
        fprintf(stderr, "error: bodySetAuxiliaries: bounds exceeded.\n");
}

/* Sets one of the body's textures. */
void bodySetTexture(bodyBody *body, int index, const texTexture *tex) {
    if (0 <= index && index < body->texNum)
        body->textures[index] = tex;
    else
        fprintf(stderr, "error: bodySetTexture: bounds exceeded.\n");
}



/*** Cylinders ***/

/* Performs the intersection calculations for a cylinder centered along the
z-axis in local coordinates. Assumes that body->auxiliaries[0] is the radius. */
rayResponse bodyIntersectCylinder(const bodyBody *body, rayQuery query) {
    rayResponse resp;
    resp.intersected = rayNONE;
    /* Transform to local coordinates. */
    double eLocal[3], dLocal[3];
    isoUntransformPoint(&(body->isometry), query.e, eLocal);
    isoUnrotateVector(&(body->isometry), query.d, dLocal);
    /* Ignore the third dimension. */
    double eE, dE, dD, rSq, disc, t;
    eE = vecDot(2, eLocal, eLocal);
    dE = vecDot(2, dLocal, eLocal);
    dD = vecDot(2, dLocal, dLocal);
    disc = dE * dE - dD * (eE - body->auxiliaries[0] * body->auxiliaries[0]);
    if (disc <= 0)
        resp.intersected = rayNONE;
    double sqrtDisc = sqrt(disc);
    t = (-dE - sqrtDisc) / dD;
    if (query.tStart <= t && t <= query.tEnd) {
        resp.intersected = rayENTER;
        resp.t = t;
    } else {
        t = (-dE + sqrtDisc) / dD;
        if (query.tStart <= t && t <= query.tEnd) {
            resp.intersected = rayEXIT;
            resp.t = t;
        }
    }
    if (resp.intersected != rayNONE) {
        /* In local coordinates, x = e + t d. */
        double xLocal[3];
        vecScale(3, resp.t, dLocal, xLocal);
        vecAdd(3, eLocal, xLocal, xLocal);
        /* Texture coordinates are essentially cylindrical coordinates. */
        resp.texCoords[0] = atan2(xLocal[1], xLocal[0]);
        if (resp.texCoords[0] < 0.0)
            resp.texCoords[0] += 2.0 * M_PI;
        resp.texCoords[0] = resp.texCoords[0] / (2.0 * M_PI);
        resp.texCoords[1] = xLocal[2];
        /* Get the unit outward-pointing normal in world coordinates. */
        double nLocal[3] = {xLocal[0], xLocal[1], 0.0};
        vecUnit(3, nLocal, nLocal);
        isoRotateVector(&(body->isometry), nLocal, resp.normal);
    }
    return resp;
}

/* auxNum must be at least 1. auxiliaries[0] is set to the radius. Returns 0 upon success and non-zero when an error occurs. Upon success, don't forget to call bodyDestroy when you are done with the body. */
int bodyInitializeCylinder(
        bodyBody *body, int auxNum, int texNum, double radius,
        rayMaterial (*getMaterial)(
            const bodyBody *body, rayQuery query, rayResponse response)) {
    int error = bodyInitialize(body, auxNum, texNum, bodyIntersectCylinder,
        getMaterial);
    if (error == 0)
        bodySetAuxiliaries(body, 0, 1, &radius);
    return error;
}



/*** Spheres ***/

/* Performs the intersection calculations for a sphere centered at the origin in
local coordinates. Assumes that body->auxiliaries[0] is the radius. */
rayResponse bodyIntersectSphere(const bodyBody *body, rayQuery query) {
    rayResponse resp;
    rayResponse result;
    double radius = body->auxiliaries[0];
    
    double d[3] = {query.d[0],query.d[1] ,query.d[2]};
    double e[3] = {query.e[0],query.e[1] ,query.e[2]};
    double center[3] = {body->isometry.translation[0],body->isometry.translation[1],body->isometry.translation[2] };
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

/* auxNum must be at least 1. auxiliaries[0] is set to the radius. Returns 0 upon success and non-zero when an error occurs. Upon success, don't forget to call bodyDestroy when you are done with the body. */
int bodyInitializeSphere(
        bodyBody *body, int auxNum, int texNum, double radius,
        rayMaterial (*getMaterial)(
            const bodyBody *body, rayQuery query, rayResponse response)) {
    int error = bodyInitialize(body, auxNum, texNum, bodyIntersectSphere,
        getMaterial);
    if (error == 0)
        bodySetAuxiliaries(body, 0, 1, &radius);
    return error;
}
