/* Jack Uchitel and Jade Kandel*/



/* Describes an isometry as a rotation followed by a translation. Can be used
to describe the position and orientation of a rigid body. If the position is
the translation, and the columns of the rotation are the local coordinate axes
in global coordinates, then the isometry takes local coordinates to global. */

/* Feel free to read from, but not write to, this struct's members. */
typedef struct isoIsometry isoIsometry;
struct isoIsometry {
    double translation[3];
    double rotation[3][3];
};

/* Sets the rotation. */
void isoSetRotation(isoIsometry *iso, const double rot[3][3]) {
    vecCopy(9, (double *)rot, (double *)(iso->rotation));
}

/* Sets the translation. */
void isoSetTranslation(isoIsometry *iso, const double transl[3]) {
    vecCopy(3, transl, iso->translation);
}

/* Applies the rotation and translation to a point. The output CANNOT safely
alias the input. */
void isoTransformPoint(const isoIsometry *iso, const double p[3], double isoP[3]) {
    mat331Multiply(iso->rotation, p, isoP);
    vecAdd(3, isoP, iso->translation, isoP);
}



/* Applies the rotation to a vector. The output CANNOT safely alias the input.
*/
void isoRotateVector(const isoIsometry *iso, const double v[3], double rotV[3]) {
    mat331Multiply(iso->rotation, v, rotV);
}

/* Applies the inverse rotation to a vector. The output CANNOT safely alias the
input. */
void isoUnrotateVector(
        const isoIsometry *iso, const double rotV[3], double v[3]) {
    mat331TransposeMultiply(iso->rotation, rotV, v);
}

/* Fills homog with the homogeneous version of the isometry. */
void isoGetHomogeneous(const isoIsometry *iso, double homog[4][4]) {
    mat44Isometry(iso->rotation, iso->translation, homog);
}

/* Fills homog with the homogeneous version of the inverse isometry. That is,
the product of this matrix and the one from isoGetHomogeneous is the identity
matrix. */
void isoGetInverseHomogeneous(const isoIsometry *iso, double homogInv[4][4]) {
   
    
    homogInv[0][0] = iso->rotation[0][0];
    homogInv[0][1] = iso->rotation[1][0];
    homogInv[0][2] = iso->rotation[2][0];
    homogInv[0][3] = -iso->rotation[0][0]*iso->translation[0] - iso->rotation[1][0]*iso->translation[1] - iso->rotation[2][0]*iso->translation[2];
    
    homogInv[1][0] = iso->rotation[0][1];
    homogInv[1][1] = iso->rotation[1][1];
    homogInv[1][2] = iso->rotation[2][1];
    homogInv[1][3] = -iso->rotation[0][1]*iso->translation[0] - iso->rotation[1][1]*iso->translation[1] - iso->rotation[2][1]*iso->translation[2];
    
    homogInv[2][0] = iso->rotation[0][2];
    homogInv[2][1] = iso->rotation[1][2];
    homogInv[2][2] = iso->rotation[2][2];
    homogInv[2][3] = -iso->rotation[0][2]*iso->translation[0] - iso->rotation[1][2]*iso->translation[1] - iso->rotation[2][2]*iso->translation[2];
    
    homogInv[3][0] = 0;
    homogInv[3][1] = 0;
    homogInv[3][2] = 0;
    homogInv[3][3] = 1;
        
    
    
    
}

/* Applies the inverse of the isometry to a point. If you transform a point and
then untransform the result, then you recover the original point. Similarly, if
you untransform a point and then transform the result, then you recover the
original point. The output CANNOT safely alias the input. */
void isoUntransformPoint(const isoIsometry *iso, const double isoP[3], double p[3]) {
    
    double untranslated[3];
    vecSubtract(3,isoP,iso->translation,untranslated);
    mat331TransposeMultiply(iso->rotation, untranslated, p);
    //printf("p: %f,%f,%f\n", p[0], p[1],p[2]);
    
}
