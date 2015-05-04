//
//  cube.c
//  
//
//  Created by Edward Nicol on 29/04/2015.
//
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

typedef struct {
    int *elem;
    double *X;
    double *Y;
    double *Z;
    int nElem;
    int nNode;
    int nLocalNode;
} femMesh;

femMesh *femMeshRead(const char *filename)
{
    femMesh *theMesh = malloc(sizeof(femMesh));
    
    int i,trash,*elem;
    
    FILE* file = fopen(filename,"r");
    if (file == NULL) printf ("No mesh file !");
    
    if (!fscanf(file, "Number of nodes %d \n", &theMesh->nNode)) printf("Cannot read the number of nodes\n");
    theMesh->X = malloc(sizeof(double)*theMesh->nNode);
    theMesh->Y = malloc(sizeof(double)*theMesh->nNode);
    theMesh->Z = malloc(sizeof(double)*theMesh->nNode);
    for (i = 0; i < theMesh->nNode; ++i) {
        if (!fscanf(file,"%d : %le %le %le \n",&trash,&theMesh->X[i],&theMesh->Y[i], &theMesh->Z[i])) printf("Cannot read node %d\n", i); }
    
    char str[256]; if (!fgets(str, sizeof(str), file)) printf("Error while reading the mesh file !\n");
    if (!strncmp(str,"Number of triangles",19))  {
        if (!sscanf(str,"Number of triangles %d \n", &theMesh->nElem)) printf("Cannot read the number of elements\n");
        theMesh->elem = malloc(sizeof(int)*3*theMesh->nElem);
        theMesh->nLocalNode = 3;
        for (i = 0; i < theMesh->nElem; ++i) {
            elem = &(theMesh->elem[i*3]);
            if(!fscanf(file,"%d : %d %d %d\n", &trash,&elem[0],&elem[1],&elem[2])) printf("Cannot read element %d\n", i); }}
    else if (!strncmp(str,"Number of quads",15))  {
        if (!sscanf(str,"Number of quads %d \n", &theMesh->nElem))printf("Cannot read the number of elements\n");
        theMesh->elem = malloc(sizeof(int)*4*theMesh->nElem);
        theMesh->nLocalNode = 4;
        for (i = 0; i < theMesh->nElem; ++i) {
            elem = &(theMesh->elem[i*4]);
            if (!fscanf(file,"%d : %d %d %d %d\n", &trash,&elem[0],&elem[1],&elem[2],&elem[3])) printf("Cannot read element %d\n", i); }}
    else if (!strncmp(str,"Number of tetrahedra",20))  {
        if (!sscanf(str,"Number of tetrahedra %d \n", &theMesh->nElem)) printf("Cannot read the number of elements\n");
        theMesh->elem = malloc(sizeof(int)*4*theMesh->nElem);
        theMesh->nLocalNode = 4;
        for (i = 0; i < theMesh->nElem; ++i) {
            elem = &(theMesh->elem[i*4]);
            if (!fscanf(file,"%d : %d %d %d %d\n", &trash,&elem[0],&elem[1],&elem[2],&elem[3])) printf("Cannot read element %d\n", i); }}
    
    fclose(file);
    return theMesh;
}

void femMeshFree(femMesh *theMesh)
{
    free(theMesh->X);
    free(theMesh->Y);
    free(theMesh->Z);
    free(theMesh->elem);
    free(theMesh);
}

typedef struct {
    int n;
    double *xsi;
    double *eta;
    double *zet;
    double *weight;
} femIntegrator;

femIntegrator *femIntegratorCreate()
{
    int i;
    int n = 5;
    double xsi[5]    = { 0.25, 0.5, 0.166666666666666666667, 0.166666666666666666667, 0.166666666666666666667};
    double eta[5]    = { 0.25, 0.166666666666666666667, 0.5, 0.166666666666666666667, 0.166666666666666666667};
    double zet[5]    = { 0.25, 0.166666666666666666667, 0.166666666666666666667, 0.5, 0.166666666666666666667};
    double weight[5] = {-0.133333333333333, 0.075, 0.075, 0.075, 0.075};
    
    
    femIntegrator *theIntegrator = malloc(sizeof(femIntegrator));
    theIntegrator->n = n;
    theIntegrator->xsi    = malloc(sizeof(double) * n);
    theIntegrator->eta    = malloc(sizeof(double) * n);
    theIntegrator->zet    = malloc(sizeof(double) * n);
    theIntegrator->weight = malloc(sizeof(double) * n);
    
    for (i = 0; i < theIntegrator->n; ++i) {
        theIntegrator->xsi[i]    = xsi[i];
        theIntegrator->eta[i]    = eta[i];
        theIntegrator->zet[i]    = zet[i];
        theIntegrator->weight[i] = weight[i]; }
    
    return theIntegrator;
}

void femIntegratorFree(femIntegrator *theIntegrator)
{
    free(theIntegrator->xsi);
    free(theIntegrator->eta);
    free(theIntegrator->zet);
    free(theIntegrator->weight);
    free(theIntegrator);
}

void tetrahedraXsi(double Xsi[4], double Eta[4], double Zet[4])
{
    Xsi[0] = 0.0;   Eta[0] = 0.0,   Zet[0] = 0.0;
    Xsi[1] = 1.0;   Eta[1] = 0.0,   Zet[1] = 0.0;
    Xsi[2] = 0.0;   Eta[2] = 1.0,   Zet[2] = 0.0;
    Xsi[3] = 0.0;   Eta[3] = 0.0,   Zet[3] = 1.0;
}


void tetrahedraPhi(double xsi, double eta, double zet, double phi[4])
{
    phi[0] = 1.0 - xsi - eta - zet;
    phi[1] = xsi;
    phi[2] = eta;
    phi[3] = zet;
}

void tetrahedraDphi(double xsi, double eta, double zet,
                    double dphidxsi[4],double dphideta[4],double dphidzet[4])
{
    
    dphidxsi[0] = -1.0; dphideta[0] = -1.0; dphidzet[0] = -1.0;
    dphidxsi[1] =  1.0; dphideta[1] =  0.0; dphidzet[1] =  0.0;
    dphidxsi[2] =  0.0; dphideta[2] =  1.0; dphidzet[2] =  0.0;
    dphidxsi[3] =  0.0; dphideta[3] =  0.0; dphidzet[3] =  1.0;
    
}

double tetrahedraJacobian(double X[4], double Y[4], double Z[4])
{
    
    
    
    int i;
    
    
    return fabs(ajac);
}


void calcResidu(double* U, double* V, double* W, double* dphidx, double* dphidy, double* dphidz)
{
    
}


void cubeCompute(double alpha, double E, double nu, const char *meshFileName, double *U, double *V, double *W)
{
    femMesh* theMesh = femMeshRead(meshFileName);
    femIntegrator* theIntegrator = femIntegratorCreate();
    int size = 12*theMesh->nElem,iElem,i,j, nLocalNode = theMesh->nLocalNode;
    double X[nLocalNode], Y[nLocalNode], Z[nLocalNode], map[nLocalNode];
    double xsi, eta, zet, weight;
    double dxdxsi, dxdeta, dxdzet, dydxsi, dydeta, dydzet, dzdxsi, dzdeta, dzdzet;
    double ajac;
    
    double dphidxsi[nLocalNode], dphideta[nLocalNode], dphidzet[nLocalNode];
    
    
    for (iElem=0; iElem<size; iElem++)
    {
        for (i=0; i<nLocalNode; i++)
        {
            X[i] = theMesh->X[iElem*nLocalNode+i];
            Y[i] = theMesh->Y[iElem*nLocalNode+i];
            Z[i] = theMesh->Z[iElem*nLocalNode+i];
            map[i] = theMesh->elem[iElem*nLocalNode+i];
        }
        
        for (i=0; i<theIntegrator->n; i++)
        {
            xsi = theIntegrator->xsi[i];
            eta = theIntegrator->eta[i];
            zet = theIntegrator->zet[i];
            weight = theIntegrator->weight[i];
            
            tetrahedraDphi(xsi,eta,zet,dphidxsi,dphideta,dphidzet);
            dxdxsi = 0.0;
            dxdeta = 0.0;
            dxdzet = 0.0;
            dydxsi = 0.0;
            dydeta = 0.0;
            dydzet = 0.0;
            dzdxsi = 0.0;
            dzdeta = 0.0;
            dzdzet = 0.0;
            
            for(j=0; j<nLocalNode; j++) {
                dxdxsi += X[j] * dphidxsi[j];
                dydxsi += Y[j] * dphidxsi[j];
                dzdxsi += Z[j] * dphidxsi[j];
                dxdeta += X[j] * dphideta[j];
                dydeta += Y[j] * dphideta[j];
                dzdeta += Z[j] * dphideta[j];
                dxdzet += X[j] * dphidzet[j];
                dydzet += Y[j] * dphidzet[j];
                dzdzet += Z[j] * dphidzet[j]; }
            
            ajac = dxdxsi * dydeta * dzdzet + dxdeta * dydzet * dzdxsi + dxdzet * dydxsi * dzdeta - dzdxsi * dydeta * dxdzet - dxdxsi * dzdeta * dydzet - dydxsi * dxdeta * dzdzet;
            ajac = fabs(ajac);
        }
    }
}



