//
//  cube.c
//  
//
//  Created by Simon Boigelot & Edward Nicol on 29/04/2015.
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





typedef struct
{
    double *R;
    double *R0;
    double *D;
    double *D0;
    double *AD;
    double *AD0;
    double error;
    int size;
    int iter;
} femIterativeSolver;


void femIterativeSolverInit(femIterativeSolver *mySolver)
{
    int i;
    mySolver->iter = 0;
    mySolver->error = 10.0e+12;
    for (i=0 ; i < mySolver->size*6 ; i++)
        mySolver->R[i] = 0;
}

femIterativeSolver *femIterativeSolverCreate(int size)
{
    femIterativeSolver *mySolver = malloc(sizeof(femIterativeSolver));
    mySolver->R = malloc(sizeof(double)*size*6);
    mySolver->R0 = mySolver->R + size;
    mySolver->D  = mySolver->R + size*2;
    mySolver->D0 = mySolver->R + size*3;
    mySolver->AD  = mySolver->R + size*4;
    mySolver->AD0 = mySolver->R + size*5;
    mySolver->size = size;
    femIterativeSolverInit(mySolver);
    return(mySolver);
}

void femIterativeSolverFree(femIterativeSolver *mySolver)
{
    free(mySolver->R);
    free(mySolver);
}

int femIterativeSolverConverged(femIterativeSolver *mySolver)
{
    int  testConvergence = 0;
    if (mySolver->iter  > 3000)     testConvergence = -1;
    if (mySolver->error < 10.0e-6)  testConvergence = 1;
    return(testConvergence);
}





typedef struct
{
    femIterativeSolver* theSolver;
    femMesh* theMesh;
    int size;
} femProblem;

femProblem* femProblemCreate(const char *meshFileName)
{
    femProblem *theProblem = malloc(sizeof(femProblem));
    
    theProblem->theMesh = femMeshRead(meshFileName);
    theProblem->size = theProblem->theMesh->nElem*12;
    theProblem->theSolver = femIterativeSolverCreate(theProblem->size);
    
    return theProblem;
}

void femProblemFree(femProblem* theProblem)
{
    femMeshFree(theProblem->theMesh);
    femIterativeSolverFree(theProblem->theSolver);
    free (theProblem);
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


void cubeEliminate (femProblem* theProblem)
{
    
}


void cubeCompute(double alpha, double E, double nu, const char *meshFileName, double *U, double *V, double *W)
{
    femProblem* theProblem = femProblemCreate(meshFileName);
    femMesh* theMesh = theProblem->theMesh;
    femIterativeSolver* theSolver = theProblem->theSolver;
    int size = theProblem->size;
    
    int iElem,i,j, nLocalNode = theMesh->nLocalNode, index;
    double X[nLocalNode], Y[nLocalNode], Z[nLocalNode], Uloc[nLocalNode], Vloc[nLocalNode], Wloc[nLocalNode];
    int map[nLocalNode];
    double xsi, eta, zet, weight;
    double dxdxsi, dxdeta, dxdzet, dydxsi, dydeta, dydzet, dzdxsi, dzdeta, dzdzet;
    double ajac;
    
    double dphidxsi[nLocalNode], dphideta[nLocalNode], dphidzet[nLocalNode];
    double dphidx[nLocalNode], dphidy[nLocalNode], dphidz[nLocalNode];
    
    /*
    double **Aloc = malloc(sizeof(double*)*12);
    double A11[nLocalNode][nLocalNode], A12[nLocalNode][nLocalNode], A13[nLocalNode][nLocalNode], A22[nLocalNode][nLocalNode], A23[nLocalNode][nLocalNode], A33[nLocalNode][nLocalNode];*/
    
    const double a = E/(2*(1+nu));
    const double b = 2*a;
    const double c = E*nu/((1+nu)*(1-2*nu));
    
    double dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz;
    double dD1dx, dD1dy, dD1dz, dD2dx, dD2dy, dD2dz, dD3dx, dD3dy, dD3dz;
    
    double *D0 = theSolver->D0;
    double *AD0 = theSolver->AD0;
    double *R0 = theSolver->R0;
    
    
    for (iElem=0; iElem<size; iElem++)
    {
        /*
        for (i=0; i<nLocalNode; i++)// remise a zero des mini-matrices pour l'element nouveau
        {
            for (j=0; j<nLocalNode; j++)
            {
                A11[i][j]=0;
                A12[i][j]=0;
                A13[i][j]=0;
                A22[i][j]=0;
                A23[i][j]=0;
                A33[i][j]=0;
            }
        }*/
        
        for (i=0; i<nLocalNode; i++)
        {
            map[i] = theMesh->elem[iElem*nLocalNode+i];
            
            X[i] = theMesh->X[map[i]];
            Y[i] = theMesh->Y[map[i]];
            Z[i] = theMesh->Z[map[i]];
            
            Uloc[i] = U[map[i]];
            Vloc[i] = V[map[i]];
            Wloc[i] = W[map[i]];
        }
        
        
        xsi = 1.0/4.0;
        eta = 1.0/4.0;
        zet = 1.0/4.0;
        weight = 1.0/6.0;
        
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
        
        for(i=0; i<nLocalNode; i++)
        {
            dxdxsi += X[i] * dphidxsi[i];
            dydxsi += Y[i] * dphidxsi[i];
            dzdxsi += Z[i] * dphidxsi[i];
            dxdeta += X[i] * dphideta[i];
            dydeta += Y[i] * dphideta[i];
            dzdeta += Z[i] * dphideta[i];
            dxdzet += X[i] * dphidzet[i];
            dydzet += Y[i] * dphidzet[i];
            dzdzet += Z[i] * dphidzet[i];
        }
        
        ajac = dxdxsi * dydeta * dzdzet + dxdeta * dydzet * dzdxsi + dxdzet * dydxsi * dzdeta - dzdxsi * dydeta * dxdzet - dxdxsi * dzdeta * dydzet - dydxsi * dxdeta * dzdzet;
        ajac = fabs(ajac);
        
        for (i=0; i<nLocalNode; i++)
        {
            dphidx[i] = (dphidxsi[i]*(dydeta*dzdzet-dydzet*dzdeta)+dphideta[i]*(dydzet*dzdxsi-dydxsi*dzdzet)+dphidzet[i]*(dydxsi*dzdeta-dydeta*dzdxsi))/ajac;
            dphidy[i] = (dphidxsi[i]*(dxdzet*dzdeta-dxdeta*dzdzet)+dphideta[i]*(dxdxsi*dzdzet-dxdzet*dzdxsi)+dphidzet[i]*(dxdeta*dzdxsi-dxdxsi*dzdeta))/ajac;
            dphidz[i] = (dphidxsi[i]*(dxdeta*dydzet-dxdzet*dydeta)+dphideta[i]*(dxdzet*dydxsi-dxdxsi*dydzet)+dphidzet[i]*(dxdxsi*dydeta-dxdeta*dydxsi))/ajac;
        }
        
        
        if (theSolver->iter==0)
        {
            dudx=0, dudy=0, dudz=0, dvdx=0, dvdy=0, dvdz=0, dwdx=0, dwdy=0, dwdz=0;
            
            for (i=0; i<nLocalNode; i++)
            {
                dudx += Uloc[i]*dphidx[i];
                dudy += Uloc[i]*dphidy[i];
                dudz += Uloc[i]*dphidz[i];
            
                dvdx += Vloc[i]*dphidx[i];
                dvdy += Vloc[i]*dphidy[i];
                dvdz += Vloc[i]*dphidz[i];
                
                dwdx += Wloc[i]*dphidx[i];
                dwdy += Wloc[i]*dphidy[i];
                dwdz += Wloc[i]*dphidz[i];
            }
        
            for (i=0; i<nLocalNode; i++)
            {
                index = map[i];
                D0[index] -= (b*dphidx[i]*dudx + c*(dphidx[i]*(dudx+dvdy+dwdz))+a*(dphidy[i]*(dudy+dvdx))+a*(dphidz[i]*(dudz+dwdx)))*ajac*weight;
            
                D0[index+size/3] -= (a*(dphidx[i]*(dudy+dvdx))+b*dphidy[i]*dvdy+c*(dphidy[i]*(dudx+dvdy+dwdz))+a*(dphidz[i]*(dvdz+dwdx)))*ajac*weight;
            
                D0[index+size*2/3] -= (a*(dphidx[i]*(dudz+dwdx))+a*(dphidy[i]*(dvdz+dwdy))+b*dphidz[i]*dwdz+c*(dphidz[i]*(dudx+dvdy+dwdz)))*ajac*weight;
            }
        
            for (i=0; i<nLocalNode; i++)
            {
                index = map[i];
                R0[index] = D0[index];
                R0[index+size/3] = D0[index+size/3];
                R0[index+size*2/3] = D0[index+size*2/3];
            }
        }
        
        dD1dx=0, dD1dy=0, dD1dz=0, dD2dx=0, dD2dy=0, dD2dz=0, dD3dx=0, dD3dy=0, dD3dz=0;
        
        for (i=0; i<nLocalNode; i++)
        {
            index = map[i];
            
            dD1dx += D0[index]*dphidx[i];
            dD1dy += D0[index]*dphidy[i];
            dD1dz += D0[index]*dphidz[i];
            
            dD2dx += D0[index+size/3]*dphidx[i];
            dD2dy += D0[index+size/3]*dphidy[i];
            dD2dz += D0[index+size/3]*dphidz[i];
            
            dD3dx += D0[index+size*2/3]*dphidx[i];
            dD3dy += D0[index+size*2/3]*dphidy[i];
            dD3dz += D0[index+size*2/3]*dphidz[i];
        }
        
        for (i=0; i<nLocalNode; i++)
        {
            index = map[i];
            AD0[index] += (b*dphidx[i]*dD1dx + c*(dphidx[i]*(dD1dx+dD2dy+dD3dz))+a*(dphidy[i]*(dD1dy+dD2dx))+a*(dphidz[i]*(dD1dz+dD3dx)))*ajac*weight;
            
            AD0[index+size/3] += (a*(dphidx[i]*(dD1dy+dD2dx))+b*dphidy[i]*dD2dy+c*(dphidy[i]*(dD1dx+dD2dy+dD3dz))+a*(dphidz[i]*(dD2dz+dD3dx)))*ajac*weight;
            
            AD0[index+size*2/3] += (a*(dphidx[i]*(dD1dz+dD3dx))+a*(dphidy[i]*(dD2dz+dD3dy))+b*dphidz[i]*dD3dz+c*(dphidz[i]*(dD1dx+dD2dy+dD3dz)))*ajac*weight;
        }
        
        
        /*
        for (i=0; i<nLocalNode; i++) // remplissage des mini-matrices :-)
        {
            for (j=0; j<nLocalNode; j++)
            {
                A11[i][j] += ((b+c)*dphidx[i]*dphidx[j] + a*dphidy[i]*dphidy[j] + a*dphidz[i]*dphidz[j]) * weight * ajac;
                A12[i][j] += (c*dphidx[i]*dphidy[j] + a*dphidy[i]*dphidx[j]) * weight * ajac;
                A13[i][j] += (c*dphidx[i]*dphidz[j] + a*dphidz[i]*dphidx[j]) * weight * ajac;
                
                A22[i][j] += ((b+c)*dphidy[i]*dphidy[j] + a*dphidx[i]*dphidx[j] + a*dphidz[i]*dphidz[j]) * weight * ajac;
                A23[i][j] += (c*dphidy[i]*dphidz[j] + a*dphidz[i]*dphidy[j]) * weight * ajac;
                
                A33[i][j] += ((b+c)*dphidz[i]*dphidz[j] + a*dphidy[i]*dphidy[j] + a*dphidx[i]*dphidx[j]) * weight * ajac;
            }
        }
        
        for (i=0; i<nLocalNode; i++)
        {
            for (j=0; j<nLocalNode; j++)
            {
                Aloc[i][j] = A11[i][j];// premiere ligne de Aloc
                Aloc[i][j+4] = A12[i][j];
                Aloc[i][j+8] = A13[i][j];
                
                Aloc[i+4][j] = A12[i][j];// deuxieme ligne de Aloc
                Aloc[i+4][j+4] = A22[i][j];
                Aloc[i+4][j+8] = A23[i][j];
                
                Aloc[i+8][j] = A13[i][j];// troisieme ligne de Aloc
                Aloc[i+8][j+4] = A23[i][j];
                Aloc[i+8][j+8] = A33[i][j];
            }
        }*/
        
    }
    femProblemFree(theProblem);
}



