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
    int elem[2];
    int node[3];
}femEdge;



typedef struct
{
    femMesh *theMesh;
    femEdge *edges;
    int nEdge;
    int nBoundary;
}femEdges;


int fmid (int n1, int n2, int n3)
{
    int diff1 = n1 - n2;
    int diff2 = n1 - n3;
    
    if ((diff1>0 || diff2>0) && diff1*diff2 <0)
    {
        return n1;
    }
    else if (diff1 < 0 && diff2 < 0)
    {
        return fmin(n2, n3);
    }
    else
    {
        return fmax(n2, n3);
    }
}

int femEdgesCompare(const void *edgeOne, const void *edgeTwo)
{
    int *nodeOne = ((femEdge*) edgeOne)->node;
    int *nodeTwo = ((femEdge*) edgeTwo)->node;
    
    int diffMin = fmin(fmin(nodeOne[0],nodeOne[1]),nodeOne[2]) - fmin(fmin(nodeTwo[0],nodeTwo[1]),nodeTwo[2]);
    int diffMax = fmax(fmax(nodeOne[0],nodeOne[1]),nodeOne[2]) - fmax(fmax(nodeTwo[0],nodeTwo[1]),nodeTwo[2]);
    int diffMid = fmid(nodeOne[0], nodeOne[1], nodeOne[2]) - fmid(nodeTwo[0], nodeTwo[1], nodeTwo[2]);
    
    if (diffMin < 0)    return  1;
    if (diffMin > 0)    return -1;
    if (diffMid < 0)    return  1;
    if (diffMid > 0)    return -1;
    if (diffMax < 0)    return  1;
    if (diffMax > 0)    return -1;
    return  0;
}

femEdges* femEdgesCreate(femMesh *theMesh)
{
    femEdges *theEdges = malloc(sizeof(femEdges));
    int nLoc = theMesh->nLocalNode;
    int i,j,n = theMesh->nElem * nLoc;
    femEdge* edges = malloc(n * sizeof(femEdge));
    theEdges->theMesh  = theMesh;
    theEdges->edges = edges;
    theEdges->nEdge = n;
    theEdges->nBoundary = n;
    
    for (i = 0; i < theMesh->nElem; i++)
    {
        int *elem = &(theMesh->elem[i*nLoc]);
        for (j = 0; j < nLoc; j++)
        {
            int id = i * nLoc + j;
            edges[id].elem[0] = i;
            edges[id].elem[1] = -1;
            
            edges[id].node[0] = elem[j];
            edges[id].node[1] = elem[(j + 1) % nLoc];
            edges[id].node[2] = elem[(j + 2) % nLoc];
        }
    }
    
    qsort(theEdges->edges, theEdges->nEdge, sizeof(femEdge), femEdgesCompare);
    
    int index = 0;
    int nBoundary = 0;
    
    for (i=0; i < theEdges->nEdge; i++)
    {
        if (i == theEdges->nEdge - 1 || femEdgesCompare(&edges[i],&edges[i+1]) != 0)
        {
            edges[index] = edges[i];
            nBoundary++;
        }
        else
        {
            edges[index] = edges[i];
            edges[index].elem[1] = edges[i+1].elem[0];
            i ++;
        }
        index++;
    }
    
    theEdges->edges = realloc(edges, index * sizeof(femEdge));
    theEdges->nEdge = index;
    theEdges->nBoundary = nBoundary;
    return theEdges;
}

void femEdgesFree(femEdges* theEdges)
{
    free(theEdges->edges);
    free(theEdges);
}

typedef struct
{
    femIterativeSolver* theSolver;
    femMesh* theMesh;
    int size;
    double* soluce;
    femEdges* theEdges;
} femProblem;

femProblem* femProblemCreate(const char *meshFileName)
{
    femProblem *theProblem = malloc(sizeof(femProblem));
    
    theProblem->theMesh = femMeshRead(meshFileName);
    theProblem->size = theProblem->theMesh->nNode*3;
    theProblem->theSolver = femIterativeSolverCreate(theProblem->size);
    theProblem->theEdges = femEdgesCreate(theProblem->theMesh);
    
    theProblem->soluce = malloc(sizeof(double)*theProblem->size);
    
    int i;
    for (i = 0; i < theProblem->size; i++)
    {
        theProblem->soluce[i] = 0;
    }
    
    return theProblem;
}

void femProblemFree(femProblem* theProblem)
{
    femMeshFree(theProblem->theMesh);
    femIterativeSolverFree(theProblem->theSolver);
    femEdgesFree(theProblem->theEdges);
    free(theProblem->soluce);
    free(theProblem);
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

void tetrahedraDphi(double xsi, double eta, double zet, double dphidxsi[4],double dphideta[4],double dphidzet[4])
{
    dphidxsi[0] = -1.0; dphideta[0] = -1.0; dphidzet[0] = -1.0;
    dphidxsi[1] =  1.0; dphideta[1] =  0.0; dphidzet[1] =  0.0;
    dphidxsi[2] =  0.0; dphideta[2] =  1.0; dphidzet[2] =  0.0;
    dphidxsi[3] =  0.0; dphideta[3] =  0.0; dphidzet[3] =  1.0;
}

/*
 * direction:
 * 0 = u
 * 1 = v
 * 2 = w
 * 3 = contrainte de colonne
 */
void cubeConstrain (femProblem* theProblem, int myNode, int myValue, int direction)
{
    femIterativeSolver* theSolver = theProblem->theSolver;
    
    
    if (direction == 3)
    {
        theSolver->R0[myNode] = 0;
        theSolver->D0[myNode] = 0;
        theSolver->AD0[myNode] = 0;
        
        theSolver->R0[myNode+theProblem->theMesh->nNode] = 0;
        theSolver->D0[myNode+theProblem->theMesh->nNode] = 0;
        theSolver->AD0[myNode+theProblem->theMesh->nNode] = 0;
        
        theSolver->R0[myNode+2*theProblem->theMesh->nNode] = myValue;
        theSolver->D0[myNode+2*theProblem->theMesh->nNode] = myValue;
        theSolver->AD0[myNode+2*theProblem->theMesh->nNode] = myValue;
    }
    else
    {
        theSolver->R0[myNode+direction*theProblem->theMesh->nNode] = myValue;
        theSolver->D0[myNode+direction*theProblem->theMesh->nNode] = myValue;
        theSolver->AD0[myNode+direction*theProblem->theMesh->nNode] = myValue;
    }
}


double MaxError(double *tab, int size)
{
    int i;
    double max = tab[0];
    for(i=1;i<size;i++)
    {
        if(tab[i] >= max) max = tab[i];
    }
    return max;
}



/*
 * Elimination suivant la methode des gradients conjugues
 * Copier-Coller de la soluce du devoir 5 :-)
 */
double *femIterativeSolverEliminate(femIterativeSolver *mySolver)
{
    mySolver->iter++;
    //double error = 0.0;
    double residus = 0.0;
    double error[mySolver->size];
    int i;
    double denAlpha = 0.0;
    
    
    for (i=0; i < mySolver->size; i++)
    {
        //error += (mySolver->R0[i])*(mySolver->R0[i]);
        residus += (mySolver->R0[i])*(mySolver->R0[i]);
        error[i] = (mySolver->R0[i])-(mySolver->R0[i]);
        denAlpha += mySolver->D0[i] * mySolver->AD0[i];
    }
    //double alpha = error/denAlpha;
    double alpha = residus/denAlpha;
    
    if (mySolver->iter == 1)
    {
        for (i=0; i < mySolver->size; i++)
        {
            mySolver->R[i] = 0.0;
        }
    }
    
    else
    {
        double numBeta = 0.0;
        for (i=0; i < mySolver->size; i++)
        {
            mySolver->R[i] = alpha * mySolver->D0[i];
            mySolver->R0[i] = mySolver->R0[i] - alpha * mySolver->AD0[i];
            numBeta += mySolver->R0[i] * mySolver->R0[i];
            mySolver->AD[i] = 0.0;
        }
        //double beta = numBeta/error;
        double beta = numBeta/residus;
        
        for (i=0; i < mySolver->size; i++)
        {
            mySolver->AD0[i] = 0.0;
            mySolver->D0[i] = mySolver->R0[i] + beta * mySolver->D0[i];
        }
    }
    //mySolver->error = sqrt(error);
    mySolver->error = MaxError(error,mySolver->size);
    return(mySolver->R);
    
}



void cubeEliminate (femProblem* theProblem, double alpha, double E, double nu, double *U, double *V, double *W)
{
    femMesh* theMesh = theProblem->theMesh;
    femIterativeSolver* theSolver = theProblem->theSolver;
    int size = theProblem->size;
    
    int iElem,i,iEdge, nLocalNode = theMesh->nLocalNode, index;
    double X[nLocalNode], Y[nLocalNode], Z[nLocalNode], Uloc[nLocalNode], Vloc[nLocalNode], Wloc[nLocalNode];
    int map[nLocalNode];
    double xsi, eta, zet, weight;
    double dxdxsi, dxdeta, dxdzet, dydxsi, dydeta, dydzet, dzdxsi, dzdeta, dzdzet;
    double ajac;
    
    double dphidxsi[nLocalNode], dphideta[nLocalNode], dphidzet[nLocalNode];
    double dphidx[nLocalNode], dphidy[nLocalNode], dphidz[nLocalNode];
    
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
        else
        {
            dD1dx=0, dD1dy=0, dD1dz=0, dD2dx=0, dD2dy=0, dD2dz=0, dD3dx=0, dD3dy=0, dD3dz=0;
            
            for (i=0; i<nLocalNode; i++)
            {
                index = map[i];
                
                dD1dx += D0[index]*dphidx[i]; //partie "en U"
                dD1dy += D0[index]*dphidy[i];
                dD1dz += D0[index]*dphidz[i];
                
                dD2dx += D0[index+size/3]*dphidx[i]; //partie "en V"
                dD2dy += D0[index+size/3]*dphidy[i];
                dD2dz += D0[index+size/3]*dphidz[i];
                
                dD3dx += D0[index+size*2/3]*dphidx[i]; //partie "en W"
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
        }
    } // Fin assemblage -> fin boucle iElem
    
    /*
     *
     * Appliquer les contraintes sur D0 R0 AD0 :-)
     *
     *
     */
    
    for (i=0; i<size; i++)
    {
        printf("R0[%d] = %f\n",i, R0[i]);
    }
    
    
    for (iEdge=0; iEdge<theProblem->theEdges->nEdge; iEdge++)
    {
        femEdge edge = theProblem->theEdges->edges[iEdge];
        if (edge.elem[1]==-1)
        {
            for (i=0; i<3; i++)
            {
                if (theMesh->X[edge.node[i] == 0])
                {
                    cubeConstrain(theProblem, edge.node[i], 0, 0);
                }
                
                if (theMesh->Y[edge.node[i] == 0])
                {
                    cubeConstrain(theProblem, edge.node[i], 0, 1);
                }
                
                if (theMesh->Z[edge.node[i] == 0])
                {
                    cubeConstrain(theProblem, edge.node[i], 0, 2);
                }
                
                if (theMesh->Z[edge.node[i]] == 1 && theMesh->X[edge.node[i]]<=4/10 && theMesh->Y[edge.node[i]]<=4/10)
                {
                    cubeConstrain(theProblem, edge.node[i], alpha, 3);
                }
            }
        }
    }
    
    /*
    for (i=0; i<size; i++)
    {
        printf("R0[%d] = %f\n",i, R0[i]);
    }*/
    
    
    /*
     * Elimination suivant la methode des gradients conjugues
     * Copier-Coller de la soluce du devoir 5 :-)
     */
    double *soluce = femIterativeSolverEliminate(theSolver);
    for (i = 0; i < theProblem->theMesh->nNode; i++)
    {
        // vecteur a redecomposer en U, V, W a la fin
        theProblem->soluce[i] += soluce[i];
    }
}


void cubeCompute(double alpha, double E, double nu, const char *meshFileName, double *U, double *V, double *W)
{
    femProblem* theProblem = femProblemCreate(meshFileName);
    int testConvergence, i, nNode = theProblem->theMesh->nNode, iter=0;
    
    do
    {
        cubeEliminate(theProblem, alpha, E, nu, U, V, W);
        testConvergence = femIterativeSolverConverged(theProblem->theSolver);
        iter++;
    }
    while (testConvergence==0);
    
    if(testConvergence==1) printf("Cool :-)\n");
    if(testConvergence == -1) printf("No convergence in %d steps",3000);
    /*
    for (i=0; i<nNode; i++)
    {
        printf("U[%d] = %f\n",i, U[i]);
    }*/
    
    for (i=0; i<nNode; i++)
    {
        U[i] = theProblem->soluce[i];
        V[i] = theProblem->soluce[i+nNode];
        W[i] = theProblem->soluce[i+2*nNode];
    }
    
    printf("Number of iterations %d\n", iter);
    /*
    for (i=0; i<nNode; i++)
    {
        printf("U[%d] = %f\n",i, U[i]);
    }*/
    
    femProblemFree(theProblem);
    
}