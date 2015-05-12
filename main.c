
#include <stdio.h>
#include <stdlib.h>

#define GLFW_INCLUDE_GLU
#include <GLFW/glfw3.h>

#include "cube.c"


int main( void )
{
    int nElem,nNode,i,j,k,index,trash,*elem;
    double *X,*Y,*Z;
    
    char *filename = "cube65.txt";
    
    const double alpha = -1/10, E=1e-3, nu = 0.3;
    
    FILE* file = fopen(filename,"r");
    if (file == NULL) {
    	printf("Error : cannot open mesh file :-) \n");
        exit(0); }
 	fscanf(file, "Number of nodes %d \n", &nNode);
  	X = malloc(sizeof(double)*nNode);
  	Y = malloc(sizeof(double)*nNode);
  	Z = malloc(sizeof(double)*nNode);
	for (i = 0; i < nNode; i++) 
    	fscanf(file,"%d : %le %le %le  \n",&trash,&X[i],&Y[i],&Z[i]); 
    fscanf(file, "Number of tetrahedra %d \n", &nElem); 
  	elem = malloc(sizeof(int)*4*nElem);
  	for (i = 0; i < nElem; i++)  {   
     	fscanf(file,"%d : %d %d %d %d \n", &trash,&elem[i*4],&elem[i*4+1],&elem[i*4+2],&elem[i*4+3]); }
  	fclose(file);
    
    double *U = malloc(sizeof(double)*nNode);
    double *V = malloc(sizeof(double)*nNode);
    double *W = malloc(sizeof(double)*nNode);
    
    for (i=0; i<nNode; i++)//initialisation des vecteurs U, V et W
    {
        U[i] = 1;
        V[i] = 1;
        W[i] = 1;
    }
    
    cubeCompute(alpha, E, nu, filename, U, V, W);
    
 /*   file = fopen("cube23000.txt","w");
    fprintf(file, "Number of nodes %d \n", nNode);
    for (i = 0; i < nNode; ++i) {
        fprintf(file,"%6d : %14.7e %14.7e %14.7e \n",i,X[i],Y[i],1.0-Z[i]); }
    
    fprintf(file, "Number of tetrahedra %d \n", nElem);  
    for (i = 0; i < nElem; ++i) {
        fprintf(file,"%6d : %6d %6d %6d %6d \n", i,elem[i*4],elem[i*4+1],elem[i*4+2],elem[i*4+3]);   }    
    fclose(file);
 */   

   	glfwInit();
	GLFWwindow* window =  glfwCreateWindow(480,480,"MECA1120 : Project",NULL,NULL);
    glfwMakeContextCurrent(window);//rajoute
    glShadeModel(GL_SMOOTH);//rajoute
    
    GLfloat mat_specular[]   = { 1.0, 1.0, 1.0, 0.0 };
    GLfloat light_position[] = { -18.0, -18.0, -18.0, 0.0 };    
	glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
    glLightfv(GL_LIGHT0, GL_POSITION, light_position);
    GLfloat light_radiance[] = {1., 1., 1., 1.};
    glLightfv(GL_LIGHT0, GL_DIFFUSE, light_radiance);
    glLightfv(GL_LIGHT0, GL_SPECULAR, light_radiance);
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glDepthFunc(GL_LEQUAL);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_COLOR_MATERIAL);
    glEnable(GL_NORMALIZE);	
    
   
    do {       
        int width,height;
        glfwGetFramebufferSize(window, &width, &height );
        height = height > 0 ? height : 1;
        glViewport( 0, 0, width, height );

        glClearColor( 0.9f, 0.9f, 0.8f, 0.0f );
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        gluPerspective(65.0f,(GLfloat)width/(GLfloat)height,1.0f,100.0f);

        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
        gluLookAt(0.0f,0.0f,15.0f,0.0f, 30.0f, 0.0f,0.0f,0.0f,1.0f);  
        glTranslatef(0.0f,14.0f,8.0f);
        int rotation = 389;
        glRotatef(0.3f*(GLfloat)rotation,0.0f,0.0f,1.0f);
              
        int map[4][3] = {{0,2,1},{0,1,3},{1,2,3},{2,0,3}};
        GLfloat colors[9], coord[9];

        for (i=0; i < nElem; ++i) {
            for (k=0; k < 4; ++k) {
            
            // Check for face to color in red :-)
            
                int test = 0;
                for (j=0; j < 3; ++j) {                       
                    index = elem[4*i+map[k][j]];
                    if ((Z[index] >= 1.0 - 1.0e-6) && (X[index]-0.4 <= 1.0e-6) && (Y[index]-0.4 <= 1.0e-6))
                        test = test + 1;} 
                        
            // Draw all faces
              
                for (j=0; j < 3; ++j) {                    
                    index = elem[4*i+map[k][j]];
                    if (test == 3) {
                        colors[j*3+0] = 1.0;
                        colors[j*3+1] = 0.0;
                        colors[j*3+2] = 0.0; }
                    else {
                        colors[j*3+0] = 1.0;
                        colors[j*3+1] = 1.0;
                        colors[j*3+2] = 1.0; }
                    coord[j*3+0] = (X[index]-0.5+U[index])*8;
                    coord[j*3+1] = (Y[index]-0.5+V[index])*8;
                    coord[j*3+2] = (Z[index]-0.5+W[index])*8;  }
  
                glEnableClientState(GL_VERTEX_ARRAY);
                glEnableClientState(GL_COLOR_ARRAY);
                glEnableClientState(GL_NORMAL_ARRAY);
                glVertexPointer(3, GL_FLOAT, 0, coord);
                glNormalPointer(GL_FLOAT, 0, coord);
                glColorPointer(3, GL_FLOAT, 0, colors);
                glDrawArrays(GL_TRIANGLES, 0, 3);
                glDisableClientState(GL_NORMAL_ARRAY);    
                glDisableClientState(GL_COLOR_ARRAY);
                glDisableClientState(GL_VERTEX_ARRAY);      
                             
                glColor3f(0.0, 0.0, 0.0);
                glEnableClientState(GL_VERTEX_ARRAY);
                for (j=0; j < 9; ++j)
                    coord[j] = coord[j] * 1.001;
                glVertexPointer(3, GL_FLOAT, 0, coord);
                glDrawArrays(GL_LINE_LOOP, 0, 3);
                glDisableClientState(GL_VERTEX_ARRAY); }}
        glfwSwapBuffers(window);
        glfwPollEvents();
   
    } while( glfwGetKey(window, GLFW_KEY_ESCAPE ) != GLFW_PRESS &&
       (!glfwWindowShouldClose(window)) );


    free(X);
    free(Y);
    free(Z);
    free(U);
    free(V);
    free(W);
    free(elem);
    
    glfwTerminate();
    exit( EXIT_SUCCESS );
}



