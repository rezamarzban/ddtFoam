#include "stdio.h"
#include "math.h"
#include "stdlib.h"

// scale to metres
#define f 1e-3
// channel length
#define L 3000.0
// number of obstacles
#define OBSTACLES 9
// obstacle spacing
#define s 300.0
// blockage ratio
#define BR 0.30
// where the first obstacle is placed:
#define firstObstacle 205.0
// channel height
#define h 60.0
// half thickness in z dimension
#define Z 1.0
// obstacle / wall thickness
#define d 12.0
// a small number to decide if cell centers coinciding with a box should be cut or not
#define bit 0.01


#define dx 2
#define dy 2 
#define zc 1



int main()
{
int vert=0; // vertex number
int i,j,k;
int xcells;
//int ycells;
int blocks=0;

//double L = (OBSTACLES+SMOOTH)*s;

FILE *bmd;
FILE *obs;

bmd = fopen("blockMeshDict","w");

fprintf(bmd,"/*--------------------------------*- C++ -*----------------------------------*\\\n");
fprintf(bmd,"| =========                 |                                                 |\n");
fprintf(bmd,"| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n");
fprintf(bmd,"|  \\    /   O peration     | Version:  1.5                                   |\n");
fprintf(bmd,"|   \\  /    A nd           | Web:      http://www.OpenFOAM.org               |\n");
fprintf(bmd,"|    \\/     M anipulation  |                                                 |\n");
fprintf(bmd,"\\*---------------------------------------------------------------------------*/\n");
fprintf(bmd,"FoamFile\n");
fprintf(bmd,"{\n");
fprintf(bmd,"    version     2.0;\n");
fprintf(bmd,"    format      ascii;\n");
fprintf(bmd,"    class       dictionary;\n");
fprintf(bmd,"    object      blockMeshDict;\n");
fprintf(bmd,"}\n");
fprintf(bmd,"// * * * * created with GraVentGen4.c * * * * * * * * * * * * * * * * * * * //\n\n");
fprintf(bmd,"convertToMeters 0.001;\n\n");


// *****************************************
fprintf(bmd,"vertices\n");
fprintf(bmd,"(\n\n");
vert=0;

fprintf(bmd,"   (%4.0f %4.0f %4.0f) // vert. %i start\n", 0.0, h/2, Z, vert++);
fprintf(bmd,"   (%4.0f %4.0f %4.0f) // vert. %i \n", 0.0, h/2,-Z, vert++);
fprintf(bmd,"   (%4.0f %4.0f %4.0f) // vert. %i \n", 0.0, -h/2, Z, vert++);
fprintf(bmd,"   (%4.0f %4.0f %4.0f) // vert. %i \n", 0.0, -h/2,-Z, vert++);
fprintf(bmd,"   (%4.0f %4.0f %4.0f) // vert. %i \n", L, -h/2, Z, vert++);
fprintf(bmd,"   (%4.0f %4.0f %4.0f) // vert. %i \n", L, -h/2,-Z, vert++);
fprintf(bmd,"   (%4.0f %4.0f %4.0f) // vert. %i start\n", L, h/2, Z, vert++);
fprintf(bmd,"   (%4.0f %4.0f %4.0f) // vert. %i \n", L, h/2,-Z, vert++);

fprintf(bmd,");\n\n\n");
// *****************************************

// *****************************************
fprintf(bmd,"blocks\n");
fprintf(bmd,"(\n\n");


	fprintf(bmd,"    hex (%4i %4i %4i %4i %4i %4i %4i %4i)  (%3i %3i %3i) simpleGrading (1 1 1) // block %i\n",3,5,7,1,2,4,6,0,lround(L/dx),lround(h/dy),zc,i,blocks++);

fprintf(bmd,");\n\n");

fprintf(bmd,"\nedges // only required for curved edges\n");
fprintf(bmd,"(\n");
fprintf(bmd,");\n\n");

fprintf(bmd,"patches\n");
fprintf(bmd,"(  // Normal vector has to point out of the domain;\n");
fprintf(bmd,"   // (right hand rule in order of enumeration)\n");

fprintf(bmd,"  wall wand\n");
fprintf(bmd,"   (\n");
  fprintf(bmd,"    (%i %i %i %i)\n", 2,0,1,3);
  fprintf(bmd,"    (%i %i %i %i)\n", 0,6,7,1);
  fprintf(bmd,"    (%i %i %i %i)\n", 5,7,6,4);
  fprintf(bmd,"    (%i %i %i %i)\n", 4,2,3,5);
fprintf(bmd,"   )\n\n");



fprintf(bmd,"    empty frontAndBack\n");
fprintf(bmd,"    (\n");
  fprintf(bmd,"      (%i %i %i %i) // front\n",2,4,6,0);
  fprintf(bmd,"      (%i %i %i %i) // back\n",3,1,7,5);
fprintf(bmd,"    )\n");

fprintf(bmd,");\n\n");

fprintf(bmd,"mergePatchPairs\n");
fprintf(bmd,"(\n");
fprintf(bmd,");\n\n");

fprintf(bmd,"// * * * * created with Vent2D * * * * * * * * * * * * * * * * * * * * * * * //\n\n");

fclose(bmd);

/*********** remove obstacles from mesh ***********************/

obs = fopen("../../setObstacles.setSet","w");

fprintf(bmd,"cellSet c0 new\n");
fprintf(bmd,"cellSet c0 invert\n\n");

for (i=0; i< OBSTACLES; i++)
{
  fprintf(bmd," cellSet c0 delete boxToCell (%2.6f %2.6f %2.6f) (%2.6f %2.6f %2.6f)\n", 
     (i*s-d+firstObstacle+bit)*f, (h*(1-BR)/2-bit)*f,-Z*f, (i*s+firstObstacle-bit)*f, (h/2+bit)*f, Z*f);
  fprintf(bmd," cellSet c0 delete boxToCell (%2.6f %2.6f %2.6f) (%2.6f %2.6f %2.6f)\n\n", 
     (i*s-d+firstObstacle+bit)*f, -(h/2+bit)*f,-Z*f, (i*s+firstObstacle-bit)*f,-(h*(1-BR)/2-bit)*f,Z*f);

}

fclose(obs);



return 0;
}
