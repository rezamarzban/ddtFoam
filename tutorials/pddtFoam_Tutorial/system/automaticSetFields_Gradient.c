#include "stdio.h"
#include "math.h"
#include "stdlib.h"

#define MO2 32.0
#define MH2 2.016
#define MN2 28.01
#define MH2O 18.016

#define T0 293.0
#define fStoich 0.028511

#define ymin -0.03
#define ymax  0.03

int main(int argc, char *argv[])
{
int i;
double h;
double y,z;
double _10H2, _20H2, _30H2, _40H2;
double Nenner, xO2,yO2, xH2,yH2, xN2,yN2, yH2mid;
double H2Ign,O2Ign,N2Ign,H2OIgn,TIgn,Tad;
double xH2max=0.0;
double xH2min=1.0;
double H2molf, cIgn, ignRadius;
int steps;

if(argc!=5)
{
printf("Error! Provide 4 arguments for automaticSetFields: H2 molar fraction in percent (0...40), cIgn (0...1), radius of ignition spot (in millimetres!) and number of steps for the gradient.\n");
return -1;
}

H2molf = ( atof(argv[1]) )/ 100;
cIgn =   atof(argv[2]);
ignRadius = atof(argv[3])*1e-3;
steps = atoi(argv[4]);

h=1.0*(ymax-ymin)/steps;

printf("/*--------------------------------*- C++ -*----------------------------------*\\\n");
printf("| =========                 |                                                 |\n");
printf("| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n");
printf("|  \\    /   O peration     | Version:  1.7.1                                 |\n");
printf("|   \\  /    A nd           | Web:      http://www.OpenFOAM.org               |\n");
printf("|    \\/     M anipulation  |                                                 |\n");
printf("\\*---------------------------------------------------------------------------*/\n");
printf("FoamFile\n");
printf("{\n");
printf("    version     2.0;\n");
printf("    format      ascii;\n");
printf("    class       dictionary;\n");
printf("    object      setFieldsDict;\n");
printf("}\n");
printf("// * * * * created with automaticSetFields_tfc_zu_Grad3s_03.c * * * * * * * *  //\n\n");


	xH2=H2molf; 
        xO2=0.21*(1.0-xH2);
	xN2=1.0-xO2-xH2;

	Nenner=xO2*MO2+xH2*MH2+xN2*MN2;
	yO2=xO2*MO2/Nenner;
	yH2=xH2*MH2/Nenner;
	yH2mid = yH2;
	yN2=1.0-yO2-yH2;	

	if((yN2<0.0)||(yH2<0.0)||(yO2<0.0)||(yN2>1.0)||(yH2>1.0)||(yO2>1.0)) 
		{printf("\n\nError: mass fractions don't sum to unity\n"); return 1;}







	    if(yH2<fStoich)
	    {
		H2Ign = yH2*(1.0-cIgn);
		/*
		O2Ign = yO2-MO2/MH2*0.5*yH2*(1-cIgn);
		H2OIgn= yH2-H2Ign + yO2-O2Ign;
		N2Ign = 1.0-H2Ign-O2Ign-H2OIgn;
		*/
		N2Ign = yN2;
		H2OIgn = (yH2-H2Ign)*MH2O/MH2;
		O2Ign = 1.0-H2Ign-N2Ign-H2OIgn;
	    }
	    else
	    {
		O2Ign = yO2*(1.0-cIgn);
		N2Ign = yN2;
		H2OIgn = (yO2-O2Ign)*MH2O/(0.5*MO2);
		H2Ign = 1.0-O2Ign-N2Ign-H2OIgn;
	    }

 	    if(yH2<0.03) // polynoms for adiabatic flame temperature from Cantera
	    {
		Tad = -1.47201e+06*pow(yH2,2) + 1.14928e+05*yH2 + 2.91241e+02 + (T0-293.0);
	    }
	    else
	    {
	    	Tad =  5.40634e+04*pow(yH2,2) - 1.88089e+04*yH2 + 2.93276e+03 + (T0-293.0);
	    }
	    TIgn = T0 + cIgn* (Tad-T0);
	


	printf("\n");




	printf("	defaultFieldValues\n");
	printf("	(\n");
	printf("            volScalarFieldValue c   %.2f\n",0.0);
	printf("            volScalarFieldValue fH %.5f\n",yH2);
	printf("            volScalarFieldValue H2  %.5f\n",yH2);
	printf("            volScalarFieldValue O2  %.5f\n",yO2);
	printf("            volScalarFieldValue H2O 0.0\n");
	printf("            volScalarFieldValue N2  %.5f\n",yN2);
	printf("            volScalarFieldValue T   %.1f\n", T0);
	printf("            volScalarFieldValue Tu   %.1f\n", T0);
	printf("        );\n");



printf("regions\n");
printf("(\n");

for (i=0; i<steps; i++)
{
	y=ymin+0.5*h+h*i;
	z=y-ymax;

// Polynome neu ausgewertet, sind flacher:
// aus der Semesterarbeit Stadlmair liegen CFX-Ergebnisse vor, z.B. /media/Elements1000/SA_Stadlmair/ANSYS-Projekte/Reaktorblock_Fertig_Fein/30Prozent_H2/Tecplot/
// diese wurden interpoliert in /nfs/home/ettner/Matlab/Gradient_Nico/Gradient_interpolate_geschlossen_02.m
// Ergebnis: /nfs/home/ettner/Matlab/Gradient_Nico/geschlossen/3s/polynoms_geschlossen_3s_order5.txt:

	_10H2 = -3.50957e+04*pow(z,5) - 3.26065e+04*pow(z,4) -5.25765e+03*pow(z,3) -2.62149e+02*pow(z,2) -3.84660e-01*z + 1.99094e-01;
	_20H2 =  1.84717e+05*pow(z,5) - 1.40719e+04*pow(z,4) -6.78365e+03*pow(z,3) -4.17966e+02*pow(z,2) -6.09813e-01*z + 3.77362e-01;
	_30H2 =  3.70457e+05*pow(z,5) + 1.89256e+04*pow(z,4) -5.74215e+03*pow(z,3) -4.60436e+02*pow(z,2) -4.98623e-01*z + 5.26445e-01;
	_40H2 =  5.69961e+05*pow(z,5) + 6.66643e+04*pow(z,4) -2.62118e+03*pow(z,3) -4.11878e+02*pow(z,2) -1.79756e-01*z + 6.48390e-01;
	

	if(H2molf < 0.10) xH2 = H2molf/0.10 * _10H2;
	else if (H2molf < 0.20) xH2 = _10H2 + (H2molf-0.10)/0.10 * (_20H2-_10H2);
	else if (H2molf < 0.30) xH2 = _20H2 + (H2molf-0.20)/0.10 * (_30H2-_20H2);
	else if (H2molf <= 0.40) xH2 = _30H2 + (H2molf-0.30)/0.10 * (_40H2-_30H2);
	else {printf("\n Error: molar H2 fraction above 0.40 \n"); return 1;}

        xO2=0.21*(1.0-xH2);
	xN2=1.0-xO2-xH2;

	xH2min=fmin(xH2,xH2min);
	xH2max=fmax(xH2,xH2max);


//	printf("%.4f, %.4f, %.4f\n",xO2,xH2,xN2);


	Nenner=xO2*MO2+xH2*MH2+xN2*MN2;
	yO2=xO2*MO2/Nenner;
	yH2=xH2*MH2/Nenner;
	yN2=1.0-yO2-yH2;	

	if((yN2<0.0)||(yH2<0.0)||(yO2<0.0)||(yN2>1.0)||(yH2>1.0)||(yO2>1.0)) 
		{printf("\n\nError: mass fractions don't sum to unity\n"); return 1;}


	printf("    boxToCell\n");
	printf("    {\n");
	printf("        box (-10 %.4f -10) (10 10 10);\n\n",y-0.5*h-1e-5);
	printf("        fieldValues\n");
	printf("        (\n");
	printf("	    // xH2 = %.3f\n",xH2);		
	printf("	    volScalarFieldValue fH  %.5f\n",yH2);		
	printf("	    volScalarFieldValue H2  %.5f\n",yH2);	
	printf("	    volScalarFieldValue O2  %.5f\n",yO2);
	printf("	    volScalarFieldValue N2  %.5f\n\n",yN2);
	printf("        );\n");
	printf("    }\n");

}
	printf("\n");


  	printf("    cylinderToCell\n");
	printf("    {\n");
	printf("        p1       (0.0 0.0 -1); // start point on cylinder axis\n");
	printf("        p2       (0.0 0.0  1);  // end point on cylinder axis\n");
	printf("        radius   %.5f;\n",ignRadius);

	printf("        fieldValues\n");
	printf("        (\n\n");

	printf("            volScalarFieldValue c   %.2f\n",cIgn);
	printf("            volScalarFieldValue fH  %.5f\n",yH2mid);
	printf("            volScalarFieldValue H2  %.5f\n",H2Ign);
	printf("            volScalarFieldValue O2  %.5f\n",O2Ign);
	printf("            volScalarFieldValue H2O %.5f\n",H2OIgn);
	printf("            volScalarFieldValue N2  %.5f\n",yN2);
	printf("            volScalarFieldValue T   %.1f\n",TIgn+cIgn*60.0); 

	printf("        );\n");


	printf("    }\n");




printf(");\n\n");


printf("// ************************************************************************* //\n");
printf("// inhomogeneous mixture for closed vessel with %.1f \% H2 generated         //\n",100*H2molf);
printf("// ************************************************************************* //\n");

return 0;
}
