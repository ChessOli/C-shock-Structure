#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*
    This codes integrates the static equations for a C shock with non-ideal MHD following Appendix A in 
    Zier, Mayer, Springel 2023. The equations are integrated backwards starting at x0 and going to  x = 0.
*/

double eta_ohm_0 = 0.0002;
double x0 = 6;
double eta_hall_0 = -0.05;
double eta_Ambipolar_0 = 0;


// functions go calculate non-ideal MHD diffusivities
double get_eta_ohm(double Bx, double By, double Bz) {
	return eta_ohm_0;
}

double get_eta_hall(double Bx, double By, double Bz) {
	double B = sqrt(Bx * Bx + By * By + Bz * Bz);
	return eta_hall_0 * B;
}

double get_eta_ambipolar(double Bx, double By, double Bz,double rho) {
	double B2 = (Bx * Bx + By * By + Bz * Bz);
	return eta_Ambipolar_0 * B2 / rho;
}


// The post shock state, can be calculated from pre-shock values using the jump conditions for an ideal MHD shock
double rho0 = 1.7942;
double velx0 = 0.9759;
double vely0 = 0.6561;
double velz0 = 0;

double Bx0 = -1;
double By0 = -1.74885;
double Bz0 = 0;

//isothermal sound speed
double cs = 0.1;

//The initial gradients of the magnetic field, are required to trigger the shock
double dBy0 = -0.0001;
double dBz0 = -0.0001;

// the number of integration steps in x-direction
#define N_STEPS 400000


//Arrays to save the values at a specific position
double* rho;
double* velx;
double* vely;
double* velz;

double* By;
double* Bz;

double* dBy;
double* dBz;


double Q, K;

// the intial values of vector M due to the inital gradients of the magnetic field
double M1_0 = 0;
double M2_0 = 0;


// the original values of the conserved quantities, used for code testing
double conserved_values_original[6];


void get_M(double Bx, double By, double Bz, double vx, double vy, double vz, double*M) {
	M[0] = vx * By - vy * Bx  + vely0 * Bx0 - velx0 * By0 + M1_0;
	M[1] = vx * Bz - vz * Bx + velz0 * Bx0 - velx0 * Bz0 + M2_0;
	return;
}

void get_R(double By, double Bz, double rho, double* R) {
	double B2 = Bx0 * Bx0 + By * By + Bz * Bz;
	double B =sqrt(B2);
	double eta_ohm = get_eta_ohm(Bx0, By, Bz);
	double eta_hall = get_eta_hall(Bx0, By, Bz);
	double eta_Ambipolar = get_eta_ambipolar(Bx0, By, Bz, rho);


	R[0] = eta_ohm + eta_Ambipolar *(1 - Bz * Bz / B2);
	R[1] = eta_hall * Bx0 / B + eta_Ambipolar * By * Bz /B2;
	R[2] = -eta_hall * Bx0 / B + eta_Ambipolar * By * Bz /B2;
	R[3] = eta_ohm + eta_Ambipolar * (1 - By * By / B2);
}

void set_inital_M() {
	double R[4];
	get_R(By0, Bz0,rho0 ,R);
	M1_0 =  R[0] * dBy0 + R[1] * dBz0;
	M2_0 =  R[2] * dBy0 + R[3] * dBz0;
}


void get_conserved_quantities(double rho, double vx, double vy, double vz, double Bx, double By, double Bz, double dBy, double dBz, double*conserved_values) {
	conserved_values[0] = rho * vx;
	conserved_values[1] = rho * vx * vx + rho * cs * cs+
	(Bx * Bx + By * By + Bz * Bz) / 2.;
	conserved_values[2] = rho * vx * vy - Bx * By;
	conserved_values[3] = rho * vx * vz - Bx * Bz;

	double R[4];
	get_R(By, Bz,rho ,R);

	conserved_values[4] = vx * By - vy * Bx - R[0] * dBy - R[1] * dBz;
	conserved_values[5] = vx * Bz - vz * Bx - R[2] * dBy - R[3] * dBz;
	return;
}

void setup_verification() {
	get_conserved_quantities(rho0, velx0, vely0, velz0, Bx0, By0, Bz0, dBy0, dBz0, conserved_values_original);
}

// verify if conservative quantities are conserved
void verify_data(int i) {
	double conserved_values_now[6];
	get_conserved_quantities(rho[i], velx[i], vely[i], velz[i], Bx0, By[i], Bz[i], dBy[i], dBz[i], conserved_values_now);
	for(int j = 0; j < 6; j++) {
		if(fabs(conserved_values_now[j] - conserved_values_original[j])/ conserved_values_original[j] > 1e-2) printf("Problem, i: %i, cons: %i %g %g\n", i,j,conserved_values_now[j],conserved_values_original[j]);
	}


}

double get_vy(double By) {
	return vely0 + Bx0 * (By - By0) / Q;
}

double get_vz(double Bz) {
	return velz0 + Bx0 * (Bz - Bz0) / Q;
}

double get_rho(double vx) {
	return Q / vx;
}

double get_vx(double By, double Bz) {
	double B2 = By * By + Bz * Bz;
	double term1 = K - B2 / 2.;
	return 1. / (2* Q) * (term1 + sqrt(term1 * term1 - 4 * Q * Q * cs *cs));
}

void update_terms(int i) {
 	vely[i+1] = get_vy(By[i+1]);
 	velz[i+1] = get_vz(Bz[i+1]);

 	velx[i+1] = get_vx(By[i+1],Bz[i+1]);

 	rho[i+1] = get_rho(velx[i+1]);
}

void get_derivatives(double* M, double* R, double* dB) {
	double term1 = R[0] * R[3] - R[1] * R[2];
	dB[0] = (M[0] * R[3] - M[1]*R[1]) /term1;
	dB[1] = (M[1] * R[0] - M[0]*R[2]) /term1;
}

int main() {
	Q = rho0 *velx0;
	K = Q *velx0 + Q / velx0 * cs * cs + By0 * By0 / 2 + Bz0 * Bz0/2.;

double x1 = 0;

double dx = (x1 - x0) / N_STEPS;

rho = (double*) malloc((N_STEPS +1) * sizeof(double));
velx = (double*) malloc((N_STEPS +1) * sizeof(double));
vely = (double*) malloc((N_STEPS +1) * sizeof(double));
velz = (double*) malloc((N_STEPS +1) * sizeof(double));

By = (double*) malloc((N_STEPS +1) * sizeof(double));
Bz = (double*) malloc((N_STEPS +1) * sizeof(double));

dBy = (double*) malloc((N_STEPS +1) * sizeof(double));
dBz = (double*) malloc((N_STEPS +1) * sizeof(double));

rho[0] = rho0;
velx[0] = velx0;
vely[0] = vely0;
velz[0] = velz0;

By[0] = By0;
Bz[0] = Bz0;

setup_verification();
dBy[0] =  dBy0;
dBz[0] =  dBz0;


FILE *fptr;
fptr = fopen("HallData.txt","w");

set_inital_M();

for(int i = 0; i < N_STEPS; i++) {
	double M[2];
	double R[4];


	get_M(Bx0, By[i], Bz[i], velx[i], vely[i], velz[i], M);
	get_R(By[i], Bz[i], rho[i], R);

	double dB1[2];
	get_derivatives(M, R, dB1);

 	By[i+1] = By[i] + dB1[0] * dx * 0.5;
 	Bz[i+1] = Bz[i] + dB1[1] * dx * 0.5;
 	update_terms(i);




 	get_M(Bx0, By[i+1], Bz[i+1], velx[i+1], vely[i+1], velz[i+1], M);
	get_R(By[i+1], Bz[i+1], rho[i+1], R);

	double dB2[2];
	get_derivatives(M, R, dB2);

 	By[i+1] = By[i] + dB2[0] / 2. * dx;
 	Bz[i+1] = Bz[i] + dB2[1] / 2. * dx;
 	update_terms(i);


 	get_M(Bx0, By[i+1], Bz[i+1], velx[i+1], vely[i+1], velz[i+1], M);
	get_R(By[i+1], Bz[i+1], rho[i+1], R);

	double dB3[2];
	get_derivatives(M, R, dB3);

 	By[i+1] = By[i] + dB3[0]  * dx;
 	Bz[i+1] = Bz[i] + dB3[1]  * dx;
 	update_terms(i);


 	get_M(Bx0, By[i+1], Bz[i+1], velx[i+1], vely[i+1], velz[i+1], M);
	get_R(By[i+1], Bz[i+1], rho[i+1], R);

	double dB4[2];
	get_derivatives(M, R, dB4);

 	By[i+1] = By[i] + (dB1[0] + 2* dB2[0] + 2* dB3[0] + dB4[0])/6.  * dx;
 	Bz[i+1] = Bz[i] + (dB1[1] + 2* dB2[1] + 2* dB3[1] + dB4[1])/6.  * dx;
 	update_terms(i);



 	get_M(Bx0, By[i+1], Bz[i+1], velx[i+1], vely[i+1], velz[i+1], M);
	get_R(By[i+1], Bz[i+1], rho[i+1], R);

	double dBfinal[2];
	get_derivatives(M, R, dBfinal);

 	dBy[i+1] =  dBfinal[0];
 	dBz[i+1] =  dBfinal[1];

 	if(velx[i+1] != velx[i+1]) break;
 	verify_data(i+1);
 	double x = i *dx + x0;
	//fprintf(fptr,"%g %g %g %g %g %g %g\n",x, rho[i], velx[i], vely[i], velz[i], By[i], Bz[i]);
}

for(int i = 0; i < 1000; i++){
	double x =  -1 + 1./1000. *i; 
	fprintf(fptr,"%g %g %g %g %g %g %g\n",x, rho[N_STEPS], velx[N_STEPS], vely[N_STEPS], velz[N_STEPS], By[N_STEPS], Bz[N_STEPS]);
}

for(int i = N_STEPS; i >= 0; i--) {
	double x = i *dx + x0;
	fprintf(fptr,"%g %g %g %g %g %g %g\n",x, rho[i], velx[i], vely[i], velz[i], By[i], Bz[i]);
}

for(int i = 1; i < 1000; i++){
	double x =  x0 + 1./1000. *i; 
	fprintf(fptr,"%g %g %g %g %g %g %g\n",x, rho[0], velx[0], vely[0], velz[0], By[0], Bz[0]);
}

fclose(fptr);

free(rho);
free(velx);
free(vely);
free(velz);

free(By);
free(Bz);
}