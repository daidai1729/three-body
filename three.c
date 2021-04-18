#pragma GCC optimize("Ofast")
#include <stdio.h>
#include <math.h>
#include <mpi.h>
//#define m1 5
//#define m2 4
//#define m3 3
#define fx1(x1,y1,x2,y2,x3,y3) (-(4*(x1-x2)/pow(pow(x2-x1,2)+pow(y2-y1,2),1.5)+3*(x1-x3)/pow(pow(x3-x1,2)+pow(y3-y1,2),1.5)))
#define fx2(x1,y1,x2,y2,x3,y3) (-(3*(x2-x3)/pow(pow(x3-x2,2)+pow(y3-y2,2),1.5)+5*(x2-x1)/pow(pow(x1-x2,2)+pow(y1-y2,2),1.5)))
#define fx3(x1,y1,x2,y2,x3,y3) (-(5*(x3-x1)/pow(pow(x1-x3,2)+pow(y1-y3,2),1.5)+4*(x3-x2)/pow(pow(x2-x3,2)+pow(y2-y3,2),1.5)))
#define fy1(x1,y1,x2,y2,x3,y3) (-(4*(y1-y2)/pow(pow(x2-x1,2)+pow(y2-y1,2),1.5)+3*(y1-y3)/pow(pow(x3-x1,2)+pow(y3-y1,2),1.5)))
#define fy2(x1,y1,x2,y2,x3,y3) (-(3*(y2-y3)/pow(pow(x3-x2,2)+pow(y3-y2,2),1.5)+5*(y2-y1)/pow(pow(x1-x2,2)+pow(y1-y2,2),1.5)))
#define fy3(x1,y1,x2,y2,x3,y3) (-(5*(y3-y1)/pow(pow(x1-x3,2)+pow(y1-y3,2),1.5)+4*(y3-y2)/pow(pow(x2-x3,2)+pow(y2-y3,2),1.5)))
void gx1(double x1, double y1, double x2, double y2, double x3, double y3, double xp1, double* xn1, double* xpn1, double h) {
	double ph;
	ph = xp1 + (h * 0.5) * fx1(x1, y1, x2, y2, x3, y3);
	*xn1 = x1 + h * ph;
	*xpn1 = ph + (h * 0.5) * fx1(*xn1, y1, x2, y2, x3, y3);
}
void gy1(double x1, double y1, double x2, double y2, double x3, double y3, double xp1, double* xn1, double* xpn1, double h) {
	double ph;
	ph = xp1 + (h * 0.5) * fy1(x1, y1, x2, y2, x3, y3);
	*xn1 = y1 + h * ph;
	*xpn1 = ph + (h * 0.5) * fy1(x1, *xn1, x2, y2, x3, y3);
}
void gx2(double x1, double y1, double x2, double y2, double x3, double y3, double xp1, double* xn1, double* xpn1, double h) {
	double ph;
	ph = xp1 + (h * 0.5) * fx2(x1, y1, x2, y2, x3, y3);
	*xn1 = x2 + h * ph;
	*xpn1 = ph + (h * 0.5) * fx2(x1, y1, *xn1, y2, x3, y3);
}
void gy2(double x1, double y1, double x2, double y2, double x3, double y3, double xp1, double* xn1, double* xpn1, double h) {
	double ph;
	ph = xp1 + (h * 0.5) * fy2(x1, y1, x2, y2, x3, y3);
	*xn1 = y2 + h * ph;
	*xpn1 = ph + (h * 0.5) * fy2(x1, y1, x2, *xn1, x3, y3);
}
void gx3(double x1, double y1, double x2, double y2, double x3, double y3, double xp1, double* xn1, double* xpn1, double h) {
	double ph;
	ph = xp1 + (h * 0.5) * fx3(x1, y1, x2, y2, x3, y3);
	*xn1 = x3 + h * ph;
	*xpn1 = ph + (h * 0.5) * fx3(x1, y1, x2, y2, *xn1, y3);
}
void gy3(double x1, double y1, double x2, double y2, double x3, double y3, double xp1, double* xn1, double* xpn1, double h) {
	double ph;
	ph = xp1 + (h * 0.5) * fy3(x1, y1, x2, y2, x3, y3);
	*xn1 = y3 + h * ph;
	*xpn1 = ph + (h * 0.5) * fy3(x1, y1, x2, y2, x3, *xn1);
}
int main(int argc,char **argv){
	MPI_Init(&argc,&argv);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	double dt1, dt2;
	double dt = 0.0000001;
	double x1, y1, x2, y2, x3, y3, xp1, yp1, xp2, yp2, xp3, yp3, xn1, yn1, xn2, yn2, xn3, yn3, xpn1, ypn1, xpn2, ypn2, xpn3, ypn3;
	double cx1, cy1, cx2, cy2, cx3, cy3;
	double d, d1, d2;
	double wide = 0.00006;
	int n;
	FILE* fp;
	/* 4th order yoshida parameters */
	d = pow(2, 1.0 / 3.0);
	d1 = 1 / (2 - d);
	d2 = -d / (2 - d);
	/* initial condition */
	x1 = 0;
	y1 = 0;
	x2 = -3;
	y2 = 0;
	x3 = 0;
	y3 = 4;
	xp1 = 0;
	yp1 = 0;
	xp2 = 0;
	yp2 = 0;
	xp3 = 0;
	yp3 = 0;
	
	cx1 = x1;
	cy1 = y1;
	cx2 = x2;
	cy2 = y2;
	cx3 = x3;
	cy3 = y3;
	/* file open */
	fp = fopen("three9.txt", "w");
	fprintf(fp, "%6.2f %6.2f %6.2f %6.2f %6.2f %6.2f\n", cx1, cy1, cx2, cy2, cx3, cy3);
	/* time development */
	for (n = 1; n < 5000000000; n++) {
		dt1 = dt * d1;
		dt2 = dt * d2;
		
		if(rank==0)
		{
			gx1(x1, cy1, cx2, cy2, cx3, cy3, xp1, &xn1, &xpn1, dt1);
			x1 = xn1;
			xp1 = xpn1;
			MPI_Bcast(&x1,1,MPI_DOUBLE,rank,MPI_COMM_WORLD);
			MPI_Bcast(&xp1,1,MPI_DOUBLE,rank,MPI_COMM_WORLD);
			gx1(x1, cy1, cx2, cy2, cx3, cy3, xp1, &xn1, &xpn1, dt2);
			x1 = xn1;
			xp1 = xpn1;
			MPI_Bcast(&x1,1,MPI_DOUBLE,rank,MPI_COMM_WORLD);
			MPI_Bcast(&xp1,1,MPI_DOUBLE,rank,MPI_COMM_WORLD);
			gx1(x1, cy1, cx2, cy2, cx3, cy3, xp1, &xn1, &xpn1, dt1);
			x1 = xn1;
			xp1 = xpn1;
			MPI_Bcast(&x1,1,MPI_DOUBLE,rank,MPI_COMM_WORLD);
			MPI_Bcast(&xp1,1,MPI_DOUBLE,rank,MPI_COMM_WORLD);
			/* write data */
			if (!(n % 200000)) {
				printf("%d\n",n);
				fprintf(fp, "%6.2f %6.2f %6.2f %6.2f %6.2f %6.2f\n", cx1, cy1, cx2, cy2, cx3, cy3);
			}
		}
		if(rank==1)
		{
			gx2(cx1, cy1, x2, cy2, cx3, cy3, xp2, &xn2, &xpn2, dt1);
			x2 = xn2;
			xp2 = xpn2;
			MPI_Bcast(&x2,1,MPI_DOUBLE,rank,MPI_COMM_WORLD);
			MPI_Bcast(&xp2,1,MPI_DOUBLE,rank,MPI_COMM_WORLD);
			gx2(cx1, cy1, x2, cy2, cx3, cy3, xp2, &xn2, &xpn2, dt2);
			x2 = xn2;
			xp2 = xpn2;
			MPI_Bcast(&x2,1,MPI_DOUBLE,rank,MPI_COMM_WORLD);
			MPI_Bcast(&xp2,1,MPI_DOUBLE,rank,MPI_COMM_WORLD);
			gx2(cx1, cy1, x2, cy2, cx3, cy3, xp2, &xn2, &xpn2, dt1);
			x2 = xn2;
			xp2 = xpn2;
			MPI_Bcast(&x2,1,MPI_DOUBLE,rank,MPI_COMM_WORLD);
			MPI_Bcast(&xp2,1,MPI_DOUBLE,rank,MPI_COMM_WORLD);
			if (!(pow(cx1 - cx2, 2) <= wide && pow(cy1 - cy2, 2) <= wide) || (pow(cx2 - cx3, 2) <= wide && pow(cy2 - cy3, 2) <= wide) || (pow(cx3 - cx1, 2) <= wide && pow(cy3 - cy1, 2) <= wide)) {
				dt = 0.0000001;
			}
			else
			{
				dt = 0.000000000001;
			}
		}
		if(rank==2)
		{
			gx3(cx1, cy1, cx2, cy2, x3, cy3, xp3, &xn3, &xpn3, dt1);
			x3 = xn3;
			xp3 = xpn3;
			MPI_Bcast(&x3,1,MPI_DOUBLE,rank,MPI_COMM_WORLD);
			MPI_Bcast(&xp3,1,MPI_DOUBLE,rank,MPI_COMM_WORLD);
			gx3(cx1, cy1, cx2, cy2, x3, cy3, xp3, &xn3, &xpn3, dt2);
			x3 = xn3;
			xp3 = xpn3;
			MPI_Bcast(&x3,1,MPI_DOUBLE,rank,MPI_COMM_WORLD);
			MPI_Bcast(&xp3,1,MPI_DOUBLE,rank,MPI_COMM_WORLD);
			gx3(cx1, cy1, cx2, cy2, x3, cy3, xp3, &xn3, &xpn3, dt1);
			x3 = xn3;
			xp3 = xpn3;
			MPI_Bcast(&x3,1,MPI_DOUBLE,rank,MPI_COMM_WORLD);
			MPI_Bcast(&xp3,1,MPI_DOUBLE,rank,MPI_COMM_WORLD);
		}
		if(rank==3)
		{
			gy1(cx1, y1, cx2, cy2, cx3, cy3, yp1, &yn1, &ypn1, dt1);
			y1 = yn1;
			yp1 = ypn1;
			MPI_Bcast(&y1,1,MPI_DOUBLE,rank,MPI_COMM_WORLD);
			MPI_Bcast(&yp1,1,MPI_DOUBLE,rank,MPI_COMM_WORLD);
			gy1(cx1, y1, cx2, cy2, cx3, cy3, yp1, &yn1, &ypn1, dt2);
			y1 = yn1;
			yp1 = ypn1;
			MPI_Bcast(&y1,1,MPI_DOUBLE,rank,MPI_COMM_WORLD);
			MPI_Bcast(&yp1,1,MPI_DOUBLE,rank,MPI_COMM_WORLD);
			gy1(cx1, y1, cx2, cy2, cx3, cy3, yp1, &yn1, &ypn1, dt1);
			y1 = yn1;
			yp1 = ypn1;
			MPI_Bcast(&y1,1,MPI_DOUBLE,rank,MPI_COMM_WORLD);
			MPI_Bcast(&yp1,1,MPI_DOUBLE,rank,MPI_COMM_WORLD);
		}
		if(rank==4)
		{
			gy2(cx1, cy1, cx2, y2, cx3, cy3, yp2, &yn2, &ypn2, dt1); 
			y2 = yn2;
			yp2 = ypn2;
			MPI_Bcast(&y2,1,MPI_DOUBLE,rank,MPI_COMM_WORLD);
			MPI_Bcast(&yp2,1,MPI_DOUBLE,rank,MPI_COMM_WORLD);
			gy2(cx1, cy1, cx2, y2, cx3, cy3, yp2, &yn2, &ypn2, dt2);
			y2 = yn2;
			yp2 = ypn2;
			MPI_Bcast(&y2,1,MPI_DOUBLE,rank,MPI_COMM_WORLD);
			MPI_Bcast(&yp2,1,MPI_DOUBLE,rank,MPI_COMM_WORLD);
			gy2(cx1, cy1, cx2, y2, cx3, cy3, yp2, &yn2, &ypn2, dt1); 
			y2 = yn2;
			yp2 = ypn2;
			MPI_Bcast(&y2,1,MPI_DOUBLE,rank,MPI_COMM_WORLD);
			MPI_Bcast(&yp2,1,MPI_DOUBLE,rank,MPI_COMM_WORLD);
		}
		if(rank==5)
		{
			gy3(cx1, cy1, cx2, cy2, cx3, y3, yp3, &yn3, &ypn3, dt1);
			y3 = yn3;
			yp3 = ypn3;
			MPI_Bcast(&y3,1,MPI_DOUBLE,rank,MPI_COMM_WORLD);
			MPI_Bcast(&yp3,1,MPI_DOUBLE,rank,MPI_COMM_WORLD);
			gy3(cx1, cy1, cx2, cy2, cx3, y3, yp3, &yn3, &ypn3, dt2);
			y3 = yn3;
			yp3 = ypn3;
			MPI_Bcast(&y3,1,MPI_DOUBLE,rank,MPI_COMM_WORLD);
			MPI_Bcast(&yp3,1,MPI_DOUBLE,rank,MPI_COMM_WORLD);
			gy3(cx1, cy1, cx2, cy2, cx3, y3, yp3, &yn3, &ypn3, dt1);
			y3 = yn3;
			yp3 = ypn3;
			MPI_Bcast(&y3,1,MPI_DOUBLE,rank,MPI_COMM_WORLD);
			MPI_Bcast(&yp3,1,MPI_DOUBLE,rank,MPI_COMM_WORLD);
		}
		
		cx1 = x1;
		cy1 = y1;
		cx2 = x2;
		cy2 = y2;
		cx3 = x3;
		cy3 = y3;

	}/* end of time development */
	MPI_Finalize();
	return 0;
}
