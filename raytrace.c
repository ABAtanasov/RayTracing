#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#define w 1000
#define h 1000
#define shadow 0.2
#define stotal 3
#define N 4
#define reflection .2
#define maxrecursion 3
typedef struct{
	double cx,cy,cz,r;
	int R,G,B;
} sphere;
typedef struct{
	double x,y,z;
	double dx,dy,dz;
} ray;
typedef struct{
	double x,y,z;
	double nx,ny,nz;
	int R,G,B;
} plane;
sphere s[stotal];
plane fl;
void init(){
	s[0]=(sphere){0.5,0.5,0.166667,0.166667,0,0,255};
	s[1]=(sphere){0.833333,0.5,0.6,0.166667,0,255,0};
	s[2]=(sphere){0.333333,0.666667,0.666667,0.333333,255,0,0};
	fl=(plane){0,0.333333,0,0,1,0,255,255,0};

}
void solvequadratic(int i, double ex, double ey, double ez,
		double dx, double dy, double dz, double* tmin, int* R,
		int* G, int* B, int* spherenum)
{
	double cx=s[i].cx;
	double cy=s[i].cy;
	double cz=s[i].cz;
	double r=s[i].r;
	double a=1.0;
	double b=2*dx*(ex-cx)+2*dy*(ey-cy)+2*dz*(ez-cz);
	double c=(ex-cx)*(ex-cx)+(ey-cy)*(ey-cy)+(ez-cz)*(ez-cz)-r*r;
	double det = b*b-4*a*c;
	if(det>0){
		double sol1=(-b-sqrt(det))/(2*a);
		double sol2=(-b+sqrt(det))/(2*a);
		if(sol1>0 && sol1<*tmin){
			*R=s[i].R; *G=s[i].G; *B=s[i].B;
			*tmin=sol1; *spherenum=i;}
		else if(sol1<0 && sol2>0 && sol2<*tmin){
			*R=s[i].R; *G=s[i].G; *B=s[i].B;
			*tmin=sol2; *spherenum=i;}
	}

}
int lightobstruction(double x0, double y0, double z0,
		double dx, double dy, double dz,double lightdistance)
{
	int i;
	for(i=0;i<stotal;i++){
	double cx=s[i].cx;
	double cy=s[i].cy;
	double cz=s[i].cz;
	double r=s[i].r;
	double a=1.0;
	double b=2*dx*(x0-cx)+2*dy*(y0-cy)+2*dz*(z0-cz);
	double c=(x0-cx)*(x0-cx)+(y0-cy)*(y0-cy)+(z0-cz)*(z0-cz)-r*r;
	double det = b*b-4*a*c;
	if(det>0){
		double sol1=(-b+sqrt(det))/(2*a);
		double sol2=(-b-sqrt(det))/(2*a);
		if(sol1>0.000001 && sol1<lightdistance) return 1;
		else if(sol2>0.000001 && sol2<lightdistance) return 1;
	}}
	return 0;
}
int mod(double a, double b,double c)
{
	while(a>c) a-=c;
	while(a<0) a+=c;
	while(b>c) b-=c;
	while(b<0) b+=c;
	if((b<0.5*c && a<0.5*c) || (b>0.5*c && a>0.5*c))
		return 1;
	else return 0;
}
void see(double ex, double ey, double ez,
		double dx, double dy, double dz,
		double lx, double ly, double lz, int* R, int* G, int* B, int level)
{
	int i;
	int rtemp,gtemp,btemp;
	double tmin,norm;
	double x0,y0,z0;
	double dot;
	ray radial,light,eye;
	int spherenum;
	int covered;
	tmin=1000; //can be higher (seeing distance)
	if(level>=maxrecursion) return;

	eye=(ray){ex,ey,ez,dx,dy,dz};

	for(i=0;i<stotal;i++){ //solving for the closest sphere
		solvequadratic(i,ex,ey,ez,dx,dy,dz,&tmin,R,G,B,&spherenum);
	}

	if(tmin==1000){ //if there is no closest sphere
		if((fl.y-ey)/dy>0){tmin=(fl.y-ey)/dy; //but the floor is seen
			x0=ex+tmin*dx; y0=ey+tmin*dy; z0=ez+tmin*dz; //establish the floor point
				if(mod(x0+.5,z0+.5,0.25)==0){*R=51; *G=51; *B=51;} //checkerboard it
				else{*R=255; *G=255; *B=255;}
				dx=fl.nx; dy=fl.ny; dz=fl.nz; //and find the normal
				norm=sqrt(dx*dx+dy*dy+dz*dz); //then normalize it
				dx/=norm; dy/=norm; dz/=norm;
				radial=(ray){fl.x,fl.y,fl.z,dx,dy,dz}; //and make it a ray
				x0+=0.000001*radial.dx; y0+=0.000001*radial.dy; z0+=0.000001*radial.dz;
				dx=lx-x0; dy=ly-y0; dz=lz-z0; //also find the ray of light to that point
				norm = sqrt(dx*dx+dy*dy+dz*dz); //and normalize it
				dx/=norm;dy/=norm;dz/=norm;
				light=(ray){x0,y0,z0,dx,dy,dz}; //and make it a ray
			}
			else{*R=0; *G=0; *B=0; //or print darkness
			return;} //and continue to the next pixel
		} //none of this will change when we add reflection

		else{i=spherenum; //if we DID hit a sphere, find its number
			x0=ex+tmin*dx; y0=ey+tmin*dy; z0=ez+tmin*dz; //and the point on it
			dx=x0-s[i].cx; dy=y0-s[i].cy; dz=z0-s[i].cz; //find the radial vector
			norm=sqrt(dx*dx+dy*dy+dz*dz); //normalize it
			dx/=norm; dy/=norm; dz/=norm;
			radial=(ray){s[i].cx,s[i].cy,s[i].cz,dx,dy,dz}; //make it a ray
			x0+=0.000001*radial.dx; y0+=0.000001*radial.dy; z0+=0.000001*radial.dz;
			dx=lx-x0; dy=ly-y0; dz=lz-z0; //find the light vector
			norm = sqrt(dx*dx+dy*dy+dz*dz); //normalize
			dx/=norm;dy/=norm;dz/=norm;
			light=(ray){x0,y0,z0,dx,dy,dz}; //make it a ray
		} //none of this will change when we add reflection
		covered=lightobstruction(x0,y0,z0,dx,dy,dz,norm); //is a sphere blocking the light?
		(*R)*=shadow; (*G)*=shadow; (*B)*=shadow; //darken everything
		if(covered==0){dot=light.dx*radial.dx+light.dy*radial.dy+light.dz*radial.dz; //dot product (<=1)
			if(dot<0) dot=0; //it is blocking itself
			*R=*R+(1-shadow)*dot*((*R)/shadow);
			*G=*G+(1-shadow)*dot*((*G)/shadow);
			*B=*B+(1-shadow)*dot*((*B)/shadow);
	}
	rtemp=*R;gtemp=*G;btemp=*B;
	dot=radial.dx*eye.dx+radial.dy*eye.dy+radial.dz*eye.dz;
	see(x0,y0,z0,eye.dx-2*dot*radial.dx,eye.dy-2*dot*radial.dy,eye.dz-2*dot*radial.dz,
			lx,ly,lz,R,G,B,level+1);
	*R=(int)((1-reflection)*rtemp+reflection*(*R)); *G=(int)((1-reflection)*gtemp+reflection*(*G)); *B=(int)((1-reflection)*btemp+reflection*(*B));
}
int main()
{
	double zscr=0;
	double ex=0.5; double ey=0.5; double ez=-1;
	double lx=0; double ly=1; double lz=-0.5;
	int i,/*j,*/frame;
	int row,col;
	//double dx,dy,dz;
	//double px,py;
	double norm;
	//int R,G,B;
	char* command=malloc(100);
	int* colors=malloc(h*w*3*sizeof(int));
	omp_set_num_threads(N);
	for(frame=0; frame<1; frame++){
	FILE* output =fopen("raycast.ppm","w");
	fprintf(output, "P3\n%d %d\n255\n",h,w);
	init();
	//#pragma omp parallel for private(zscr,ex,ey,ez,lx,ly,lz,i,j,frame,row,col,dx,dy,dz,px,py,norm,R,G,B,command)
	for(row=0; row<3*h;row+=3){ //looping over each pixel
		#pragma omp parallel for private(i)//,j,dx,dy,dz,px,py,norm,R,G,B)
		for(col=0; col<3*w; col+=3)
		{
			int colorblock[9][3];
			int j;
			double dx,dy,dz;
			double px,py;
	//		double norm;
			int R,G,B;
			for(i=0;i<3;i++)
				for(j=0;j<3;j++){
					px=(double)(col+j)/(3*w);
					py=1-(double)(row+i)/(3*h);
					dx=(px-ex);dy=(py-ey);dz=(zscr-ez); //direction of the eye ray
					norm=sqrt(dx*dx+dy*dy+dz*dz); //normalizing it
					dx/=norm; dy/=norm; dz/=norm;
					R=0;G=0;B=0;
					see(ex,ey,ez,dx,dy,dz,lx,ly,lz,&R,&G,&B,0);
					colorblock[3*i+j][0]=R;colorblock[3*i+j][1]=G;colorblock[3*i+j][2]=B;
				}
			R=0; G=0; B=0;
			for(i=0; i<9;i++){
				R+=colorblock[i][0];G+=colorblock[i][1];B+=colorblock[i][2];}
			R/=9; G/=9; B/=9;
			//fprintf(output,"%d %d %d\n",R,G,B);
			colors[row*h+col]=R;
			colors[row*h+col+1]=G;
			colors[row*h+col+2]=B;
		}}
		printf("here\n");
	for(i=0;i<w*h;i++)
	  fprintf(output,"%d %d %d\n",colors[3*i],colors[3*i+1],colors[3*i+2]);
	fclose(output);
	ez=-1+0.05*frame;
	zscr+=0.05*frame;
	sprintf(command,"convert raycast.ppm frame/frame%03d.png",frame);
	system(command);
	}
	return 0;
}
