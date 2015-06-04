//Mip-Map section by Momo

int coo (int i,int j,int w, int pd){
	return (i+j*w)*pd;
}


//on périodise par défaut
int coord_mipmap (int u,int v,int w,int pd,int d,int l){
	int x,y;
	int z = w/pow(2,d-1);
	x=good_modulus(u,z);
	y=good_modulus(v,z);
	return 2*w*w*(1-1/pow(2,d-1))*pd+(x+y*w)*pd+l;
}

//pour plus d'option, on déclare ce filtre qui est la moyenne
//float filter_mipmap(int i,int j,int w,int pd,int d,int l,float *r){
//	return (r[coord_mipmap(2*i,2*j,w,pd,d-1,l)]+r[coord_mipmap(2*i+1,2*j,w,pd,d-1,l)]+r[coord_mipmap(2*i,2*j+1,w,pd,d-1,l)]+r[coord_mipmap(2*i+1,2*j+1,w,pd,d-1,l)])/4;
//}

#define SIGC 0.36 //paramètre du filtre gaussien au carre
#define TAPSM 5 //racine du nombre de coefficient interpolé, doit être impaire

float filter_mipmap_aux(int u,int v,int w,int pd,int d,int l,float *r,double *g){
	float total=0;
	int i,j;
	for (i=0;i<TAPSM;i++){
		for(j=0;j<TAPSM;j++){
			total = total + g[i+TAPSM*j]*r[coord_mipmap(u-(TAPSM-1)/2+i,v-(TAPSM-1)/2+j,w,pd,d-1,l)];
		}
	}
	return total;
}

float filter_mipmap(int u,int v,int w,int pd,int d,int l,float *r,double *g){
	float a[4];
	a[0] = filter_mipmap_aux(u,v,w,pd,d,l,r,g);
	a[1] = filter_mipmap_aux(u+1,v,w,pd,d,l,r,g);
	a[2] = filter_mipmap_aux(u,v+1,w,pd,d,l,r,g);
	a[3] = filter_mipmap_aux(u+1,v+1,w,pd,d,l,r,g);
	return (a[0]+a[1]+a[2]+a[3])/4;
}

// on admet pour l'instant que l'image et de taille une puissance de 2, cette fonction construit le mip_map d'une image
int build_mipmap(float *img, float *r,int w,int h, int pd){
int i,j,d,l,ll;
double gauss[TAPSM*TAPSM];
for (i=0;i<TAPSM;i++){
	for(j=0;j<TAPSM;j++){
		gauss[i+TAPSM*j]=exp(-(pow(j-(TAPSM-1)/2,2)+pow(i-(TAPSM-1)/2,2))/(2*SIGC))/(2*M_PI*SIGC);   //on a prit sigma = 0,6 (ptet le prendre variable)
	}
}
float norm = 0;
for(i=0;i<pow(TAPSM,2);i++){norm = norm + gauss[i];}
for(i=0;i<pow(TAPSM,2);i++){gauss[i]=gauss[i]/norm;}
for(ll=0;ll<3;ll++){
	if(ll>pd-1){l=pd-1;}else{l=ll;}
	for(i=0;i<w;i++){
		for(j=0;j<h;j++){
			r[(i+j*w)*3+ll]=img[(i+j*w)*pd+l];
		}
	}
}
for(l=0;l<3;l++){
	for(d=2;pow(2,d-1)<=w;d++){
		for(j=0;j<w/pow(2,d-1);j++){
			for(i=0;i<w/pow(2,d-1);i++){
				r[coord_mipmap(i,j,w,3,d,l)]= filter_mipmap(2*i,2*j,w,3,d,l,r,gauss);			
			}
		}
	}
}
return 0;
}



// on a utilisé cela comme ref pour écrire les fonctions
// du/dx = D[1]y+D[2]
// du/dy = D[5]x+D[6]
// dv/dx = D[3]y+D[4]
// dv/dy = D[7]x+D[8]

void precal_D_vsing(double *D){
	//et pour E le coeff dans le polynôme
	D[12]=D[5]*D[7];
	D[13]=D[6]*D[7]+D[5]*D[8];
	D[14]=D[1]*D[3];
	D[15]=D[1]*D[4]+D[3]*D[2];
	D[16]=D[6]*D[8]+D[2]*D[4];
}


//D(x,y)=sqrt(12xy+13x+14y+15)/c^2
void precal_D_det(double *D){
	D[12]=D[1]*D[7]+D[5]*D[3];
	D[13]=D[2]*D[7]+D[5]*D[4];
	D[14]=D[3]*D[6]+D[1]*D[8];
	D[15]=D[2]*D[8]+D[6]*D[4];
}


void precal_D_sqr(double H[3][3],double *D){
	// 0 1 2
	// 3 4 5
	// 6 7 8
	double *a = H[0];
	D[1]=a[0]*a[7]-a[6]*a[1];
	D[2]=a[0]*a[8]-a[6]*a[2];
	D[3]=a[3]*a[7]-a[6]*a[4];
	D[4]=a[3]*a[8]-a[6]*a[5];
	D[5]=-D[1];
	D[6]=a[1]*a[8]-a[7]*a[2];
	D[7]=-D[3];
	D[8]=a[4]*a[8]-a[7]*a[5];
	D[9]=a[6];
	D[10]=a[7];
	D[11]=a[8];
}

void precal_D_prgm(double H[3][3],double *D){
	double *a = H[0];
	D[1]=a[0]*a[7]-a[6]*a[1]-(a[3]*a[7]-a[6]*a[4]);
	D[2]=a[0]*a[8]-a[6]*a[2]-(a[3]*a[8]-a[6]*a[5]);
	D[3]=a[3]*a[7]-a[6]*a[4]-(a[0]*a[7]-a[6]*a[1]);
	D[4]=a[1]*a[8]-a[7]*a[2]-(a[4]*a[8]-a[7]*a[5]);
	D[5]=a[0]*a[7]-a[6]*a[1]+(a[3]*a[7]-a[6]*a[4]);
	D[6]=a[0]*a[8]-a[6]*a[2]+(a[3]*a[8]-a[6]*a[5]);
	D[7]=a[3]*a[7]-a[6]*a[4]+(a[0]*a[7]-a[6]*a[1]);
	D[8]=a[1]*a[8]-a[7]*a[2]+(a[4]*a[8]-a[7]*a[5]);	
	D[9]=a[6];
	D[10]=a[7];
	D[11]=a[8];
}


int precal_D_mipmap(double H[3][3],double *D,int fun_dist){
	if(fun_dist==1 || fun_dist==2){precal_D_prgm(H,D); return 0;}
	if(fun_dist==3 || fun_dist==4){precal_D_sqr(H,D); precal_D_vsing(D); return 0;}
	if(fun_dist==5){precal_D_sqr(H,D); precal_D_det(D); return 0;}
	precal_D_sqr(H,D); return 0;
}



// on note D = (1,2,3,4,5,6,7,8,9,10,11) avec 
// D = max(sqrt( (1y+2)^2 + (3y+4)^2 )),sqrt( (5x+6)^2 + (7x+8)^2 ))/(9x + 10y +11)^2

//lien 1 = 1-3
//2 = 2-4
//3 = 5-7
//4 = 6-8

//et de même avec plus

//Distance variables par defaut max_sqr
//1: prgm diag    max ( sqrt ( (1y+2)^2 + (3x+4)^2 ) , sqrt ( (5y+6)^2 + (7x+8)^2 )/(9x + 10y +11)^2

#define D_BIAS 0.

//double absd(double p){if(p>0){return p;}{return -p;}}

double cal_D_det(int x,int y,double *D){
	double p[2]={x,y};
	double a,b,c;
	c = pow(D[9]*p[0]+D[10]*p[1]+D[11],2);
	a = sqrt(absd(D[12]*p[0]*p[1]+D[13]*p[0]+D[14]*p[1]+D[15]));
	return(a/c+D_BIAS);	
}

// elle prend des int est retourne un double
double cal_D_sqr(int x,int y,double *D){
	double p[2]={x,y};
	double a,b,c;
	c = pow(D[9]*p[0]+D[10]*p[1]+D[11],2);
	c = c;
	a = sqrt (pow(D[1]*p[1]+D[2],2) + pow(D[3]*p[1]+D[4],2));
	b = sqrt (pow(D[5]*p[0]+D[6],2) + pow(D[7]*p[0]+D[8],2));
	if(a>b){return a/c+D_BIAS;}{return b/c+D_BIAS;}
}

//diff subtil, c'est du code pour rien
double cal_D_prgm(int x,int y,double *D,int opt){
	double p[2]={x,y};
	double a,b,c;
	c = pow(D[9]*p[0]+D[10]*p[1]+D[11],2);
	c = c*sqrt(2);
	a = sqrt (pow(D[1]*p[1]+D[2],2) + pow(D[3]*p[0]+D[4],2));
	b = sqrt (pow(D[5]*p[1]+D[6],2) + pow(D[7]*p[0]+D[8],2));
	if(opt==1){if(a>b){return a/c+D_BIAS;}{return b/c+D_BIAS;}}{if(a<b){return a/c+D_BIAS;}{return b/c+D_BIAS;}}
}

double cal_D_vsing(int x,int y,double *D,int opt){
	double p[2]={x,y};
	double A,B,C,E;
	C = pow(D[9]*p[0]+D[10]*p[1]+D[11],2)*sqrt(2);  //(car le /2 dans le pol de deg 2)
	A = pow(D[1]*p[1]+D[2],2) + pow(D[3]*p[1]+D[4],2);
	B = pow(D[5]*p[0]+D[6],2) + pow(D[7]*p[0]+D[8],2);
	E = pow(D[12]*pow(p[0],2)+D[13]*p[0]+D[14]*pow(p[1],2)+D[15]*p[1]+D[16],2);    // E = du/dx dv/dx + du/dy dv/dy, E = (12x^2+13x+14y^2+15y+16)^2
	if(opt==3){
		return sqrt(A+B+sqrt(pow(A+B,2)-4*(A*B-E)))/C + D_BIAS;
	}{
		return sqrt(A+B-sqrt(pow(A+B,2)-4*(A*B-E)))/C + D_BIAS;
	} //et le moins pour avoir lambda au lieu de t //sqrt(A+B-sqrt(pow(A+B,2)-4*(AB-E)))/(2*c)
}


double cal_D_mipmap(int x,int y,double *D,int fun_dist){
	if(fun_dist==1 || fun_dist==2){return cal_D_prgm(x,y,D,fun_dist);}
	if(fun_dist==3 || fun_dist==4){return cal_D_vsing(x,y,D,fun_dist);}
	if(fun_dist==5){return cal_D_det(x,y,D);}
	return cal_D_sqr(x,y,D);
}


//int coord_mipmap (int u,int v,int w, int pd,int d,int l)
// for bilinear interpolation during mipmap
static float bilinear_mipmap(float *x, int w, int h,
		float p, float q, int l, int d)
{
	int ip = floor(p);
	int iq = floor(q);
	float a = x[coord_mipmap(ip,iq,w,3,d,l)];
	float b = x[coord_mipmap(ip+1,iq,w,3,d,l)];
	float c = x[coord_mipmap(ip,iq+1,w,3,d,l)];
	float dd = x[coord_mipmap(ip+1,iq+1,w,3,d,l)];
	return evaluate_bilinear_cell(a, b, c, dd, p-ip, q-iq);
}



//mipmapping interpolator, img is a mipmapping
static float mip_mapping_interpolator_at(float *img, int w, int h, int pd,
		float x, float y, int l,double D,int opt){
	int d ;
	float a,b;
	for(d=0;(pow(2,d))<=D && (pow(2,d))<=w;d++){;}
	if (d==1||d==0){return bilinear_mipmap(img,w,h,x,y,l,1);}
	if ((pow(2,d))>w){return img[6*w*(w-1)+l];} 
	a = bilinear_mipmap(img,w,h,x/pow(2,d-1),y/pow(2,d-1),l,d);
	b = bilinear_mipmap(img,w,h,x/pow(2,d-2),y/pow(2,d-2),l,d-1);
//	if(opt == 7){return ((D-pow(2,d-1))/(pow(2,d-1)))*a + ((pow(2,d)-D)/pow(2,d-1))*b;}
	return (log2(D)-d+1)*a + (d-log2(D))*b;  //((D-pow(2,d-1))/(pow(2,d-1)))*a + ((pow(2,d)-D)/pow(2,d-1))*b; //logarithmique ? (log(D)-d+1)*a + (d-log(D))*b;
}


// on a l'impression que img[i][j][l] avec l la couleur, pd le nombre de couleur = img[(i+j*w)*pd + l]
// et dans le mip-map ?
// r[u][v][d] = h*(w+..+w/(2^(d-1)))*pd + (u+v*w/d)*pd+l

