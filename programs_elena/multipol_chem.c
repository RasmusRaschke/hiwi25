#include <complex.h>  /* This also includes iostream.h. */
#include <string.h>
#include <stdio.h>
#include <fcntl.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <conio.h>

//#define ERROR -1
#define OK 1
#define DEGREE 0.017453292519943296
#define PI 3.1415926535897932

//typedefs
typedef struct _vec3D{
	double x;
	double y;
	double z;
}VEC3D;

typedef struct _ql{
	  int l;
	  int m;
	  double _Complex qlm;
	  struct _ql *next;
}QL;

typedef struct _multipoles{
	int l;
	QL *moments;
	struct _multipoles *next;
}MULTIPOLES;


/*#################################################################*/
/*######################### prototypes ############################*/
/*#################################################################*/

int readMaster(void);
int readCoordinates(void);
int writeCoordinates(char *name);
int writeEnergies(char *name);

double _Complex wignerD(int l, int m,int n,VEC3D *winkel);
double wignerd(int j,int m,int n,double beta,int k);
double fac(long i);
double sqrtfac(long i);
double plgndr(int l, int m, double x);
double _Complex ylm(int l, int m, double theta, double phi);

void createMoments(QL **mHndl, int l, double* array);
void addMultipol(MULTIPOLES **mHndl, int l, double* array);

void readOutMultipol(MULTIPOLES **mHndl);
void freeMultipol(MULTIPOLES **mHndl);
void copyRotateMultipol(MULTIPOLES **dHndl,VEC3D *winkel);

double pairEnergy(MULTIPOLES **aHndl, MULTIPOLES **bHndl,VEC3D *rVec);
double scaleFun(int l1, int l2, int m1, int m2);

void calcEnergies(int j, double *energy);//j ist der Geänderte in der Zählweise 1...N
unsigned long int nu(int n, int k);// summe von i=1 bis k über (N-i)

void newAngle3D(VEC3D *newAngle);
void newAngleIsing(VEC3D *newAngle);
void newAngleXY(VEC3D *newAngle);
void replaceEnergies(int j);
void montecarlo(int i, int frequence ,int offset,double kT1,double kT2,double kappa);

void xyz_To_RThetaPhi( VEC3D *rVec );
void printVec3D( VEC3D *vec);
void fprintVec3D( VEC3D *vec );
void rDiff( VEC3D *R1, VEC3D *R2,VEC3D *R3);
void copyVec3D( VEC3D *s, VEC3D *d);


/*#################################################################*/
/*########################## globals ##############################*/
/*#################################################################*/

//complex gI(0.,1.);

MULTIPOLES *master=NULL;
MULTIPOLES *mp1=NULL;
MULTIPOLES *mp2=NULL;

VEC3D **coordinates;
unsigned long int eSize;
double *energies; //is double huge *energies; in borland
double *energiesTmp;
double eGes;

int cLength;

FILE *gFilePntr;

/*#################################################################*/
/*##########################   main  ##############################*/
/*#################################################################*/

int main(void){
	double _Complex x;
	double kT1,kT2,kappa;
	FILE *inData;
	int i,f,offset;
	unsigned long int uli;

	gFilePntr=fopen("energy.bin","ab");

	if(!(inData=fopen("indata.bin","rb")))return(ERROR);

	fread(&i,sizeof(i),1,inData);
	fread(&f,sizeof(f),1,inData);
	fread(&offset,sizeof(offset),1,inData);
	fread(&kT1,sizeof(double),1,inData);
	fread(&kT2,sizeof(double),1,inData);
	fread(&kappa,sizeof(double),1,inData);
	fclose(inData);


	readMaster();
	readOutMultipol(&master);

	readCoordinates();

	eSize=(((unsigned long int)(cLength)*(unsigned long int)(cLength-1))/2);
	energies=(double *)calloc(eSize,sizeof(double));//Array für alle Paarenergien is huge & farcalloc in borland
	energiesTmp=(double*)calloc(cLength,sizeof(double));//Array für die neu zu berechnenden Energien

	printf("cLength=%d\n",cLength);
	printf("eSize=%ld\n",eSize);
	printf("=>%gMb\n",8.*(double)eSize/(1024*1024));

	for(uli=0;uli<eSize;uli++)energies[uli]=0.0;
	for(uli=0;uli<cLength;uli++)energiesTmp[uli]=0.0;

	eGes=0;
	for(uli=1;uli<=cLength;uli++)
	{
		printf("\r");
		printf("k(%d)=%ld\tE=%5.5e",cLength,uli,eGes);
		calcEnergies(uli,&eGes);
	}
	//mit der Eingabe haben sich virtuell alle Multipole geändert.

	printf("\n eGes(Start)=%g\n",eGes);

	printf("i=%d  f=%d   kT1=%g  kT2=%g  kappa=%g\n",i,f,kT1,kT2,kappa);

	montecarlo(i,f,offset,kT1,kT2,kappa);

	fclose(gFilePntr);


	free(energies);// is farfree in borland
	free(energiesTmp);
	freeMultipol(&master);
	freeMultipol(&mp1);
	freeMultipol(&mp2);
	return OK;
}


/*#########################################################
#####################   file operations   #################
#########################################################*/

int readMaster(void)
{
	FILE *filePntr;
	int l,i;
	double *myArray;

	if(!(filePntr=fopen("master.bin","rb")))return(ERROR);

	fread(&l,sizeof(l),1,filePntr);
	do
	{
		myArray= (double*)calloc(2*(2*l+1),sizeof(double));//2*, da real und imag
		for(i=0;i<(2*(2*l+1));i++)
		{
			fread(&(myArray[i]),8,1,filePntr);
			//printf("%d->%g\t",i,myArray[i]);
		}
		//printf("\n");

		addMultipol(&master,l,myArray);
		addMultipol(&mp1,l,myArray);//damit sie überhaupt existieren!
		addMultipol(&mp2,l,myArray);

		free(myArray);
		fread(&l,sizeof(l),1,filePntr);//read next

	}while(!feof(filePntr));

	fclose(filePntr);

	return(OK);
}


/*############################################*/


int readCoordinates(void)
{
	FILE *filePntr;
	int i;

	if(!(filePntr=fopen("lattice.bin","rb")))
	{
		printf("could not open lattice.bin");
		return(ERROR);
	}

	fread(&cLength,sizeof(cLength),1,filePntr);//länge des arrays, ist global
	coordinates= calloc(cLength,sizeof(VEC3D*));
	for(i=0;i<cLength;i++)coordinates[i]= calloc(2,sizeof(VEC3D));//(x,y,z)&(aplha,beta,gamma)

	for(i=0;i<cLength;i++)
	{
		fread(&(coordinates[i][0]),sizeof(VEC3D),1,filePntr);
		fread(&(coordinates[i][1]),sizeof(VEC3D),1,filePntr);
	}

	fclose(filePntr);

	return(OK);
}

/*############################################*/


int writeCoordinates(char *name)
{
	FILE *filePntr;
	int i;

	if(!(filePntr=fopen(name,"wb")))
	{
		printf("could not open %s\n",name);
		return(ERROR);
	}

	fwrite(&cLength,sizeof(cLength),1,filePntr);//länge des arrays, ist global

	for(i=0;i<cLength;i++)
	{
		fwrite(&(coordinates[i][0]),sizeof(VEC3D),1,filePntr);
		fwrite(&(coordinates[i][1]),sizeof(VEC3D),1,filePntr);
	}

	fclose(filePntr);

	return(OK);
}

/*############################################*/

int writeEnergies(char *name)
{
	FILE *filePntr;
	int i;

	if(!(filePntr=fopen(name,"wb")))
	{
		printf("could not open %s\n",name);
		return(ERROR);
	}

	fwrite(&eSize,sizeof(eSize),1,filePntr);//länge der liste, ist global

	for(i=0;i<eSize;i++)
	{
		fwrite(&(energies[i]),sizeof(double),1,filePntr);
	}

	fclose(filePntr);

	return(OK);
}

/*#########################################################
#####################   wigner function   #################
#########################################################*/

double _Complex wignerD(int l, int m,int n,VEC3D *winkel)
{
  double _Complex myOut=0;
  double mySum=0;
  int i;

  	myOut.im-=(winkel->x*(double)m+winkel->z*(double)n);//damit ist es äquivalent zu (-I)*(...)
	myOut=cexp(myOut);

	for(i=0;i<=(l+n);i++) mySum+=wignerd(l,m,n,winkel->y,i);
	return(myOut*mySum);
}

double wignerd(int j,int m,int n,double beta,int k)
{
	long sig;
	double N[5];
	double D[5];
	long test, cpow, spow;
	double num_den=1;
	int i;

		sig = j + m + k;
		N[1] = sqrtfac((long)(j + n));
		N[2] = sqrtfac((long)(j - n));
		N[3] = sqrtfac((long)(j + m));
		N[4] = sqrtfac((long)(j - m));
		D[1] = j + m - k;
		D[2] = k;
		D[3] = j + n - k;
		D[4] = k - n - m;

		if((D[1] < 0.) || (D[2] < 0.) || (D[3] < 0.) || (D[4] < 0.))test=1;
		else test=0;

		//d1 checks whether kmax=min(j+m,j+n)*)
		//d2 checks whether k>=0 should be redundant*)
		//d3 checks whether k<=j+n should be redundant*)
		//d4 checks whether kmin=max(0,m+n)*)
		if(test==1)return(0);

		D[1] = fac(D[1]);
		D[2] = fac(D[2]);
		D[3] = fac(D[3]);
		D[4] = fac(D[4]);

		for(i=1;i<5;i++)num_den *= (N[i]/D[i]);
		spow = 2*j - 2*k + n + m;
		cpow = 2*k - n - m;

		if( (fabs(beta)<1.e-15)&&(spow==0))if(n==m)return(1);else return(0);
		if( (fabs(beta-180.*DEGREE)<1.e-15)&&(cpow==0))
			if(n+m==0) return(pow(-1.,j+m)); else return(0);

	  // one has to calculate <jm'|exp(beta Jy)|jm>*)
	  // if beta=0 that is <jm'|1|jm> = <jm'|jm> = delta(m,m')*)
	  // if beta=180°,djmn(180°)=(-1)^(j+m)delta(m,-n)*)
	  // this must be checked to avoid 0^0 ! *)
	  // if n != m the algorithm would work.*)
	  //cout << "{erg=" <<	pow(-1,sig)*num_den*pow(cos(beta/2.),cpow)*pow( sin(beta/2.),spow)<< "}\n";
	  return(
		pow(-1,sig)*num_den*pow(cos(beta/2.),cpow)*pow( sin(beta/2.),spow)
		);

}



/*#########################################################
##########################   basic math   #################
#########################################################*/

//###############fakultät

double fac(long i)
{
	long j;
	double out=1;

	for(j=1;j<=i;j++)out*=(double)j;
 	return(out);
}

//##############sqrtfakultät
double sqrtfac(long i)
{
	long j;
	double out=1;

 	for(j=1;j<=i;j++)out*=sqrt((double)j);
 	return(out);
}

//##################legendre

double plgndr(int l, int m, double x)
{

	double fact,pll,pmm,pmmp1,somx2;
	int i,ll;


	if ( (m > l) || fabs(x) > 1.0) return(ERROR);


	pmm=1.0;
	if (m > 0)
	{
		somx2=sqrt((1.0-x)*(1.0+x));
		fact=1.0;
		for (i=1;i<=m;i++)
		{
			pmm *= -fact*somx2;
			fact += 2.0;
		}
	}
	if (l == m) return(pmm);
	else
	{
		pmmp1=x*(2*m+1)*pmm;
		if (l == (m+1)) return(pmmp1);
		else
		{
			for (ll=m+2;ll<=l;ll++)
			{
				pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m);
				pmm=pmmp1;
				pmmp1=pll;
			}
			return (pll);
		}
	}
}

//sperical harmonic

double _Complex ylm(int l, int m, double theta, double phi)
{
	double out=0, sign=0;
	double _Complex outC=0;
	int M;

	M=abs(m);

	if(m>0)
	{
		sign=1.;
		outC.im=phi*(double)M;
		outC=cexp(outC);//is also exp(I m phi)
	}
	else
	{
		sign=pow(-1.,m);
		outC.im-=phi*(double)M;
		outC=cexp(outC);//(-1^m*conjugate(...))
	}

	out=sqrtfac( (long)(l+M));
	out=sqrtfac((long)(l-M))/out;
	out*=sqrt((2.*l+1.) / (4.*PI));

	//cout<<out<<"..."<<plgndr( l,M,cos(theta) )<<"..."<<exp(gI*m*phi)<<"\n";

	out*=plgndr( l,M,cos(theta) );

	outC*=(out*sign);

	return(outC);
}
/*#########################################################
########################   list operations   ##############
#########################################################*/

void createMoments(QL **mHndl,int l, double* array)
{
		QL *newPntr;
		QL *tmpPntr=NULL;
		int m;
		double _Complex qlm=0;
		double c;
		double r;

		for(m=l;m>=-l;m--)//schleife läuft rückwärts->pointer auf next ist easy!
		{
			newPntr= calloc(1,sizeof(QL));
			newPntr->l=l;
			newPntr->m=m;

			r=array[2*(l+m)];
			c=array[2*(l+m)+1];
			qlm.re=r;
			qlm.im=c;

			newPntr->qlm=qlm;
			newPntr->next=tmpPntr;
			tmpPntr=newPntr;
		}

		*mHndl=newPntr;
		return;
}

//##########################

void addMultipol(MULTIPOLES **mHndl, int l, double* array)
{
	 MULTIPOLES *newPntr;
	 MULTIPOLES *tmpPntr;


	 newPntr=*mHndl;
	 //falls schon einträge existieren, suche den letzten
	 if((*mHndl)!=NULL)
	 {
		tmpPntr=(*mHndl)->next;

		while(tmpPntr!=NULL)
		{
			newPntr=tmpPntr;
			tmpPntr=tmpPntr->next;
		}
	 }
	 tmpPntr=newPntr;

	 //erzeuge neuen
		newPntr= calloc(1,sizeof( MULTIPOLES ));
		newPntr->l=l;
		//mit den momenten:
		createMoments(&(newPntr->moments),l, array);
		newPntr->next=NULL;


		if((*mHndl)==NULL)//falls noch keiner existiert zeigt mHndl auf den erzeugten.
		{
			(*mHndl)=newPntr;
		}
		else //der letzte muß auf den neuen verweisen
		{
			tmpPntr->next=newPntr;
		}
		return;
}

//################################

void readOutMultipol(MULTIPOLES **mHndl)
{
	MULTIPOLES *tmpPntr;
	QL *qlPntr;

	tmpPntr=(*mHndl);
	while(tmpPntr!=NULL)
	{
		qlPntr=tmpPntr->moments;
		printf("\nL=%d\n",tmpPntr->l);
		while(qlPntr!=NULL)
		{
			printf("q(%d,%d)=(%g+i*%g)",qlPntr->l,qlPntr->m,creal(qlPntr->qlm),cimag(qlPntr->qlm));
			qlPntr=qlPntr->next;
		}
		tmpPntr=tmpPntr->next;
		printf("\n");
	}
	return;

}

//################################

void freeMultipol(MULTIPOLES **mHndl)
{
	MULTIPOLES *tmpPntr,*nextPntr;
	QL *qlPntr,*nextqlPntr;

	tmpPntr=(*mHndl);
	while(tmpPntr!=NULL)
	{
		qlPntr=tmpPntr->moments;
		while(qlPntr!=NULL)
		{
			nextqlPntr=qlPntr->next;
			free(qlPntr);
			qlPntr=nextqlPntr;
		}
		nextPntr=tmpPntr->next;
		free(tmpPntr);
		tmpPntr=nextPntr;

	}
	return;

}

//###########################

void copyRotateMultipol(MULTIPOLES **dHndl,VEC3D *winkel) //source is always master
{
	MULTIPOLES *sTmpPntr,*dTmpPntr;
 	QL *fixedQl,*sTmpQl,*dTmpQl;
 	int l,m,n;
 	double _Complex qValue,wD;

 	sTmpPntr=master;
 	dTmpPntr=*dHndl;

 	while(sTmpPntr!=NULL)//für alle L
 	{
		l=sTmpPntr->l;
		dTmpQl=dTmpPntr->moments;
		fixedQl=sTmpPntr->moments;

		while(dTmpQl!=NULL)//für alle n
		{
			qValue=0;
			//fprintf(gFilePntr,"new m\n");
			n=dTmpQl->m;
			sTmpQl=fixedQl;
			while(sTmpQl!=NULL)//für alle m
			{
				m=sTmpQl->m;
				wD=wignerD(l,m,n,winkel);
				wD*=(sTmpQl->qlm);
				qValue=(qValue+wD);
				sTmpQl=sTmpQl->next;
			}
			dTmpQl->qlm=qValue;
			dTmpQl=dTmpQl->next;
		}
		sTmpPntr=sTmpPntr->next;
		dTmpPntr=dTmpPntr->next;
 	}
 	return;
}

/*#########################################################
#######################   energy operations   #############
#########################################################*/

double pairEnergy(MULTIPOLES **aHndl, MULTIPOLES **bHndl,VEC3D *rVec)
{

	double energy=0;
	double _Complex eSum=0;
	MULTIPOLES *aTmpPntr,*bTmpPntr;
	QL *aTmpQ,*bTmpQ;
	int l1,l2,m1,m2;
	double _Complex aQ,bQ;
	double scale,rPow;
	double _Complex Y;

	aTmpPntr=*aHndl;

	xyz_To_RThetaPhi(rVec);

 	while(aTmpPntr!=NULL)
 	{
		aTmpQ=aTmpPntr->moments;
		while(aTmpQ!=NULL)
		{
			l1=aTmpQ->l;
			m1=aTmpQ->m;
			aQ=conj(aTmpQ->qlm);

			bTmpPntr=*bHndl;
			while(bTmpPntr!=NULL)
			{
				bTmpQ=bTmpPntr->moments;
				while(bTmpQ!=NULL)
				{
					l2=bTmpQ->l;
					m2=bTmpQ->m;
					bQ=conj(bTmpQ->qlm);
					scale=scaleFun(l1,l2,m1,m2);
					rPow=pow(rVec->x,l1+l2+1);
					Y=ylm(l1+l2,-m1-m2, rVec->y,rVec->z);
					Y/=rPow;
					Y*=scale;
					eSum=(eSum+(Y*aQ*bQ));
					bTmpQ=bTmpQ->next;
				}
				bTmpPntr=bTmpPntr->next;
			}
			aTmpQ=aTmpQ->next;
		}
		aTmpPntr=aTmpPntr->next;
  	}
	energy=creal(eSum);
	return energy;
}


//###########################

double scaleFun(int l1, int l2, int m1, int m2)
{

	double num1,num2,den1,den2,den3,den4;

  num1=sqrtfac(l1+l2+m1+m2);
  num2=sqrtfac(l1+l2-m1-m2);
  den1=sqrtfac(l1+m1);
  den2=sqrtfac(l1-m1);
  den3=sqrtfac(l2+m2);
  den4=sqrtfac(l2-m2);
  num1/=den1;
  num1/=den2;
  num2/=den3;
  num2/=den4;

  num1*=num2;
  /*this factor if the multipole moments are defined in the usual way in physics...
    so they are q(l,m) in contrast to Q(l,m) in chemestry.
   */
  /*num1*=(sqrt(
				pow(4.*PI,3)/(
				(2.*(l1+l2)+1)*(2.*l1+1)*(2.*l2+1)
				)
				)
			);
			*/
  /*this factor if the multipole moments are defined in the usual way in chemestry...
    so they are Q(l,m) in contrast to q(l,m) in physics.
   */
  num1*=(sqrt(
				(4.*PI)/( (2.*(l1+l2)+1) )
			 )
		);
	return(num1*pow(-1.,-l2+m1+m2));
}

//###########################

void calcEnergies(int j, double *energy)//geändert wurde j(in der  Nummerierung 1...N)
{
	//*energy ist ein pointer auf locE in monteCarlo
	int k;
	unsigned long int position,nuJ;
	double e=0;
	VEC3D *r1,*r2,deltaR;
	r1=&(coordinates[j-1][0]);
	copyRotateMultipol(&mp1,&(coordinates[j-1][1]));

	//zeile: alle i mit i<j
	for(k=0;k<(j-1);k++)
	{
		position= nu(cLength,k)+(unsigned long int)(j-k-2);//wo in der paarwechselwirkungsliste?
		energiesTmp[k]=energies[position];
		//merke die alte Energie
		(*energy)-=energiesTmp[k];
		//entferne aus der alten energie die alte paarwechselwirkung;
		r2=&(coordinates[k][0]);
		copyRotateMultipol(&mp2,&(coordinates[k][1]));
		rDiff(r1, r2, &deltaR);
		e=pairEnergy(&mp1, &mp2, &deltaR);
		*energy+=e;
		//addiere die neue paarwechselwirkung
		energies[position]=e;
		//schreibe die neue energie in die paarwechselwirkungsliste.

	}
	nuJ=nu(cLength,j-1);
	//spalte: alle i in der j-ten spalte mit i>j (nie j-j, immer i!=j)
	for(k=(j-1);k<(cLength-1);k++)
	{
		position= nuJ+(unsigned long int)(k-j+1);

		energiesTmp[k]=energies[position];
		*energy-=energiesTmp[k];
		r2=&(coordinates[k+1][0]);
		copyRotateMultipol(&mp2,&(coordinates[k+1][1]));
		rDiff(r1, r2, &deltaR);
		e=pairEnergy(&mp1, &mp2, &deltaR);
		*energy+=e;
		energies[position]=e;

	}

	return;
}

//###########################

unsigned long int nu(int n, int k)
//damit bloß kein variablentypüberlauf stattfindet auf nummer sicher!
{
	unsigned long int out;
	out=(2*n-k-1);
	out*=(unsigned long int)k;
	out/=2;
	return(out);
}


/*#########################################################
################### monte carlo operations ################
#########################################################*/

 void newAngle3D(VEC3D *newAngle)
 {
	double a,b,g;

			a=((double)rand())/((double)RAND_MAX)*2*PI;
			g=((double)rand())/((double)RAND_MAX)*2*PI;
			b=((double)rand())/((double)RAND_MAX)*2-1;
			b=acos(b);
			newAngle->x=a;
			newAngle->y=b;
			newAngle->z=g;

	return;
 }

//###################

 void newAngleIsing(VEC3D *newAngle)
 {
	double a,b,g;

			a=((double)rand())/((double)RAND_MAX)*2*PI;
			g=0.0;
			b=2*( ( (double)( rand()%2 ) )-.5);
			b=acos(b);
			newAngle->x=a;
			newAngle->y=b;
			newAngle->z=g;

	return;
 }

//###################

 void newAngleXY(VEC3D *newAngle)
 {
	double a,b,g;

			a=((double)rand())/((double)RAND_MAX)*2*PI;
			g=0.0;
			b=0.0;
			newAngle->x=a;
			newAngle->y=b;
			newAngle->z=g;

	return;
 }

//###################

void replaceEnergies(int j)
{
	int k;
	unsigned long int position, nuJ;

	for(k=0;k<(j-1);k++)
	{
		position= nu(cLength,k)+(unsigned long int)(j-k-2);//wo in der paarwechselwirkungsliste?
		energies[position]=energiesTmp[k];
		//schreibe die ursprüngliche energie in die paarwechselwirkungsliste.
	}
	nuJ=nu(cLength,j-1);
	for(k=(j-1);k<(cLength-1);k++)
	{
		position= nuJ+(unsigned long int)(k-j+1);
		energies[position]=energiesTmp[k];
	}
	return;
}

//###################


void montecarlo(int i, int frequence, int offset ,double kT1,double kT2,double kappa)
{
	int j,k;
	int l=0;
	int sig=30;
	VEC3D tmpAngle,newAngle;
	double a,eLoc,myExp;
	time_t t;
	char name[100];
	char numb[100];
	char *begin="coor";
	char *beginE="enrg";
	char *endung="bin";//der punkt für *.bin kommt in diesem compiler _gvct
	unsigned long int excepted=0;
	char c;
	double kT,x0,y0,locExp;

	locExp=(-1.)*kappa*(double)(i-1);

	y0=(kT2-kT1*exp( locExp ))/(1-exp(locExp));
	x0=(kT1-kT2)/(1-exp( locExp ));


	srand((unsigned) time(&t));

	for(j=(0+offset);j<(i+offset);j++)
	{


		kT=y0+x0*exp((-1.)*kappa*(double)(j-offset));
		if(kT<0)kT=0;

		for(k=0;k<cLength;k++)
		{


		printf("\r");
		printf("j=%d\tk(%d)=%d\tE=%5.5e\tkT=%g\tt-Up=%ld",j,cLength,k+1,eGes,kT,excepted);
		fflush(stdout);

			//newAngleXY(&newAngle);
			newAngle3D(&newAngle);
			//newAngleIsing(&newAngle);

			copyVec3D(&(coordinates[k][1]),&tmpAngle);
			copyVec3D(&newAngle,&(coordinates[k][1]));

			eLoc=eGes;
			calcEnergies(k+1,&eLoc);


			if(eLoc<eGes)//ist die neue energie kleiner
			{
				eGes=eLoc;
			}
			else
			{
				a=((double)rand())/((double)RAND_MAX);
				myExp=(eLoc-eGes)/kT;
				if(myExp>705.)myExp=705.;//nur dann macht e^... überhaupt sinn
				if( a<exp( (-1)*myExp) )//vielleicht wegen der temperatur trotzdem?
				{
					excepted++;
					eGes=eLoc;
				}
				else
				{
					copyVec3D(&tmpAngle,&(coordinates[k][1]));
					//setze den winkel zurück
					replaceEnergies(k+1);
					//setze die energien zurück

				}

			}
		}


		fwrite(&j,sizeof(int),1,gFilePntr);
		fwrite(&kT,sizeof(double),1,gFilePntr);
		fwrite(&eGes,sizeof(double),1,gFilePntr);

		if(l>=frequence)
		{
			strcpy(name,begin);
			_gcvt(j,sig,numb);
			strcat(name,numb);
			strcat(name,endung);
			writeCoordinates(name);
			/*
			eTest=0;
			eSum=0;
			for(uli=0;uli<eSize;uli++)eSum+=energies[uli];
			for(uli=0;uli<eSize;uli++)energies[uli]=0.0;
			for(uli=0;uli<cLength;uli++)energiesTmp[uli]=0.0;
			for(uli=1;uli<=cLength;uli++)calcEnergies(uli,&eTest);
			printf("\nj=%d\teGes=%g\teTest=%g\teSum=%g\n",j,eGes,eTest,eSum);
			*/
			strcpy(name,beginE);
			_gcvt(j,sig,numb);
			strcat(name,numb);
			strcat(name,endung);
			writeEnergies(name);
			l=0;
			printf("\n");
		}
		l++;
	}
			strcpy(name,begin);
			_gcvt(j,sig,numb);
			strcat(name,numb);
			strcat(name,endung);
			writeCoordinates(name);
			printf("\nj=%d\teGes=%g\n",j,eGes);
			strcpy(name,beginE);
			_gcvt(j,sig,numb);
			strcat(name,numb);
			strcat(name,endung);
			writeEnergies(name);
			l=0;
	return;
}

/*#########################################################
#######################   struct operations   #############
#########################################################*/

void xyz_To_RThetaPhi(VEC3D *rVec)
{
	double r,theta;

	r=sqrt(rVec->x*rVec->x+rVec->y*rVec->y+rVec->z*rVec->z);
	theta=acos(rVec->z/r);
	if((fabs(rVec->x)<1.e-15)&&(fabs(rVec->y)<1.e-15))	(rVec->z)=0;
	else (rVec->z)=atan2(rVec->y,rVec->x);//Phi
	(rVec->y)=theta;
	(rVec->x)=r;
	return;
}

//###########################

void rDiff(VEC3D *R1, VEC3D *R2,VEC3D *R3)
{
	R3->x=R2->x-R1->x;
	R3->y=R2->y-R1->y;
	R3->z=R2->z-R1->z;
	return;
}

//###########################

void copyVec3D(VEC3D *s, VEC3D *d)
{
	d->x=s->x;
	d->y=s->y;
	d->z=s->z;
	return;
}

//###########################

void printVec3D(VEC3D *vec)
{
 printf("(%g,%g,%g)",vec->x,vec->y,vec->z);
 return;
}

void fprintVec3D(VEC3D *vec)
{
 fprintf(gFilePntr,"(%g,%g,%g)",vec->x,vec->y,vec->z);
 return;
}
