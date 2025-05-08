#include <math.h>
#include <complex.h>
#include <stdlib.h>
#include <stdio.h>

#define OK 1
#define ERROR -1
#define DEGREE 0.017453292519943296
#define PI 3.1415926535897932
#define ylm_Table_Max 3

/*#########################################################
####################  typrdefinitions  ################
#########################################################*/

typedef complex<double> (*ylm_Pntr)(int, int ,double, double);

/*#########################################################
#################### prototypes ################
#########################################################*/

double fac(long i);
double facsqrt(long i);
double plgndr(int l, int m, double x);

complex<double> chooseYlm(int l, int m, double theta, double phi);
//*** l>3
complex<double> ylm(int l, int m, double theta, double phi);
//*** l=0
complex<double> y00(int l, int m, double theta, double phi);
//*** l=1
complex<double> y1m1(int l, int m, double theta, double phi);
complex<double> y1(int l, int m, double theta, double phi);
complex<double> y1p1(int l, int m, double theta, double phi);
complex<double> y2m2(int l, int m, double theta, double phi);
//*** l=2
complex<double> y2m2(int l, int m, double theta, double phi);
complex<double> y2m1(int l, int m, double theta, double phi);
complex<double> y2(int l, int m, double theta, double phi);
complex<double> y2p1(int l, int m, double theta, double phi);
complex<double> y2p2(int l, int m, double theta, double phi);
//*** l=3
complex<double> y3m3(int l, int m, double theta, double phi);
complex<double> y3m2(int l, int m, double theta, double phi);
complex<double> y3m1(int l, int m, double theta, double phi);
complex<double> y3(int l, int m, double theta, double phi);
complex<double> y3p1(int l, int m, double theta, double phi);
complex<double> y3p2(int l, int m, double theta, double phi);
complex<double> y3p3(int l, int m, double theta, double phi);


/*#########################################################
######################### globals  #######################
#########################################################*/

ylm_Pntr ylm_Array[(ylm_Table_Max+1)*(ylm_Table_Max+1)]={
y00,
y1m1,y1,y1p1,
y2m2,y2m1,y2,y2p1,y2p2,
y3m3,y3m2,y3m1,y3,y3p1,y3p2,y3p3
};


/*#########################################################
##########################   main  #####33#################
#########################################################*/

void main(void)
{
        double phi,theta;
        double div=PI/5.;
        complex<double> Y;
        for(theta=0.;theta<=PI;theta+=div)
        {
               for(phi=0.;phi <=(2*PI);phi+=div)
               {
                 Y=ylm(2,-1,theta,phi);
                 cout<<"{"<<real(Y)<<","<<imag(Y)<<"}_";
               }
               cout<<"\n";
        }
        return;

}
/*#########################################################
######################## hardcoded faculty ################
#########################################################*/

static double fac_Table[61]={1.,1., 2., 6., 24., 120., 720.,
 5040., 40320., 362880.,
 3.6288e6, 3.99168e7, 4.790016e8, 6.2270208e9,
 8.71782912e10, 1.307674368e12, 2.0922789888e13,
 3.55687428096e14, 6.402373705728e15, 1.21645100408832e17,
  2.43290200817664e18, 5.109094217170944e19, 1.124000727777608e21,
  2.585201673888498e22, 6.204484017332394e23, 1.551121004333099e25,
  4.032914611266056e26, 1.088886945041835e28, 3.048883446117139e29,
  8.841761993739702e30, 2.652528598121911e32, 8.222838654177923e33,
  2.631308369336935e35, 8.683317618811886e36, 2.952327990396041e38,
  1.033314796638614e40, 3.719933267899012e41, 1.376375309122635e43,
  5.230226174666011e44, 2.039788208119744e46, 8.159152832478977e47,
  3.345252661316381e49, 1.40500611775288e51, 6.041526306337384e52,
  2.658271574788449e54, 1.196222208654802e56, 5.502622159812089e57,
  2.586232415111682e59, 1.241391559253607e61, 6.082818640342676e62,
  3.041409320171338e64, 1.551118753287382e66, 8.065817517094388e67,
  4.274883284060026e69, 2.308436973392414e71, 1.269640335365828e73,
  7.109985878048635e74, 4.052691950487722e76, 2.350561331282879e78,
  1.386831185456898e80, 8.32098711274139e81};

 static double fac_Sqrt_Table[101]={1., 1., 1.414213562373095, 2.449489742783178,
4.898979485566356, 10.95445115010332, 26.83281572999748, 70.99295739719539,
200.7984063681781, 602.3952191045344, 1904.940943966505, 6317.974358922328,
21886.10518114176, 78911.47445080468, 295259.7012800765, 1.143535905863913e6,
  4.574143623455652e6, 1.885967730625315e7, 8.001483428544984e7,
  3.487765766344294e8, 1.559776268628498e9, 7.147792818185866e9,
  3.352612008237171e10, 1.607856235454059e11, 7.876854713229383e11,
  3.938427356614691e12, 2.008211794424596e13, 1.04349745809074e14,
  5.521669535672285e14, 2.973510046012911e15, 1.628658527169496e16,
  9.067986906793549e16, 5.129628026803635e17, 2.946746955341073e18,
  1.718233974287565e19, 1.016520927791757e20, 6.099125566750542e20,
  3.709953246501409e21, 2.28696877430935e22, 1.428211541796153e23,
  9.032802905233224e23, 5.783815921445271e24, 3.748341123420973e25,
  2.457951648494613e26, 1.630420674178431e27, 1.093719437815202e28,
  7.417966136220958e28, 5.085501366740237e29, 3.523338699662023e30,
  2.466337089763416e31, 1.743963680863606e32, 1.245439180886559e33,
  8.980989654316716e33, 6.538259159791714e34, 4.804619624270389e35,
  3.563201278858419e36, 2.666455677120592e37, 2.013129889124823e38,
  1.533154046820762e39, 1.177637968756484e40, 9.121944481710788e40,
  7.124466393192018e41, 5.609810447812648e42, 4.452649004137245e43,
  3.562119203309796e44, 2.871872314724746e45, 2.333120097803461e46,
  1.909741105966688e47, 1.574812859496909e48, 1.308137807832727e49,
  1.094466613011557e50, 9.222139602976428e50, 7.825244940376377e51,
  6.685892207860283e52, 5.751421947239992e53, 4.980877514193197e54,
  4.342228346904444e55, 3.810289910601106e56, 3.365156932181068e57,
  2.991016905800262e58, 2.675246849288189e59, 2.40772216435937e60,
  2.180285150390389e61, 1.986334304622628e62, 1.820505461284133e63,
  1.678423103505356e64, 1.556505553593457e65, 1.451811729660402e66,
  1.361920123419132e67, 1.284832874770429e68, 1.218899489080934e69,
  1.162756005221389e70, 1.115276380752381e71, 1.075533591796017e72,
  1.042768505784838e73, 1.016365017512855e74, 9.958302741285533e74,
  9.807790764615756e75, 9.709217501366033e76, 9.66054943799493e77,
  9.66054943799493e78};

/*#########################################################
########################## hardcoded Ylm ####################
#########################################################*/

#pragma option push
#pragma option -w-par
//### l=0
complex<double> y00(int l, int m, double theta, double phi) {return (0.28209479177387814,0);}
//### l=1
complex<double> y1m1(int l, int m, double theta, double phi)
{
   complex<double> out=(cos(phi),-sin(phi));
   return (sqrt(3./(8.*PI))*sin(theta)*out );
};
complex<double> y1(int l, int m, double theta, double phi){return (sqrt(3/(4*PI))*cos(theta),0);}
complex<double> y1p1(int l, int m, double theta, double phi){return ( (-1.)*conj( y1m1(l,m,theta,phi) ) );};
//### l=2
complex<double> y2m2(int l, int m, double theta, double phi)
{
        complex<double> out=(cos(2.*phi),-sin(2.*phi));
        return (sqrt(15./(32.*PI))*sin(theta)*sin(theta)*out);
}
complex<double> y2m1(int l, int m, double theta, double phi)
{
        complex<double> out=(cos(phi),-sin(phi));
        return (sqrt(15./(8.*PI))*cos(theta)*sin(theta)*out);
}
complex<double> y2(int l, int m, double theta, double phi){ return( sqrt(5/(16*PI))*(3*cos(theta)*cos(theta)-1.),0);}
complex<double> y2p1(int l, int m, double theta, double phi){return ( (-1.)*conj( y2m2(l,m,theta,phi) ) );}
complex<double> y2p2(int l, int m, double theta, double phi){return ( conj( y2m2(l,m,theta,phi) ) );}  //*(-1)^2
//### l=3
complex<double> y3m3(int l, int m, double theta, double phi)
{
        complex<double> out=(cos(3.*phi),-sin(3.*phi));
        return (sqrt(35./(64.*PI))*sin(theta)*sin(theta)*sin(theta)*out);
}
complex<double> y3m2(int l, int m, double theta, double phi)
{
        complex<double> out=(cos(2.*phi),-sin(2.*phi));
        return (sqrt(105./(32.*PI))*cos(theta)*sin(theta)*sin(theta)*out);
}
complex<double> y3m1(int l, int m, double theta, double phi)
{
        complex<double> out=(cos(phi),-sin(phi));
        return (sqrt(21./(64.*PI))*(5*cos(theta)*cos(theta)-1.)*sin(theta)*out);
}
complex<double> y3(int l, int m, double theta, double phi)
{
        return (sqrt(7./(64.*PI))*(5*cos(2.*theta)*cos(2*theta)-1.)*cos(theta),0);
}
complex<double> y3p1(int l, int m, double theta, double phi){return ( (-1.)*conj( y3m1(l,m,theta,phi) ) );}
complex<double> y3p2(int l, int m, double theta, double phi){return ( conj( y3m2(l,m,theta,phi) ) );}
complex<double> y3p3(int l, int m, double theta, double phi){return ( (-1.)*conj( y2m2(l,m,theta,phi) ) );}

#pragma option pop
  /*#########################################################
##########################   basic math   #################
#########################################################*/

//###############fakultät

double fac(long i)
{

        if(i<=60)return fac_Table[i];
        else
        {
                long j;
	        double out=1;
                for(j=1;j<=i;j++)out*=(double)j;
 	        return(out);
        }
}

//##############sqrtfakultät

double sqrtfac(long i)
{
        if(i<=100)return fac_Sqrt_Table[i];
        else
        {
                long j;
	        double out=1;
                for(j=1;j<=i;j++)out*=sqrt((double)j);
 	        return(out);
        }
}

//##################legendre  (from numerical recepies)

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

complex<double> ylm(int l, int m, double theta, double phi)
{
	double out, sign;
	complex<double> outC(0,0);
	int M;

	M=abs(m);

	if(m>0)
	{
		sign=1.;
		outC=complex<double>(0.,phi*(double)M);
		outC=exp(outC);//is also exp(I m phi)
	}
	else
	{
		sign=pow(-1.,m);
		outC=complex<double>(0.,(-1.)*phi*(double)M);
		outC=exp(outC);//(-1^m*conjugate(...))
	}

	out=sqrtfac( (long)(l+M));
	out=sqrtfac((long)(l-M))/out;
	out*=sqrt((2.*l+1.) / (4.*PI));

	out*=plgndr( l,M,cos(theta) );

	outC*=(out*sign);

	return(outC);
}

//########################### select hardcoding

complex<double> chooseYlm(int l, int m, double theta, double phi)
{
        if(l>ylm_Table_Max) return ylm(l,m,theta,phi);
        else  return ylm_Array[l*(l+1)+m](l,m,theta,phi);

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
  num1*=(sqrt(
        (4.*PI)/( (2.*(l1+l2)+1) )
			 )
		);
  return(num1*pow(-1.,-l2+m1+m2));
}



/*#########################################################
    #######################   end   #############
#########################################################*/



