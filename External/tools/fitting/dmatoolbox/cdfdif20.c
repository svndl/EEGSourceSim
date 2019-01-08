/*
  Name: cdfdif.c
  Copyright: Academic use
  Author: Joachim Vandekerckhove
          Department of Psychology
          Katholieke Universiteit Leuven
          Belgium
          email: joachim.vandekerckhove@psy.kuleuven.be
  Date: The generating algorithm was written on 09/01/06
  Description: Computes Cumulative Distribution Function for the Diffusion
               model with random trial to trial mean drift (normal), starting
               point and non-decision (ter) time (both rectangular).
               Uses 20 quadrature points for drift and 20 for the others.
               Based on methods described in:
                     Tuerlinckx, F. (2004). The efficient computation of the
                     cumulative distribution and probability density functions
                     in the diffusion model, Behavior Research Methods,
                     Instruments, & Computers, 36 (4), 702-716.
 
 */

#include <stdlib.h>
#include <math.h>
#include <mex.h>
#include <matrix.h>
#include <tmwtypes.h>

#define PI 3.1415926535897932384626433832795028841971693993751


double cdfdif(double t, int x, double *par, double *prob);
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
double absol(double a){ return a>0 ? a : -a; }

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
    double *p, *x, *t, *y, *pr, nonzero = 1e-10;
    int nt,i;
    
    t = mxGetPr(prhs[0]);
    x = mxGetPr(prhs[1]);
    p = mxGetPr(prhs[2]);
    nt = mxGetN(prhs[0]);
    
    plhs[0] = mxCreateDoubleMatrix(1,nt, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1,1, mxREAL);
    y = mxGetPr(plhs[0]);
    pr = mxGetPr(plhs[1]);
    
    if (p[2]<nonzero) { p[2] = nonzero; }
    if (p[4]<nonzero) { p[4] = nonzero; }
    if (p[5]<nonzero) { p[5] = nonzero; }
    for (i=0;i<nt;i++) y[i] = cdfdif(t[i],*x,p,pr);
}

/* The main function for the seven-parameter diffusion model */
double cdfdif(double t, int x, double *par, double *prob)
{
    double a = par[0], Ter = par[1], eta = par[2], z = par[3], sZ = par[4],
    st = par[5], nu = par[6], a2 = a*a,
    Z_U = (1-x)*z+x*(a-z)+sZ/2, /* upper boundary of z distribution */
    Z_L = (1-x)*z+x*(a-z)-sZ/2, /* lower boundary of z distribution */
    lower_t = Ter-st/2, /* lower boundary of Ter distribution */
    upper_t, /* upper boundary of Ter distribution */
    delta = 1e-29, /* convergence values for terms added to partial sum */
    epsilon = 1e-7, /* value to check deviation from zero */
    min_RT=0.001; /* realistic minimum rt to complete decision process */
    
    int v_max = 5000; /* maximum number of terms in a partial sum */
    /* approximating infinite series */
    
    double Fnew, sum_z=0, sum_nu=0, p1, p0, sum_hist[3]={0,0,0},
    denom, sifa, upp, low, fact, exdif, su, sl, zzz, ser;
    double nr_nu = 20, nr_z = 20;
    
    double gk[20]={-5.3874808900112327592,-4.6036824495507442379,-3.9447640401156252032,-3.3478545673832162954,-2.7888060584281304521,-2.2549740020892756753,-1.7385377121165861425,-1.2340762153953230840,-.73747372854539439135,-.24534070830090126680,.24534070830090126680,.73747372854539439135,1.2340762153953230840,1.7385377121165861425,2.2549740020892756753,2.7888060584281304521,3.3478545673832162954,3.9447640401156252032,4.6036824495507442379,5.3874808900112327592},
           w_gh[20]={.22293936455341646528e-12,.43993409922731809293e-9,.10860693707692821796e-6,.78025564785320632605e-5,.22833863601635307687e-3,.32437733422378566862e-2,.24810520887340880430e-1,.10901720602002162863,.28667550536283409324,.46224366960061014087,.46224366960061014087,.28667550536283409324,.10901720602002162863,.24810520887340880430e-1,.32437733422378566862e-2,.22833863601635307687e-3,.78025564785320632605e-5,.10860693707692821796e-6,.43993409922731809293e-9,.22293936455341646528e-12},
           gz[20]={-.99312859919393192687,-.96397192725643798816,-.91223442827299427993,-.83911697180954192277,-.74633190646438019034,-.63605368072610402042,-.51086700195056844453,-.37370608871551552754,-.22778585114163685255,-.76526521133497449334e-1,.76526521133497338312e-1,.22778585114163671377,.37370608871551153074,.51086700195060519292,.63605368072587842310,.74633190646520863876,.83911697180775823846,.91223442827524447996,.96397192725477587327,.99312859919448925883},
           w_g[20]={.17614007115893910022e-1,.40601429824919738065e-1,.62672048327592072559e-1,.83276741580984969815e-1,.10193011981663034626,.11819453196160842334,.13168863844919270756,.14209610931837018954,.14917298647260451849,.15275338713072578178,.15275338713072586505,.14917298647260440747,.14209610931836397230,.13168863844922096273,.11819453196171304798,.10193011981627890516,.83276741581628954680e-1,.62672048318034509484e-1,.40601429821751071347e-1,.17614007115704332501e-1};
    int i,m,v;
    
    for(i=0; i<nr_nu; i++)
    {
        gk[i] = 1.41421356237309505*gk[i]*eta+nu;
        /* appropriate scaling of GH quadrature nodes based on normal kernel */
        w_gh[i] = w_gh[i]/1.772453850905515882;
        /* appropriate weighting of GH quadrature weights based on normal kernel */
    }
    for(i=0; i<nr_z; i++)
        gz[i] = (.5*sZ*gz[i])+z;
    /* appropriate scaling of Gaussian quadrature nodes */
    
    /* numerical integration with respect to z_0 */
    for (i=0; i<nr_z; i++)
    {
        sum_nu=0;
        /* numerical integration with respect to xi */
        for (m=0; m<nr_nu; m++)
        {
            if (absol(gk[m])>epsilon) /* test for gk(m) being close to zero */
            {
                sum_nu+=(exp(-200*gz[i]*gk[m])-1)/(exp(-200*a*gk[m])-1)*w_gh[m];
            }
            else
            {
                sum_nu+=gz[i]/a*w_gh[m];
            }
        }
        sum_z+=sum_nu*w_g[i]/2;
    }
    *prob=sum_z;
    /* end of prob */
    
    if (t-Ter+st/2>min_RT)/* is t larger than lower boundary Ter distribution? */
    {
        upper_t = t<Ter+st/2 ? t : Ter+st/2;
        p1=*prob*(upper_t-lower_t)/st; /* integrate probability with respect to t */
        p0=(1-*prob)*(upper_t-lower_t)/st;
        if (t>Ter+st/2) /* is t larger than upper boundary Ter distribution? */
        {
            sum_hist[0] = 0;
            sum_hist[1] = 0;
            sum_hist[2] = 0;
            for (v=0;v<v_max;v++) /* infinite series */
            {
                sum_hist[0]=sum_hist[1];
                sum_hist[1]=sum_hist[2];
                sum_nu=0;
                sifa=PI*v/a;
                for (m=0;m<nr_nu;m++) /* numerical integration with respect to xi */
                {
                    denom=(100*gk[m]*gk[m]+(PI*PI)*(v*v)/(100*a2));
                    upp=exp((2*x-1)*Z_U*gk[m]*100-3*log(denom)+log(w_gh[m])-2*log(100));
                    low=exp((2*x-1)*Z_L*gk[m]*100-3*log(denom)+log(w_gh[m])-2*log(100));
                    fact=upp*((2*x-1)*gk[m]*sin(sifa*Z_U)*100-sifa*cos(sifa*Z_U))-
                    low*((2*x-1)*gk[m]*sin(sifa*Z_L)*100-sifa*cos(sifa*Z_L));
                    exdif=exp((-.5*denom*(t-upper_t))+
                    log(1-exp(-.5*denom*(upper_t-lower_t))));
                    sum_nu+=fact*exdif;
                }
                sum_hist[2]=sum_hist[1]+v*sum_nu;
                if ((absol(sum_hist[0]-sum_hist[1])<delta) &&
                (absol(sum_hist[1]-sum_hist[2])<delta) && (sum_hist[2]>0))
                    break;
            }
            Fnew=(p0*(1-x)+p1*x)-sum_hist[2]*4*PI/(a2*sZ*st);
            /* cumulative distribution function for t and x */
        }
        else if (t<=Ter+st/2) /* is t lower than upper boundary Ter distribution? */
        {
            sum_nu=0;
            for (m=0;m<nr_nu;m++)
            {
                if (absol(gk[m])>epsilon)
                {
                    sum_z=0;
                    for (i=0;i<nr_z;i++)
                    {
                        zzz=(a-gz[i])*x+gz[i]*(1-x);
                        ser=-((a*a2)/((1-2*x)*gk[m]*PI*.01))*sinh(zzz*(1-2*x)*gk[m]/.01)/
                        (sinh((1-2*x)*gk[m]*a/.01)*sinh((1-2*x)*gk[m]*a/.01))
                        +(zzz*a2)/((1-2*x)*gk[m]*PI*.01)*cosh((a-zzz)*(1-2*x)
                        *gk[m]/.01)/sinh((1-2*x)*gk[m]*a/.01);
                        sum_hist[0] = 0;
                        sum_hist[1] = 0;
                        sum_hist[2] = 0;
                        for (v=0;v<v_max;v++)
                        {
                            sum_hist[0]=sum_hist[1];
                            sum_hist[1]=sum_hist[2];
                            sifa=PI*v/a;
                            denom=(gk[m]*gk[m]*100+(PI*v)*(PI*v)/(a2*100));
                            sum_hist[2]=sum_hist[1]+v*sin(sifa*zzz)*
                            exp(-.5*denom*(t-lower_t)-2*log(denom));
                            if ((absol(sum_hist[0]-sum_hist[1])<delta) &&
                            (absol(sum_hist[1]-sum_hist[2])<delta) && (sum_hist[2]>0))
                                break;
                        }
                        sum_z+=.5*w_g[i]*(ser-4*sum_hist[2])*(PI/100)/
                        (a2*st)*exp((2*x-1)*zzz*gk[m]*100);
                    }
                }
                else
                {
                    sum_hist[0] = 0;
                    sum_hist[1] = 0;
                    sum_hist[2] = 0;
                    su=-(Z_U*Z_U)/(12*a2)+(Z_U*Z_U*Z_U)/
                    (12*a*a2)-(Z_U*Z_U*Z_U*Z_U)/(48*a2*a2);
                    sl=-(Z_L*Z_L)/(12*a2)+(Z_L*Z_L*Z_L)/
                    (12*a*a2)-(Z_L*Z_L*Z_L*Z_L)/(48*a2*a2);
                    for (v=1;v<v_max;v++)
                    {
                        sum_hist[0]=sum_hist[1];
                        sum_hist[1]=sum_hist[2];
                        sifa=PI*v/a;
                        denom=(PI*v)*(PI*v)/(a2*100);
                        sum_hist[2]=sum_hist[1]+1/(PI*PI*PI*PI*v*v*v*v)*(cos(sifa*Z_L)-
                        cos(sifa*Z_U))*exp(-.5*denom*(t-lower_t));
                        if ((absol(sum_hist[0]-sum_hist[1])<delta) &&
                        (absol(sum_hist[1]-sum_hist[2])<delta) && (sum_hist[2]>0))
                            break;
                    }
                    sum_z=400*a2*a*(sl-su-sum_hist[2])/(st*sZ);
                }
                sum_nu+=sum_z*w_gh[m];
            }
            Fnew=(p0*(1-x)+p1*x)-sum_nu;
        }
    }
    else if (t-Ter+st/2<=min_RT) /* is t lower than lower boundary Ter distr? */
    {
        Fnew=0;
    }
    
    Fnew = Fnew>delta ? Fnew : 0;
    
    return Fnew;
}
