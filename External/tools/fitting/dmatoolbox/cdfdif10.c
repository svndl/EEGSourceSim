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
               Uses 10 quadrature points for drift and 10 for the others.
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
    double nr_nu = 10, nr_z = 10;
    
    double gk[10]={-3.4361591188377378359,-2.5327316742327896648,-1.7566836492998818553,-1.0366108297895135770,-.34290132722370464391,.34290132722370464391,1.0366108297895135770,1.7566836492998818553,2.5327316742327896648,3.4361591188377378359},
           w_gh[10]={.76404328552326426869e-5,.13436457467812345894e-2,.33874394455481036947e-1,.24013861108231340791,.61086263373532512233,.61086263373532512233,.24013861108231340791,.33874394455481036947e-1,.13436457467812345894e-2,.76404328552326426869e-5},
           gz[10]={-.97390652851716952298,-.86506336668898531350,-.67940956829902787728,-.43339539412924654727,-.14887433898163121571,.14887433898163096591,.43339539412924582562,.67940956829902887648,.86506336668898320408,.97390652851716963401},
           w_g[10]={.66671344308693342162e-1,.14945134915057770031,.21908636251598398448,.26926671930999629412,.29552422471475281451,.29552422471475292554,.26926671930999596105,.21908636251597693456,.14945134915057806113,.66671344308693036851e-1};
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
                if (v==v_max) mexPrintf("V_MAX reached!\n");
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
