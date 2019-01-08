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
               Uses 30 quadrature points for drift and 30 for the others.
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
    double nr_nu = 30, nr_z = 30;
    
    double gk[30]={-6.8633452935298917552,-6.1382792201239348984,-5.5331471515674959250,-4.9889189685899442139,-4.4830553570925184559,-4.0039086038612285989,-3.5444438731553500332,-3.0999705295864417032,-2.6671321245356174323,-2.2433914677615041100,-1.8267411436036880001,-1.4155278001981885794,-1.0083382710467234666,-.60392105862555234275,-.20112857654887147940,.20112857654887147940,.60392105862555234275,1.0083382710467234666,1.4155278001981885794,1.8267411436036880001,2.2433914677615041100,2.6671321245356174323,3.0999705295864417032,3.5444438731553500332,4.0039086038612285989,4.4830553570925184559,4.9889189685899442139,5.5331471515674959250,6.1382792201239348984,6.8633452935298917552},
           w_gh[30]={.29082547001312413853e-20,.28103336027509134270e-16,.28786070805487111582e-13,.81061862974186541567e-11,.91785804241723095244e-9,.51085224507759163527e-7,.15790948873247266304e-5,.29387252289230019294e-4,.34831012431687488661e-3,.27379224730672033043e-2,.14703829704826700334e-1,.55144176870234540289e-1,.14673584754089091797,.28013093083921419835,.38639488954179251889,.38639488954179251889,.28013093083921419835,.14673584754089091797,.55144176870234540289e-1,.14703829704826700334e-1,.27379224730672033043e-2,.34831012431687488661e-3,.29387252289230019294e-4,.15790948873247266304e-5,.51085224507759163527e-7,.91785804241723095244e-9,.81061862974186541567e-11,.28786070805487111582e-13,.28103336027509134270e-16,.29082547001312413853e-20},
           gz[30]={-.99689347569913733249,-.98366814687687087471,-.96002183275848396171,-.92620007922136493583,-.88256051021502268661,-.82956577967703815091,-.76777742238051593926,-.69785049925031195084,-.62052618134050430143,-.53662414864355034716,-.44703376940617506330,-.35270472556011439602,-.25463692616386679468,-.15386991360870561074,-.51471842555319571866e-1,.51471842555319474721e-1,.15386991360870905243,.25463692616367417099,.35270472556286180943,.44703376938919159311,.53662414870169872216,.62052618121801539353,.69785049941597576684,.76777742219215738828,.82956578014466575421,.88256050862388146783,.92620008294608824340,.96002182717677508883,.98366815193013146246,.99689347372784720136},
           w_g[30]={.79682123113510149731e-2,.18466444117175715572e-1,.28784731884966344162e-1,.38799170621555065241e-1,.48402693019156332876e-1,.57493147731047591908e-1,.65974233309089441724e-1,.73755973060031954081e-1,.80755895614425199369e-1,.86899787130305980454e-1,.92122522251870925247e-1,.96368737172339208330e-1,.99593420587013981038e-1,.10176238974840204343,.10285265289355881302,.10285265289355881302,.10176238974840158547,.99593420587027581270e-1,.96368737172055199403e-1,.92122522252858718428e-1,.86899787119621638154e-1,.80755895576700001404e-1,.73755973129618832007e-1,.65974233055547407134e-1,.57493142718566606075e-1,.48402690665996082886e-1,.38799184808631557997e-1,.28784732358834302923e-1,.18466444868141827784e-1,.79682185831451077251e-2};
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
