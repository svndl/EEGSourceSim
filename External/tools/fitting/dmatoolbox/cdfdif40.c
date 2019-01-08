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
               Uses 40 quadrature points for drift and 40 for the others.
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
    double nr_nu = 40, nr_z = 40;
    
    double gk[40]={-8.0987611392508505048,-7.4115825314854690831,-6.8402373052493556926,-6.3282553512200818702,-5.8540950560303999239,-5.4066542479701276136,-4.9792609785452555116,-4.5675020728443946894,-4.1682570668325000796,-3.7792067534352233871,-3.3985582658596285022,-3.0248798839012844830,-2.6569959984428956901,-2.2939171418750832210,-1.9347914722822958655,-1.5788698949316137821,-1.2254801090462890123,-.87400661235708809738,-.52387471383227712796,-.17453721459758239631,.17453721459758239631,.52387471383227712796,.87400661235708809738,1.2254801090462890123,1.5788698949316137821,1.9347914722822958655,2.2939171418750832210,2.6569959984428956901,3.0248798839012844830,3.3985582658596285022,3.7792067534352233871,4.1682570668325000796,4.5675020728443946894,4.9792609785452555116,5.4066542479701276136,5.8540950560303999239,6.3282553512200818702,6.8402373052493556926,7.4115825314854690831,8.0987611392508505048},
           w_gh[40]={.25910437138470341194e-28,.85440569637754419606e-24,.25675933654112184211e-20,.19891810121165100450e-17,.60083587894908528971e-15,.88057076452083577402e-13,.71565280526903646397e-11,.35256207913654222548e-9,.11212360832275848610e-7,.24111441636019625888e-6,.36315761506840953354e-5,.39369339810916718697e-4,.31385359454132611551e-3,.18714968295979511528e-2,.84608880082581404414e-2,.29312565536172417724e-1,.78474605865404459260e-1,.16337873271327166269,.26572825187708215555,.33864327742558919532,.33864327742558919532,.26572825187708215555,.16337873271327166269,.78474605865404459260e-1,.29312565536172417724e-1,.84608880082581404414e-2,.18714968295979511528e-2,.31385359454132611551e-3,.39369339810916718697e-4,.36315761506840953354e-5,.24111441636019625888e-6,.11212360832275848610e-7,.35256207913654222548e-9,.71565280526903646397e-11,.88057076452083577402e-13,.60083587894908528971e-15,.19891810121165100450e-17,.25675933654112184211e-20,.85440569637754419606e-24,.25910437138470341194e-28},
           gz[40]={-.99823007668897911770,-.99074456319832715501,-.97724421201887445854,-.95791741396179397316,-.93282786927060368232,-.90207721750936598060,-.86597830737194936290,-.82460009839214354344,-.77831187481202868117,-.72731565096148653726,-.67195756216544155759,-.61255367214190326042,-.54946715082724151280,-.48307580931786880951,-.41377919892381465061,-.34199409250705214980,-.26815218467656593004,-.19269758074344939258,-.11608407067285775316,-.38772417506061480907e-1,.38772417506061751524e-1,.11608407067286004299,.19269758074297566042,.26815218467978524375,.34199409252852924768,.41377919860509299044,.48307581089111756301,.54946714666631024659,.61255367783130543202,.67195756257359817720,.72731562526758453124,.77831199205257750595,.82459962105729867066,.86597995148433870582,.90207281795888860643,.93283689405635383807,.95790339574755600793,.97726019847530665174,.99073231193021737973,.99823447294005429598},
           w_g[40]={.45428882756907047746e-2,.10453455320676924636e-1,.16482989657810115380e-1,.22164991930851658114e-1,.28017827927605223892e-1,.33423196334136061336e-1,.38786353334588977160e-1,.43864985993934020592e-1,.48695136389870545546e-1,.53227948696192871336e-1,.57439638454765934439e-1,.61306281142436344633e-1,.64804008321460299102e-1,.67912045499704906670e-1,.70611647767459290170e-1,.72886582299559465881e-1,.74723169072356279696e-1,.76110361899336689828e-1,.77039818164291340441e-1,.77505947978424735711e-1,.77505947978424735711e-1,.77039818164291395952e-1,.76110361899354522786e-1,.74723169072191286677e-1,.72886582298574711936e-1,.70611647782002948115e-1,.67912045343886007220e-1,.64804009416781313546e-1,.61306272961870023064e-1,.57439688629177282353e-1,.53228532984586046650e-1,.48695390131306269532e-1,.43865719785909869366e-1,.38813149008100872317e-1,.33454189176501537839e-1,.27958786172135020875e-1,.22208422002312246113e-1,.16506448409281091610e-1,.10467435840146142520e-1,.45336236211998759069e-2};
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
