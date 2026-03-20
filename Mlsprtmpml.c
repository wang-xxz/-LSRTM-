/* least-squares prestack RTM in 2-D */
/*
  Copyright (C) 
  - 2014  Xi'an Jiaotong University, UT Austin (Pengliang Yang)
  - 2017, Uinversity of Calgary, Daniel Trad:
  Operator modified to pass the adjoint. Added adjoint tests as well. 

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published byｚ
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
#include <rsf.h>
#include <rsfpwd.h>
#include "conjgradp.h"
#include "smoothpwd.h"
#include "prtm2d.h"
#include "prertm2d.h"
const float TOLERANCE = 1.e-12;
static int np, nxs, nr, nd;
static float *r, *sp, *sx, *sr, *gp, *gx, *gr,*xtmp,*gc,*gil,*gtmp,*gz,*sz,*gm,*x2;
static float eps, tol;
static bool verb, hasp0;
static float *illum;
float sum_alpha12(float alpha1, float alpha2,float epsil ,float *dcaltmp, float *dobs, float *derr, int nr)
/*< calculate numerator and denominator of alpha >*/
{
    int ig;
    float a, b, c;
    for(ig=0; ig<nr; ig++){
	c=derr[ig];
	a=dobs[ig]+c;/* since f(mk)-dobs[id]=derr[id], thus f(mk)=b+c; */
	b=dcaltmp[ig]-a;/* f(mk+epsil*cg)-f(mk) */
	alpha1-=b*c; 
	alpha2+=b*b;
    }
	return (alpha1*epsil/(alpha2+SF_EPS));
}
float cal_alpha(float *alpha1, float *alpha2, float epsil, int ng)
/*< calculate alpha >*/
{
    int ig;
    float a,b;

    a=b=0;
    for(ig=0; ig<ng; ig++){
	a+=alpha1[ig];
	b+=alpha2[ig];
    }

    return (a*epsil/(b+SF_EPS));
}
float cal_epsilon(float *vv, float *cg, int nx)
/*< calculate epsilcon >*/
{
    int i;
    float vvmax, cgmax;
    vvmax=cgmax=0.0;

    for(i=0; i<nx; i++){
	
	    vvmax=SF_MAX(vvmax, fabsf(vv[i]));
	    cgmax=SF_MAX(cgmax, fabsf(cg[i]));
	
    }

    return 0.01*vvmax/(cgmax+SF_EPS);
}
void grad_init(int np1     /* preconditioned size */, 
		      int nx1     /* model size */, 
		      int nd1     /* data size */, 
		      int nr1     /* residual size */, 
		      float eps1  /* scaling */,
		      float tol1  /* tolerance */, 
		      bool verb1  /* verbosity flag */, 
		      bool hasp01 /* if has initial model */) 
/*< solver constructor >*/
{
    np = np1; 
    nxs = nx1;
    nr = nr1;
    nd = nd1;
	
    eps = eps1*eps1;
    tol = tol1;
    verb = verb1;
    hasp0 = hasp01;

    r = sf_floatalloc(nr);  
    sp = sf_floatalloc(nxs);
    gp = sf_floatalloc(nxs);
    sx = sf_floatalloc(nxs);
    gx = sf_floatalloc(nxs);
    sr = sf_floatalloc(nr);
    gr = sf_floatalloc(nr);
    gtmp = sf_floatalloc(nr);
    xtmp = sf_floatalloc(nxs);
    gc = sf_floatalloc(nxs);
    gil = sf_floatalloc(nxs);
    illum = sf_floatalloc(nxs);
    gz = sf_floatalloc(nxs);
    gm = sf_floatalloc(nxs);
    sz = sf_floatalloc(nxs);
    x2 = sf_floatalloc(nxs);

}

void grad_close(void) 
/*< Free allocated space >*/
{
    free (r);
    free (sp);
    free (gp);
    free (sx);
    free (gx);
    free (sr);
    free (gr);
    free(xtmp);
    free(gc);
    free(gtmp);
    free(gil);
    free(illum);
    free(gz);
    free(sz);
}

void grad(
		 float* x          /* estimated model */, 
		 float* dat        /* data */, 
		 int niter         /* number of iterations */) 
/*< Conjugate gradient solver with shaping >*/
{
    double gn, gnp, alpha, beta, g0, dg, r0,ttt,epsil;
    float *d=NULL,alpha1,alpha2;
    int i, iter;
	for (i=0; i < nxs; i++) {
	    xtmp[i] =0.0;
	} 	

	for (i=0; i < nr; i++) {
	    r[i] =dat[i];
	}
    	illumination(nxs,nr,illum,r);
    r0 = cblas_dsdot(nr,r,1,r,1);
    if (r0 == 0.) {
	if (verb) sf_warning("zero residual: r0=%g",r0);
	return;
    }
	float xmax=0.0;
    for (iter=0; iter < niter; iter++) {
	
	
	prtm2d_lop(false,false,nxs,nr,x,gr);
	
	for (i=0; i < nr; i++) {
		 gr[i]=r[i]-gr[i];
	    }
	if (verb) sf_warning("iteration %d res: %f", iter+1,cblas_dsdot(nr,gr,1,gr,1)/r0);
	
	
	if(iter<0)
	{
		prtm2d_lop(true,false,nxs,nr,gx,gr);
		for (i=0; i < nxs; i++) {
		       gz[i]=gx[i]/(sqrtf(illum[i])+SF_EPS);
			//gz[i]=gx[i];
        	    }
		//smoothpwd(20, 2, NULL, gm,gx,false,eps);
		/*prtm2d_lop(true,false,nxs,nr,gm,gr);*/
		/*if(iter<5){
			pwsmooth_lop(false,false,nxs,nxs,gm,gx);
		}*/
	}
	else
	{
		prtm2d_lop(true,false,nxs,nr,gx,gr);
		/*for (i=0; i < nxs; i++) {
		      //gz[i]=gx[i]/sqrtf(illum[i]);
			gz[i]=gx[i];
        	    }*/
		

	}
	if (iter==0) {
	   for (i=0; i < nxs; i++){
	    gp[i]=-gx[i];
          }
	}
#if 0
	else	
	{

		for (i=0; i < nxs; i++) {
	     		 gc[i]=gx[i]-sx[i];
		
           	 }
		ttt=1+fabs(cblas_dsdot(nxs,gz,1,sx,1))/cblas_dsdot(nxs,gz,1,gx,1);
	    	sf_warning("ttt=%f",ttt);

	    	beta=cblas_dsdot(nxs,gz,1,gc,1)/(cblas_dsdot(nxs,sz,1,sx,1)*ttt);
		for (i=0; i < nxs; i++){
	        	gp[i]=-(1+beta*cblas_dsdot(nxs,gz,1,sp,1)/cblas_dsdot(nxs,gz,1,gx,1))*gx[i]+beta*sp[i];
		}
        }
#else
	
 	else {
          	for (i=0; i < nxs; i++) {
	     		 gc[i]=gx[i]-sx[i];
		
           	 }
	    beta=cblas_dsdot(nxs,gx,1,gx,1)/cblas_dsdot(nxs,sx,1,sx,1);
            
	    for (i=0; i < nxs; i++){
	        gp[i]=-gx[i]+beta*sp[i];
            }
	}
#endif
	// alpha1=alpha2=0.;
	//if(iter==0){
	  prtm2d_lop(false,false,nxs,nr,gp,gtmp);
	  alpha=cblas_dsdot(nxs,gx,1,gp,1)/cblas_dsdot(nr,gtmp,1,gtmp,1);
	/*}
	else{
	  epsil=cal_epsilon(x, gp, nxs);
	  //sf_warning("epsil=%e",epsil);
	   for(i=0; i<nxs; i++){
	    xtmp[i]=x[i]-epsil*gp[i];
	  }
	    prtm2d_lop(false,false,nxs,nr,xtmp,gtmp);
	   
	    alpha=sum_alpha12(alpha1, alpha2,epsil, gtmp, r, gr, nr);
		
	}*/
	//prtm2d_lop(false,false,nxs,nr,gp,gtmp);
	//alpha=cblas_dsdot(nxs,gx,1,gp,1)/cblas_dsdot(nr,gtmp,1,gtmp,1);
	//sf_warning("alpha=%e",alpha);
	//t2=(1+sqrtf(1+4*t1*t1))/2;


	/*if(iter==0)
	{
		for (i=0; i < nxs; i++) {
			x[i]=xtmp[i]-((t1-1)/t2)*xtmp[i]- alpha*gp[i];
	 	   }
	}
	else if(iter<30){
		for (i=0; i < nxs; i++) {
			x[i]=xtmp[i]-((t1-1)/t2)*(xtmp[i]-x2[i])- alpha*gp[i];
	 	   }	
	
	}
	else{*/
		for (i=0; i < nxs; i++) {
			x[i]=x[i]+ alpha*gp[i];
	 	   }


	//}
		
	
	
	//t1=t2;
	for (i=0; i < nxs; i++) {
		sx[i] = gx[i];
		sp[i] = gp[i];
		sz[i] = gz[i];
		//x2[i] = xtmp[i];
		//xtmp[i]=x[i];
	    }
    }
		
		
      		
	
}


int main(int argc, char* argv[])
{   
    bool verb, fromBoundary;
    int nb, nz, nx, nt, ns,nss, ng, niter, csd, sxbeg, szbeg, jsx, jsz, gxbeg, gzbeg, jgx, jgz;
    float dz, dx, dt, fm, o1, o2, amp;
    float **v0, *mod, *dat; 
    int testadj;

    sf_file shots, imglsm, imgrtm, velo;/* I/O files */

    /* initialize Madagascar */
    sf_init(argc,argv);

    shots = sf_input ("in"); /* shot records, data 	*/
    velo = sf_input ("vel"); /* velocity */
    imglsm = sf_output("out"); /* output LSRTM imglsme, model */
    imgrtm = sf_output("imgrtm"); /* output RTM imglsme */
    sf_file Fdip = sf_input("dip"); // structure dip file
    if (!sf_histint(velo,"n1",&nz)) sf_error("n1");
    /* 1st dimension size */
    if (!sf_histint(velo,"n2",&nx)) sf_error("n2");
    /* 2nd dimension size */
    if (!sf_histfloat(velo,"d1",&dz)) sf_error("d1");
    /* d1 */
    if (!sf_histfloat(velo,"d2",&dx)) sf_error("d2");
    /* d2 */
    if (!sf_histfloat(velo,"o1",&o1)) sf_error("o1");
    /* o1 */
    if (!sf_histfloat(velo,"o2",&o2)) sf_error("o2");
    /* o2 */
    if (!sf_getbool("verb",&verb)) verb=true;
    /* verbosity */
    if (!sf_getint("niter",&niter)) niter=10;
    /* totol number of least-squares iteration*/
    if (!sf_getint("nb",&nb)) nb=30;
    /* number (thickness) of ABC on each side */
    if (!sf_getbool("fromBoundary",&fromBoundary)) fromBoundary=false;
    /* if fromBoundary=true, reconstruct source wavefield from stored boundary */
    if (!sf_getint("testadj",&testadj)) testadj=0;
    /* if testadj = 1 then program only testadj without calculating */

    if (!sf_histint(shots,"n1",&nt)) sf_error("no nt");
    /* total modeling time steps */
    if (!sf_histint(shots,"n2",&ng)) sf_error("no ng");
    /* total receivers in each shot */
    if (!sf_histint(shots,"n3",&ns)) sf_error("no ns");
    /* number of shots */
    if (!sf_histfloat(shots,"d1",&dt)) sf_error("no dt");
    /* time sampling interval */
    if (!sf_histfloat(shots,"amp",&amp)) sf_error("no amp");
    /* maximum amplitude of ricker */
    if (!sf_histfloat(shots,"fm",&fm)) sf_error("no fm");
    /* dominant freq of ricker */
    if (!sf_histint(shots,"sxbeg",&sxbeg)) sf_error("no sxbeg");
    /* x-begining index of sources, starting from 0 */
    if (!sf_histint(shots,"szbeg",&szbeg)) sf_error("no szbeg");
    /* x-begining index of sources, starting from 0 */
    if (!sf_histint(shots,"gxbeg",&gxbeg)) sf_error("no gxbeg");
    /* x-begining index of receivers, starting from 0 */
    if (!sf_histint(shots,"gzbeg",&gzbeg)) sf_error("no gzbeg");
    /* x-begining index of receivers, starting from 0 */
    if (!sf_histint(shots,"jsx",&jsx)) sf_error("no jsx");
    /* source x-axis  jump interval  */
    if (!sf_histint(shots,"jsz",&jsz)) sf_error("no jsz");
    /* source z-axis jump interval  */
    if (!sf_histint(shots,"jgx",&jgx)) sf_error("no jgx");
    /* receiver x-axis jump interval  */
    if (!sf_histint(shots,"jgz",&jgz)) sf_error("no jgz");
    /* receiver z-axis jump interval  */
    if (!sf_histint(shots,"csdgather",&csd)) sf_error("csdgather or not required");
    /* default, common shot-gather; if n, record at every point*/
    if (!sf_getint("nss", &nss)) nss = ns;
     /* nss: supergather number in simulatneous-source or
           shot number(ns) in conventional acquisition */
    sf_putint(imglsm,"n1",nz);
    sf_putint(imglsm,"n2",nx);
    sf_putint(imglsm,"n3",1);
    sf_putfloat(imglsm,"d1",dz);
    sf_putfloat(imglsm,"d2",dx);
    sf_putfloat(imglsm,"o1",o1);
    sf_putfloat(imglsm,"o2",o2);
    sf_putstring(imglsm,"label1","Depth");
    sf_putstring(imglsm,"label2","Distance");

    sf_putint(imgrtm,"n1",nz);
    sf_putint(imgrtm,"n2",nx);
    sf_putint(imgrtm,"n3",1);
    sf_putfloat(imgrtm,"d1",dz);
    sf_putfloat(imgrtm,"d2",dx);
    sf_putfloat(imgrtm,"o1",o1);
    sf_putfloat(imgrtm,"o2",o2);
    sf_putstring(imgrtm,"label1","Depth");
    sf_putstring(imgrtm,"label2","Distance");

    /* In rtm, vv is the velocity model [modl], which is input parameter; 
       mod is the imglsme/reflectivity [imglsm]; dat is seismogram [data]! */
    v0=sf_floatalloc2(nz,nx);
    mod=sf_floatalloc(nz*nx);
    dat=sf_floatalloc(nt*ng*ns);
float **slope = sf_floatalloc2(nz, nx);
float *pp = sf_floatalloc(nx*nz);
sf_floatread(slope[0], nx*nz, Fdip);
    /* initialize velocity, model and data */
    sf_floatread(v0[0], nz*nx, velo);
    memset(mod, 0, nz*nx*sizeof(float));
    sf_floatread(dat, nt*ng*ns, shots);

    pt2d *src2d = pt2dalloc1(ns);
    pt2d *rec2d = pt2dalloc1(ng);
	int is,ig;
   for(is=0; is<ns; is++)
    {
	src2d[is].z=szbeg+is*jsz;
	src2d[is].x=sxbeg+is*jsx;
    }
    for(ig=0; ig<ng; ig++)
    {
	rec2d[ig].z=gzbeg+ig*jgz;
	rec2d[ig].x=gxbeg+ig*jgx;
    }

    prtm2d_init(verb, csd, fromBoundary, dz, dx, dt, amp, fm, nz, nx, nb, nt, ns, ng, 
		sxbeg, szbeg, jsx, jsz, gxbeg, gzbeg, jgx, jgz, v0, mod, dat);
    //prertm2d_init(verb, nx, nz, nb, nt, ns, nss, ng, dx, dz, dt,amp,fm ,v0, src2d, rec2d);
    // run adjoint test first.
    /*sf_warning("adjoint test \n");
    if (testadj){
	prtm2d_adjtest();
	sf_warning("exiting after testadj\n");
	mod=NULL;
	//sf_floatwrite(mod,nz*nx,imgrtm);
	//sf_floatwrite(mod,nz*nx,imglsm);
	prtm2d_close();
	free(*v0); free(v0);
	free(mod);
	free(dat); 	  
	exit (0);
    }*/
	
    //sf_warning("start migration\n");
    /* original RTM is simply applying adjoint of prtm2d_lop once!*/
    //prtm2d_lop(true, false, nz*nx, nt*ng*ns, mod, dat); 
    //sf_floatwrite(mod, nz*nx, imgrtm);/* output RTM image */
    //sf_warning("migration ends");

    /*if (niter>0){
	memset(mod, 0, nz*nx*sizeof(float));//very important

	sf_warning("inside LS loop niter= %d\n",niter);
	sf_solver(prtm2d_lop, sf_cgstep, nz*nx, nt*ng*ns, mod, dat, niter, "verb", verb, "end");
	sf_floatwrite(mod, nz*nx, imglsm); 
	sf_cgstep_close();

    }*/
	
	pwsmooth_init(10, nz, nx, 2, 0.01);
        pwsmooth_set(slope);
	if (niter>0){
	memset(mod, 0, nz*nx*sizeof(float));//very important

	sf_warning("inside LS loop niter= %d\n",niter);
	grad_init(nz*nx, nz*nx, nt*ng*ns,nt*ng*ns, eps, 1.e-6, true, false);
    	grad(mod, dat, niter);
	sf_floatwrite(mod, nz*nx, imglsm); 
    	grad_close();

    }
	//pwsmooth_init(3, nz, nx, 3, 0.01);
    	//pwsmooth_set(slope);
	/*if (niter>0){
	memset(mod, 0, nz*nx*sizeof(float));//very important

	sf_warning("inside LS loop niter= %d\n",niter);
        sf_conjgradp_init(nx*nz, nx*nz, nt*ng*ns, nt*ng*ns, eps, 1.e-6, true, false);
        sf_conjgradp(NULL, prtm2d_lop, pwsmooth_lop, pp, mod, dat, niter);
	sf_floatwrite(mod, nz*nx, imglsm); 
        sf_conjgradp_close();


    }*/

    pwsmooth_close();
    //pwsmooth_close();
    //prertm2d_close();
    prtm2d_close();
    free(*v0); free(v0);
    free(mod);
    free(dat); 

    exit(0);
}

