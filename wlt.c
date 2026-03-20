/* 2-D prestack forward modeling using sponge ABC using 4-th order FD
   NB: prepare high quality prestack seismic data for LSRTM and FWI
   Top boundary is free surface (no ABC applied)!
*/
/*
  Copyright (C) 2014  Xi'an Jiaotong University, UT Austin (Pengliang Yang)

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
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
#include <time.h>


#define NN 4
static int nt,fm;
static float dt;
static float *wlt;


int main(int argc, char* argv[])
{
    int it;
    float tmp, amp
    sf_file wlt_;

    /* initialize Madagascar */
    sf_init(argc,argv);

    wlt_ = sf_output("out"); /* 输出扰动模型 */
    /*< set up I/O files >*/

    sf_putint(wlt_,"n1",1);
    sf_putint(wlt_,"n2",nt);

    sf_putfloat(wlt_,"d1",dt);
    sf_putfloat(wlt_,"d2",dt);
    sf_putfloat(wlt_,"o1",0);
    sf_putfloat(wlt_,"o2",0);
    sf_putstring(wlt_,"label1","t");
    sf_putstring(wlt_,"label2","s");  /* output image with correlation imaging condition */ 



    if (!sf_getfloat("amp",&amp)) amp=1000;
    /* maximum amplitude of ricker */
    if (!sf_getfloat("fm",&fm)) fm=10;	
    /* dominant freq of ricker */
    /* thickness of sponge ABC  */
    if (!sf_getfloat("dt",&dt)) sf_error("no dt");	
    /* time interval */
    if (!sf_getint("nt",&nt))   sf_error("no nt");	
    /* total modeling time steps */
 
    wlt=sf_floatalloc(nt);
    for(it=0; it<nt; it++){
	tmp=SF_PI*fm*(it*dt-1/fm);tmp=tmp*tmp;
	wlt[it]=amp*(1.0-2.0*tmp)*expf(-tmp);
    }
   
    sf_floatwrite(wlt, nt, wlt_);

    
    free(wlt);
    

    exit(0);
}

