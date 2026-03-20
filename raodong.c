
#include <rsf.h>

int main(int argc, char* argv[])
{   
 
    int  nz, nx,i,j;
    float dz, dx, o1, o2;
    float **v0, **v1, *vr;

    sf_file velcon, rao, velo;/* I/O files */

    /* initialize Madagascar */
    sf_init(argc,argv);

    velo = sf_input ("in"); /* velocity */
    
    velcon = sf_input("velcon"); /* 输入背景速度*/
    rao = sf_output("out"); /* 输出扰动模型 */

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
    

    sf_putint(rao,"n1",nz);
    sf_putint(rao,"n2",nx);
    sf_putint(rao,"n3",1);
    sf_putfloat(rao,"d1",dz);
    sf_putfloat(rao,"d2",dx);
    sf_putfloat(rao,"o1",o1);
    sf_putfloat(rao,"o2",o2);
    sf_putstring(rao,"label1","Depth");
    sf_putstring(rao,"label2","Distance");


    v0=sf_floatalloc2(nz,nx);
    v1=sf_floatalloc2(nz,nx);
    vr=sf_floatalloc(nz*nx);

    sf_floatread(v0[0], nz*nx, velo);
    sf_floatread(v1[0], nz*nx, velcon);

	for(i=0;i<nx;i++)
		for(j=0;j<nz;j++)
		{
			//vr[i*nz+j]=-(1/(v0[i][j]*v0[i][j])-1/(v1[i][j]*v1[i][j]));
			vr[i*nz+j]=2*(v0[i][j]-v1[i][j])/v0[i][j];
		}
	

	sf_floatwrite(vr, nz*nx, rao);  /* output inverted image */


    free(*v0); free(v0);
	free(*v1); free(v1);
    free(vr);

    exit(0);
}

