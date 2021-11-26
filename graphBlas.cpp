
#include <stdio.h>
#include <assert.h>
#include <GraphBLAS.h>
#include "tutorial_utils.h"

double PI = 3.141592,time1,time1max;
double RR = 8.3145;


struct i3D{
	int x,y,z;
};

int get1DIndex(int x, int y, int z, int GridSizeX, int GridSizeY, int GridSizeZ){
	return (z*GridSizeX*GridSizeY) + (y*GridSizeY) + x;
}

struct i3D get3DIndex(int idx,int GridSizeX, int GridSizeY, int GridSizeZ){
	struct i3D tmp;
	tmp.z = idx/ (GridSizeX*GridSizeY);
	idx -= tmp.z*(GridSizeX*GridSizeY);
	tmp.y = idx/(GridSizeX);
	tmp.x = idx % GridSizeX;
	return tmp;
}

//****************************************************************************
int main(int argc, char** argv)
{
	 GrB_Index const GridSize = 3;
	 int nx,ny ;
	 nx = GridSize;
	 ny = nx;
    GrB_Index const nv = GridSize*GridSize*GridSize;
    GrB_Index const no = 5;
    GrB_Index const nom = no-1;
    GrB_Matrix graph,aij,wij,tij,eij,Laplacian,ph, ph2, sph2, lap_ph;

    // Initialize a GraphBLAS context
    GrB_init(GrB_BLOCKING);

    GrB_Matrix_new(&graph, GrB_UINT64, nv, nv);
    
    
    GrB_Matrix_new(&ph, GrB_FP32,  no, nv);   // ph field
    GrB_Matrix_new(&ph2, GrB_FP32,  no, nv);  // auxilliary array
    GrB_Matrix_new(&sph2, GrB_FP32,  no, nv);
    GrB_Matrix_new(&lap_ph, GrB_FP32,  no, nv);
    GrB_Matrix_new(&aij, GrB_UINT64, no, no);  // Gradient energy coefficient
    GrB_Matrix_new(&wij, GrB_UINT64, no, no);  // Penalty term coefficient
    GrB_Matrix_new(&tij, GrB_UINT64, no, no);  // Mobility of grain boundary    
    GrB_Matrix_new(&eij, GrB_UINT64, no, no);  // Driving force of grain boundary migration   
    
    // Declaring Laplacian matrix
	 GrB_Matrix_new(&Laplacian, GrB_UINT64, nv, nv);
        
    //int delt,temp,L;g
    //float RR,vm0,gama,K1,W1,amobi,M1,E1;
    	  double L, b1, temp;
        double  M1, W1, K1, E1;



        double gamma, delta, amobi;
        double aaa, vm0;

        //****** reg data ****************************************

        //printf("delt(2.0)=  "); scanf(" %lf",&delt);
        delta = 2.0;

        temp = 1000.0;
        L = 500.0;
        b1 = L / (double)GridSize * 1.0e-9;

        vm0 = 7.0e-6;
        gamma = 0.5 * vm0 / RR / temp / b1;
        delta = 3.0;

        aaa = 2.0 / PI * sqrt(2.0 * delta * gamma);
        K1 = aaa * aaa;
        W1 = 4.0 * gamma / delta;
        amobi = 1.;
        M1 = amobi * PI * PI / (8.0 * delta);
        E1 = 500.0 / RR / temp;
        //E1=0.0;

        time1 = 0.;  time1max = 1.0e+20 + 1.0;

        for (int i = 0; i < nom; i++) {
                for (int j = 0; j < nom; j++) {
                		GrB_Matrix_setElement(wij, W1, i, j);   //wij[i][j] = W1;
                     GrB_Matrix_setElement(aij, K1, i, j);   //aij[i][j] = K1;
                     GrB_Matrix_setElement(tij, M1, i, j);   //tij[i][j] = M1;
                     GrB_Matrix_setElement(eij, 0, i, j);    //eij[i][j] = 0.0;
                        if ((i == nom-1) || (j == nom-1)) { GrB_Matrix_setElement(eij, E1, i, j);}
                        if (i > j) { GrB_Matrix_setElement(eij, -E1, i, j); }
                        if (i == j) { 
                        	GrB_Matrix_setElement(wij, 0, i, j);   
                     		GrB_Matrix_setElement(aij, 0, i, j);   
                     		GrB_Matrix_setElement(tij, 0, i, j);   
                     		GrB_Matrix_setElement(eij, 0, i, j);    
                        	//wij[i][j] = 0.0;  aij[i][j] = 0.0;  tij[i][j] = 0.0; eij[i][j] = 0.0; 
                        	}
                }
        }
    
    
    
	// building laplacian matrix
	for(int i=0;i<nv;i++){    
	// setting ph(nom,:) = 1
    for(int i=0;i<nv;i++){
    	    GrB_Matrix_setElement(ph, 1, nom-1, i);
    }
    
    
    GrB_Matrix_setElement(graph, 4, 5, 2);  // set 1,2 element to 4
	 //GrB_Matrix_setElement(g2, 4, 2, 1);
	 //graph = GrB_eWiseAdd(graph,g2);
    pretty_print_matrix_UINT64(graph, "GRAPH");
    //printf("Vande Mataram %d %d %d\n",test.x,test.y,test.z);

    GrB_Index nvals;
    GrB_Matrix_nvals(&nvals, graph);
    assert(nvals == 1);
    
    GrB_Vector sphi = NULL;
    //GrB_Descriptor d0 = NULL;
	 GrB_Vector_new(&sphi,GrB_FP32,no);   
	 GrB_Vector lap_phi;
	 GrB_Vector_new(&lap_phi,GrB_FP32,no);  
	 GrB_Vector Lp;
    
    for (time1;time1<time1max;time1++){
    	GrB_Matrix_new(&sph2, GrB_FP32,  no, nv);
    	//sph2 = sparse(no,nv))
    	
    	//lap_ph = ph*laplacian;
    	for (int i = 0;i<nv;i++){
    		GrB_extract(sphi,GrB_NULL,GrB_NULL,ph,GrB_ALL,no,i,1,GrB_NULL);
    		GrB_extract(lap_phi,GrB_NULL,GrB_NULL,lap_ph,GrB_ALL,no,i,1,GrB_NULL);
    		GxB_select(Lp,NULL,NULL,GxB_NONZERO,sphi,NULL,NULL);
			// get column of sph at vertex i
			// get col o laplacian of sph
			// ni =  len(row_active)
			// get len of row active
			
			
			int ni;
			for(int n1 = 0 ; n1<ni ; n1++){
				//int ii = row_active(n1));
				//pddtt = 0.0;
				for(int n2=0;n2<ni; n2++){
					//jj = row_active(n2);
					//sum1 = 0.5*				
				}			
			}
    	}
    
    }







    // Cleanup
    GrB_free(&graph);
    GrB_finalize();
}
