
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <GraphBLAS.h>
#include "tutorial_utils.h"

double PI = 3.141592;
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

int randint(int n) {
  int r = rand() % 20;
  return r;
}


void dataVtk(int fNum,int ndxm, int ndym, int ndzm, GrB_Vector field )
{
	FILE* stream;

	int 		i, j, k, kk;

	double 	col;
	
	char str[5];
	sprintf(str, "%d", fNum);
	char fname[120];
	sprintf(fname,"/home/dkumar/GraphBLAS/src/vtkFilesGrid30/file%s.vtk",str);
   	printf("crashed here file");

	stream = fopen(fname, "a");
	//fprintf(stream, "%e  \n", time1);
	//int pts = ndxm * ndym * ndzm;
	int pts = ndxm * ndym * ndzm;
	fprintf(stream, "# vtk DataFile Version 2.0\n");
	fprintf(stream, "Mesh Field\n");
	fprintf(stream, "ASCII\n");
	fprintf(stream, "DATASET STRUCTURED_GRID\n");
	fprintf(stream, "DIMENSIONS %d %d %d\n", ndxm, ndym, ndzm);
	fprintf(stream, "POINTS %d int\n", pts);

	for (i = 0; i < ndxm; i++) {
		for (j = 0; j < ndym; j++) {
			for (k = 0; k < ndzm; k++) {
				fprintf(stream, "%d %d %d\n", i, j, k);
				/*
				fprintf(stream, "%d  \n", n00[i][j][k]);
				for (kk = 1; kk <= n00[i][j][k]; kk++) {
					fprintf(stream, "%d , %e  ", qh[kk][i][j][k], ph[kk][i][j][k]);
				}*/
				//fprintf(stream, "\n");



			}
		}
	}
	printf("Here before acessing data");
	fprintf(stream, "POINT_DATA %d\n", pts);
	fprintf(stream, "FIELD Field_Data 1\n");
	fprintf(stream, "field 1 %d double\n", pts);
	int id=0;
	for (i = 0; i < ndxm; i++) {
		for (j = 0; j < ndym; j++) {
			for (k = 0; k < ndzm; k++) {
				float val = 0;
				GrB_Vector_extractElement_FP32(&val,field,id);
				fprintf(stream, "%f\n", val);
				id += 1;
			}
		}
	}
	//fprintf(stream, "---------------------------------------------------------------------------\n");
	fclose(stream);
	printf("\n End Printing file");
}



//****************************************************************************
int main(int argc, char** argv)
{
	  printf("Beginning");
	  
    srand(time(NULL));  
	 GrB_Index const GridSize = 30;
	 GrB_Index const GridSizeSq = GridSize*GridSize;
	 int nx,ny ;
	 nx = GridSize;
	 ny = nx;
    GrB_Index const nv = GridSize*GridSize*GridSize;
    GrB_Index const no = 6;
    GrB_Index const nom = no-1;
    int time1,time1max;
    time1 = 0.;  time1max = 100;
    GrB_Matrix graph,aij,wij,tij,eij,ph, ph2, sph2, lap_ph, testMult, test2, finalPh, divideMatrix;

    // Initialize a GraphBLAS context
    GrB_init(GrB_BLOCKING);

    GrB_Matrix_new(&graph, GrB_UINT64, nv, nv);
    
    
    GrB_Matrix_new(&ph, GrB_FP32,  no, nv);   // ph field
    GrB_Matrix_new(&ph2, GrB_FP32,  no, nv);  // auxilliary array
    GrB_Matrix_new(&sph2, GrB_FP32,  no, nv);
    GrB_Matrix_new(&lap_ph, GrB_FP32,  no, nv);
    GrB_Matrix_new(&aij, GrB_FP32, no, no);  // Gradient energy coefficient
    GrB_Matrix_new(&wij, GrB_FP32, no, no);  // Penalty term coefficient
    GrB_Matrix_new(&tij, GrB_FP32, no, no);  // Mobility of grain boundary    
    GrB_Matrix_new(&eij, GrB_FP32, no, no);  // Driving force of grain boundary migration   
    GrB_Matrix_new(&testMult, GrB_FP32, no, nv); 
    GrB_Matrix_new(&test2, GrB_FP32, nv, no); 
	 GrB_Matrix_new(&finalPh, GrB_FP32, nv, no); 
	 GrB_Matrix_new(&divideMatrix, GrB_FP32,  no, nv); 
	// lets create seeds and set them
	  printf("Beginning");

  printf("Setup done initialized matrix");
  
  
  
  
  // lets create a 2D Kronecker first
  GrB_Matrix diagonal, lap1D, lap2Da,lap2Db, lap3D, lap3Da, lap3Db;
  // create two matrices
  
  GrB_Matrix_new(&diagonal, GrB_FP32, GridSize, GridSize);
  GrB_Matrix_new(&lap1D, GrB_FP32, GridSize, GridSize);  
  GrB_Matrix_new(&lap2Da, GrB_FP32, GridSizeSq, GridSizeSq);  

  GrB_Matrix_new(&lap2Db, GrB_FP32, GridSizeSq, GridSizeSq);  
  GrB_Matrix_new(&lap3D, GrB_FP32, nv, nv);  
  GrB_Matrix_new(&lap3Da, GrB_FP32, nv, nv);  
  GrB_Matrix_new(&lap3Db, GrB_FP32, nv, nv);  
  
  for (int i=0;i<GridSize; i++){
		GrB_Matrix_setElement(diagonal, 1, i, i);		
		GrB_Matrix_setElement(lap1D, -2, i, i);  
		GrB_Matrix_setElement(lap1D, 1, i, i+1);  
		GrB_Matrix_setElement(lap1D, 1, i, i-1);  
  }
  
  //pretty_print_matrix_FP64(diagonal,"Diagonal matrix");
  //pretty_print_matrix_FP64(lap1D,"Laplacian 1D matrix");
  
  // Lets do the kronecker now
  
  
  GxB_kron(lap2Da, GrB_NULL, GrB_NULL,GrB_TIMES_FP64 ,diagonal, lap1D, GrB_NULL);
  //pretty_print_matrix_FP64(lap2Da,"Laplacian 2D a matrix");
  
  GxB_kron(lap2Db, GrB_NULL, GrB_NULL,GrB_TIMES_FP64 ,lap1D, diagonal, GrB_NULL);
  //pretty_print_matrix_FP64(lap2Db,"Laplacian 2D b matrix");
  // time to do the 3D kronecker
  GxB_kron(lap3Da, GrB_NULL, GrB_NULL, GrB_TIMES_FP64, diagonal, lap2Db, GrB_NULL);
  GxB_kron(lap3Db, GrB_NULL, GrB_NULL, GrB_TIMES_FP64, diagonal, lap2Da, GrB_NULL);
  GxB_kron(lap3D, GrB_NULL, GrB_NULL, GrB_TIMES_FP64, lap2Db, diagonal, GrB_NULL);
    
  GrB_eWiseAdd_Matrix_BinaryOp(lap3D, GrB_NULL, GrB_NULL, GrB_PLUS_FP32, lap3D, lap3Da, GrB_NULL);
  GrB_eWiseAdd_Matrix_BinaryOp(lap3D, GrB_NULL, GrB_NULL, GrB_PLUS_FP32, lap3D, lap3Db, GrB_NULL);
  //pretty_print_matrix_FP64(lap3D,"Laplacian 3D matrix");
  // Laplacian is correct now and stored in lap3D

  GrB_free(&lap3Da);
  GrB_free(&lap3Db);
  GrB_free(&lap2Da);
  GrB_free(&lap2Db);



        
    //int delt,temp,L;g
    //float RR,vm0,gama,K1,W1,amobi,M1,E1;
    	  double L, b1, temp;
        double  M1, W1, K1, E1;



        float gamma, delta, amobi;
        float aaa, vm0;

        //****** reg data ****************************************

        //printf("delt(2.0)=  "); scanf(" %lf",&delt);
        delta = 5.0;

        temp = 2000.0;
        L = 2000.0;
        b1 = L / (float)100 * 1.0e-9;

        vm0 = 7.0e-6;
        gamma = 0.5 * vm0 / RR / temp / b1;
        delta = 5.0;

        aaa = 2.0 / PI * sqrt(2.0 * delta * gamma);
        K1 = aaa * aaa;
        W1 = 4.0 * gamma / delta;
        amobi = 1.;
        M1 = amobi * PI * PI / (8.0 * delta);
        E1 = 500.0 / RR / temp;
        //E1=0.0;
			


        for (int i = 0; i < nom; i++) {
                for (int j = 0; j < nom; j++) {

                		GrB_Matrix_setElement(wij, W1, i, j);   //wij[i][j] = W1;
                     GrB_Matrix_setElement(aij, K1, i, j);   //aij[i][j] = K1;
							if (i < nom-1){
								GrB_Matrix_setElement(tij, M1, i, nom-1);
								GrB_Matrix_setElement(eij, E1, i, nom-1);
								GrB_Matrix_setElement(tij, M1, nom-1, i);		
								GrB_Matrix_setElement(eij,-E1, nom-1, i);					
							}                     
                     
                     //GrB_Matrix_setElement(eij, 0, i, j);    //eij[i][j] = 0.0;
                        
                        if (i == j) { 
                        	GrB_Matrix_setElement(wij, 0, i, j);   
                     		GrB_Matrix_setElement(aij, 0, i, j);   
                     		GrB_Matrix_setElement(tij, 0, i, j);   
                     		GrB_Matrix_setElement(eij, 0, i, j);    
                        	//wij[i][j] = 0.0;  aij[i][j] = 0.0;  tij[i][j] = 0.0; eij[i][j] = 0.0; 
                        	}
                }
                GrB_Matrix_setElement(tij, 0, i, nom);
								GrB_Matrix_setElement(eij, 0, i, nom);
								GrB_Matrix_setElement(tij, 0, i, nom);		
								GrB_Matrix_setElement(wij, 0, i, nom);	
								GrB_Matrix_setElement(aij, 0, i, nom);	
								GrB_Matrix_setElement(eij, 0, nom, i);
								GrB_Matrix_setElement(tij, 0, nom, i);		
								GrB_Matrix_setElement(wij, 0, nom, i);	
								GrB_Matrix_setElement(aij, 0, nom, i);	
        }
    printf("no is %d nom is %d\n",no,nom);
    printf("Values of W1 %f K1 %f  M1 %f E1 %f \n",W1,K1,M1,E1);
    printf("Printing wij \n");
    pretty_print_matrix_FP64(wij,"Wij Values");
    printf("Printing aij \n");
    pretty_print_matrix_FP64(aij,"Aij Values");
    printf("Printing eij \n");
    pretty_print_matrix_FP64(eij,"Eij Values");
    printf("Printing tij \n");
    pretty_print_matrix_FP64(tij,"Tij Values");
    


    
    
	// building laplacian matrix
	/*
	GrB_Matrix_new(&Laplacian, GrB_FP32, nv, nv);
	printf("Initializing start Laplacian matrix\n");
	for(int i=0;i<nv;i++){
			// set value of (i,i) to -6)
			GrB_Matrix_setElement(Laplacian, -6, i, i);   
			if (i+1 < nv){
				GrB_Matrix_setElement(Laplacian, 1, i, i+1);  
			}
			if (i+nx < nv){
				GrB_Matrix_setElement(Laplacian, 1, i, i+nx);  
			}
			if(i+nx*ny < nv){
				GrB_Matrix_setElement(Laplacian, 1, i, i+nx*ny);			
			}
			if (i-1 > 0){
				GrB_Matrix_setElement(Laplacian, 1, i, i-1);  
			}
			if (i-nx >0){
				GrB_Matrix_setElement(Laplacian, 1, i, i-nx);  
			}
			if(i-nx*ny > 0){
				GrB_Matrix_setElement(Laplacian, 1, i, i-nx*ny);			
			}
			
	}    
	    
	 printf("Finished Laplacian matrix \n");
	 GxB_Matrix_fprint(Laplacian,"laplacian",3,stdout);
    */
    
    
	// setting ph(nom,:) = 1
 
 int cxA[] = {15,8,20,5,25,7};
 int cyA[] = {15,5,5,15,5,20};
 int czA[] = {15,5,10,10,20,14};
 int cx,cy,cz;
    for(int i=0;i<nv;i++){
    	    GrB_Matrix_setElement(ph, 1, nom-1, i);
    }
    
    int r0 = 3;

    for (int layer=0;layer<nom-1;layer++){
    	/*cx = 5 + rand()%(GridSize-5) ;
    	cy = 5 + rand()%(GridSize-5) ;
    	cz = 5 + rand()%(GridSize-5) ;*/
    	cx = cxA[layer];
    	cy = cyA[layer];
    	cz = czA[layer];
    	
		for(int i = cx-r0; i<=(cx+r0);i++){
			for(int j = cy-r0; j <= cy+r0; j++){
				for (int k = cz-r0; k <= cz + r0; k++){
					if (sqrt((i-cx)*(i-cx) + (j-cy)*(j-cy) + (k-cz)*(k-cz)) <= r0){
						int ind = (k*nx*ny)+(j*nx) + i;
						GrB_Matrix_setElement(ph, 1.0, layer, ind-1);
						GrB_Matrix_setElement(ph, 0.0, nom-1, ind-1);  					
					} 				
				}			
			}		
		}		

	}    

	pretty_print_matrix_FP64(ph,"ph Values");   	
	


	 
	 GrB_Index M,N;
    GrB_Matrix_nrows(&M, ph);
    GrB_Matrix_ncols(&N, ph);
    printf("Rows in Field %lu cols %lu\n",M,N);
    //pretty_print_matrix_FP64(Laplacian,"Laplacian Values");    
    
    GrB_Matrix_setElement(graph, 4, 5, 2);  // set 1,2 element to 4
    //-------------------------------------------------------------------------------------------------



	 //GrB_Matrix_setElement(g2, 4, 2, 1);
	 //graph = GrB_eWiseAdd(graph,g2);
    //pretty_print_matrix_UINT64(graph, "GRAPH")
    //pretty_print_matrix_FP64(testMult,"lap_ph = ph*laplcaian Values after laplacian multiplication");   	


    int va;
    GrB_Matrix_extractElement(&va,graph,4,2);


    GrB_Index nvals;
    GrB_Matrix_nvals(&nvals, graph);
    //assert(nvals == 1);
    
    GrB_Vector sphi = NULL;
    //GrB_Descriptor d0 = NULL;
	 GrB_Vector_new(&sphi,GrB_FP32,no);   
	 GrB_Vector lap_phi;
	 GrB_Vector_new(&lap_phi,GrB_FP32,no);  

    GrB_Vector field;
    GrB_Vector_new(&field,GrB_FP32,nv);  
    time1 = 1;
    for (time1;time1<time1max;time1++){
    	
    	//sph2 = sparse(no,nv))
    	GrB_Matrix_clear(sph2);
    	printf("Iteration %d\n",time1);
    	//lap_ph = ph*laplacian;
    	GrB_mxm(lap_ph,GrB_NULL,GrB_NULL,GxB_PLUS_TIMES_FP32,ph,lap3D,GrB_NULL);
    	
    	

    	for (int i = 0;i< nv;i++){
    		//printf("i is %d",i)
    		GrB_extract(sphi,GrB_NULL,GrB_NULL,ph,GrB_ALL,nv,i,GrB_NULL);
    		//pretty_print_vector_FP64(sphi,"sphi values");   	
    		GrB_extract(lap_phi,GrB_NULL,GrB_NULL,lap_ph,GrB_ALL,no,i,GrB_NULL);
    		//pretty_print_vector_FP64(lap_phi,"lap_phi values");
    		

    		//GxB_select(Lp,NULL,NULL,GxB_NONZERO,sphi,NULL,NULL);
			// get column of sph at vertex i
			// get col o laplacian of sph
			// ni =  len(row_active)
			// get len of row active
			int ni;
			ni = 0;
			int row_active[no];

			for(int n1 = 0 ; n1<no ; n1++){
				float sphi_val,lap_phi_val ;
				sphi_val = 0;
				lap_phi_val = 0;
				GrB_Vector_extractElement_FP32(&sphi_val,sphi,n1);
				GrB_Vector_extractElement_FP32(&lap_phi_val,lap_phi,n1);
				//printf("\n sphi_val is %f and lap_phi_val is %f", sphi_val,lap_phi_val);
    			if ((sphi_val>0) || (sphi_val == 0 && lap_phi_val > 0)){
					// here we track number of elements and put them in a list we ll use later
					
					row_active[ni] = n1;
					ni++;
					    			
    			}
			}
			//printf("Printing row active here of length %d\n",ni);

			

			for(int n1 = 0 ; n1<ni ; n1++){
				int ii = row_active[n1];
				float pddtt = 0.0;
				float sph1_val, sph2_val;
				GrB_Matrix_extractElement(&sph1_val,ph,ii,i);
				//printf("ph[%d][%d] is %f", ii,i,sph1_val);

				for(int n2=0;n2<ni; n2++){
					int jj = row_active[n2];
					float sum1 = 0;
					float t1,e1;
					
					for (int n3=0;n3<ni;n3++){
						int kk = row_active[n2];
						float a1,a2,w1,w2;
						float lap_phi_scalar, sph_scalar;
						GrB_Matrix_extractElement(&lap_phi_scalar,lap_ph,kk,i);
						//printf("lap_ph[%d][%d] is %f ", ii,i,lap_phi_scalar);
						GrB_Matrix_extractElement(&sph_scalar,ph,kk,i);
						//printf("    ph[%d][%d] is %f ", kk,i,sph_scalar);
						GrB_Matrix_extractElement(&a1,aij,ii,kk);
						GrB_Matrix_extractElement(&a2,aij,jj,kk);
						GrB_Matrix_extractElement(&w1,wij,ii,kk);
						GrB_Matrix_extractElement(&w2,wij,jj,kk);
						//printf("\nfor [%d][%d] a1 is %f w1 is %f  for [%d][%d] a2 is %f w2 is %f",ii,kk,a1,w1,jj,kk,a2,w2);
						
												
						sum1 += 0.5*(a1 - a2)*lap_phi_scalar + (w1 - w2)*sph_scalar;		
						//printf("Sum1 is %f jj is %d\n",sum1,jj);
									
					}				
					GrB_Matrix_extractElement(&t1,tij,ii,jj);
					GrB_Matrix_extractElement(&e1,eij,ii,jj);
					
					GrB_Matrix_extractElement(&sph2_val,ph,jj,i);
					pddtt = pddtt - 2.0*(t1/ni)*(sum1-(8.0/PI)*e1*sqrt(sph1_val*sph2_val));
									
				}		
				float value;
				value =sph1_val + pddtt*delta;
				//printf("\nValue is %f \n",value);
				GrB_Matrix_setElement(sph2, value, ii, i);	
				
				if(value>1.0){
					GrB_Matrix_setElement(sph2, 1.0, ii, i);								
				}
				if (value < 0.0){
					GrB_Matrix_setElement(sph2, 0.0, ii, i);					
				}
				
			}
			
    	}
    
    	GrB_Vector_clear(field);
    	//pretty_print_vector_FP64(field,"Cleaned Field Values");
    	//pretty_print_vector_FP64(ph,"Cleaned ph Values");
    	//pretty_print_vector_FP64(lap_ph,"Cleaned lap_ph Values");
    	GrB_Matrix_dup(&ph,sph2);
    	
    	
    //GrB_Matrix_nrows(&M, ph);
    //GrB_Matrix_ncols(&N, ph);
    //printf("Rows in Field ph %lu cols %lu\n",M,N);
    	GrB_transpose(finalPh,GrB_NULL,GrB_NULL,ph,GrB_NULL);
    	
    	//need to do two things divide matrix by sum along columns, then multiply matrix by itself ewise and again do col wise add
    	// 
    	//GrB_eWiseMult_Matrix_Semiring(finalPh,GrB_NULL,GrB_NULL, , finalPh, finalPh, GrB_NULL);
		GrB_Matrix_reduce_Monoid(field,GrB_NULL,GrB_NULL, GxB_PLUS_FP32_MONOID, finalPh,GrB_NULL);
		
		// create a divide matrix of sum of field values along column
		for(int it = 0;it <nv; it++){
			float divider = 0;
			GrB_Vector_extractElement_FP32(&divider, field, it);
			for (int it2=0;it2<no;it2++){
				GrB_Matrix_setElement_FP32(divideMatrix,divider,it2,it);
			}			
		}
		printf("Before printing field values");
		//pretty_print_vector_FP64(field,"Field Values");
		
		GrB_Matrix_nrows(&M, divideMatrix);
      GrB_Matrix_ncols(&N, divideMatrix);
      printf("Rows in divideMatrix are %lu cols %lu\n",M,N);		
		
		
		//pretty_print_matrix_FP64(divideMatrix,"Field Values");
		// Now lets do the normalization of values
		GrB_eWiseMult_Matrix_BinaryOp(ph,GrB_NULL,GrB_NULL,GrB_DIV_FP32 , ph,divideMatrix, GrB_NULL);
		printf("ewise Mult done");
		GrB_transpose(finalPh,GrB_NULL,GrB_NULL,ph,GrB_NULL);
		//GrB_Matrix_reduce_Monoid(field,GrB_NULL,GrB_NULL, GxB_PLUS_FP32_MONOID, finalPh,GrB_NULL);
		//pretty_print_vector_FP64(field,"Field Values");
		
		// here lets multiply finalPh.^2 and then add along col into field
		GrB_eWiseMult_Matrix_BinaryOp(finalPh,GrB_NULL,GrB_NULL,GrB_TIMES_FP32 , finalPh,finalPh, GrB_NULL);
		GrB_Matrix_reduce_Monoid(field,GrB_NULL,GrB_NULL, GxB_PLUS_FP32_MONOID, finalPh,GrB_NULL);
		//printf("Values after normalizing and then adding along column");
		//pretty_print_matrix_FP64(finalPh,"Full field Values");		
		//pretty_print_vector_FP64(field,"Field Values");

		// Saving data as VTK
		dataVtk( time1,GridSize,GridSize,GridSize,  field );

		int nonZero = 0;
		//GrB_Vector_nvals(&nonZero,field);
		//printf("Non zero entries in vector : %d \n",nonZero);
		 
    
    }

    // Cleanup
    GrB_free(&graph);
    GrB_finalize();
}
