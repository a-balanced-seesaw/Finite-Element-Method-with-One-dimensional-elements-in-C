	#include <stdio.h>
	#include <stdlib.h>
	#include <string.h>
	#include <math.h>

	
	struct Coordinates
	{
		long double x;
		long double y;
		long coord_no;
	};

	struct Connectivity
	{
	 	long conn_no;
	 	long conn_i;
	 	long conn_j;
	};

	struct Ele_props
	{
		long elementno;
		long double E_ele;
		long double A_ele;
		long double I_ele;
		long double L_ele;
		long double l_truss;
		long double m_truss;
		long double a_frame;
		long double b_frame;
		long double c_frame;
		long double d_frame;
		long double l_frame;
		long double m_frame;
	};	

	struct pointForce
	{
		long forceno;
		long double force_mag;
	};

	long N;

// Functions
// Function to get cofactor of A[p][q] in temp[][]. n is current dimension of A[][]
void getCofactor(long double A[N][N], long double temp[N][N], int p, int q, int n);

/* Recursive function for finding determinant of matrix.
   n is current dimension of A[][]. */
long double determinant(long double A[N][N], int n);

// Function to get adjoint of A[N][N] in adj[N][N].
void adjoint(long double A[N][N],long double adj[N][N]);

// Function to calculate and store inverse, returns false if matrix is singular
int inverse(long double A[N][N], long double inverse[N][N]);

// Generic function to display the matrix.  We use it to display both adjoin and inverse. adjoin is integer matrix and inverse is a float.
void display(long double A[N][N]);


	//Main Function

	long main()
	{
		printf("\n WELCOME TO THIS 1D-ELEMENTS F.E.M. SOLVER --- A COURSE PROJECT MADE BY BONTHA MUKUND SAI (MDM15B036) AS PART OF THE 'ADVANCED MECHANICS OF MATERIALS' COURSE\n");

		
		//ENTERING COORDINATES OF NODES...
		long n_nodes;
		printf("Enter the number of points: ");
		scanf("%ld",&n_nodes);

		struct Coordinates c[n_nodes];

		printf("Enter the %ld points with their node numbers and coordinates in the form 'Node no : x-coordinate, y-coordinate':\n", n_nodes);

		for (long i = 0; i < n_nodes; i++)
		  {
			scanf("%ld:%Lf,%Lf",&c[i].coord_no,&c[i].x,&c[i].y);
		  }
		

	/*
		//TESTING.....
		for (long i=0; i< n_nodes; i++)
		  {
			printf("The node number and coordinates are:\n");		
			printf("\n Node %ld - (%Lf,%Lf) \n",c[i].coord_no,c[i].x,c[i].y);
		  }

	*/  


	//!!!!!!!!!!!!!!OK till here!!!!!!!!!!!!!!!


		long n_elements;	
	

		printf("Enter the number of elements:");
		scanf("%ld",&n_elements);

		printf("\nI have assumed three degrees of freedom (dofs) at each node\n Horizontal DoF - 3*nodenumber-2 \nVertical DoF - 3*nodenumber-1\nRotational DoF - 3*nodenumber\n I have considered only those degrees of freedom relevant for each 1D element type in their respective formulations.\n");

		struct Connectivity conn[n_elements];
									
		long n_dofs=3*n_nodes;
							
		long conn_no[n_elements];
		long conn_i[n_elements];
		long conn_j[n_elements];


		struct Ele_props elep[n_elements];
		
		for(long k=0;k<n_elements;k++) 		//INITIALIZING MATERIAL PROPERTIES TO ZERO....
			{
			elep[k].elementno=0;
			elep[k].E_ele=0;
			elep[k].A_ele=0;
			elep[k].I_ele=0;
			elep[k].L_ele=0;
			elep[k].l_truss=0;
			elep[k].m_truss=0;
			elep[k].a_frame=0;
			elep[k].b_frame=0;
			elep[k].c_frame=0;
			elep[k].d_frame=0;
			elep[k].l_frame=0;
			elep[k].m_frame=0;
			}
	
		
		long double K_ele[n_elements][n_dofs][n_dofs];  	//ELEMENT STIFFNESS MATRIX

		for (long i = 0; i < n_elements; i++) 		//INITIALIZING ELEMENT STIFFNESS MATRIX TO ZERO....
		{
			for (long j = 0; j < n_dofs; j++)
			{
				for (long k = 0; k < n_dofs; k++)
				{
					K_ele[i][j][k]=0.0000000000000000000000000000;
				}
			}
		}


		long ele_type;
		
		for (long i = 0; i < n_elements; i++)
		{
			elep[i].elementno=i+1;

			printf("\n *************Element Number %ld************* \n", i+1);

			printf("\nEnter the nodes that the element connects. Type this in the form -> 'ith node - jth node'\n");

			scanf("%ld-%ld",&conn[i].conn_i,&conn[i].conn_j);

			printf("What 1D element is it? Enter the corresponding number: \n 1 - Bar \t 2 - Beam \t 3 - Frame \t 4 - Truss \n");
			
			scanf("%ld", &ele_type);
			
			switch (ele_type)
				{
					case 1:
						//**********BAR FORMULATION***********

						printf("\nEnter the value of Young's Modulus E: ");
						scanf("%Lf", &elep[i].E_ele);
						printf("\nEnter the value of Cross Sectional Area A: ");
						scanf("%Lf", &elep[i].A_ele);
					
						//element length
						elep[i].L_ele=sqrt(pow((c[conn[i].conn_j-1].x-c[conn[i].conn_i-1].x),2)+pow((c[conn[i].conn_j-1].y-c[conn[i].conn_i-1].y),2));

						
						//CONSTRUCTING ELEMENT STIFFNESS MATRIX

						K_ele[i][3*conn[i].conn_i-1-2][3*conn[i].conn_i-1-2]=(elep[i].E_ele*elep[i].A_ele)/elep[i].L_ele;
						K_ele[i][3*conn[i].conn_i-1-2][3*conn[i].conn_j-1-2]=-(elep[i].E_ele*elep[i].A_ele)/elep[i].L_ele;
						K_ele[i][3*conn[i].conn_j-1-2][3*conn[i].conn_i-1-2]=-(elep[i].E_ele*elep[i].A_ele)/elep[i].L_ele;
						K_ele[i][3*conn[i].conn_j-1-2][3*conn[i].conn_j-1-2]=(elep[i].E_ele*elep[i].A_ele)/elep[i].L_ele;			
						
						
						//PRINTING ELEMENT STIFFNESS MATRIX

						printf("\n The stiffness matrix of this bar element is as follows:\n");
						
						for (long j = 0; j < n_dofs; j++)
						{
							for (long k = 0; k < n_dofs; k++)
							{
								printf("%Lf\t", K_ele[i][j][k]);
							}
							printf("\n");
						}
						
						break;
					
					case 2:
						//**********BEAM FORMULATION***********
						
						printf("\nEnter the value of Young's Modulus E: ");
						scanf("%Lf", &elep[i].E_ele);
						printf("\nEnter the value of Moment of Inertia I: ");
						scanf("%Lf", &elep[i].I_ele);

						//Element length
						elep[i].L_ele=sqrt(pow((c[conn[i].conn_j-1].x-c[conn[i].conn_i-1].x),2)+pow((c[conn[i].conn_j-1].y-c[conn[i].conn_i-1].y),2));


						//CONSTRUCTING ELEMENT STIFFNESS MATRIX

						K_ele[i][3*conn[i].conn_i-1-1][3*conn[i].conn_i-1-1]=(elep[i].E_ele*elep[i].I_ele)/(elep[i].L_ele*elep[i].L_ele*elep[i].L_ele)*12;
						K_ele[i][3*conn[i].conn_i-1-1][3*conn[i].conn_i-1]=(elep[i].E_ele*elep[i].I_ele)/(elep[i].L_ele*elep[i].L_ele*elep[i].L_ele)*6*elep[i].L_ele;
						K_ele[i][3*conn[i].conn_i-1-1][3*conn[i].conn_j-1-1]=(elep[i].E_ele*elep[i].I_ele)/(elep[i].L_ele*elep[i].L_ele*elep[i].L_ele)*(-12);
						K_ele[i][3*conn[i].conn_i-1-1][3*conn[i].conn_j-1]=(elep[i].E_ele*elep[i].I_ele)/(elep[i].L_ele*elep[i].L_ele*elep[i].L_ele)*6*elep[i].L_ele;

						K_ele[i][3*conn[i].conn_i-1][3*conn[i].conn_i-1-1]=(elep[i].E_ele*elep[i].I_ele)/(elep[i].L_ele*elep[i].L_ele*elep[i].L_ele)*6*elep[i].L_ele;
						K_ele[i][3*conn[i].conn_i-1][3*conn[i].conn_i-1]=(elep[i].E_ele*elep[i].I_ele)/(elep[i].L_ele*elep[i].L_ele*elep[i].L_ele)*4*(pow((elep[i].L_ele),2));
						K_ele[i][3*conn[i].conn_i-1][3*conn[i].conn_j-1-1]=(elep[i].E_ele*elep[i].I_ele)/(elep[i].L_ele*elep[i].L_ele*elep[i].L_ele)*(-6)*(elep[i].L_ele);
						K_ele[i][3*conn[i].conn_i-1][3*conn[i].conn_j-1]=(elep[i].E_ele*elep[i].I_ele)/(elep[i].L_ele*elep[i].L_ele*elep[i].L_ele)*2*(pow((elep[i].L_ele),2));
						
						K_ele[i][3*conn[i].conn_j-1-1][3*conn[i].conn_i-1-1]=(elep[i].E_ele*elep[i].I_ele)/(elep[i].L_ele*elep[i].L_ele*elep[i].L_ele)*(-12);
						K_ele[i][3*conn[i].conn_j-1-1][3*conn[i].conn_i-1]=(elep[i].E_ele*elep[i].I_ele)/(elep[i].L_ele*elep[i].L_ele*elep[i].L_ele)*(-6)*(elep[i].L_ele);
						K_ele[i][3*conn[i].conn_j-1-1][3*conn[i].conn_j-1-1]=(elep[i].E_ele*elep[i].I_ele)/(elep[i].L_ele*elep[i].L_ele*elep[i].L_ele)*12;
						K_ele[i][3*conn[i].conn_j-1-1][3*conn[i].conn_j-1]=(elep[i].E_ele*elep[i].I_ele)/(elep[i].L_ele*elep[i].L_ele*elep[i].L_ele)*(-6)*(elep[i].L_ele);
						
						K_ele[i][3*conn[i].conn_j-1][3*conn[i].conn_i-1-1]=(elep[i].E_ele*elep[i].I_ele)/(elep[i].L_ele*elep[i].L_ele*elep[i].L_ele)*(6)*(elep[i].L_ele);
						K_ele[i][3*conn[i].conn_j-1][3*conn[i].conn_i-1]=(elep[i].E_ele*elep[i].I_ele)/(elep[i].L_ele*elep[i].L_ele*elep[i].L_ele)*2*(pow((elep[i].L_ele),2));
						K_ele[i][3*conn[i].conn_j-1][3*conn[i].conn_j-1-1]=(elep[i].E_ele*elep[i].I_ele)/(elep[i].L_ele*elep[i].L_ele*elep[i].L_ele)*(-6)*(elep[i].L_ele);
						K_ele[i][3*conn[i].conn_j-1][3*conn[i].conn_j-1]=(elep[i].E_ele*elep[i].I_ele)/(elep[i].L_ele*elep[i].L_ele*elep[i].L_ele)*4*(pow((elep[i].L_ele),2));

						
						//PRINTING ELEMENT STIFFNESS MATRIX

						printf("\n The stiffness matrix of this beam element is as follows:\n");

						for (long j = 0; j < n_dofs; j++)
						{
							for (long k = 0; k < n_dofs; k++)
							{
								printf("%Lf\t", K_ele[i][j][k]);
							}
							printf("\n");
						}
						
						break;
						
					case 3:
						//**********FRAME FORMULATION***********

						printf("\nEnter the value of Young's Modulus E: ");
						scanf("%Lf", &elep[i].E_ele);
						printf("\nEnter the value of Cross Sectional Area A: ");
						scanf("%Lf", &elep[i].A_ele);
						printf("\nEnter the value of Moment of Inertia I: ");
						scanf("%Lf", &elep[i].I_ele);

						//Element length
						elep[i].L_ele=sqrt(pow((c[conn[i].conn_j-1].x-c[conn[i].conn_i-1].x),2)+pow((c[conn[i].conn_j-1].y-c[conn[i].conn_i-1].y),2));

						//I OBTAINED THE GENERAL STIFFNESS MATRIX FOR FRAME AS A SINGLE MATRIX

						//CONSTRUCTING ELEMENT STIFFNESS MATRIX

						elep[i].a_frame=((elep[i].E_ele)*(elep[i].A_ele))/elep[i].L_ele;
						elep[i].b_frame=(12*elep[i].E_ele*elep[i].I_ele)/(elep[i].L_ele*elep[i].L_ele*elep[i].L_ele);
						elep[i].c_frame=(6*elep[i].E_ele*elep[i].I_ele)/(elep[i].L_ele*elep[i].L_ele);
						elep[i].d_frame=(2*elep[i].E_ele*elep[i].I_ele)/(elep[i].L_ele);
						elep[i].l_frame=(c[conn[i].conn_j-1].x-c[conn[i].conn_i-1].x)/elep[i].L_ele;
						elep[i].m_frame=(c[conn[i].conn_j-1].y-c[conn[i].conn_i-1].y)/elep[i].L_ele;


						K_ele[i][3*conn[i].conn_i-1-2][3*conn[i].conn_i-1-2]=elep[i].a_frame*elep[i].l_frame*elep[i].l_frame+elep[i].b_frame*elep[i].m_frame*elep[i].m_frame;
						K_ele[i][3*conn[i].conn_i-1-2][3*conn[i].conn_i-1-1]=(elep[i].a_frame-elep[i].b_frame)*elep[i].l_frame*elep[i].m_frame;
						K_ele[i][3*conn[i].conn_i-1-2][3*conn[i].conn_i-1-0]=-elep[i].c_frame*elep[i].m_frame;
						K_ele[i][3*conn[i].conn_i-1-2][3*conn[i].conn_j-1-2]=-(elep[i].a_frame*elep[i].l_frame*elep[i].l_frame+elep[i].b_frame*elep[i].m_frame*elep[i].m_frame);
						K_ele[i][3*conn[i].conn_i-1-2][3*conn[i].conn_j-1-1]=-((elep[i].a_frame-elep[i].b_frame)*elep[i].l_frame*elep[i].m_frame);
						K_ele[i][3*conn[i].conn_i-1-2][3*conn[i].conn_j-1-0]=-elep[i].c_frame*elep[i].m_frame;

						K_ele[i][3*conn[i].conn_i-1-1][3*conn[i].conn_i-1-2]=(elep[i].a_frame-elep[i].b_frame)*elep[i].l_frame*elep[i].m_frame;
						K_ele[i][3*conn[i].conn_i-1-1][3*conn[i].conn_i-1-1]=elep[i].a_frame*elep[i].m_frame*elep[i].m_frame+elep[i].b_frame*elep[i].l_frame*elep[i].l_frame;
						K_ele[i][3*conn[i].conn_i-1-1][3*conn[i].conn_i-1-0]=elep[i].c_frame*elep[i].l_frame;
						K_ele[i][3*conn[i].conn_i-1-1][3*conn[i].conn_j-1-2]=-((elep[i].a_frame-elep[i].b_frame)*elep[i].l_frame*elep[i].m_frame);
						K_ele[i][3*conn[i].conn_i-1-1][3*conn[i].conn_j-1-1]=-(elep[i].a_frame*elep[i].m_frame*elep[i].m_frame+elep[i].b_frame*elep[i].l_frame*elep[i].l_frame);
						K_ele[i][3*conn[i].conn_i-1-1][3*conn[i].conn_j-1-0]=elep[i].c_frame*elep[i].l_frame;

						K_ele[i][3*conn[i].conn_i-1-0][3*conn[i].conn_i-1-2]=-elep[i].c_frame*elep[i].m_frame;
						K_ele[i][3*conn[i].conn_i-1-0][3*conn[i].conn_i-1-1]=elep[i].c_frame*elep[i].l_frame;
						K_ele[i][3*conn[i].conn_i-1-0][3*conn[i].conn_i-1-0]=2*elep[i].d_frame;
						K_ele[i][3*conn[i].conn_i-1-0][3*conn[i].conn_j-1-2]=elep[i].c_frame*elep[i].m_frame;
						K_ele[i][3*conn[i].conn_i-1-0][3*conn[i].conn_j-1-1]=-elep[i].c_frame*elep[i].l_frame;
						K_ele[i][3*conn[i].conn_i-1-0][3*conn[i].conn_j-1-0]=elep[i].d_frame;

						K_ele[i][3*conn[i].conn_j-1-2][3*conn[i].conn_i-1-2]=-(elep[i].a_frame*elep[i].l_frame*elep[i].l_frame+elep[i].b_frame*elep[i].m_frame*elep[i].m_frame);
						K_ele[i][3*conn[i].conn_j-1-2][3*conn[i].conn_i-1-1]=-((elep[i].a_frame-elep[i].b_frame)*elep[i].l_frame*elep[i].m_frame);
						K_ele[i][3*conn[i].conn_j-1-2][3*conn[i].conn_i-1-0]=elep[i].c_frame*elep[i].m_frame;
						K_ele[i][3*conn[i].conn_j-1-2][3*conn[i].conn_j-1-2]=elep[i].a_frame*elep[i].l_frame*elep[i].l_frame+elep[i].b_frame*elep[i].m_frame*elep[i].m_frame;
						K_ele[i][3*conn[i].conn_j-1-2][3*conn[i].conn_j-1-1]=(elep[i].a_frame-elep[i].b_frame)*elep[i].l_frame*elep[i].m_frame;
						K_ele[i][3*conn[i].conn_j-1-2][3*conn[i].conn_j-1-0]=elep[i].c_frame*elep[i].m_frame;

						K_ele[i][3*conn[i].conn_j-1-1][3*conn[i].conn_i-1-2]=-((elep[i].a_frame-elep[i].b_frame)*elep[i].l_frame*elep[i].m_frame);
						K_ele[i][3*conn[i].conn_j-1-1][3*conn[i].conn_i-1-1]=-(elep[i].a_frame*elep[i].m_frame*elep[i].m_frame+elep[i].b_frame*elep[i].l_frame*elep[i].l_frame);
						K_ele[i][3*conn[i].conn_j-1-1][3*conn[i].conn_i-1-0]=-elep[i].c_frame*elep[i].l_frame;
						K_ele[i][3*conn[i].conn_j-1-1][3*conn[i].conn_j-1-2]=(elep[i].a_frame-elep[i].b_frame)*elep[i].l_frame*elep[i].m_frame;
						K_ele[i][3*conn[i].conn_j-1-1][3*conn[i].conn_j-1-1]=elep[i].a_frame*elep[i].m_frame*elep[i].m_frame+elep[i].b_frame*elep[i].l_frame*elep[i].l_frame;
						K_ele[i][3*conn[i].conn_j-1-1][3*conn[i].conn_j-1-0]=-elep[i].c_frame*elep[i].l_frame;

						K_ele[i][3*conn[i].conn_j-1-0][3*conn[i].conn_i-1-2]=-elep[i].c_frame*elep[i].m_frame;
						K_ele[i][3*conn[i].conn_j-1-0][3*conn[i].conn_i-1-1]=elep[i].c_frame*elep[i].l_frame;
						K_ele[i][3*conn[i].conn_j-1-0][3*conn[i].conn_i-1-0]=elep[i].d_frame;
						K_ele[i][3*conn[i].conn_j-1-0][3*conn[i].conn_j-1-2]=elep[i].c_frame*elep[i].m_frame;
						K_ele[i][3*conn[i].conn_j-1-0][3*conn[i].conn_j-1-1]=-elep[i].c_frame*elep[i].l_frame;
						K_ele[i][3*conn[i].conn_j-1-0][3*conn[i].conn_j-1-0]=2*elep[i].d_frame;

						//PRINTING ELEMENT STIFFNESS MATRIX

						printf("\n The stiffness matrix of this frame element is as follows:\n");

						for (long j = 0; j < n_dofs; j++)
						{
							for (long k = 0; k < n_dofs; k++)
							{
								printf("%Lf", K_ele[i][j][k]);
							}
							printf("\n");
						}
						
						break;
						
					case 4:
					//**********TRUSS FORMULATION***********

						printf("\nEnter the value of Young's Modulus E: ");
						scanf("%Lf", &elep[i].E_ele);
						printf("\nEnter the value of Cross Sectional Area A: ");
						scanf("%Lf", &elep[i].A_ele);
					

						elep[i].L_ele=sqrt(pow((c[conn[i].conn_j-1].x-c[conn[i].conn_i-1].x),2)+pow((c[conn[i].conn_j-1].y-c[conn[i].conn_i-1].y),2));
						elep[i].l_truss=(c[conn[i].conn_j-1].x-c[conn[i].conn_i-1].x)/elep[i].L_ele;
						elep[i].m_truss=(c[conn[i].conn_j-1].y-c[conn[i].conn_i-1].y)/elep[i].L_ele;

						//CONSTRUCTING ELEMENT STIFFNESS MATRIX

						K_ele[i][3*conn[i].conn_i-1-2][3*conn[i].conn_i-1-2]=((elep[i].E_ele*elep[i].A_ele)/elep[i].L_ele)*elep[i].l_truss*elep[i].l_truss;
						K_ele[i][3*conn[i].conn_i-1-2][3*conn[i].conn_i-1-1]=((elep[i].E_ele*elep[i].A_ele)/elep[i].L_ele)*elep[i].l_truss*elep[i].m_truss;
						K_ele[i][3*conn[i].conn_i-1-2][3*conn[i].conn_j-1-2]=((elep[i].E_ele*elep[i].A_ele)/elep[i].L_ele)*(-elep[i].l_truss*elep[i].l_truss);
						K_ele[i][3*conn[i].conn_i-1-2][3*conn[i].conn_j-1-1]=((elep[i].E_ele*elep[i].A_ele)/elep[i].L_ele)*(-elep[i].l_truss*elep[i].m_truss);

						K_ele[i][3*conn[i].conn_i-1-1][3*conn[i].conn_i-1-2]=((elep[i].E_ele*elep[i].A_ele)/elep[i].L_ele)*elep[i].l_truss*elep[i].m_truss;
						K_ele[i][3*conn[i].conn_i-1-1][3*conn[i].conn_i-1-1]=((elep[i].E_ele*elep[i].A_ele)/elep[i].L_ele)*elep[i].m_truss*elep[i].m_truss;
						K_ele[i][3*conn[i].conn_i-1-1][3*conn[i].conn_j-1-2]=((elep[i].E_ele*elep[i].A_ele)/elep[i].L_ele)*(-elep[i].l_truss*elep[i].m_truss);
						K_ele[i][3*conn[i].conn_i-1-1][3*conn[i].conn_j-1-1]=((elep[i].E_ele*elep[i].A_ele)/elep[i].L_ele)*(-elep[i].m_truss*elep[i].m_truss);

						K_ele[i][3*conn[i].conn_j-1-2][3*conn[i].conn_i-1-2]=((elep[i].E_ele*elep[i].A_ele)/elep[i].L_ele)*(-elep[i].l_truss*elep[i].l_truss);
						K_ele[i][3*conn[i].conn_j-1-2][3*conn[i].conn_i-1-1]=((elep[i].E_ele*elep[i].A_ele)/elep[i].L_ele)*(-elep[i].l_truss*elep[i].m_truss);
						K_ele[i][3*conn[i].conn_j-1-2][3*conn[i].conn_j-1-2]=((elep[i].E_ele*elep[i].A_ele)/elep[i].L_ele)*elep[i].l_truss*elep[i].l_truss;
						K_ele[i][3*conn[i].conn_j-1-2][3*conn[i].conn_j-1-1]=((elep[i].E_ele*elep[i].A_ele)/elep[i].L_ele)*elep[i].l_truss*elep[i].m_truss;

						K_ele[i][3*conn[i].conn_j-1-1][3*conn[i].conn_i-1-2]=((elep[i].E_ele*elep[i].A_ele)/elep[i].L_ele)*(-elep[i].l_truss*elep[i].m_truss);
						K_ele[i][3*conn[i].conn_j-1-1][3*conn[i].conn_i-1]=((elep[i].E_ele*elep[i].A_ele)/elep[i].L_ele)*(-elep[i].m_truss*elep[i].m_truss);
						K_ele[i][3*conn[i].conn_j-1-1][3*conn[i].conn_j-1-2]=((elep[i].E_ele*elep[i].A_ele)/elep[i].L_ele)*elep[i].l_truss*elep[i].m_truss;
						K_ele[i][3*conn[i].conn_j-1-1][3*conn[i].conn_j-1-1]=((elep[i].E_ele*elep[i].A_ele)/elep[i].L_ele)*elep[i].m_truss*elep[i].m_truss;

						
						printf("\n The stiffness matrix of this truss element is as follows:\n");

						for (long j = 0; j < n_dofs; j++)
						{
							for (long k = 0; k < n_dofs; k++)
							{
								printf("%Lf\t", K_ele[i][j][k]);
							}
							printf("\n");
						}
						
						break;
						
					default:
						printf("\n Please enter one of 1, 2, 3 or 4.\n");
						break;
					
			}
					
		}


		long double K_global[n_dofs][n_dofs];
		
		for (long i = 0; i < n_dofs; i++) 		//INITIALIZING GLOBAL STIFFNESS MATRIX TO ZERO....
		{
			for (long j = 0; j < n_dofs; j++)
			{
					K_global[i][j]=0;
			}
		}

		//ASSEMBLY OF ELEMENT STIFFNESS MATRICES TO FORM THE GLOBAL STIFFNESS MATRIX
		for (long i = 0; i < n_dofs; i++)
		{
			for (long j = 0; j < n_dofs; j++)
			{
				for (long k = 0; k < n_elements; k++)
				{
					K_global[i][j] = K_global[i][j] + K_ele[k][i][j];
				}
			}
		}

		//PRINTING GLOBAL STIFFNESS MATRIX

		printf("\n The stiffness matrix of the structure is as follows:\n");

						
						for (long i = 0; i < n_dofs; i++)
						{
							for (long j = 0; j < n_dofs; j++)
							{
								printf("%Lf\t", K_global[i][j]);
							}
							printf("\n");
						}


	//REMOVING ZERO ROWS AND COLUMNS (TO BE ABLE TO TAKE INVERSE WHILE SOLVING)

	long flag;
    long ctr = 0;
    long zero_rows[n_dofs];
    long num_zerodisp;
    long temp_row;

    for (long i=0;i<n_dofs;i++)
    	zero_rows[i] = 0;

	printf("\nEnter the number of dofs along which displacement is zero: ");
	scanf("%ld",&num_zerodisp);

	printf("\nMention the %ld dofs along which displacement is zero:\n",num_zerodisp);

	for(long i=0;i<num_zerodisp;i++)
	{
		scanf("%ld",&temp_row);
		zero_rows[temp_row-1] = 1;
	}
    
    printf("The zero rows and columns (to be removed) are: ");
    // count the number of rows that have zeros
    for( long i = 0; i < n_dofs; i++)
    {
        flag = 0;
        for(long j=0; j<n_dofs; j++)
        {
            if(K_global[i][j] != 0)
                flag = 1;
        }

        if(flag == 0)
        {
            zero_rows[i] = 1;
            printf("%ld ", i+1);
            ctr++;
        }	
    }

    // count the number of rows that will become zero are
    printf("\nTotally the rows and columns that are going to be removed are: ");
    	ctr = 0;         						// reinitialize the counter
    	for( long i =0; i <n_dofs; i++)
    	if(zero_rows[i] == 1)
    	{
    		ctr++;
    		printf("%ld ", i+1);
    	}

    printf("\nNumber of rows that need to be removed is %ld\n", ctr);
   
    printf("Totally the rows and columns that are going to remain are: ");
    for( long i =0; i <n_dofs; i++)
    	if(zero_rows[i] == 0)
    		{
    		printf("%ld ", i+1);
       		}

    printf("\n");

    long double K_global_reduced[n_dofs-ctr][n_dofs-ctr];
    long m=0,n=0;
    for ( long i =0; i < n_dofs; i++)
    {
        for ( long j=0; j< n_dofs; j++)
        {
            if(zero_rows[i] != 1 && zero_rows[j]!=1)
            {
                // copy in b row
                K_global_reduced[m][n] = K_global[i][j];
                n = (n+1)%(n_dofs-ctr);
                if(n == 0)
                    m++;
            }
        }
    }


    for( long i = 0; i < n_dofs-ctr; i++)
    {
        for(long j=0; j<n_dofs-ctr; j++)
        {
            printf("%Lf\t\t", K_global_reduced[i][j]);
        }
        printf("\n" );

	}		

//FORMULATING FORCE VECTORS

	long no_forces;

	printf("\nEnter the number of point forces acting on structure: ");
	scanf("%ld",&no_forces);

	struct pointForce frc[n_dofs];

	for(long i=0; i<n_dofs;i++)   //Initializing force
	{
		frc[i].forceno=0;
		frc[i].force_mag=0.0;
	}

	long double forcevec[n_dofs][1];

	for (long i = 0; i < n_dofs; i++)
		forcevec[i][1]=0;

	printf("\nEnter the dofs and corresponding force in the form: 'dof:force'. Again...\n\tHorizontal DoFs are of the form '3*nodenumber-2'.\n\tVertical DoFs are of the form '3*nodenumber-1'.\n\tRotational DoFs are of the form '3*nodenumber'.\n");
		
	long temp_dof; 
	long double temp_force_mag;
	for (long i = 0; i < no_forces; i++)
	{
		scanf("%ld:%Lf",&temp_dof,&temp_force_mag);
		frc[temp_dof-1].force_mag = temp_force_mag;
	}

	printf("\nThe point forces act along the following dofs:\n");

	for (long i = 0; i < n_dofs; i++)
	{
		if(frc[i].force_mag != 0.0)
			printf("\t%Lf along DoF %ld.\n",frc[i].force_mag,i+1);
	}

	for (long i=0;i<n_dofs;i++)
	{
		forcevec[i][0]=frc[i].force_mag;
	}

	for (long i = 0; i < n_dofs; i++)
	{
		printf("%Lf\t", forcevec[i][0]);
	}
	
	printf("\n\n\n");

	printf("Required Modified Force Vector, F_modified is:\n");

	long double reduced_forcevec[n_dofs-ctr][1];

	long temp_index_ctr=0;
	for(long i=0; i<n_dofs; i++)
	{
		if(zero_rows[i]==0)
		{
			reduced_forcevec[temp_index_ctr][0] = forcevec[i][0];
			temp_index_ctr++;
		}
	}	

	for(long i=0; i<n_dofs-ctr; i++)
		printf("%Lf\t", reduced_forcevec[i][0]);

	printf("\n");



//SOLUTION: Kmod*Q=Fmod:

// matrix of size K_global_reduced is going to be processed using the functions defined at the beginning of this code

N = n_dofs-ctr;

long double adj[N][N];
long double inv[N][N];

for(int i=0; i<N; i++)
	for(int j=0; j<N; j++)
		adj[i][j] = 0;

for(int i=0; i<N; i++)
	for(int j=0; j<N; j++)
		inv[i][j] = 0;

printf("\n'K_modified' matrix is :\n");;
display(K_global_reduced);

printf("\nThe Adjoint is :\n");
adjoint(K_global_reduced, adj);
display(adj);

printf("\nThe Inverse of K_modified is :\n");
if (inverse(K_global_reduced, inv))
    display(inv);

printf("\nMultiplying inverse with force vector.\n");
long double answer[N][1];

for (int i =0; i<N; i++)
{
	answer[i][1] = 0;
}


long double temp1=0;

for( int i =0; i< N; i++)
{
	for(int j=0; j<1;j++)
	{
		for(int k=0; k<N; k++ )
		{
			temp1 = temp1 + inv[i][k] * reduced_forcevec[k][j];
		}

		answer[i][j] = temp1;
		
		temp1 = 0;
	}
}

// PRINT DISPLACEMENTS
printf("The displacements along DoFs ");
    for( long i =0; i <n_dofs; i++)
    	if(zero_rows[i] == 0)
    		{
    		printf("%ld, ", i+1);
       		}
     printf("is ");
for ( int i =0; i<N; i++)
{
	printf("%Lf, ",answer[i][0]);
}

printf("respectively.\n");


	return 0;

	}


// functions definitions
// Function to get cofactor of A[p][q] in temp[][]. n is current dimension of A[][]
void getCofactor(long double A[N][N], long double temp[N][N], int p, int q, int n)
{
    int i = 0, j = 0;
 
    // Looping for each element of the matrix
    for (int row = 0; row < n; row++)
    {
        for (int col = 0; col < n; col++)
        {
            //  Copying into temporary matrix only those element
            //  which are not in given row and column
            if (row != p && col != q)
            {
                temp[i][j++] = A[row][col];
 
                // Row is filled, so increase row index and
                // reset col index
                if (j == n - 1)
                {
                    j = 0;
                    i++;
                }
            }
        }
    }
}
 
/* Recursive function for finding determinant of matrix.
   n is current dimension of A[][]. */
long double determinant(long double A[N][N], int n)
{
    long double D = 0; // Initialize result
 
    //  Base case : if matrix contains single element
    if (n == 1)
        return A[0][0];
 
    long double temp[N][N]; // To store cofactors
 
    long double sign = 1;  // To store sign multiplier
 
     // Iterate for each element of first row
    for (int f = 0; f < n; f++)
    {
        // Getting Cofactor of A[0][f]
        getCofactor(A, temp, 0, f, n);
        D += sign * A[0][f] * determinant(temp, n - 1);
 
        // terms are to be added with alternate sign
        sign = -sign;
    }
 
    return D;
}
 
// Function to get adjoint of A[N][N] in adj[N][N].
void adjoint(long double A[N][N],long double adj[N][N])
{
    if (N == 1)
    {
        adj[0][0] = 1;
        return;
    }
 
    // temp is used to store cofactors of A[][]
    int sign = 1;
    long double temp[N][N];
 
    for (int i=0; i<N; i++)
    {
        for (int j=0; j<N; j++)
        {
            // Get cofactor of A[i][j]
            getCofactor(A, temp, i, j, N);
 
            // sign of adj[j][i] positive if sum of row
            // and column indexes is even.
            sign = ((i+j)%2==0)? 1: -1;
 
            // Interchanging rows and columns to get the
            // transpose of the cofactor matrix
            adj[j][i] = (sign)*(determinant(temp, N-1));
        }
    }
}
 
// Function to calculate and store inverse, returns false if
// matrix is singular
int inverse(long double A[N][N], long double inverse[N][N])
{
    // Find determinant of A[][]
    long double det = determinant(A, N);
    if (det == 0)
    {
        printf("Singular matrix, can't find its inverse.\n");
        return 0;
    }
 
    // Find adjoint
    long double adj[N][N];
    adjoint(A, adj);
 
    // Find Inverse using formula "inverse(A) = adj(A)/det(A)"
    for (int i=0; i<N; i++)
        for (int j=0; j<N; j++)
            inverse[i][j] = adj[i][j]/det;
 
    return 1;
}

// Generic function to display the matrix.  We use it to display
// both adjoint and inverse. adjoint is integer matrix and inverse
// is a float.
void display(long double A[N][N])
{
    for (int i=0; i<N; i++)
    {
        for (int j=0; j<N; j++)
            printf("%.10Lf ",A[i][j]);
        printf("\n");
    }
}