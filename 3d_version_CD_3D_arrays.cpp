#include<iostream>
#include<fstream>
#include<cmath>
#include<string.h>
#include<iterator>
#include<array>

using namespace std;

const double rey= 400.0; // length of container= 1unit, rho= 10^3 SI units, viscosity= 8.9x 10^-4
const double length= 1.0;// so that reynolds number is equal to 100 and peclet number less than two
const double dt= 0.0005;// from peclet no. less than 1
const int N= 100;
const double rho= 1.0;
double tolerance=pow(10,-4);
const double NT=100000;
double max_err=1.0;

// defining variables
double u_infinity= 1.0;

double p[N][N][N],p_error[N][N][N]={0.0};
double stream_function[N][N], vorticity[N][N][N]={0.0};
double u[N][N][N],u_n_1[N][N][N]={0.0};
double v[N][N][N],v_n_1[N][N][N]={0.0};
double w[N][N][N],w_n_1[N][N][N]={0.0};

double a(int ix, int iy,int iz, int n);
double a_w(int ix, int iy,int iz, int n);
double a_e(int ix, int iy,int iz, int n);
double a_s(int ix, int iy,int iz, int n);
double a_n(int ix, int iy,int iz, int n);
double a_b(int ix, int iy,int iz, int n);
double a_t(int ix, int iy,int iz, int n);
double p_e(int ix, int iy,int iz);
double p_s(int ix, int iy,int iz);
double p_n(int ix, int iy,int iz);
double p_w(int ix, int iy,int iz);
double p_t(int ix, int iy,int iz);
double p_b(int ix, int iy,int iz);
double boundary(int ix, int iy,int iz);

double dx= length/(N-1);
double viscosity= (rho*u_infinity*length)/rey ;
const double D= viscosity/dx;
//under relaxation factors for u,v and p
const double ur_u=1;
const double ur_v=1;
const double ur_w=1;
const double ur_p=1;

void gs_vel(void);
void gs_p(void);


int main(){

    // WRITING VALUES TO TOLERANCE FILE

    ofstream fout15;
    fout15.open("tolerance.dat");


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 //                       SIMPLE algorithm
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    for(int it=NT;it>0;it--){

        cout<<"Iteration number: "<<NT- it<<endl;

            // copying values of velocities from previous iteration
        std::copy(std::begin(u),std::end(u),std::begin(u_n_1));
        std::copy(std::begin(v),std::end(v),std::begin(v_n_1));
        std::copy(std::begin(w),std::end(w),std::begin(w_n_1));

        // solving the system of linear equations for u and v momentum equations and pressure correction equation
        gs_vel();

        memset(p_error,0.0,sizeof(p_error));


        gs_p();



       // updating pressure with under-relaxation factor
        for(int iz=1;iz<N-1;iz++){
            for(int ix=1;ix<N-1;ix++){
                for(int iy=1;iy<N-1;iy++){
                    p[ix][iy][iz]= p[ix][iy][iz]+ ur_p*p_error[ix][iy][iz];// here '0.5' is the under relaxation  factor
               }
            }
        }

      // update the boundary conditions

        //side of box
        for(int iy=0;iy<N;iy++){
            for(int ix=0;ix<N;ix++){
                p[N-1][ix][iy]=p[N-2][ix][iy];
                p[ix][N-1][iy]=p[ix][N-2][iy];
                p[iy][ix][N-1]=p[iy][ix][N-2];
                p[0][ix][iy]=p[1][ix][iy];
                p[ix][0][iy]=p[ix][1][iy];
                p[ix][iy][0]=p[ix][iy][1];
            }
        }



        // update velocity values using correction equations
        for(int iz=N-2;iz>0;iz--){
               for(int ix=N-2;ix>0;ix--){
                    for(int iy=N-2;iy>0;iy--){
                        u[ix][iy][iz]= u[ix][iy][iz]+ (dx*dx/a(ix,iy,iz,1))*(p_error[ix-1][iy][iz]-p_error[ix][iy][iz]);//updating velocity values without using under relaxation factor
                        u[ix][iy][iz]=ur_u*u[ix][iy][iz]+ (1-ur_u)*u_n_1[ix][iy][iz];// updating velocity values using under relaxation factor


                    }
                }
        }

        // update velocity values using correction equations
        for(int iz=N-2;iz>0;iz--){
               for(int ix=N-2;ix>0;ix--){
                    for(int iy=N-2;iy>0;iy--){
                        v[ix][iy][iz]= v[ix][iy][iz]+ (dx*dx/a(ix,iy,iz,2))*(p_error[ix][iy-1][iz]-p_error[ix][iy][iz]);//updating velocity values without using under relaxation factor
                        v[ix][iy][iz]=ur_v*v[ix][iy][iz]+ (1-ur_v)*v_n_1[ix][iy][iz];// updating velocity values using under relaxation factor
                    }
                }
        }

        // update velocity values using correction equations
        for(int iz=N-2;iz>0;iz--){
               for(int ix=N-2;ix>0;ix--){
                    for(int iy=N-2;iy>0;iy--){
                        w[ix][iy][iz]= w[ix][iy][iz]+ (dx*dx/a(ix,iy,iz,3))*(p_error[ix][iy][iz-1]-p_error[ix][iy][iz]);//updating velocity values without using under relaxation factor
                        w[ix][iy][iz]=ur_w*w[ix][iy][iz]+ (1-ur_w)*w_n_1[ix][iy][iz];// updating velocity values using under relaxation factor
                    }
                }
        }

         //boundary conditions

        for(int iy=0;iy<N;iy++){
            for(int ix=0;ix<N;ix++){
                // u-velocity
                u[ix][N-1][iy]=u_infinity;//lid
                u[ix][0][iy]=0.0;//bottom
                u[ix][iy][0]=0.0;//side
                u[ix][iy][N-1]=0.0;//side

                //v-velocity
                v[0][ix][iy]=0.0;//side
                v[N-1][ix][iy]=0.0;//side
                v[iy][ix][0]=0.0;//side
                v[iy][ix][N-1]=0.0;//side

                //w-velocity
                w[0][iy][ix]=0.0;//side
                w[N-1][iy][ix]=0.0;//side
                w[iy][0][ix]=0.0;//bottom
                w[iy][N-1][ix]=0.0;//lid
            }
        }

        for(int iy=0;iy<N;iy++){
            for(int ix=0;ix<N;ix++){
                // u-velocity
                u[0][ix][iy]=0.0;//side
                u[N-1][ix][iy]=0.0;//side

                //v-velocity
                v[ix][0][iy]=0.0;//bottom
                v[ix][N-1][iy]=0.0;//lid


                //w-velocity
                w[ix][iy][0]=0.0;//side
                w[ix][iy][N-1]=0.0;//side
            }
        }

        //checking convergence

        double maximum1=0.0;
        double maximum2=0.0;
        double maximum3=0.0;
        double maximum4=0.0;
        double maximum=0.0;

        //error in x-momentum equation
        for(int iz=0;iz<N;iz++){
            for(int ix=0;ix<N;ix++){
                for(int iy=0;iy<N;iy++){
                    if((abs(u_n_1[ix][iy][iz]-u[ix][iy][iz])/dt)>maximum1){
                        maximum1= abs(u_n_1[ix][iy][iz]-u[ix][iy][iz])/dt;
                    }
                }
            }
        }

        //error in y-momentum equation
        for(int iz=0;iz<N;iz++){
            for(int ix=0;ix<N;ix++){
                for(int iy=0;iy<N;iy++){
                    if(abs((v_n_1[ix][iy][iz]-v[ix][iy][iz])/dt)>maximum2){
                        maximum2= abs(v_n_1[ix][iy][iz]-v[ix][iy][iz])/dt;
                    }
                }
            }
        }

        //error in z-momentum equation
        for(int iz=0;iz<N;iz++){
            for(int ix=0;ix<N;ix++){
                for(int iy=0;iy<N;iy++){
                    if(abs((w_n_1[ix][iy][iz]-w[ix][iy][iz])/dt)>maximum3){
                        maximum3= abs(w_n_1[ix][iy][iz]-w[ix][iy][iz])/dt;
                    }
                    if(abs(p[ix][iy][iz])>maximum4)
                        maximum4=abs(p[ix][iy][iz]);
                }
            }
        }

        maximum= maximum1+maximum2+maximum3;


        if(maximum<tolerance*tolerance*max_err)
            break;
	if(it==NT)
	   max_err=maximum;

        cout<<"Value of error, u: "<<maximum1<<",v: "<<maximum2<<",w: "<<maximum3<<",p: "<<maximum4<<endl;

        fout15<<NT-it<<" "<<maximum<<endl;


}




 cout<<"Writing values to file."<<endl;
            //write values to a file

         ofstream fout12;
         ofstream fout13;
         ofstream fout14;
         ofstream fout17;
         ofstream fout18;
         ofstream fout22;



            fout12.open("pressure.dat");
            fout13.open("mid_line_x_Velocity.dat");
            fout14.open("mid_line_y_velocity.dat");
            fout17.open("x_vel.dat");
            fout18.open("y_vel.dat");
            fout22.open("z_vel.dat");

            for(int iz=0;iz<N;iz++){
                for(int iy=0;iy<N;iy++){
                    for(int ix=0;ix<N;ix++){
                        fout12<<ix*dx<<","<<iy*dx<<","<<iz*dx<<","<<p[ix][iy][iz]<<endl;// pressure field
                    }
                }
            }

            for(int i=0;i<N;i++){
                fout13<<u[16][i][16]/u_infinity<<" "<<(dx*i)/length<<endl;// mid line u- velocity
            }

            for(int i=0;i<N;i++){
               fout14<<dx*i/length<<" "<<v[i][16][16]/u_infinity<<endl;// mid line v- velocity

            }


            for(int iz=0;iz<N;iz++){
                for(int iy=0;iy<N;iy++){
                    for(int ix=0;ix<N;ix++){
                        fout17<<dx*ix/length<<","<<dx*iy/length<<","<<dx*iz/length<<","<<u[ix][iy][iz]<<endl;
                    }
                }
            }

            for(int iz=0;iz<N;iz++){
                for(int iy=0;iy<N;iy++){
                    for(int ix=0;ix<N;ix++){
                        fout18<<dx*ix/length<<","<<dx*iy/length<<","<<dx*iz/length<<","<<v[ix][iy][iz]<<endl;
                    }
                }
            }

            for(int iz=0;iz<N;iz++){
                for(int iy=0;iy<N;iy++){
                    for(int ix=0;ix<N;ix++){
                        fout22<<dx*ix/length<<","<<dx*iy/length<<","<<dx*iz/length<<","<<w[ix][iy][iz]<<endl;

                    }
                }
            }

            fout12.close();
            fout13.close();
            fout14.close();
            fout17.close();
            fout18.close();
            fout22.close();
            fout15.close();
}

double boundary(int ix, int iy,int iz){

    double ans=0.0;

    return ans;
}

// momentum equation coefficients

double a_n(int ix, int iy,int iz,int n){

    double ans=0.0;

    if(n==1)
     ans= dx*dx*(D -rho*(v_n_1[ix-1][iy+1][iz]+v_n_1[ix][iy+1][iz])/4.0);
     else if(n==2)
        ans= dx*dx*(D -rho*(v_n_1[ix][iy][iz]+v_n_1[ix][iy+1][iz])/4.0);
        else if(n==3)
            ans= dx*dx*(D -rho*(v_n_1[ix][iy+1][iz]+v_n_1[ix][iy+1][iz-1])/4.0);

    return ans;
}


double a_s(int ix, int iy,int iz,int n){

    double ans=0.0;

    if(n==1)
     ans= dx*dx*(D + rho*(v_n_1[ix-1][iy][iz]+v_n_1[ix][iy][iz])/4.0);
     else if(n==2)
        ans= dx*dx*(D + rho*(v_n_1[ix][iy-1][iz]+v_n_1[ix][iy][iz])/4.0);
        else if(n==3)
            ans= dx*dx*(D + rho*(v_n_1[ix][iy][iz]+v_n_1[ix][iy][iz-1])/4.0);

    return ans;
}


double a_e(int ix, int iy,int iz,int n){

    double ans=0.0;

    if(n==1)
     ans= dx*dx*(D -rho*(u_n_1[ix][iy][iz]+u_n_1[ix+1][iy][iz])/4.0);
     else if(n==2)
        ans= dx*dx*(D -rho*(u_n_1[ix+1][iy][iz]+u_n_1[ix+1][iy-1][iz])/4.0);
        else if(n==3)
            ans= dx*dx*(D -rho*(u_n_1[ix+1][iy][iz]+u_n_1[ix+1][iy][iz-1])/4.0);

    return ans;
}

double a_w(int ix, int iy,int iz,int n){

    double ans=0.0;

    if(n==1)
     ans= dx*dx*(D + rho*(u_n_1[ix-1][iy][iz]+u_n_1[ix][iy][iz])/4.0);
     else if(n==2)
        ans= dx*dx*(D + rho*(u_n_1[ix][iy][iz]+u_n_1[ix][iy-1][iz])/4.0 );
        else if(n==3)
            ans= dx*dx*(D + rho*(u_n_1[ix][iy][iz]+u_n_1[ix][iy][iz-1])/4.0);

    return ans;
}

double a_t(int ix, int iy,int iz,int n){

    double ans=0.0;

    if(n==1)
     ans= dx*dx*(D -rho*(w_n_1[ix][iy][iz+1]+w_n_1[ix-1][iy][iz+1])/4.0);
     else if(n==2)
        ans= dx*dx*(D -rho*(w_n_1[ix][iy][iz+1]+w_n_1[ix][iy-1][iz+1])/4.0);
        else if(n==3)
            ans= dx*dx*(D -rho*(w_n_1[ix][iy][iz]+w_n_1[ix][iy][iz+1])/4.0);

    return ans;
}

double a_b(int ix, int iy,int iz,int n){

   double ans=0.0;

    if(n==1)
     ans= dx*dx*(D + rho*(w_n_1[ix-1][iy][iz]+w_n_1[ix][iy][iz])/4.0);
     else if(n==2)
        ans= dx*dx*(D +rho*(w_n_1[ix][iy][iz]+w_n_1[ix][iy-1][iz])/4.0);
        else if(n==3)
            ans= dx*dx*(D + rho*(w_n_1[ix][iy][iz-1]+w_n_1[ix][iy][iz])/4.0);

    return ans;
}

double a(int ix, int iy,int iz,int n){

    double ans=0.0;

    if(n==1)
        ans= 12.0*dx*dx*D -(a_w(ix,iy,iz,1)+a_e(ix,iy,iz,1)+a_n(ix,iy,iz,1)+a_s(ix,iy,iz,1)+a_t(ix,iy,iz,1)+a_b(ix,iy,iz,1));
        else if(n==2)
            ans= 12.0*dx*dx*D -(a_w(ix,iy,iz,2)+a_e(ix,iy,iz,2)+a_n(ix,iy,iz,2)+a_s(ix,iy,iz,2)+a_t(ix,iy,iz,2)+a_b(ix,iy,iz,2));
            else if(n==3)
            ans= 12.0*dx*dx*D -(a_w(ix,iy,iz,3)+a_e(ix,iy,iz,3)+a_n(ix,iy,iz,3)+a_s(ix,iy,iz,3)+a_t(ix,iy,iz,3)+a_b(ix,iy,iz,3));

    return ans;
}

// pressure correction equation coefficients

double p_e(int ix, int iy,int iz){

    double answer=0.0;

    answer= rho*dx*dx*dx*dx/a(ix+1,iy,iz,1);

    return answer;

}

double p_n(int ix, int iy,int iz){

    double answer=0.0;


    answer= rho*dx*dx*dx*dx/a(ix,iy+1,iz,2);


    return answer;

}

double p_w(int ix, int iy,int iz){

    double answer=0.0;

    answer= rho*dx*dx*dx*dx/a(ix,iy,iz,1);

    return answer;

}

double p_s(int ix, int iy,int iz){

    double answer=0.0;

    answer= rho*dx*dx*dx*dx/a(ix,iy,iz,2);


    return answer;

}

double p_t(int ix, int iy,int iz){

    double answer=0.0;

    answer= rho*dx*dx*dx*dx/a(ix,iy,iz+1,3);


    return answer;

}

double p_b(int ix, int iy,int iz){

    double answer=0.0;

    answer= rho*dx*dx*dx*dx/a(ix,iy,iz,3);


    return answer;

}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



void gs_vel(void){

    for(int iz=1;iz<N-1;iz++){
        for(int ix=1;ix<N-1;ix++){
            for(int iy=1;iy<N-1;iy++){
                u[ix][iy][iz]= (dx*dx*((p[ix-1][iy][iz]-p[ix][iy][iz])) + u[ix][iy][iz+1]*(a_t(ix,iy,iz,1))+ u[ix][iy][iz-1]*(a_b(ix,iy,iz,1))+ u[ix+1][iy][iz]*(a_e(ix,iy,iz,1)) + u[ix][iy+1][iz]*(a_n(ix,iy,iz,1)) + u[ix-1][iy][iz]*(a_w(ix,iy,iz,1)) + u[ix][iy-1][iz]*(a_s(ix,iy,iz,1)))/(a(ix,iy,iz,1));
                v[ix][iy][iz]= (dx*dx*((p[ix][iy-1][iz]-p[ix][iy][iz])) + v[ix+1][iy][iz]*(a_e(ix,iy,iz,2)) + v[ix][iy+1][iz]*(a_n(ix,iy,iz,2)) + v[ix-1][iy][iz]*(a_w(ix,iy,iz,2)) + v[ix][iy-1][iz]*(a_s(ix,iy,iz,2))+ v[ix][iy][iz+1]*(a_t(ix,iy,iz,2))+ v[ix][iy][iz-1]*(a_b(ix,iy,iz,2)))/(a(ix,iy,iz,2));
                w[ix][iy][iz]= (dx*dx*((p[ix][iy][iz-1]-p[ix][iy][iz])) + w[ix+1][iy][iz]*(a_e(ix,iy,iz,3)) + w[ix][iy+1][iz]*(a_n(ix,iy,iz,3)) + w[ix-1][iy][iz]*(a_w(ix,iy,iz,3)) + w[ix][iy-1][iz]*(a_s(ix,iy,iz,3))+ w[ix][iy][iz+1]*(a_t(ix,iy,iz,3))+ w[ix][iy][iz-1]*(a_b(ix,iy,iz,3)))/(a(ix,iy,iz,3));
                   if(u[ix][iy][iz]!=u[ix][iy][iz]){

                      cout<<"Check u momentum equation"<<endl;

                   }
            }
        }
    }
}


void gs_p(void){

    for(int iz=1;iz<N-1;iz++){
        for(int ix=1;ix<N-1;ix++){
            for(int iy=1;iy<N-1;iy++){
                p_error[ix][iy][iz]= ((rho*dx*dx*(u[ix][iy][iz]-u[ix+1][iy][iz]+v[ix][iy][iz]-v[ix][iy+1][iz]+w[ix][iy][iz]-w[ix][iy][iz+1]))+(p_error[ix+1][iy][iz]*p_e(ix,iy,iz) + p_error[ix-1][iy][iz]*p_w(ix,iy,iz) + p_error[ix][iy-1][iz]*p_s(ix,iy,iz) + p_error[ix][iy+1][iz]*p_n(ix,iy,iz)+ p_error[ix][iy][iz+1]*p_t(ix,iy,iz)+ p_error[ix][iy][iz-1]*p_b(ix,iy,iz)) + boundary(ix,iy,iz))/(p_n(ix,iy,iz)+p_w(ix,iy,iz)+p_e(ix,iy,iz)+p_s(ix,iy,iz)+p_t(ix,iy,iz)+p_b(ix,iy,iz));

            }
        }
    }
}



