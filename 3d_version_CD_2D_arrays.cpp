#include<iostream>
#include<fstream>
#include<cmath>
#include<string.h>
#include<iterator>
#include<array>

using namespace std;

const double rey= 100.0; // length of container= 1unit, rho= 10^3 SI units, viscosity= 8.9x 10^-4
const double length= 1.0;// so that reynolds number is equal to 100 and peclet number less than two
const double dt= 0.005;// from peclet no. less than 1
const int N= 100;
const double rho= 1.0;
double tolerance=pow(10,-4);
const double NT=1000000;
double max_err=1.0;

// defining variables
double u_inf= 1.0;

double p[N*N][N],p_error[N*N][N]={0.0};
double u[N*N][N],u_n_1[N*N][N]={0.0};
double v[N*N][N],v_n_1[N*N][N]={0.0};
double w[N*N][N],w_n_1[N*N][N]={0.0};

double a(int ix, int iy, int n);
double a_w(int ix, int iy, int n);
double a_e(int ix, int iy, int n);
double a_s(int ix, int iy, int n);
double a_n(int ix, int iy, int n);
double a_b(int ix, int iy, int n);
double a_t(int ix, int iy, int n);
double p_e(int ix, int iy);
double p_s(int ix, int iy);
double p_n(int ix, int iy);
double p_w(int ix, int iy);
double p_t(int ix, int iy);
double p_b(int ix, int iy);
double boundary(int ix, int iy);

double dx= length/(N-1);
double viscosity= (rho*u_inf*length)/rey ;
const double D= viscosity/dx;
//under relaxation factors for u,v and p
const double ur_u=0.5;
const double ur_v=0.5;
const double ur_w=0.5;
const double ur_p=0.5;

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

        for(int ix=0;ix<N*N;ix++){
            p[ix][N-1]=0.0;
        }

       // updating pressure with under-relaxation factor
            for(int ix=0;ix<N*N;ix++){
                for(int iy=0;iy<N;iy++){
                    p[ix][iy]= p[ix][iy]+ ur_p*p_error[ix][iy];// here '0.5' is the under relaxation  factor
               }
            }


      // update the boundary conditions

        //side of box
        for(int iy=1;iy<N+1;iy++){
            for(int ix=0;ix<N;ix++){
                p[(N)*iy-1][ix]=p[(N)*iy-2][ix];//north side of cavity w.r.t to velocity in +ve x- direction
                p[ix+(iy-1)*N][N-1]=p[ix+(iy-1)*N][N-2];// lid of cavity i.e. top of cavity
                p[(iy-1)+(N-1)*N][ix]=p[(iy-1)+(N-2)*N][ix];// east side of cavity
                p[(iy-1)*N][ix]=p[(iy-1)*N+1][ix];//south side of cavity
                p[ix+(iy-1)*N][0]=p[ix+(iy-1)*N][1];// bottom of cavity
                p[iy-1][ix]=p[(iy-1)+N][ix];// west side of cavity
            }
        }



        // update velocity values using correction equations
        // in the correction equations below, values are not updated in boundary planes: ix=0:N-1 and ix=N*N-N:N*N-1
         for(int ix=N;ix<N*N-N;ix++){
             for(int iy=1;iy<N;iy++){
                u[ix][iy]= u[ix][iy]+ (dx*dx/a(ix,iy,1))*(p_error[ix-1][iy]-p_error[ix][iy]);//updating velocity values without using under relaxation factor
                u[ix][iy]= ur_u*u[ix][iy]+ (1-ur_u)*u_n_1[ix][iy];// updating velocity values using under relaxation factor
            }
         }

        // update velocity values using correction equations
         for(int ix=N;ix<N*N-N;ix++){
             for(int iy=1;iy<N;iy++){
                v[ix][iy]= v[ix][iy]+ (dx*dx/a(ix,iy,2))*(p_error[ix][iy-1]-p_error[ix][iy]);//updating velocity values without using under relaxation factor
                v[ix][iy]=ur_v*v[ix][iy]+ (1-ur_v)*v_n_1[ix][iy];// updating velocity values using under relaxation factor
            }
         }

        // update velocity values using correction equations
        for(int ix=N;ix<N*N-N;ix++){
            for(int iy=1;iy<N;iy++){
                w[ix][iy]= w[ix][iy]+ (dx*dx/a(ix,iy,3))*(p_error[ix-N][iy]-p_error[ix][iy]);//updating velocity values without using under relaxation factor
                w[ix][iy]=ur_w*w[ix][iy]+ (1-ur_w)*w_n_1[ix][iy];// updating velocity values using under relaxation factor
            }
         }

         //boundary conditions

        for(int ix=1;ix<N+1;ix++){
            for(int iy=0;iy<N;iy++){
                // u-velocity
                u[N*(ix-1)+iy][N-1]=u_inf;//lid
                u[N*(ix-1)+iy][0]=0.0;//bottom
                u[ix-1][iy]=0.0;//west side of cavity
                u[(ix-1)+N*N-N][iy]=0.0;//east side of cavity
                u[(ix-1)*N][iy]=0.0;//south side of cavity
                u[N*ix-1][iy]=0.0;//north side of cavity


                //v-velocity
                v[N*(ix-1)+iy][N-1]=0.0;//lid
                v[N*(ix-1)+iy][0]=0.0;//bottom
                v[ix-1][iy]=0.0;//west side of cavity
                v[(ix-1)+N*N-N][iy]=0.0;//east side of cavity
                v[(ix-1)*N][iy]=0.0;//south side of cavity
                v[N*ix-1][iy]=0.0;//north side of cavity

                //w-velocity
                w[N*(ix-1)+iy][N-1]=0.0;//lid
                w[N*(ix-1)+iy][0]=0.0;//bottom
                w[ix-1][iy]=0.0;//west side of cavity
                w[(ix-1)+N*N-N][iy]=0.0;//east side of cavity
                w[(ix-1)*N][iy]=0.0;//south side of cavity
                w[N*ix-1][iy]=0.0;//north side of cavity

            }
        }


        //checking convergence

        double maximum1=0.0;
        double maximum2=0.0;
        double maximum3=0.0;
        double maximum4=0.0;
        double maximum=0.0;

        //error in x-momentum equation
            for(int ix=0;ix<N*N;ix++){
                for(int iy=0;iy<N;iy++){
                    if((abs(u_n_1[ix][iy]-u[ix][iy])/dt)>maximum1){
                        maximum1= abs(u_n_1[ix][iy]-u[ix][iy])/dt;
                    }
                }
            }


        //error in y-momentum equation

            for(int ix=0;ix<N*N;ix++){
                for(int iy=0;iy<N;iy++){
                    if(abs((v_n_1[ix][iy]-v[ix][iy])/dt)>maximum2){
                        maximum2= abs(v_n_1[ix][iy]-v[ix][iy])/dt;
                    }
                }
            }

        //error in z-momentum equation
            for(int ix=0;ix<N*N;ix++){
                for(int iy=0;iy<N;iy++){
                    if(abs((w_n_1[ix][iy]-w[ix][iy])/dt)>maximum3){
                        maximum3= abs(w_n_1[ix][iy]-w[ix][iy])/dt;
                    }
                    if(abs(p[ix][iy])>maximum4)
                        maximum4=abs(p[ix][iy]);
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
         ofstream fout17;
         ofstream fout18;
         ofstream fout22;



            fout12.open("pressure.dat");
            fout17.open("x_vel.dat");
            fout18.open("y_vel.dat");
            fout22.open("z_vel.dat");

            for(int iz=0;iz<N;iz++){
                for(int iy=0;iy<N;iy++){
                    for(int ix=0;ix<N;ix++){
                            fout12<<p[ix+N*iz][iy]<<endl;// pressure field
                    }
                }
            }



            for(int iz=0;iz<N;iz++){
                for(int iy=0;iy<N;iy++){
                    for(int ix=0;ix<N;ix++){
                        fout17<<u[ix+N*iz][iy]<<endl;
                    }
                }
            }

            for(int iz=0;iz<N;iz++){
                for(int iy=0;iy<N;iy++){
                    for(int ix=0;ix<N;ix++){
                        fout18<<v[ix+N*iz][iy]<<endl;
                    }
                }
            }

            for(int iz=0;iz<N;iz++){
                for(int iy=0;iy<N;iy++){
                    for(int ix=0;ix<N;ix++){
                        fout22<<w[ix+N*iz][iy]<<endl;

                    }
                }
            }

            fout12.close();
            fout17.close();
            fout18.close();
            fout22.close();
            fout15.close();
}

double boundary(int ix, int iy){

    double ans=0.0;

    return ans;
}

// momentum equation coefficients

double a_n(int ix, int iy,int n){

    double ans=0.0;

    if(n==1)
     ans= dx*dx*(D -rho*(v_n_1[ix-1][iy+1]+v_n_1[ix][iy+1])/4.0);
     else if(n==2)
        ans= dx*dx*(D -rho*(v_n_1[ix][iy]+v_n_1[ix][iy+1])/4.0);
        else if(n==3)
            ans= dx*dx*(D -rho*(v_n_1[ix][iy+1]+v_n_1[ix-N][iy+1])/4.0);

    return ans;
}


double a_s(int ix, int iy,int n){

    double ans=0.0;

    if(n==1)
     ans= dx*dx*(D + rho*(v_n_1[ix-1][iy]+v_n_1[ix][iy])/4.0);
     else if(n==2)
        ans= dx*dx*(D + rho*(v_n_1[ix][iy-1]+v_n_1[ix][iy])/4.0);
        else if(n==3)
            ans= dx*dx*(D + rho*(v_n_1[ix][iy]+v_n_1[ix-N][iy])/4.0);

    return ans;
}


double a_e(int ix, int iy,int n){

    double ans=0.0;

    if(n==1)
     ans= dx*dx*(D -rho*(u_n_1[ix][iy]+u_n_1[ix+1][iy])/4.0);
     else if(n==2)
        ans= dx*dx*(D -rho*(u_n_1[ix+1][iy]+u_n_1[ix+1][iy-1])/4.0);
        else if(n==3)
            ans= dx*dx*(D -rho*(u_n_1[ix+1][iy]+u_n_1[ix+1-N][iy])/4.0);

    return ans;
}

double a_w(int ix, int iy,int n){

    double ans=0.0;

    if(n==1)
     ans= dx*dx*(D + rho*(u_n_1[ix-1][iy]+u_n_1[ix][iy])/4.0);
     else if(n==2)
        ans= dx*dx*(D + rho*(u_n_1[ix][iy]+u_n_1[ix][iy-1])/4.0 );
        else if(n==3)
            ans= dx*dx*(D + rho*(u_n_1[ix][iy]+u_n_1[ix-N][iy])/4.0);

    return ans;
}

double a_t(int ix, int iy,int n){

    double ans=0.0;

    if(ix+N>N*N-1){
        if(n==1)
         ans= dx*dx*(D -rho*(0.0+0.0)/4.0);
         else if(n==2)
            ans= dx*dx*(D -rho*(0.0+0.0)/4.0);
            else if(n==3)
                ans= dx*dx*(D -rho*(w_n_1[ix][iy]+0.0)/4.0);
    }else{
        if(n==1)
         ans= dx*dx*(D -rho*(w_n_1[ix+N][iy]+w_n_1[ix-1+N][iy])/4.0);
         else if(n==2)
            ans= dx*dx*(D -rho*(w_n_1[ix+N][iy]+w_n_1[ix+N][iy-1])/4.0);
            else if(n==3)
                ans= dx*dx*(D -rho*(w_n_1[ix][iy]+w_n_1[ix+N][iy])/4.0);
    }
    return ans;
}

double a_b(int ix, int iy,int n){

   double ans=0.0;

    if(n==1)
     ans= dx*dx*(D + rho*(w_n_1[ix-1][iy]+w_n_1[ix][iy])/4.0);
     else if(n==2)
        ans= dx*dx*(D +rho*(w_n_1[ix][iy]+w_n_1[ix][iy-1])/4.0);
        else if(n==3)
            ans= dx*dx*(D + rho*(w_n_1[ix-N][iy]+w_n_1[ix][iy])/4.0);

    return ans;
}

double a(int ix, int iy,int n){

    double ans=0.0;

    if(n==1)
        ans= 12.0*dx*dx*D -(a_w(ix,iy,1)+a_e(ix,iy,1)+a_n(ix,iy,1)+a_s(ix,iy,1)+a_t(ix,iy,1)+a_b(ix,iy,1));
        else if(n==2)
            ans= 12.0*dx*dx*D -(a_w(ix,iy,2)+a_e(ix,iy,2)+a_n(ix,iy,2)+a_s(ix,iy,2)+a_t(ix,iy,2)+a_b(ix,iy,2));
            else if(n==3)
            ans= 12.0*dx*dx*D -(a_w(ix,iy,3)+a_e(ix,iy,3)+a_n(ix,iy,3)+a_s(ix,iy,3)+a_t(ix,iy,3)+a_b(ix,iy,3));

    return ans;
}

// pressure correction equation coefficients

double p_e(int ix, int iy){

    double answer=0.0;

    answer= rho*dx*dx*dx*dx/a(ix+1,iy,1);

    return answer;

}

double p_n(int ix, int iy){

    double answer=0.0;


    answer= rho*dx*dx*dx*dx/a(ix,iy+1,2);


    return answer;

}

double p_w(int ix, int iy){

    double answer=0.0;

    answer= rho*dx*dx*dx*dx/a(ix,iy,1);

    return answer;

}

double p_s(int ix, int iy){

    double answer=0.0;

    answer= rho*dx*dx*dx*dx/a(ix,iy,2);


    return answer;

}

double p_t(int ix, int iy){

    double answer=0.0;

    answer= rho*dx*dx*dx*dx/a(ix+N,iy,3);


    return answer;

}

double p_b(int ix, int iy){

    double answer=0.0;

    answer= rho*dx*dx*dx*dx/a(ix,iy,3);


    return answer;

}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



void gs_vel(void){

        for(int ix=N;ix<N*N-N;ix++){
            for(int iy=1;iy<N-1;iy++){

                 if(ix%N==0||(ix+1)%(N)==0){
                    continue;
                 }else{
                    u[ix][iy]= (dx*dx*(p[ix-1][iy]-p[ix][iy]) + u[ix+N][iy]*(a_t(ix,iy,1))+ u[ix-N][iy]*(a_b(ix,iy,1))+ u[ix+1][iy]*(a_e(ix,iy,1)) + u[ix][iy+1]*(a_n(ix,iy,1)) + u[ix-1][iy]*(a_w(ix,iy,1)) + u[ix][iy-1]*(a_s(ix,iy,1)))/(a(ix,iy,1));
                    v[ix][iy]= (dx*dx*(p[ix][iy-1]-p[ix][iy]) + v[ix+1][iy]*(a_e(ix,iy,2))+ v[ix][iy+1]*(a_n(ix,iy,2))+ v[ix-1][iy]*(a_w(ix,iy,2)) + v[ix][iy-1]*(a_s(ix,iy,2)) + v[ix+N][iy]*(a_t(ix,iy,2)) + v[ix-N][iy]*(a_b(ix,iy,2)))/(a(ix,iy,2));
                    w[ix][iy]= (dx*dx*(p[ix-N][iy]-p[ix][iy]) + w[ix+1][iy]*(a_e(ix,iy,3))+ w[ix][iy+1]*(a_n(ix,iy,3))+ w[ix-1][iy]*(a_w(ix,iy,3)) + w[ix][iy-1]*(a_s(ix,iy,3)) + w[ix+N][iy]*(a_t(ix,iy,3)) + w[ix-N][iy]*(a_b(ix,iy,3)))/(a(ix,iy,3));
                 }
                if(u[ix][iy]!=u[ix][iy]){
                   cout<<"Pressure grad: "<<p[ix-1][iy]-p[ix][iy]<<endl;
                   cout<<"U_top: "<< u[ix+N][iy]<<" U_bottom: "<<u[ix-N][iy]<<" U_north: "<<u[ix][iy+1]<<" U_south: "<<u[ix][iy-1]<<" U_east: "<<u[ix+1][iy]<<" U_west: "<<u[ix-1][iy]<<endl;
                   cout<<" Top: "<<a_t(ix,iy,1)<<" Bottom: "<<a_b(ix,iy,1)<<" North: "<<a_n(ix,iy,1)<<" South: "<<a_s(ix,iy,1)<<" East: "<<a_e(ix,iy,1)<<" West: "<<a_w(ix,iy,1)<<endl;
                   cout<<" Total: "<<a(ix,iy,1)<<endl;
                   cout<<" Location: "<<ix<<","<<iy<<endl;
                   break;
                }
            }
        }

}


void gs_p(void){

        for(int ix=N;ix<N*N-N;ix++){
            for(int iy=1;iy<N-1;iy++){

                if(ix%N==0||(ix+1)%(N)==0){
                   continue;
                }else{
                p_error[ix][iy]= ((rho*dx*dx*(u[ix][iy]-u[ix+1][iy]+v[ix][iy]-v[ix][iy+1]+w[ix][iy]-w[ix+N][iy]))+(p_error[ix+1][iy]*p_e(ix,iy) + p_error[ix-1][iy]*p_w(ix,iy) + p_error[ix][iy-1]*p_s(ix,iy) + p_error[ix][iy+1]*p_n(ix,iy)+ p_error[ix+N][iy]*p_t(ix,iy)+ p_error[ix-N][iy]*p_b(ix,iy)) + boundary(ix,iy))/(p_n(ix,iy)+p_w(ix,iy)+p_e(ix,iy)+p_s(ix,iy)+p_t(ix,iy)+p_b(ix,iy));
                }
                if(p_error[ix][iy]!=p_error[ix][iy]){
                   cout<<"Vel. correction: "<<u[ix][iy]-u[ix+1][iy]+v[ix][iy]-v[ix][iy+1]+w[ix][iy]-w[ix+N][iy]<<endl;
                   cout<<"P_top: "<< p_error[ix+N][iy]<<" P_bottom: "<<p_error[ix-N][iy]<<" P_north: "<<p_error[ix][iy+1]<<" P_south: "<<p_error[ix][iy-1]<<" P_east: "<<p_error[ix+1][iy]<<" P_west: "<<p_error[ix-1][iy]<<endl;
                   cout<<" Top: "<<p_t(ix,iy)<<" Bottom: "<<p_b(ix,iy)<<" North: "<<p_n(ix,iy)<<" South: "<<p_s(ix,iy)<<" East: "<<p_e(ix,iy)<<" West: "<<p_w(ix,iy)<<endl;
                   cout<<" Momentum eq coefficients: west,east,north,south,top,bottom "<<a_w(ix+N,iy,3)<<","<<a_e(ix+N,iy,3)<<","<<a_n(ix+N,iy,3)<<","<<a_s(ix+N,iy,3)<<","<<a_t(ix+N,iy,3)<<","<<a_b(ix+N,iy,3)<<endl;
                   cout<<" Location: "<<ix<<","<<iy<<endl;
                   break;
                }
            }
        }

}



