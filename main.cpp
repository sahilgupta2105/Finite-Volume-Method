#include<iostream>
#include<fstream>
#include<cmath>

using namespace std;

const double rey= 100.0; // length of container= 1unit, rho= 10^3 SI units, viscosity= 8.9x 10^-4
const double length= 1.0;// so that reynolds number is equal to 100 and peclet number less than two
const double dt= 0.0001;// from peclet no. less than 1
const int N= 33;
const double rho= 1.0;
double tolerance=pow(10,-6);
const double NT=300000;

// defining variables
double u_infinity= 1;

double p[N][N],p_error[N][N]={0.0};
double stream_function[N][N], vorticity[N][N]={0.0};
double u[N-1][N], un[N-1][N],u_n_1[N-1][N]={0.0};
double v[N][N-1], vn[N][N-1],v_n_1[N][N-1]={0.0};

double a(int ix, int iy,int z);
double a_w(int ix, int iy,int z);
double a_e(int ix, int iy,int z);
double a_s(int ix, int iy,int z);
double a_n(int ix, int iy,int z);
double p_e(int ix, int iy);
double p_s(int ix, int iy);
double p_n(int ix, int iy);
double p_w(int ix, int iy);
double boundary(int ix, int iy);

double dx= length/(N-1);
double viscosity= (rho*u_infinity*length)/rey ;
double D= viscosity/dx;

void gs_u(void);
void gs_v(void);
void gs_p(void);


int main(){



    // WRITING VALUES TO TOLERANCE FILE

    ofstream fout15;
    fout15.open("tolerance.dat");

    for(int i=0;i<N;i++){
      u[i][N-1]=u_infinity;// this value from calculation from Reynold's number =100.0
    }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 //                       SIMPLE algorithm
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    for(int it=NT;it>0;it--){

        cout<<"Iteration number: "<<NT- it<<endl;

            // copying values of U-velocity
        for(int ix=0;ix<N-1;ix++){
            for(int iy=0;iy<N;iy++){
                u_n_1[ix][iy]=u[ix][iy];
            }
       }

       // copying values of V-velocity

       for(int ix=0;ix<N;ix++){
        for(int iy=0;iy<N-1;iy++){
            v_n_1[ix][iy]=v[ix][iy];
        }
       }

        // solving the system of linear equations for u and v momentum equations and pressure correction equation

        gs_u();
        gs_v();

        for(int ix=0;ix<N;ix++){
            for(int iy=0;iy<N;iy++){
                p_error[ix][iy]=0.0;
            }
        }


        gs_p();

          // updating pressure with under-relaxation factor

        for(int ix=1;ix<N-1;ix++){
            for(int iy=1;iy<N-1;iy++){

                p[ix][iy]= p[ix][iy]+ 0.5*p_error[ix][iy];// here '0.5' is the under relaxation  factor

               }
        }
      // update the boundary conditions


        for(int ix=0;ix<N;ix++){
            p[N-1][ix]=p[N-2][ix];
        }
        for(int ix=0;ix<N;ix++){
            p[ix][N-1]=p[ix][N-2];
        }
         for(int ix=0;ix<N;ix++){
            p[0][ix]=p[1][ix];
        }
        for(int ix=0;ix<N;ix++){
            p[ix][0]=p[ix][1];
        }

        p[0][0]=0.0;


        // update velocity values using correction equations

        //updating u velocity





        for(int ix=N-2;ix>0;ix--){
            for(int iy=N-2;iy>0;iy--){

                    if(ix==65&&iy==127){
                        cout<<u[ix][iy]<<".."<<(dx*1.0/a(ix,iy,1))<<".."<<(p_error[ix-1][iy]-p_error[ix][iy])<<endl;
                    }
                    u[ix][iy]= u[ix][iy]+ (dx*1.0/a(ix,iy,1))*(p_error[ix-1][iy]-p_error[ix][iy]);//updating velocity values without using under relaxation factor
                    v[ix][iy]= v[ix][iy]+ (dx*1.0/a(ix,iy,2))*(p_error[ix][iy-1]-p_error[ix][iy]);//updating velocity values without using under relaxation factor
                    }
            }

        for(int ix=1;ix<N-1;ix++){
            for(int iy =1;iy<N-1;iy++){

                    u[ix][iy]=0.5*u[ix][iy]+ 0.5*u_n_1[ix][iy];// updating velocity values using under relaxation factor
                    v[ix][iy]=0.5*v[ix][iy]+ 0.5*v_n_1[ix][iy];// updating velocity values using under relaxation factor
            }
        }


         for(int ix=0;ix<N;ix++){
            u[ix][N-1]=u_infinity;
        }
        for(int ix=0;ix<N-1;ix++){
            u[0][ix]=0.0;
            u[N-1][ix]=0.0;
            u[ix][0]=0.0;
        }



        for(int iy=0; iy<N; iy++){
            v[iy][0]=0.0;
            v[iy][N-1]=0.0;
        }
        for(int iy=0; iy<N; iy++){
            v[N-1][iy]=0.0;
            v[0][iy]=0.0;
        }

        cout<<u[65][1]<<endl;
        cout<<u[65][127]<<endl;

        //checking convergence

        double maximum1=0.0;
        double maximum2=0.0;
        double maximum=0.0;

        for(int ix=0;ix<N-1;ix++){
            for(int iy=0;iy<N;iy++){

                if((abs(u_n_1[ix][iy]-u[ix][iy])/dt)>maximum1){
                    maximum1= abs(u_n_1[ix][iy]-u[ix][iy])/dt;
                   // cout<<ix<<","<<iy<<endl;
                }
            }
        }
        for(int ix=0;ix<N;ix++){
            for(int iy=0;iy<N-1;iy++){

           if(abs((v_n_1[ix][iy]-v[ix][iy])/dt)>maximum2){
                    maximum2= abs(v_n_1[ix][iy]-v[ix][iy])/dt;
                }
            }
        }

        if(maximum1<maximum2){
            maximum = maximum2;
        }
        else if(maximum2<maximum1){
            maximum = maximum1;
        }
        if(maximum<tolerance)
            break;

        cout<<"Value of error:"<<maximum<<endl;

        fout15<<NT-it<<" "<<maximum<<endl;


}



for(int ix=1;ix<N;ix++){
    for(int iy=1;iy<N;iy++){
            if(ix<=N/2){
        stream_function[ix][iy]=  -(v[ix-1][iy]*dx) + stream_function[ix-1][iy];
            }

            else{
        stream_function[ix][iy]= (u[ix][iy-1]*dx) + stream_function[ix][iy-1];
            }
    }
}



for(int ix=1;ix<N-1;ix++){
    for(int iy=1;iy<N-1;iy++){
        vorticity[ix][iy]= (v[ix+1][iy] -v[ix][iy])/dx - (u[ix][iy+1] -u[ix][iy])/dx;
    }
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
            fout17.open("streamfunction.dat");
            fout18.open("velocity.dat");
            fout22.open("vorticity.dat");

            for(int ix=0;ix<N;ix++){
                for(int iy=0;iy<N;iy++){
                    fout12<<ix*dx<<","<<iy*dx<<","<<p[ix][iy]<<endl;// pressure field
                }
            }
            for(int i=0;i<N;i++){
                fout13<<u[65][i]/u_infinity<<" "<<(dx*i)/length<<endl;// mid line u- velocity
            }

            for(int i=0;i<N;i++){

               fout14<<dx*i/length<<" "<<v[i][65]/u_infinity<<endl;// mid line v- velocity

            }

            for(int ix=0;ix<N;ix++){
                for(int iy=0;iy<N;iy++){
                    fout17<<dx*ix/length<<" "<<dx*iy/length<<" "<<stream_function[ix][iy]<<endl;// stream function values
                }
            }


            for(int ix=0;ix<N-1;ix++){
                for(int iy=0;iy<N-1;iy++){
                    fout18<<dx*ix/length<<" "<<dx*iy/length<<" "<<sqrt(v[ix][iy]*v[ix][iy]+u[ix][iy]*u[ix][iy])/u_infinity<<endl;// total velocity field
                }
            }

            for(int ix=0;ix<N;ix++){
                for(int iy=0;iy<N;iy++){
                    fout22<<dx*ix/length<<" "<<dx*iy/length<<" "<<vorticity[ix][iy]<<endl;// vorticity field
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

double boundary(int ix, int iy){

    double ans=0.0;
    /*if(iy==0){
        ans= -rho*dx*1*v[ix][iy]+ans;
    }

    if(ix==0){
        ans= -rho*dx*1*un[ix][iy]+ans;
    }*/
    ///////// the mass flux is zero at ix=0 and iy=0 so the value of boundary contribution is anyway zero
    return ans;
}

// momentum equation coefficients

double a_n(int ix, int iy,int z){

    double ans=0.0;

    if(z==1)
     ans= dx*1.0*(D -rho*(v_n_1[ix-1][iy+1]+v_n_1[ix][iy+1])/4.0);
     else if(z==2)
     ans= dx*1.0*(D -rho*(v_n_1[ix][iy]+v_n_1[ix][iy+1])/4.0);

    return ans;
}

double a_s(int ix, int iy,int z){

    double ans=0.0;

        if(z==1)
        ans= dx*1.0*(D + rho*(v_n_1[ix][iy]+v_n_1[ix-1][iy])/4.0);
        else if(z==2)
        ans= dx*1.0*(D + rho*(v_n_1[ix][iy]+v_n_1[ix][iy-1])/4.0);

    return ans;
}


double a_e(int ix, int iy,int z){

    double ans=0.0;

    if(z==1)
      ans= dx*1.0*(D - rho*(u_n_1[ix][iy]+u_n_1[ix+1][iy])/4.0);
      else if(z==2)
        ans= dx*1.0*(D - rho*(u_n_1[ix+1][iy]+u_n_1[ix+1][iy-1])/4.0);

    return ans;
}

double a_w(int ix, int iy,int z){

    double ans=0.0;

    if(z==1)
    ans= dx*1.0*(D + rho*(u_n_1[ix][iy]+u_n_1[ix-1][iy])/4.0);
    else if(z==2)
    ans= dx*1.0*(D + rho*(u_n_1[ix][iy]+u_n_1[ix][iy-1])/4.0);


    return ans;
}


double a(int ix, int iy,int z){

    double ans=0.0;

    if(z==1)
     ans= 8.0*dx*1.0*D -(a_w(ix,iy,1)+a_e(ix,iy,1)+a_n(ix,iy,1)+a_s(ix,iy,1));
     else if(z==2)
     ans= 8.0*dx*1.0*D -(a_w(ix,iy,2)+a_e(ix,iy,2)+a_n(ix,iy,2)+a_s(ix,iy,2));

    return ans;
}

// pressure correction equation coefficients

double p_e(int ix, int iy){

    double answer=0.0;

    answer= rho*dx*dx*1.0*1.0/a(ix+1,iy,1);

    return answer;

}

double p_n(int ix, int iy){

    double answer=0.0;


    answer= rho*dx*dx*1.0*1.0/a(ix,iy+1,2);


    return answer;

}

double p_w(int ix, int iy){

    double answer=0.0;

    answer= rho*dx*dx*1.0*1.0/a(ix,iy,1);

    return answer;

}

double p_s(int ix, int iy){

    double answer=0.0;

    answer= rho*dx*dx*1.0*1.0/a(ix,iy,2);


    return answer;

}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



void gs_u(void){

for(int ix=1;ix<N-2;ix++){

    for(int iy=1;iy<N-1;iy++){

        u[ix][iy]= (dx*1.0*((p[ix-1][iy]-p[ix][iy])) + u[ix+1][iy]*(a_e(ix,iy,1)) + u[ix][iy+1]*(a_n(ix,iy,1)) + u[ix-1][iy]*(a_w(ix,iy,1)) + u[ix][iy-1]*(a_s(ix,iy,1)))/(a(ix,iy,1));
        if(u[ix][iy]!=u[ix][iy])
            cout<<ix<<"..."<<iy<<"....././"<<u_n_1[ix][iy]<<endl;

    }
}

}

void gs_p(void){

for(int ix=1;ix<N-1;ix++){

    for(int iy=1;iy<N-1;iy++){

        p_error[ix][iy]= (rho*dx*1.0*(u[ix][iy] -u[ix+1][iy] +v[ix][iy] -v[ix][iy+1]) + (p_error[ix+1][iy]*p_e(ix,iy) + p_error[ix-1][iy]*p_w(ix,iy) + p_error[ix][iy-1]*p_s(ix,iy) + p_error[ix][iy+1]*p_n(ix,iy)) + boundary(ix,iy))/(p_n(ix,iy)+p_w(ix,iy)+p_e(ix,iy)+p_s(ix,iy));

//            if(ix==16&&iy==31){
//                    cout<<"Coefficient: "<<dx/(p_n(ix,iy)+p_w(ix,iy)+p_e(ix,iy)+p_s(ix,iy))<<endl;
//                    cout<<"Boundary term: "<<boundary(ix,iy)<<endl;
//                    cout<<"North(lid): "<<p_error[ix][iy+1]<<","<<p_n(ix,iy)<<endl;
//                    cout<<"South: "<<p_error[ix][iy-1]<<","<<p_s(ix,iy)<<endl;
//                    cout<<"East: "<<p_error[ix+1][iy]<<","<<p_e(ix,iy)<<endl;
//                    cout<<"West: "<<p_error[ix-1][iy]<<","<<p_w(ix,iy)<<endl;
//                   }
    }
}

}


void gs_v(void){

for(int ix=1;ix<N-1;ix++){

    for(int iy=1;iy<N-2;iy++){

        v[ix][iy]= (dx*1.0*((p[ix][iy-1]-p[ix][iy])) + v[ix+1][iy]*(a_e(ix,iy,2)) + v[ix][iy+1]*(a_n(ix,iy,2)) + v[ix-1][iy]*(a_w(ix,iy,2)) + v[ix][iy-1]*(a_s(ix,iy,2)))/(a(ix,iy,2));

    }
}

}




