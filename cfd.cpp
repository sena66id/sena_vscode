#include "iostream"
#include "fstream"
#include "Eigen/Core"
#include "Eigen/LU"

using namespace Eigen;
using namespace std;

typedef Matrix<float,12,12> Matrix12f;
typedef Matrix<float,14,14> Matrix14f;
typedef Matrix<float,12*12+1,12*12+1> Matrix145f;
typedef Matrix<float,12*12,12*12> Matrix144f;
typedef Matrix<float,12*12,1> Vector144f;

const int nx = 12, //mesh number of x direction
          ny = 12, //mesh number of y direction 
          nt = 100,
          imin = 1,
          jmin = 1,
          imax = imin+nx,
          jmax = jmin+ny;

const float Re = 10,
            rho = 1000, //[kg/m^3]
            Lx = 1, //[m]
            Ly = 1, //[m]
            dt = 0.0001,
            mu = 0.001, // [Pa*s]
            nu = mu/rho,
            dx = Lx/(nx),
            dy = Ly/(ny),
            dxi = 1/dx,
            dyi = 1/dy,
            u_top = mu*Re/Lx;

void location_x(float[]);
void location_y(float[]);
void location_xm(float[], float[]);
void location_ym(float[], float[]);
Matrix144f Laplacian(Matrix145f,Matrix144f);
Matrix14f discretization_u(Matrix14f, Matrix14f, Matrix14f);
Matrix14f discretization_v(Matrix14f, Matrix14f, Matrix14f);
Matrix12f Poisson_equation(Matrix144f, Matrix14f, Matrix14f);
Matrix14f nextstep_u(Matrix14f, Matrix14f, Matrix12f);
Matrix14f nextstep_v(Matrix14f, Matrix14f, Matrix12f);
void print(Matrix14f);

int main()
{
    // create mesh
    float x[imax], y[jmax], xm[imax], ym[jmax];

    location_x(x);
  
    location_y(y);
   
    location_xm(x,xm);
   
    location_ym(y,ym);
    
    //initial condition
    Matrix14f u,v;
    u = Matrix14f::Zero();
    v = Matrix14f::Zero();
    for (int i = 1; i<imax+1; i++){
        u(i,jmax) = u(i,imax-1)-2*(u(i,jmax-1)-u_top);
    }

    Matrix14f un, vn, us, vs;
    Matrix145f L;
    Matrix144f Lap;
    Matrix12f p;

    Lap = Laplacian(L,Lap);

    for (int n=0; n<nt; n++){

        
        vn = v;un = u;
        
        us = discretization_u(us,un,vn);

        vs = discretization_v(vs,un,vn);

        p = Poisson_equation(Lap,us,vs);

        u = nextstep_u(u,us,p);

        v = nextstep_v(v,vs,p);

        //boundary conditions
        for (int i=0; i<imax; i++){
            u(imin,i)= 0; //left wall velocity u
            u(imax,i) = 0; // right wall velocity u
            u(i,jmin-1) = u(i,jmin)-2*(u(i,jmin)-0); //bottom wall velocity u
            v(i,jmin) = 0; // bottom wall velocity v
            v(i,jmax) = 0; // top velocity v
            v(imin-1,i) = v(imin,i)-2*(v(imin,i)-0); //left wall velocity v
            v(imax,i) = v(imax-1,i)-2*(v(imax-1,i)-0); //right wall velocity v        
            }
            for (int i = 1; i<imax+1; i++){
            u(i,jmax) = u(i,imax-1)-2*(u(i,jmax-1)-u_top);
            }

    }
    print(v);

   return 0;
}

void location_x(float x[imax]){
    for (int i = 0; i<imax; i++){
        x[i] = 0.0 + i*dx;
    }
}

void location_y(float y[jmax]){
    for ( int i = 0; i<jmax; i++){
        y[i] = 0.0 + i*dy;
    }
}

void location_xm(float x[imax], float xm[imax]){
    for (int i = 1; i<imax; i++){
        xm[i] = 0.5*(x[i-1]+x[i]);
    }
}

void location_ym(float y[jmax], float ym[jmax]){
    for (int i = 1; i<jmax; i++){
        ym[i] = 0.5*(y[i-1]+y[i]);
    }
}

Matrix14f discretization_u(Matrix14f us, Matrix14f un, Matrix14f vn){
    float v_here;
        
    for (int j = jmin; j<jmax; j++){
        for (int i = imin+1; i<imax; i++ ){
            v_here = 0.25*(vn(i-1,j)+vn(i-1,j+1)+vn(i,j)+vn(i,j+1));
            us(i,j) = un(i,j)+dt*(nu*(un(i-1,j)-2*un(i,j)+un(i+1,j))*pow(dxi,2.0)+nu*(un(i,j-1)-2*un(i,j)+un(i,j+1))*pow(dyi,2.0)-un(i,j)*(un(i+1,j)-un(i-1,j))*0.5*dxi-v_here*(un(i,j+1)-un(i,j-1))*0.5*dyi);
        }
    }
    return us;
}

Matrix14f discretization_v(Matrix14f vs, Matrix14f un, Matrix14f vn){
    float u_here;
        
    for (int j = jmin; j<jmax; j++){
        for (int i = imin+1; i<imax; i++ ){
            u_here = 0.25*(un(i,j-1)+un(i,j)+un(i+1,j-1)+un(i+1,j));
            vs(i,j) = vn(i,j)+dt*(nu*(vn(i-1,j)-2*vn(i,j)+vn(i+1,j))*pow(dxi,2.0)+nu*(vn(i,j-1)-2*vn(i,j)+vn(i,j+1))*pow(dyi,2.0)-vn(i,j)*(vn(i,j+1)-vn(i,j-1))*0.5*dyi-u_here*(vn(i+1,j)-vn(i-1,j))*0.5*dxi);
        }
    }
    return vs;
}

Matrix144f Laplacian(Matrix145f L, Matrix144f Lap){
    L = Matrix145f::Zero();
    for(int j=1; j<=ny; j++){
        for(int i=1; i<=nx; i++){
            L(i+(j-1)*nx, i+(j-1)*nx) = 2*pow(dxi,2.0)+2*pow(dyi,2.0);
            for (int k=i-1; k<=i+1; k = k+2){
                if(k>0 && k<=nx){
                    L(i+(j-1)*nx, k+(j-1)*nx) = -pow(dxi,2.0);
                }
                else{
                    L(i+(j-1)*nx,i+(j-1)*nx) = L(i+(j-1)*nx,i+(j-1)*nx)-pow(dxi,2.0);
                }
            }
            for(int l=j-1; l<=j+1; l = l+2){
                if(l>0 && l<=ny){
                    L(i+(j-1)*nx,i+(l-1)*nx) = -pow(dyi,2.0);
                }
                else{
                    L(i+(j-1)*nx,i+(j-1)*nx) = L(i+(j-1)*nx,i+(j-1)*nx)-pow(dyi,2.0);
                }
            }
        }
    }
    for( int i=0; i<=nx*nx; i++){
        L(1,i) = 0;
    }
    L(1,1) = 1;

    for (int j=0; j<144; j++){
        for (int i=0; i<144; i++){
            Lap(i,j) = L(i+1,j+1);
        }
    }
    return Lap;
}

Matrix12f Poisson_equation(Matrix144f Lap, Matrix14f us, Matrix14f vs){
    Vector144f R, pv;
    Matrix12f p;
    int n = 0;
    for (int j = jmin; j<jmax; j++){
        for (int i=imin; i<imax; i++){
            n = n+1;
            R(n-1,0) = -rho/dt*((us(i+1,j)-us(i,j))*dxi+(vs(i,j+1)-vs(i,j))*dyi);
        }
    }
    
    FullPivLU<Matrix144f> lu(Lap);
    pv = lu.solve(R);
    
    n = 0;
    for (int j = 0; j<jmax-1; j++){
        for ( int i = 0; i<imax-1; i++){
            n = n+1;
            p(i,j) = pv(n-1,0);
        }
    }

    return p;
}

Matrix14f nextstep_u(Matrix14f u, Matrix14f us, Matrix12f p){
    for (int j = jmin; j<jmax; j++){
        for ( int i = imin+1; i<imax; i++){
            u(i,j) = us(i,j)-dt/rho*(p(i-1,j-1)-p(i-2,j-1))*dxi;
        }
    }
    return u;
}

Matrix14f nextstep_v(Matrix14f v, Matrix14f vs, Matrix12f p){
    for (int j = jmin+1; j<jmax; j++){
        for ( int i = imin; i<imax; i++){
            v(i,j) = vs(i,j)-dt/rho*(p(i-1,j-1)-p(i-1,j-2))*dyi;
        }
    }
    return v;
}

void print(Matrix14f v){
    ofstream fout ("out.dat");

    if(fout.fail()){  
        cout << "出力ファイルをオープンできません" << endl;
    }

    fout << "v = \n" << v << endl;
}