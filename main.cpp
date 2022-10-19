#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>
using namespace std;

int main()
{
  int nx = 32;
  int ny = 32;
  int nu = 1;
  int rho = 1;
  double err;
  double u[nx + 2][ny + 1], v[nx + 1][ny + 2], p[nx + 1][ny + 1];
  double us[nx][ny];
  double vs[nx][ny];
  double u_top = 10;
  double dx = 0.1;
  double dy = 0.1;
  double dxi = 1 / dx;
  double dyi = 1 / dy;
  double dt = 0.001;
  int i, j, k;
  double C1 = dy * dy / (2.0 * (dx * dx + dy * dy));
  double C2 = dx * dx / (2.0 * (dx * dx + dy * dy));
  double C3 = dx * dx * dy * dy / (2.0 * (dx * dx + dy * dy)) / dt;
  double tmp;
  int km = 250;
  int m, n = 0;


  for (i = 0; i < nx + 1; i++)
  {
    for (j = 0; j < ny + 1; j++)
    {
      p[i][j] = 0.;
    }
  }

  for (i = 0; i < nx + 2; i++)
  {
    for (j = 0; j < ny + 1; j++)
    {
      u[i][j] = 0.;
    }
  }

  for (i = 0; i < nx + 1; i++)
  {
    for (j = 0; j < ny + 2; j++)
    {
      v[i][j] = 0.0;
    }
  }

  //time step//
  for (k = 1; k <= km; k++)
  {

   //BC for right and left//
    for (j = 0; j < ny + 1; j++)
    {
      u[1][j] = 0.;
      v[0][j] = -v[1][j];

      u[nx][j] = 0.;
      v[nx][j] = -v[nx - 1][j];
    }
    v[0][ny + 1] = -v[1][ny + 1];
    v[nx][ny + 1] = -v[nx - 1][ny + 1];

    // BC for bottom and top//
    for (i = 0; i < nx + 1; i++)
    {
      v[i][1] = 0.;
      u[i][0] = -u[i][1];

      v[i][ny] = 0.;
      u[i][ny] = 2.0 * u_top - u[i][ny - 1]; // for wall
    }
    u[nx + 1][0] = -u[nx][1];
    u[nx + 1][ny] = 2.0 * u_top - u[nx+1][ny - 1];

    // Apply us//
    for (i = 2; i < nx; i++)
    {
      for (j = 1; j < ny; j++)
      {
        double v_here = 0.25 * (v[i - 1][j] + v[i - 1][j + 1] + v[i][j] + v[i][j + 1]);
        us[i][j] = u[i][j] + dt * (nu * (u[i - 1][j] - 2 * u[i][j] + u[i][j + 1]) * (dxi * dxi) + nu * (u[i][j - 1] - 2 * u[i][j] + u[i][j + 1]) * (dyi * dyi) - u[i][j] * (u[i + 1][j] - u[i - 1][j]) * 0.5 * dxi - v_here * (u[i][j + 1] - u[i][j - 1]) * 0.5 * dyi);
      }
    }

    // Apply vs//
    for (i = 1; i < nx; i++)
    {
      for (j = 2; j < ny; j++)
      {
        double u_here = 0.25 * (v[i][j - 1] + v[i][j] + v[i + 1][j - 1] + v[i + 1][j]);
        vs[i][j] = v[i][j] + dt * (nu * (v[i - 1][j] - 2 * v[i][j] + v[i][j + 1]) * (dxi * dxi) + nu * (v[i][j - 1] - 2 * v[i][j] + v[i][j + 1]) * (dyi * dyi) - v[i][j] * (v[i + 1][j] - v[i - 1][j]) * 0.5 * dxi - u_here * (v[i][j + 1] - v[i][j - 1]) * 0.5 * dyi);
      }
    }

    while (1)
    {
      err = 0;
      for (j = 0; j < ny + 1; j++)
      {
        p[0][j] = p[1][j];       // left//
        p[nx][j] = p[nx - 1][j]; // right//
      }

      for (i = 0; j < nx + 1; i++)
      {
        p[i][0] = p[i][1];       // bottom//
        p[i][ny] = p[i][ny - 1]; // top//
      }

      for (i = 1; i < nx; i++)
      {
        for (j = 1; j < ny; j++)
        {
          tmp = p[i][j];
          p[i][j] = C1 * (p[i + 1][j] + p[i - 1][j]) + C2 * (p[i][j + 1] + p[i][j - 1]) - C3 * ((u[i + 1][j] - u[i][j]) / dx + (v[i][j + 1] - v[i][j]) / dy);
          err += (p[i][j] - tmp) * (p[i][j] - tmp);
        }
      }

      if (err <= 0.001)
        break;
    }

    for (i = 2; i < nx; i++)
    {
      for (j = 1; j < ny; j++)
      {
        u[i][j] = us[i][j] - dt / rho * (p[i][j] - p[i - 1][j]) / dx;
      }
    }

    for (i = 1; i < nx; i++)
    {
      for (j = 2; j < ny; j++)
      {
        v[i][j] = vs[i][j] - dt / rho * (p[i][j] - p[i][j - 1]) / dy;
      }
    }
  
    if(k == 1 || k % 10 == 0){

      ofstream fk, ff;
      fk.open("Data/vel"+to_string(m)+".csv");
      ff.open("Data/pre"+to_string(n)+".csv");
      m++;
      n++;

      int z = 0;
      int w = 0;
      fk << "x" << " " << "y" << " " << "z" << " " << "u" << " " << "v" << " " << "w" <<endl;
      ff << "x" << " " << "y" << " " << "z" << " " << "p" << endl;
      for (i = 1; i < nx; i++)
      {
        for (j = 1; j < ny; j++)
      {
          fk << double((i - 0.5) * dx) << " " << double((j - 0.5) * dy) << " " << z << " " << (u[i][j] + u[i + 1][j]) / 2.0 << " " << (v[i][j] + v[i][j + 1]) / 2.0 << " " << w << endl;
          ff << double((i - 0.5) * dx) << " " << double((j - 0.5) * dy) << " " << z << " " << p[i][j] << endl;
        }
      }
    }

  }

  return 0;
}