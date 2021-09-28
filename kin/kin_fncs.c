#include <stdio.h>
#include <math.h>

#ifndef M_PI
#define M_PI 3.1415927
#endif

//link lengths
const double l0 = 0.25;
const double l1 = 0.2;
const double l2 = 0.2;
const double l3 = 0.15;
//offsets
const double d1 = -0.04;
const double d2 = 0.04;
const double d3 = -0.04;
const double d4 = -0.04;


void multiply(double first[4][4], double second[4][4], double result[4][4], double temp[4][4]){
    for (int i = 0; i < 4; ++i) {
      for (int j = 0; j < 4; ++j) {
         result[i][j] = 0;
      }
   }

    for (int i = 0; i < 4; ++i) {
      for (int j = 0; j < 4; ++j) {
         for (int k = 0; k < 4; ++k) {
            result[i][j] += first[i][k] * second[k][j];
         }
      }
   }
    copy(result, temp);
}

void copy(double result[4][4], double temp[4][4]){
    for (int i = 0; i < 4; ++i) {
      for (int j = 0; j < 4; ++j) {
          temp[i][j] = result[i][j];
      }
   }
}

//Rz(Θ0) Dz(l0) Ry(Θ1) Dy(d1) Dx(l1) Ry1(Θ2) Dy1(d2)
//Dx1(l2) Ry2(Θ3) Dy2(d3) Dz1(d4) Dx2(l3)

// | x j j i |  j, x ,y ,z represents area where rotations happen
// | j y j i |  i represents area where position of new origin
// | j j z i |
// | 0 0 0 1 |

fwd_kin(theta, x)
double theta[6];
double x[3];
{
    double Rz[4][4] =  {{cos(theta[0]), -sin(theta[0]),0,0},
                        {sin(theta[0]),cos(theta[0]),0,0},
                        {0,0,1,0},
                        {0,0,0,1}};

    double Dz[4][4] =  {{1,0,0,0},
                        {0,1,0,0},
                        {0,0,1,l0},
                        {0,0,0,1}};

    double Ry[4][4] =  {{cos(theta[1]),0,sin(theta[1]),0},
                        {0,1,0,0},
                        {-sin(theta[1]),0,cos(theta[1]),0},
                        {0,0,0,1}};

    double Dy[4][4] =  {{1,0,0,0},
                        {0,1,0,d1},
                        {0,0,1,0},
                        {0,0,0,1}};

    double Dx[4][4] =  {{1,0,0,l1},
                        {0,1,0,0},
                        {0,0,1,0},
                        {0,0,0,1}};

    double Ry1[4][4] = {{cos(theta[2]),0,sin(theta[2]),0},
                        {0,1,0,0},
                        {-sin(theta[2]),0,cos(theta[2]),0},
                        {0,0,0,1}};

    double Dy1[4][4] = {{1,0,0,0},
                        {0,1,0,d2},
                        {0,0,1,0},
                        {0,0,0,1}};

    double Dx1[4][4] = {{1,0,0,l2},
                        {0,1,0,0},
                        {0,0,1,0},
                        {0,0,0,1}};

    double Ry2[4][4] = {{cos(theta[3]),0,sin(theta[3]),0},
                        {0,1,0,0},
                        {-sin(theta[3]),0,cos(theta[3]),0},
                        {0,0,0,1}};

    double Dy2[4][4] = {{1,0,0,0},
                        {0,1,0,d3},
                        {0,0,1,0},
                        {0,0,0,1}};

    double Dz1[4][4] = {{1,0,0,0},
                        {0,1,0,0},
                        {0,0,1,d4},
                        {0,0,0,1}};
    
    double Dx2[4][4] = {{1,0,0,l3},
                        {0,1,0,0},
                        {0,0,1,0},
                        {0,0,0,1}};

    double result[4][4] = {{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}};
    double temp[4][4] = {{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}};
    
    multiply(Rz, Dz, result, temp);
    multiply(temp, Ry, result, temp);
    multiply(temp, Dy, result, temp);
    multiply(temp, Dx, result, temp);
    multiply(temp, Ry1, result, temp);
    multiply(temp, Dy1, result, temp);
    multiply(temp, Dx1, result, temp);
    multiply(temp, Ry2, result, temp);
    multiply(temp, Dy2, result, temp);
    multiply(temp, Dz1, result, temp);
    multiply(temp, Dx2, result, temp);

    x[0] = result[0][3];
    x[1] = result[1][3];
    x[2] = result[2][3];
}

// void printMatrix(double matrix[4][4]){
//     for(int i = 0; i < 4; i++){
//         for(int j = 0; j < 4; j++){
//             printf("%f ", matrix[i][j]);
//         }
//         printf("\n");
//     }
//     printf("\n");
// }


inv_kin(x, theta)
double x[3];
double theta[6];
{

}