#include <stdio.h>
#include <stdlib.h>
#include "iostream"
#include <math.h>
#include <iomanip>

using namespace std;

#define VEC_SIZE 3
#define MATRIX_SIZE 3

//---------------------------------------------------------------
//出力用
//---------------------------------------------------------------
void printVector(double vec[VEC_SIZE]){
    for(int i = 0; i< MATRIX_SIZE; i++){
        cout << vec[i] << " ";
    }
    cout << endl;
}

void printMatrix(double mat[MATRIX_SIZE][MATRIX_SIZE]){
    for(int i = 0; i < MATRIX_SIZE; i++){
        printVector(mat[i]);
    }
}

double vec_mlt(double v0[MATRIX_SIZE], double v1[MATRIX_SIZE]){
    double result = 0;
    for(int i = 0; i < MATRIX_SIZE; i++){
        result += v0[i] * v1[i];
    }
    return result;
}

void mat_mlt(double a[MATRIX_SIZE][MATRIX_SIZE], double b[MATRIX_SIZE][MATRIX_SIZE]){
    
    for(int i = 0; i < MATRIX_SIZE; i++){
        cout << endl;
        for(int j = 0; j < MATRIX_SIZE; j++){
            //double vec[MATRIX_SIZE] = {b[0][j], b[1][j], b[2][j]};
            double vec[MATRIX_SIZE];
            for(int k = 0; k < MATRIX_SIZE; k++){
                vec[k] = b[k][j];
            }
            cout << vec_mlt(a[i], vec) << " ";
        }
    }
}

//前進代入計算で数値解を求める
void forwardSubstitution(const double mat[VEC_SIZE][VEC_SIZE], const double vec[VEC_SIZE], double resultVec[VEC_SIZE]){
    for(int i = 0; i < VEC_SIZE; i++){
        double result = 0;
        result += vec[i];
        if(i == 0){
            resultVec[i] = result;
            continue;
        }
        for(int k = 0; k < i; k++){
            result -= mat[i][k] * resultVec[k];
        }
        resultVec[i] = result;
    }
    return;
}

//後進代入計算で数値解を求める
void backSubstitution(const double mat[VEC_SIZE][VEC_SIZE], const double vec[VEC_SIZE], double resultVec[VEC_SIZE]){
    for(int i = VEC_SIZE - 1; i >= 0; i--){
        double result = 0;
        result += vec[i];
        if(i == VEC_SIZE - 1){
            resultVec[i] = result/mat[i][i];
            continue;
        }
        for(int k = i + 1; k < VEC_SIZE; k++){
            result -= mat[i][k] * resultVec[k];
        }
        resultVec[i] = result/mat[i][i];
    }
    return;
}

//LU分解
void LUdecomposition(const double a[MATRIX_SIZE][MATRIX_SIZE], double l[MATRIX_SIZE][MATRIX_SIZE], double u[MATRIX_SIZE][MATRIX_SIZE]){
    //U0j = a0j から縦に解いていく
    for(int w = 0; w < MATRIX_SIZE; w++){
        for (int h = 0; h < MATRIX_SIZE; h++) {
            //L
            if(h == w){
                l[h][w] = 1.0;
            }else if(h > w){
                l[h][w] = a[h][w];
                for(int k = 0; k < w; k++ ){
                    l[h][w] -= l[h][k] * u[k][w];
                }
                l[h][w] /= u[w][w];
            }
            //about U
            if(w >= h){
                u[h][w] = a[h][w];
                for(int k = 0; k < h; k++){
                    u[h][w] -= l[h][k] * u[k][w];
                }
            }
        }
    }
    return;
}

//連立 Ax = b
void solveEquation(const double a[MATRIX_SIZE][MATRIX_SIZE], const double b[MATRIX_SIZE], double result[MATRIX_SIZE]){
    double x[VEC_SIZE];
    double x1[VEC_SIZE];
    double lMat[VEC_SIZE][VEC_SIZE];
    double uMat[VEC_SIZE][VEC_SIZE];
    
    LUdecomposition(a,lMat,uMat);
    forwardSubstitution(lMat, b, x1);
    backSubstitution(uMat, x1, x);
    
    for(int i = 0; i < MATRIX_SIZE; i++){
        result[i] = x[i];
    }
    printVector(result);
}

void solveEquation(const double a[MATRIX_SIZE][MATRIX_SIZE], const double b[MATRIX_SIZE]){
    double result[MATRIX_SIZE];
    solveEquation(a, b, result);
}

void getReverseMatrix(const double A[MATRIX_SIZE][MATRIX_SIZE],double reverseA[MATRIX_SIZE][MATRIX_SIZE]){
    double result[MATRIX_SIZE];
    for(int i = 0; i < MATRIX_SIZE; i++){
        //i番目の要素のみ１のベクトルを利用して解く
        //A * A^(-1)[0~MATRIX_SIZE][i] = {0,0,0,....1(i番目), 0,0,0 ,,,}を利用
        double simpleVec[MATRIX_SIZE];
        simpleVec[i] = 1;
        solveEquation(A,simpleVec,result);
        
        //方程式の解をreverseAに適用
        for(int j = 0; j < MATRIX_SIZE; j++){
            reverseA[j][i] = result[j];
            simpleVec[j] = 0;
        }
    }
    printMatrix(reverseA);
}

void getReverseMatrix(const double A[MATRIX_SIZE][MATRIX_SIZE]){
    double reverseA[MATRIX_SIZE][MATRIX_SIZE];
    getReverseMatrix(A, reverseA);
};

int main(int argc, char* argv[]){
    
    //-------------------------------------------------
    //課題１、２
    //-------------------------------------------------
    double fMat[VEC_SIZE][VEC_SIZE] = {
        {1.0, 0.0, 0.0},
        {3.0, 1.0, 0.0},
        {-2.0, 2.0, 1.0}
    };
    double fVec[VEC_SIZE] = {2.0, 3.0, -1.0};
    double bMat[VEC_SIZE][VEC_SIZE] = {
        {2.0, 1.0, -1.0},
        {0.0, 3.0, 2.0},
        {0.0, 0.0, -3.0}
    };
    double bVec[VEC_SIZE] = {2.0, -3.0, 9.0};
    double forwardResult[VEC_SIZE] = {0,0,0};
    double backResult[VEC_SIZE] = {0,0,0};
    cout << "Kadai1" << endl;
    forwardSubstitution(fMat, fVec, forwardResult);
    cout << forwardResult[0] << endl;
    cout << forwardResult[1] << endl;
    cout << forwardResult[2] << endl;
    cout << "Kadai2" << endl;
    backSubstitution(bMat, bVec, backResult);
    cout << backResult[0] << endl;
    cout << backResult[1] << endl;
    cout << backResult[2] << endl;
    cout << endl;
    //-------------------------------------------------
    //課題4
    //-------------------------------------------------
    cout << "Kadai4" << endl;
    double a[MATRIX_SIZE][MATRIX_SIZE] = {
        {2.0, 1.0, -1.0},
        {6.0, 6.0, -1.0},
        {-4.0,4.0, 3.0}
    };
    
    double l[MATRIX_SIZE][MATRIX_SIZE]; //LMatrix
    double u[MATRIX_SIZE][MATRIX_SIZE]; //UMatrix
    
    LUdecomposition(a,l,u);
    
    for(int h = 0; h < MATRIX_SIZE; h++){
        for(int w = 0; w < MATRIX_SIZE; w++){
            cout << l[h][w] << "  ";
        }
        cout << endl;
    }
    cout << endl;
    //-------------------------------------------------
    //課題５ 済み
    //-------------------------------------------------
    cout << "Kadai5" << endl;
    //まずLU分解 L * x' = b の形にして
    //U * x = x' を前進代入により求める
    //その後 xを後進代入により求める
    //定義してメソッド読んでくだけ
    double b[VEC_SIZE] = {2.0 , 3.0 , -1.0};
    solveEquation(a,b);
    cout << endl;
    //-------------------------------------------------
    //課題６　課題５の確認 SIZEを3と7確認済み
    //-------------------------------------------------
    cout << "Kadai6" << endl;
    double H6[MATRIX_SIZE][MATRIX_SIZE];
    double b6[MATRIX_SIZE];
    for(int i = 0; i < MATRIX_SIZE; i++){
        for(int j = 0; j < MATRIX_SIZE; j++){
            H6[i][j] = pow(0.5, abs(i - j));
        }
        b6[i] = 3 - pow(2.0, i - MATRIX_SIZE + 1) - pow(2, -i);
    }
    solveEquation(H6,b6);
    cout <<  endl;
    //-------------------------------------------------
    //課題７　課題５の応用逆行列
    //-------------------------------------------------
    cout << "Kadai7" << endl;
    //AA^-1 = I より
    double exampleA[MATRIX_SIZE][MATRIX_SIZE] = {
        {1, 0, 0},
        {0, 1, 0},
        {0, 0, 1}
    };
    getReverseMatrix(exampleA);
    cout << endl;
    //-------------------------------------------------
    //課題８　課題７の確認
    //-------------------------------------------------
    cout << "Kadai8" << endl;
    double H[MATRIX_SIZE][MATRIX_SIZE];
    double reverseH[MATRIX_SIZE][MATRIX_SIZE];
    for(int i = 0; i < MATRIX_SIZE; i++){
        for(int j = 0; j < MATRIX_SIZE; j++){
            H[i][j] = pow(0.5, abs(i - j));
        }
    }
    getReverseMatrix(H, reverseH);
    mat_mlt(H, reverseH);
    cout << endl;
}

