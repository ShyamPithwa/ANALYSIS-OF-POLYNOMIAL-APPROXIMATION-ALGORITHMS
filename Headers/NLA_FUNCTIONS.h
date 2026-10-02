#pragma once
#ifndef GRAMSCHMIDT_H
#define GRAMSCHMIDT_H

#include "matrix.h"

Matrix MG_Schmidt(Matrix &A);

Matrix vandermonde(const vector<double>& x, int degree);
Matrix identity(int n);

Matrix vectorToMatrix(const vector<double> &v);
Matrix Sub_Matrix(const Matrix &A, int i1, int i2, int j3, int j4);
void InsertSubMatrix(Matrix& A,const Matrix& sub,unsigned startRow,unsigned startCol);

double tnorm(const Matrix& M);
double tnorm(const vector<double>& a);
double tnorm(Matrix &A, int k);
double frobenius_norm(const Matrix &A);

vector<double> access_colunm_vector(const Matrix &A, int k);
Matrix access_colunm(const Matrix &A,int k);
void set_column(Matrix &A, int k, const vector<double>& col);
void set_column(Matrix &A, int k, const Matrix &col);



Matrix MV_multiplication(const Matrix &A, const vector<double> &v);
vector<double> VS_Multiplication(const vector<double>& v, double k);
double Dot_product(const vector<double>& Q,const vector<double>& P);
Matrix VM_multiplication(const vector<double>& v,const Matrix& A);
Matrix OuterProduct(const vector<double>& v, const Matrix& row);
Matrix OuterProduct(const vector<double> &a,const vector<double> &b);
double inner_product(const Matrix &A,const Matrix &B);

vector<double> Vec_Sub(const vector<double>& A,const vector<double>& B);
vector<double> Vec_Add(const vector<double>& A,const vector<double>& B);

#endif