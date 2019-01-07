/*
author: cclplus
date:2018/12/09
if you think it is necessary to reward me,
my alipay account number is 707101557@qq.com
*/

#ifndef __MATRIX_CLL_H__
#define __MATRIX_CCL_H__
#include "pch.h"

class Matrix {
private:
	int rows_num, cols_num;
	double **p;
	void initialize();//初始化矩阵

public:
	Matrix(int, int);
	Matrix(int, int, double);//预配分空间
	virtual ~Matrix();//析构函数应当是虚函数，除非此类不用做基类
	Matrix& operator=(const Matrix&);//矩阵的复制
	Matrix& operator=(double *);//将数组的值传给矩阵
	Matrix& operator+=(const Matrix&);//矩阵的+=操作
	Matrix& operator-=(const Matrix&);//-=
	Matrix& operator*=(const Matrix&);//*=
	Matrix operator*(const Matrix & m)const;
	static Matrix Hadamard(const Matrix & m, const Matrix & n);
	static Matrix Solve(const Matrix&, const Matrix&);//求解线性方程组Ax=b
	void Show() const;//矩阵显示
	void swapRows(int, int);
	double det();//求矩阵的行列式
	double Point(int i, int j) const;
	static Matrix inv(Matrix);//求矩阵的逆矩阵
	static Matrix eye(int);//制造一个单位矩阵
	int row() const;
	int col() const;
	static Matrix T(const Matrix & m);//矩阵转置的实现,且不改变矩阵
	Matrix gaussianEliminate();//高斯消元法
	friend std::istream& operator>>(std::istream&, Matrix&);//实现矩阵的输入
};


#endif
