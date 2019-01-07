
#include "pch.h"
#include "matrix.h"
using namespace std;

int main()
{
	/*
	Matrix A = Matrix(3, 3);
	cin >> A;
	Matrix b = Matrix(3, 1);
	cin >> b;
	Matrix x = Matrix::Solve(A, b);
	x.Show();
	*/
	Matrix C= Matrix(3, 3,2);
	Matrix D = Matrix(3, 3, 4);
	Matrix S = Matrix::Hadamard(C, D);
	S.Show();
	return 0;
}

