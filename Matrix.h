/* Copyright by zyt 2016.2 | All rights reserved */
class Matrix
{
private:

public:
	int row; int col;
	double** M;
	Matrix() {};
	Matrix(Matrix&);
	Matrix(int i, int j);
	Matrix Matrix::operator+(const Matrix & m2);
	Matrix Matrix::operator-(const Matrix & m2);
	Matrix Matrix::operator*(const Matrix & m2);
	Matrix Inverse();
	Matrix Transpose();
	Matrix Deleterow(int r);
	Matrix Deletecol(int c);
	void DisplayCross();
	void Display();
};
