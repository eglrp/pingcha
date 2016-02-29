/* Copyright by zyt 2016.2 | All rights reserved */
#include "stdafx.h"
#include<iostream>
#include<math.h>
#include"Matrix.h"
using namespace std;

Matrix::Matrix(Matrix & m)
{
	M = new double*[m.row + 1];
	row = m.row;
	col = m.col;
	int p, q;
	for (p = 1; p <= m.row; p++)
	{
		M[p] = new double[m.col + 1];
	}
	for (p = 1; p <= row; p++)
	{
		for (q = 1; q <= col; q++)
		{
			M[p][q] = m.M[p][q];
		}
	}
}
Matrix::Matrix(int i, int j)
{
	if (i <= 0 || j <= 0)
	{
		cerr << "维数要大于〇";
		return;
	}
	row = i;
	col = j;
	M = new double*[row + 1];
	int p, q;
	for (p = 1; p <= row; p++)
	{
		M[p] = new double[col + 1];
	}
	for (p = 1; p <= row; p++)
	{
		for (q = 1; q <= col; q++)
		{
			M[p][q] = 0;
		}
	}
}

Matrix Matrix::operator+(const Matrix & m2)
{
	if (this->row != m2.row || this->col != m2.col)
	{
		cerr << "维数不匹配，无法相加" << endl;
		exit(1);
	}
	Matrix res = Matrix(row, col);
	for (int p = 1; p <= row; p++)
	{
		for (int q = 1; q <= col; q++)
		{
			res.M[p][q] = this->M[p][q] + m2.M[p][q];
		}
	}
	return res;
}

Matrix Matrix::operator-(const Matrix & m2)
{
	if (this->row != m2.row || this->col != m2.col)
	{
		cerr << "维数不匹配，无法相减" << endl;
		exit(1);
	}
	Matrix res = Matrix(row, col);
	for (int p = 1; p <= row; p++)
	{
		for (int q = 1; q <= col; q++)
		{
			res.M[p][q] = this->M[p][q] - m2.M[p][q];
		}
	}
	return res;
}

Matrix Matrix::operator*(const Matrix & m2)
{
	if (this->col != m2.row)
	{
		cerr << "维数不匹配，无法相乘" << endl;
		exit(1);
	}
	Matrix res = Matrix(row, m2.col);
	for (int p = 1; p <= row; p++)
	{
		for (int q = 1; q <= m2.col; q++)
		{
			double total = 0;
			for (int r = 1; r <= col; r++)
			{
				total += this->M[p][r] * m2.M[r][q];
			}
			res.M[p][q] = total;
		}
	}
	return res;
}

Matrix Matrix::Inverse()
{
	if (this->row != this->col)
	{
		cerr << "行列不相同，无法求逆" << endl;
		exit(1);
	}
	Matrix res = Matrix(row, col);
	Matrix temp = Matrix(*this);
	for (int k = 1; k <= row; k++)
	{
		res.M[k][k] = 1;
	}

	//寻找主元

	for (int p = 1; p <= row; p++)
	{
		double max_value = -1;
		int max_row;
		for (int i = p; i <= row; i++)
		{
			if (max_value < abs((temp.M[i][p])))
			{
				max_value = abs(temp.M[i][p]);
				max_row = i;
			}
		}
		if (max_value == 0)
		{
			cerr << "该矩阵无法求逆" << endl;
			exit(1);
		}

		//交换主元所在行与当前起始行
		if (max_row != p)
		{
			for (int q = 1; q <= col; q++)
			{
				double t;
				t = temp.M[p][q]; temp.M[p][q] = temp.M[max_row][q]; temp.M[max_row][q] = t;
				t = res.M[p][q]; res.M[p][q] = res.M[max_row][q]; res.M[max_row][q] = t;
			}
		}

		// 首元化为1
		for (int q = p + 1; q <= col; q++)
		{
			temp.M[p][q] /= temp.M[p][p];
		}
		for (int q = 1; q <= col; q++)
		{
			res.M[p][q] /= temp.M[p][p];
		}
		temp.M[p][p] = 1;


		// 当前列化为0

		for (int q = p + 1; q <= col; q++)
		{
			for (int i = 1; i < p; i++)
			{
				temp.M[i][q] -= temp.M[i][p] * temp.M[p][q];
			}
			for (int i = p + 1; i <= row; i++)
			{
				temp.M[i][q] -= temp.M[i][p] * temp.M[p][q];
			}
		}

		for (int q = 1; q <= col; q++)
		{
			for (int i = 1; i < p; i++)
			{
				res.M[i][q] -= temp.M[i][p] * res.M[p][q];
			}
			for (int i = p + 1; i <= row; i++)
			{
				res.M[i][q] -= temp.M[i][p] * res.M[p][q];
			}
		}

		for (int i = 1; i <= row; i++)
		{
			temp.M[i][p] = 0;
		}
		temp.M[p][p] = 1;
	}

	return res;
}

Matrix Matrix::Transpose()
{
	Matrix res = Matrix(col, row);
	for (int p = 1; p <= row; p++)
	{
		for (int q = 1; q <= col; q++)
		{
			res.M[q][p] = this->M[p][q];
		}
	}
	return res;
}

void Matrix::Display()
{
	for (int i = 1; i <= row; i++)
	{
		for (int j = 1; j <= col; j++)
		{
			cout << this->M[i][j] << "  ";
		}
		cout << endl;
	}
}

void Matrix::DisplayCross()
{
	if (row == col)
	{
		for (int i = 1; i <= row; i++)
		{
			cout << this->M[i][i] << endl;
		}
	}
	else
	{
		cout << "不是方阵，不能输出对角线";
	}
}

Matrix Matrix::Deleterow(int r)
{
	if (row >= 1 && col >= 1)
	{
		Matrix res = Matrix(row-1, col);
		for (int i = 1; i < r; i++)
		{
			for (int j = 1; j <= col; j++)
			{
				res.M[i][j] = this->M[i][j];
			}
		}
		for (int i = r; i < row; i++)
		{
			for (int j = 1; j <= col; j++)
			{
				res.M[i][j] = this->M[i+1][j];
			}
		}
		return res;
	}
}

Matrix Matrix::Deletecol(int c)
{
	if (row >= 1 && col >= 1)
	{
		Matrix res = Matrix(row, col-1);
		for (int i = 1; i < c; i++)
		{
			for (int j = 1; j <= row; j++)
			{
				res.M[j][i] = this->M[j][i];
			}
		}
		for (int i = c; i < col; i++)
		{
			for (int j = 1; j <= row; j++)
			{
				res.M[j][i] = this->M[j][i+1];
			}
		}
		return res;
	}
}
/* Copyright by zyt 2016.2 | All rights reserved */
