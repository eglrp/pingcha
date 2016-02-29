/* Copyright by zyt 2016.2 | All rights reserved */
#include "stdafx.h"
#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include"Matrix.h"
using namespace std;

struct edge
{
	int id;
	double pre_height;
	double after_height;
	string frompoint;
	string topoint;
	double length;
};
struct point
{
	bool iskonwn;
	double height;
	string id;
};

/* 从文件中取得边，点数据，并保存在两个vector中 */
void getvalue(vector<edge> & edges, vector<point> & known_points, int& all_edges,
	int& known_point_num, int& unknown_point_num);

/* 计算B矩阵和l矩阵，还有权阵 */
void calcuB_P_l(Matrix & B, Matrix & l, Matrix & P, const vector<edge>& edges,
	const vector<point>& known_points);

/* 计算参数的协因数阵，参数平差值 */
void calcuQXX_X(Matrix & QXX, Matrix & X, Matrix & B, const Matrix & l, const Matrix & P);

/* 计算改正数，观测值的平差值 */
void calcuV_L(Matrix & V, Matrix & L, Matrix & B, const Matrix & X, const Matrix & l,
	vector<edge> & edges);

/* 计算单位权中误差 */
double calcuDelta0(Matrix & V, Matrix & P, int r);

/* 计算参数中误差 */
void calcuDelta_X(Matrix & D_X, Matrix & QXX, double delta0);

/* Copyright by zyt 2016.2 | All rights reserved */
int main()
{
	vector<edge> edges;
	vector<point> known_points;
	int edges_num = 0;
	int known_point_num = 0;
	int unknowm_point_num = 0;
	getvalue(edges, known_points, edges_num, known_point_num, unknowm_point_num);
	Matrix B = Matrix(edges_num, unknowm_point_num);
	Matrix l = Matrix(edges_num, 1);
	Matrix P = Matrix(edges_num, edges_num);
	Matrix Q;
	Matrix res;  //保存中间结果
	Matrix QXX;  //参数的协因数阵
	Matrix X;    //参数平差值矩阵
	Matrix V;    //改正数矩阵
	Matrix L;    //观测值平差值矩阵
	Matrix D_X;  //参数中误差
	Matrix D_L;  //观测值中误差
	double delta0 = 0;  //单位权中误差

	calcuB_P_l(B, l, P, edges, known_points);
	Q = P.Inverse();
	calcuQXX_X(QXX, X, B, l, P);
	calcuV_L(V, L, B, X, l, edges);
	delta0 = calcuDelta0(V, P, edges_num - unknowm_point_num);
	calcuDelta_X(D_X, QXX, delta0);

	//cout << endl;
	QXX.DisplayCross();
	//cout << endl;
	//V.Display();
	//cout << endl;
	//L.Display();
	X.Display();
	//cout<<delta0;
	//D_X.Display();
	//QXX.Display();
	D_X.Display();
	
    return 0;
}
void getvalue(vector<edge> & edges, vector<point> & known_points, int& all_edges,
int& known_point_num, int& unknown_point_num)
{
	string gc_filename = "D://学习//平差//平差课程设计//gc_data.txt";
	ifstream in;
	char buf[1024];

	in.open(gc_filename);

	in >> all_edges >> known_point_num >> unknown_point_num;
	in.getline(buf, 1024);
	in.getline(buf, 1024);
	int i = 1;
	for (i = 1; i < all_edges + 1; i++)
	{
		int id;
		double height;
		double length;
		string frompoint;
		string topoint;
		try
		{
			in >> id >> height >> length >> frompoint >> topoint;
			edge n;
			n.id = id; n.pre_height = height; n.length = length; n.frompoint = frompoint; n.topoint = topoint;
			edges.push_back(n);
		}
		catch (exception e)
		{
			cout << "input error";
			break;
		}
	}
	in.clear();
	in.getline(buf, 1024);
	in.getline(buf, 1024);
	i = 1;
	while (in.good() && !in.eof())
	{
		string id;
		double height;
		in >> id >> height;
		point p;
		p.iskonwn = true;
		p.height = height;
		p.id = id;
		known_points.push_back(p);
	}
}

void calcuB_P_l(Matrix & B, Matrix & l, Matrix & P, const vector<edge>& edges,
	const vector<point>& known_points)
{
	for (int i = 0; i < edges.size(); i++)
	{
		char a = edges[i].frompoint[0];
		char b = edges[i].topoint[0];
		if (a == 'P')
		{
			a = edges[i].frompoint[1];
			int fp = int(a) - '0';
			B.M[edges[i].id][fp] = -1;
		}
		else
		{
			a = edges[i].frompoint[1];
			int fp = int(a) - '0';
			l.M[edges[i].id][1] += known_points[fp - 1].height;
		}
		if (b == 'P')
		{
			b = edges[i].topoint[1];
			int tp = int(b) - '0';
			B.M[edges[i].id][tp] = 1;
		}
		else
		{
			b = edges[i].topoint[1];
			int tp = int(b) - '0';
			l.M[edges[i].id][1] -= known_points[tp - 1].height;
		}

		l.M[edges[i].id][1] += edges[i].pre_height;
		P.M[edges[i].id][edges[i].id] = 1 / edges[i].length;
	}
}
void calcuQXX_X(Matrix & QXX, Matrix & X, Matrix & B, const Matrix & l, const Matrix & P)
{
	Matrix res;
	res = (B.Transpose() * P * B).Inverse();
	QXX = res;
	X = res * B.Transpose() * P * l;
}

void calcuV_L(Matrix & V, Matrix & L, Matrix & B, const Matrix & X, const Matrix & l,
	vector<edge> & edges)
{
	V = B * X - l;
	Matrix L0(edges.size(), 1);
	for (int i = 0; i < edges.size(); i++)
	{
		L0.M[edges[i].id][1] = edges[i].pre_height;
	}
	L = V + L0;
	for (int i = 0; i < edges.size(); i++)
	{
		edges[i].after_height = L.M[edges[i].id][1];
	}

}

double calcuDelta0(Matrix & V, Matrix & P, int r)
{
	double res;
	res = (V.Transpose() * P * V).M[1][1];
	res = sqrt(res / r);

	return res;
}

void calcuDelta_X(Matrix & D_X, Matrix & QXX, double delta0)
{
	D_X = Matrix(QXX.row, 1);
	double delta2 = delta0 * delta0;
	for (int i = 1; i <= QXX.row; i++)
	{
		D_X.M[i][1] = sqrt(QXX.M[i][i] * delta2);
	}
}
/* Copyright by zyt 2016.2 | All rights reserved */
