/* Copyright by zyt 2016.2 | All rights reserved */
#include "stdafx.h"
#include<fstream>
#include<vector>
#include<iostream>
#include<string>
#include<math.h>
#include<map>
#include<sstream>
#include"Matrix.h"
using namespace std;

class angleBase
{
public:
	int degree;
	int minute;
	double second;

	angleBase()
	{
		degree = 0;
		minute = 0;
		second = 0.0;
	}
	angleBase(int d, int m, double s)
	{
		degree = d;
		minute = m;
		second = s;
	}
	angleBase(const angleBase& a)
	{
		degree = a.degree;
		minute = a.minute;
		second = a.second;
	}
	angleBase(double da)
	{
		degree = int(da);
		minute = int((da - degree) * 60);
		second = (da - degree - minute / 60.0) * 3600;
	}
	angleBase angleBase::operator+(angleBase a)
	{
		this->second += a.second;
		while (this->second >= 60)
		{
			this->second -= 60;
			this->minute++;
		}
		this->minute += a.minute;
		while (this->minute >= 60)
		{
			this->minute -= 60;
			this->degree++;
		}
		this->degree += a.degree;
		this->degree %= 360;

		return *this;
	}
	angleBase angleBase::operator-(angleBase a)
	{
		this->second -= a.second;
		while (this->second < 0)
		{
			this->second += 60;
			this->minute--;
		}
		this->minute -= a.minute;
		while (this->minute < 0)
		{
			this->minute += 60;
			this->degree--;
		}
		this->degree -= a.degree;
		while (this->degree < 0)
		{
			this->degree += 360;
		}
		this->degree %= 360;
		return *this;
	}

	double toReg()
	{
		double reg = (degree + minute / 60.0 + second / 3600.0) * 3.1415926 / 180.0;
		return reg;
	}

	double toDeg()
	{
		double deg = (degree + minute / 60.0 + second / 3600.0);
		return deg;
	}

	void DisplayOrigin()
	{
		if(second < 0 || minute <0 || degree <0)
		{
			cout << "-";
		}
		cout << abs(degree) << "°\t" << abs(minute) << "'\t" << abs(second) << "''\t";
	}
};

class angle
{
public:
	int id;
	angleBase pre;
	angleBase after;
	int startP;
	int middleP;
	int endP;

	angle(int sP, int mP, int eP, int d, int m, double s, int id)
	{
		pre = angleBase(d, m, s);
		after = angleBase(0, 0, 0);
		startP = sP;
		middleP = mP;
		endP = eP;
		this->id = id;
	}
	angle(const angle& a)
	{
		pre = a.pre;
		after = a.after;
		startP = a.startP;
		middleP = a.middleP;
		endP = a.endP;
		this->id = a.id;
	}
};

struct point
{
	int id;
	bool isKnown;
	bool isCalcu;
	double pre_x;
	double pre_y;
	double after_x;
	double after_y;
};

struct edge
{
	int id;
	int Pid1;
	int Pid2;
	bool isknown;
	double pre_value;
	double after_value;
};

/* 常量 */
const int ruo = 206265;
const int all_point_num = 18;
const int known_point_num = 3;

/* 全局变量 */
vector<angle> angles;
vector<point> points;
vector<point> points1;
vector<point> points2;
vector<edge> edges;
map<string, int> index_edge;

/******************************* 函数声明 ********************************/

/* 读取文件 */
void readdata(string filename,int& start1, int & start2, int & the3, vector<int>& known_point, int all_point_num , int known_point_num);
//void readdata2(string filename, int& start1, int & start2, int & the3, vector<int>& known_point, int all_point_num , int known_point_num);
/* 计算近似坐标 */
void approZB(int P1, int P2);

/* 按顺序给出两点，计算方位角 */
double azimuth(int P1, int P2);

/* 根据两点编号，在边的索引表中找到该边的编号 */
int getIndexOfEdge(int P1, int P2);

/* 角度转弧度 */
double DegtoReg(double degree)
{
	double reg = degree * 3.1415926 / 180.0;
	return reg;
}

/* 弧度转角度 */
double RegtoDeg(double reg)
{
	double degree = reg * 180 / 3.1415926;
	return degree;
}

/* 计算角度误差方程中，系数ab的值 */
void calcuAB(int P1, int P2, double & a, double & b);

/* 计算边长误差方程中，系数ab的值 */
void calcuAB_forS(int P1, int P2, double & a, double & b, double & S);

/* 给出3点，计算夹角（右角） */
double interAngle(int P1, int P2, int P3);

/**********与矩阵计算有关的函数***********/
/* 计算B、l、P矩阵 */
void calcuB_l_P(Matrix & B, Matrix & l, Matrix & P, Matrix & L0, int & known_edge);

/* 计算参数的协因数阵，参数平差值 */
void calcuQXX_x(Matrix & QXX, Matrix & x, Matrix & B, const Matrix & l, const Matrix & P);

/* 计算观测值平差值的协因数阵 */
void calcuQLL(Matrix & QLL, const Matrix & QXX, Matrix & B);

/* 计算X的平差值*/
void calcuX(Matrix & X, const Matrix & x, vector<int>& known_point);

/* 计算改正数，观测值的平差值 */
void calcuV_L(Matrix & V, Matrix & L, Matrix & B, const Matrix & L0, const Matrix & x, const Matrix & l);

/* 计算单位权中误差 */
double calcuDelta0(Matrix & V, Matrix & P, int r);

/* 计算参数中误差 */
void calcuDelta_X(Matrix & D_X, Matrix & QXX, double delta0);

/* 计算观测值中误差 */
void calcuDelta_L(Matrix & D_L, Matrix & QLL, double delta0);

/* 计算极大值、极小值，极值方向 */
void calcuE_F_faiE_faiF(int i,const double & delta0, const Matrix & QXX, double & faiE, double & faiF, double & E, double & F);

int main()
{
	int P1, P2, the3;
	vector<int> known_point;
	readdata("D:\\学习\\平差\\平差课程设计\\观测值表(1).txt", P1, P2, the3, known_point, all_point_num, known_point_num);
	approZB(P1, P2);
	vector<int> unknowm_point;
	vector<int>::iterator itt;
	for (int i = 0; i < points.size(); i++)
	{
		itt = find(known_point.begin(), known_point.end(), i);
		if (itt == known_point.end())
		{
			unknowm_point.push_back(i);
		}
	}
	//for (int i = 0; i < points.size(); i++)
	//{
	//	point p ;
	//	p.pre_x = points[i].pre_x;
	//	p.pre_y = points[i].pre_y;
	//	points1.push_back(p);
	//}
	//for (int i = 0; i < points.size(); i++)
	//{
	//	if (i != P1 && i != P2 & i != the3)
	//	{
	//		points[i].isCalcu = false;
	//	}
	//}
	//approZB(P2, P1);

	//for (int i = 0; i < points.size(); i++)
	//{
	//	point p;
	//	p.pre_x = points[i].pre_x;
	//	p.pre_y = points[i].pre_y;
	//	points2.push_back(p);
	//}
	//for (int i = 0; i < points.size(); i++)
	//{
	//	points[i].pre_x = (points1[i].pre_x + points2[i].pre_x) / 2;
	//	points[i].pre_y = (points1[i].pre_y + points2[i].pre_y) / 2;
	//}

	cout << "坐标近似值为:" << endl;
	for (int i = 0; i < points.size(); i++)
	{
		cout << i<< ":"<<points[i].pre_x << " " << points[i].pre_y << endl;
	}

	int all_l_value = edges.size() + angles.size();
	int all_x_value = 2 * (all_point_num);
	Matrix B_big(all_l_value, all_x_value);
	Matrix l_big(all_l_value, 1);
	Matrix L0_big(all_l_value, 1);
	Matrix P_big(all_l_value, all_l_value);
	Matrix v(all_l_value, 1);
	Matrix P;
	Matrix l;
	Matrix L0;
	Matrix B;
	int known_edge;
	calcuB_l_P(B_big, l_big, P_big, L0_big, known_edge);
	P = P_big.Deletecol(known_edge).Deleterow(known_edge);
	l = l_big.Deleterow(known_edge);
	L0 = L0_big.Deleterow(known_edge);

	// 去掉B_big中的全0的列
	B = B_big;
	for (int i = known_point_num -1 ; i >=0 ; i--)
	{
		B = B.Deletecol(2 * known_point[i] + 1).Deletecol(2 * known_point[i] + 1);
	}
	B = B.Deleterow(known_edge);
	Matrix QXX;  //参数的协因数阵
	Matrix x;	//参数改正数矩阵
	Matrix X;    //参数平差值矩阵
	Matrix V;    //改正数矩阵
	Matrix L;    //观测值平差值矩阵
	Matrix D_X;  //参数中误差
	Matrix D_L;  //观测值中误差
	Matrix QLL;  //观测值协因数
	double delta0 = 0;  //单位权中误差

	calcuQXX_x(QXX, x, B, l, P);
	calcuV_L(V, L, B, L0, x, l);
	calcuX(X, x, known_point);
	calcuQLL(QLL, QXX, B);
	delta0 = calcuDelta0(V, P, x.row);
	calcuDelta_X(D_X, QXX, delta0);
	calcuDelta_L(D_L, QLL, delta0);

	// 输出单位权中误差
	cout << "单位权中误差为:" << delta0 << endl;

	// 输出坐标平差值
	cout << "坐标的平差值为：" << endl;
	for (int i = 0; i < points.size(); i++)
	{
		cout << i << ":" << points[i].after_x << " , " << points[i].after_y << endl;
	}

	// 输出坐标中误差和改正数
	cout << "坐标的中误差和改正数（mm）为：" << endl;
	for (int i = 0; i < D_X.row; i++)
	{
		cout << unknowm_point[i / 2];
		if (i % 2 == 0)
		{
			cout << "X:" << D_X.M[i+1][1];
		}
		else
		{
			cout << "Y:" << D_X.M[i+1][1];
		}
		cout << endl;
		if (i % 2 == 0)
		{
			cout << "改正数:" << x.M[i + 1][1];
		}
		else
		{
			cout << "改正数:" << x.M[i + 1][1];
		}
		cout << endl;
	}

	// 点位中误差最大的点
	double max_DW = -1;
	int max_i;
	for (int i = 0; i < D_X.row / 2; i++)
	{
		double temp = D_X.M[2 * i + 1][1] * D_X.M[2 * i + 1][1] + D_X.M[2 * i + 2][1] * D_X.M[2 * i + 2][1];
		if (max_DW < temp)
		{
			max_DW = temp;
			max_i = i;
		}
	}
	cout << "最弱点是 " << unknowm_point[max_i] <<" 号点, 点位中误差为:"<< sqrt(max_DW) << endl;

	// 输出角度平差值
	int angles_size = angles.size();
	cout << "角度的平差值和改正数为：" << endl;
	for (int i = 0; i < angles_size; i++)
	{
		angleBase a = angleBase(L.M[i + 1][1]);
		cout << i << ":" << a.degree << "°" << a.minute << "'" << a.second << "''" << endl;
		a = angleBase(V.M[i + 1][1] / 3600);
		cout << "改正数" << ":" << a.degree << "\t°" << a.minute << "\t'" << a.second << "\t''" << endl;
	}

	// 输出边长平差值
	cout << "边长的平差值和改正数(除去了已知边)" << endl;
	for (int i = 0; i < edges.size()-1; i++)  //减一是因为有一条边已知，而在D_L没有它
	{
		cout << i+1 << ":" << L.M[angles_size + i + 1][1]<< "m" << endl;
		cout << "改正数" << ":" << V.M[angles_size + i + 1][1] << "m" << endl;
	}

	// 寻找中误差最大的 边
	int max_edge = 0;
	double max = -1;
	for (int i = 0; i < edges.size()-1; i++)  //减一是因为有一条边已知，而在D_L没有它
	{
		if (max < D_L.M[angles_size + i + 1][1])
		{
			max = D_L.M[angles_size + i + 1][1];
			max_edge =  i;
		}
	}
	// 输出最弱边和它的相对中误差
	cout << "中误差最大的边是【未知边】中的第 " << max_edge + 1 << " 条，" <<"中误差为： "<< max<<", 相对中误差为：" << max / 1000 / L.M[max_edge + angles_size + 1][1]<<endl;

	// 计算E/F,极大方向
	double faiE, faiF;
	double E, F;
	for (int i = 0; i < all_point_num - known_point_num ; i++)
	{
		calcuE_F_faiE_faiF(i + 1, delta0, QXX, faiE, faiF, E, F);
		cout << "第" << unknowm_point[i] << "号点的位差极大值是：" << E << ",极小值是：" << F << "极大值方向是:";
		angleBase(faiE).DisplayOrigin();
		cout<< endl;
	}

	B.Display();
    return 0;
}

/************************函数定义*************************/
void readdata(string filename, int& start1, int & start2, int & the3, vector<int>& known_point, int all_point_num ,int known_point_num)
{
	ifstream myfile;
	myfile.open(filename);
	char buf[1024];
	myfile.getline(buf, 1024);
	myfile.getline(buf, 1024);
	myfile.getline(buf, 1024);
	int PL1; 
	int PL2;
	int PL3;
	double Length;
	int d;
	int m;
	double s;
	string blank;
	char extra[3];
	int i = 0;
	
	while (true)
	{
		myfile >> PL1;
		myfile >> blank;
		myfile >> PL2 >> Length;
		if (!myfile.good())
		{
			break;
		}
		edge e;
		e.Pid1 = PL1;
		e.Pid2 = PL2;
		e.id = i;
		e.isknown = false;
		e.pre_value = Length;
		edges.push_back(e);

		//建立边的索引表
		stringstream ss1, ss2;
		string s1,s2,s;
		ss1 << PL1;
		ss1 >> s1;
		ss2 << PL2;
		ss2 >> s2;
		s = (s1 + "-" + s2);
		index_edge[s] = i;
		s = (s2 + "-" + s1);
		index_edge[s] = i;

		i++;
	}
	myfile.clear();
	myfile.getline(buf, 1024);
	myfile.getline(buf, 1024);
	myfile.getline(buf, 1024);
	i = 0;
	while (true)
	{
		myfile.get(extra, 3); 
		myfile >> PL1 >> PL2 >> PL3;
		myfile >> d;
		myfile.get(extra, 3);
		myfile >> m;
		myfile.get(extra, 3); 
		myfile >> s;
		myfile.get(extra, 3);
		myfile.get();
		if (!myfile.good())
		{
			break;
		}
		angle myangle = angle(PL1, PL2, PL3, d, m, s, i);
		angles.push_back(myangle);
		//cout << PL1 <<" " << PL2 << " " << PL3 << " " << d << " " << m << " " << s << " " << endl;
		i++;
	}
	double X, Y;
	char dot;
	myfile.clear();
	myfile.getline(buf, 19);
	myfile.get();
	myfile.clear();

	for (int m = 0; m < all_point_num; m++)
	{
		point p;
		p.id = m;
		p.pre_x = 0;
		p.pre_y = 0;
		p.after_x = 0;
		p.after_y = 0;
		p.isKnown = false;
		p.isCalcu = false;
		points.push_back(p);
	}
	for (int m = 0; m < known_point_num; m++)
	{
		myfile >> PL1;;
		myfile.get(extra, 3);
		myfile >> X >> dot >> Y >> dot;
		points[PL1].isCalcu = true;
		points[PL1].isKnown = true;
		points[PL1].pre_x = X;
		points[PL1].pre_y = Y;
		points[PL1].after_x = X;
		points[PL1].after_y = Y;
		known_point.push_back(PL1);
		//cout << PL1 <<" " << " "<< X <<" "<<Y <<endl;
	}
	myfile >> PL1 >> PL2 >> the3;
	start1 = PL1;
	start2 = PL2;
}

void approZB(int P1, int P2)
{
	for (int i = 0; i < angles.size(); i++)
	{
		if (P2 == angles[i].middleP)
		{
			if (P1 == angles[i].startP)
			{
				int P3 = angles[i].endP;
				if (points[P3].isCalcu == true)
				{
					if (points[P3].isKnown == true && points[P2].isKnown == false)
					{
						int t;
						t = P1; P1 = P2; P2 = P3;
						approZB(P1, P2);
						P2 = P1;
						P1 = t;
						continue;
					}
					else
					{
						continue;  //若此点已算过，则从下一个角开始寻找
					}
				}
				else
				{
					// 由于题目给出坐标的xy值相反，故这里，左角右角也要相反
					double F_angle = azimuth(P1, P2);
					double old = F_angle - angles[i].pre.toDeg() + 180;
					double S = edges[getIndexOfEdge(P2, P3)].pre_value;
					double delta_x = S * cos(DegtoReg(old));
					double delta_y = S * sin(DegtoReg(old));

					points[P3].pre_x = points[P2].pre_x + delta_x;
					points[P3].pre_y = points[P2].pre_y + delta_y;
					points[P3].isCalcu = true;

					int t;
					t = P1; P1 = P2; P2 = P3;
					approZB(P1, P2);
					P2 = P1;
					P1 = t;
					continue;
				}
			}
			else if(P1 == angles[i].endP)
			{
				int P3 = angles[i].startP;
				if (points[P3].isCalcu == true)
				{
					if (points[P3].isKnown == true && points[P2].isKnown == false)
					{
						int t;
						t = P1; P1 = P2; P2 = P3;
						approZB(P1, P2);
						P2 = P1;
						P1 = t;
						continue;
					}
					else
					{
						continue;  //若此点已算过，则从下一个角开始寻找
					}
					
				}
				else
				{
					// 由于题目给出坐标的xy值相反，故这里，左角右角也要相反
					double F_angle = azimuth(P1, P2);
					double old = F_angle + angles[i].pre.toDeg() - 180;
					double S = edges[getIndexOfEdge(P2, P3)].pre_value;
					double delta_x = S * cos(DegtoReg(old));
					double delta_y = S * sin(DegtoReg(old));

					points[P3].pre_x = points[P2].pre_x + delta_x;
					points[P3].pre_y = points[P2].pre_y + delta_y;
					points[P3].isCalcu = true;

					int t;
					t = P1; P1 = P2; P2 = P3;
					approZB(P1, P2);
					P2 = P1;
					P1 = t;
					continue;
				}
			}
		}
	}
}

double azimuth(int P1, int P2)
{
	double x1 = points[P1].pre_x;
	double y1 = points[P1].pre_y;
	double x2 = points[P2].pre_x;
	double y2 = points[P2].pre_y;
	double theta = 0;
	if (x1 == x2)
	{
		if (y1 > y2)
		{
			theta = 270;
		}
		else
		{
			theta = 90;
		}
		return theta;
	}
	theta = atan(abs(y2 - y1) / abs(x2 - x1));
	theta = RegtoDeg(theta);
	if (x1 < x2 && y1 > y2)
	{
		theta = 360 - theta;
	}
	else if(x1 > x2 && y1 < y2)
	{
		theta = 180 - theta;
	}
	else if (x1 > x2 && y1 > y2)
	{
		theta = 180 + theta;
	}
	return theta;
}

int getIndexOfEdge(int P1, int P2)
{
	stringstream ss1, ss2;
	string s1, s2, s;
	ss1 << P1;
	ss1 >> s1;
	ss2 << P2;
	ss2 >> s2;
	
	map<string, int>::iterator l_it;
	
	s = (s1 + "-" + s2);
	l_it = index_edge.find(s);
	if (l_it == index_edge.end())
	{
		s = (s2 + "-" + s1);
		l_it = index_edge.find(s);
		if (l_it == index_edge.end())
		{
			return -1;
		}
		else
		{
			return index_edge[s];
		}
	}
	else
	{
		return index_edge[s];
	}
}

void calcuB_l_P(Matrix & B, Matrix & l, Matrix & P, Matrix & L0, int & known_edge)
{
	for (int i = 0; i < angles.size(); i++)
	{
		double a1, a2, b1, b2;
		int start = angles[i].startP;
		int middle = angles[i].middleP;
		int end = angles[i].endP;
		calcuAB(middle, start, a1, b1);
		calcuAB(middle, end, a2, b2);

		if (!points[middle].isKnown )
		{
			B.M[i + 1][2 * middle + 1] = a1 - a2;
			B.M[i + 1][2 * middle + 2] = b1 - b2;
		}
		if (!points[start].isKnown)
		{
			B.M[i + 1][2 * start + 1] = -a1;
			B.M[i + 1][2 * start + 2] = -b1;
		}
		if (!points[end].isKnown)
		{
			B.M[i + 1][2 * end + 1] = a2;
			B.M[i + 1][2 * end + 2] = b2;
		}
		//角度以秒为单位
		double x = interAngle(start, middle, end);
		//l.M[i + 1][1] = (angles[i].pre.toDeg() - interAngle(start, middle, end)) * 3600;
		l.M[i + 1][1] = (angles[i].pre.toDeg() - x) * 3600;
		L0.M[i + 1][1] = angles[i].pre.toDeg();
		P.M[i + 1][i + 1] = 1;
	}

	for (int i = 0; i < edges.size(); i++)
	{
		int P1, P2;
		double a, b, S;
		P1 = edges[i].Pid1;
		P2 = edges[i].Pid2;
		calcuAB_forS(P1, P2, a, b, S);
		int from = angles.size() + 1;
		if (points[P1].isKnown && points[P2].isKnown)
		{
			known_edge = from + i;
			continue;
		}
		if (!points[P1].isKnown)
		{
			B.M[from + i][2 * P1 + 1] = -b;
			B.M[from + i][2 * P1 + 2] = -a;
		}
		if (!points[P2].isKnown)
		{
			B.M[from + i][2 * P2 + 1] = b;
			B.M[from + i][2 * P2 + 2] = a;
		}
		//边以m为单位
		l.M[from + i][1] = (edges[i].pre_value - S) * 1000;
		L0.M[from + i][1] = edges[i].pre_value;
		P.M[from + i][from + i] =  57600 / (S * S);  /********!!!!!!!!!!!*****/
	}
}

void calcuAB(int P1, int P2, double & a, double & b)
{
	double delta_x = points[P2].pre_x - points[P1].pre_x;
	double delta_y = points[P2].pre_y - points[P1].pre_y;
	double S2 = delta_x * delta_x + delta_y * delta_y;
	a = (ruo * delta_y) / S2  / 1000;
	b = -(ruo * delta_x) / S2 / 1000;
}

void calcuAB_forS(int P1, int P2, double & a, double & b, double & S)
{
	double delta_x = points[P2].pre_x - points[P1].pre_x;
	double delta_y = points[P2].pre_y - points[P1].pre_y;
	S = sqrt(delta_x * delta_x + delta_y * delta_y);
	a = delta_y / S;
	b = delta_x / S;
}
double interAngle(int P1, int P2, int P3)
{
	double a1 = azimuth(P1, P2);
	double a2 = azimuth(P2, P3);
	double res = (180 + (a1 - a2));
	
	while(res < 0)
	{
		res += 360;
	}
	while(res >= 360)
	{
		res -= 360;
	}
	return res;
}

/*********************************************/
void calcuQXX_x(Matrix & QXX, Matrix & x, Matrix & B, const Matrix & l, const Matrix & P)
{
	Matrix res;
	//(B.Transpose() * P * B).Display();
	res = (B.Transpose() * P * B).Inverse();
	QXX = res;
	x = res * B.Transpose() * P * l;
}

void calcuX(Matrix & X, const Matrix & x, vector<int>& known_point)
{
	X = Matrix(x.row, 1);
	vector<int>::iterator itt;
	for (int i = 0, j = 0; i < points.size(); i++)
	{
		itt = find(known_point.begin(), known_point.end(), i);
		if (itt == known_point.end())
		{
			points[i].after_x = points[i].pre_x + x.M[2 * j + 1][1] / 1000;
			points[i].after_y = points[i].pre_y + x.M[2 * j + 2][1] / 1000;
			X.M[2 * j + 1][1] = points[i].after_x;
			X.M[2 * j + 2][1] = points[i].after_y;
			j++;
		}
	}
}
void calcuQLL(Matrix & QLL, const Matrix & QXX, Matrix & B)
{
	QLL = B * QXX * B.Transpose();
}
void calcuV_L(Matrix & V, Matrix & L, Matrix & B,const Matrix & L0, const Matrix & x, const Matrix & l)
{
	V = B * x - l;
	for(int i = 1; i <= angles.size(); i++)
	{
		V.M[i][1] = V.M[i][1] / 3600;
	}
	for(int i = angles.size() + 1; i <= V.row; i++)
	{
		V.M[i][1] = V.M[i][1] / 1000;
	}
	L = V + L0;
	for (int i = 1; i <= angles.size(); i++)
	{
		V.M[i][1] = V.M[i][1] * 3600;
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
	//double delta2 = delta0 * delta0;
	for (int i = 1; i <= QXX.row; i++)
	{
		D_X.M[i][1] = sqrt(QXX.M[i][i]) * delta0;
	}
}

void calcuDelta_L(Matrix & D_L, Matrix & QLL, double delta0)
{
	D_L = Matrix(QLL.row, 1);
	//double delta2 = delta0 * delta0;
	for (int i = 1; i <= QLL .row; i++)
	{
		D_L.M[i][1] = sqrt(QLL.M[i][i]) * delta0;
	}
}

void calcuE_F_faiE_faiF(int i, const double & delta0, const Matrix & QXX, double & faiE, double & faiF, double & E, double & F)
{
	double Qxx = QXX.M[2 * i - 1][2 * i - 1];
	double Qyy = QXX.M[2 * i][2 * i];
	double Qxy = QXX.M[2 * i - 1][2 * i];
	double K = sqrt((Qxx - Qyy)*(Qxx - Qyy) + 4 * Qxy * Qxy);
	double QEE = (Qxx + Qyy + K) / 2;
	double QFF = (Qxx + Qyy - K) / 2;

	E = delta0 * sqrt(QEE);
	F = delta0 * sqrt(QFF);

	if (Qxy != 0)
	{
		faiE = atan((QEE - Qxx) / Qxy);
		faiE = RegtoDeg(faiE);
		if (faiE < 0)
		{
			faiE = 180 + faiE;
		}
		faiF = atan((QFF - Qxx) / Qxy);
		faiF = RegtoDeg(faiF);
		if (faiF < 0)
		{
			faiF = 180 + faiF;
		}
	}
	else
	{
		if (Qxx > Qyy)
		{
			faiE = 0;
			faiF = 90;
		}
		else if (Qxx < Qyy)
		{
			faiE = 90;
			faiF = 0;
		}
		else
		{
			faiE = -1;
			faiF = -1;
		}
	}
}
/* Copyright by zyt 2016.2 | All rights reserved */
