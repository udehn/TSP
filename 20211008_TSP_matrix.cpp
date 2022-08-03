#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <time.h>
#include <string.h>
#include <queue>
#include <chrono>
#include <stack>
#include <windows.h>
#include <thread>

using namespace std;
using chrono::duration;
using chrono::duration_cast;
using chrono::high_resolution_clock;
// upper bound (PS: min prize )
double h_upper_bound = INT_MAX;

//funtion for getting prize now
void get_h_now(double T) {
	Sleep(0.01*T);
	cout << h_upper_bound << endl;	
	Sleep(0.03*T);
	cout << h_upper_bound << endl;
	Sleep(0.1*T);
	cout << h_upper_bound << endl;
	Sleep(0.3*T);
	cout << h_upper_bound << endl;
	Sleep(0.5*T);
	cout << h_upper_bound << endl;
}

//Nodes used for analysis 
class matrix {
private:
	int row;
	int col; 
	double h; // prize now
	double* matr; // prize matrix
	int* path; // paths
	int num; // numbers of paths
public:
	matrix();
	matrix(int n, int m);
	matrix(const matrix& m1);
	void InitPlace();
	int rows() const;
	int columns() const;
	void InitRnd();
	void Init0();
	void Set(int nX, int nY, double r);
	double Get(int nX, int nY) const;
	double Geth() const;
	void Seth(double r);
	void Set_path(int index, int data);
	int Get_path_point(int i)const;
	void Set_num(int i);
	int Get_num()const;
	~matrix() {
		delete[] matr;
		delete[] path;
	};
	friend ostream& operator<<(ostream& os, matrix& m);
	void simplify(); // simplify matrix that every row and col have prize equal 0
	double find_min_sum(int i, int j);  // find the min sum prize from row i and col j
	matrix cut_matr_left(int oldI, int oldJ); // use i and j to cut matrix, get matrix don't have path(i.j)
	matrix cut_matr_right(int oldI, int oldJ); // use i and j to cut matrix, get matrix have path(i.j)
};

matrix::matrix() {
	row = 1;
	col = 1;
	h = 0;
	num = 0;
}

matrix::matrix(int n, int m) {
	row = n;
	col = m;
	h = 0;
	num = 0;
	matr = new double[row * col];
	path = new int[row * 2];
}

matrix::matrix(const matrix& m1) {
	row = m1.rows();
	col = m1.columns();
	h = m1.Geth();
	matr = new double[row * col];
	for (int i = 1; i <= row; i++) for (int j = 1; j <= col; j++)
		matr[(i - 1) * col + (j - 1)] = m1.Get(i, j);
	path = new int[row * 2];
	for (int i = 0; i < row * 2; i++) {
		path[i] = m1.Get_path_point(i);
	}
	num = m1.Get_num();
}

void matrix::InitPlace() {
	matr = new double[row * col];
	path = new int[row * 2];
}

int matrix::rows() const { return row; }

int matrix::columns() const { return col; }

void matrix::Set(int nX, int nY, double r) { matr[col * (nX - 1) + nY - 1] = r; }

double matrix::Get(int nX, int nY)const { return matr[col * (nX - 1) + nY - 1]; }

void matrix::Seth(double r) { h = r; }

double matrix::Geth()const { return h; }

void matrix::Set_path(int index, int data) { path[index] = data; }

int matrix::Get_path_point(int i)const { return path[i]; }

void matrix::Set_num(int i) { num = i; }

int matrix::Get_num()const { return num; }

ostream& operator<<(ostream& os, matrix& m) {
	for (int i = 0; i < m.row; i++) {
		for (int j = 0; j < m.col; j++) {
			cout.width(10);
			cout.precision(3);
			os << m.matr[i * m.col + j] << " ";
		}
		os << endl;
	}
	return os;
}

double GetRnd() {
	double x1, x2;
	while (1) {
		x1 = rand() % 100;
		x2 = rand() % 100;
		if (x2 <= (3 / 5.0) * (5 / 6.0) * (1 + pow(x1, 4)))
			return x1;
	}
}

void matrix::InitRnd() {
	for (int i = 0; i < row; i++)
		for (int j = 0; j < col; j++)
			if (i == j)
				matr[i * col + j] = -1;
			else
				matr[i * col + j] = GetRnd();
	for (int i = 0; i < row * 2; i++)
		path[i] = -1;
}

void matrix::Init0() {
	for (int i = 0; i < row; i++)
		for (int j = 0; j < col; j++)
				matr[i * col + j] = -1;
	for (int i = 0; i < row * 2; i++)
		path[i] = -1;
}

// simplify matrix that every row and col have prize equal 0
void matrix::simplify() {
	double pmin;
	int flag;
	for (int i = 1; i <= row; i++) {
		pmin = INT_MAX;
		flag = 0;
		for (int j = 1; j <= col; j++) {
			if (Get(i, j) == 0) {
				flag = 1;
				break;
			}
		}
		if (flag == 0) {
			for (int j = 1; j <= col; j++) {
				if ((Get(i, j) != -1) && (Get(i, j) < pmin))
					pmin = Get(i, j);
			}
			if (pmin != INT_MAX) {
				h += pmin;
				for (int j = 1; j <= col; j++) {
					if (Get(i, j) != -1)
						Set(i, j, Get(i, j) - pmin);
				}
			}
		}
	}
	for (int j = 1; j <= col; j++) {
		pmin = INT_MAX;
		flag = 0;
		for (int i = 1; i <= row; i++) {
			if (Get(i, j) == 0) {
				flag = 1;
				break;
			}
		}
		if (flag == 0) {
			for (int i = 1; i <= row; i++) {
				if ((Get(i, j) != -1) && (Get(i, j) < pmin))
					pmin = Get(i, j);
			}
			if (pmin != INT_MAX) {
				h += pmin;
				for (int i = 1; i <= row; i++) {
					if (Get(i, j) != -1)
						Set(i, j, Get(i, j) - pmin);
				}
			}
		}
	}
}

// find the min sum prize of two numbers; one is from row i and another one is from col j
double matrix::find_min_sum(int i, int j) {
	double minI = INT_MAX, minJ = INT_MAX;
	for (int k = 1; k <= row; k++) {
		if ((k != j) && (minI > Get(i, k)) && (Get(i, k) != -1)) {
			minI = Get(i, k);
		}
		if ((k != i) && (minJ > Get(k, j)) && (Get(k, j) != -1)) {
			minJ = Get(k, j);
		}
	}
	if (minI == INT_MAX)
		minI = 0.0;
	if (minJ == INT_MAX)
		minJ = 0.0;
	return minI + minJ;
}

// use i and j to cut matrix, get matrix don't have path(i.j)
matrix matrix::cut_matr_left(int oldI, int oldJ) {
	matrix m(row, col);
	for (int i = 1; i <= row; i++) {
		for (int j = 1; j <= col; j++)
			m.Set(i, j, Get(i, j));
	}
	m.Seth(Geth());
	m.Set_num(num);
	for (int i = 0; i < row * 2; i++) {
		m.Set_path(i, Get_path_point(i));
	}
	m.Set(oldI, oldJ, -1);
	m.simplify();
	return m;
}


// use i and j to cut matrix, get matrix have path(i.j)
matrix matrix::cut_matr_right(int oldI, int oldJ) {
	matrix m(row, col);
	//m.InitPlace();
	for (int i = 1; i <= row; i++) {
		for (int j = 1; j <= col; j++)
			m.Set(i, j, Get(i, j));
	}
	for (int i = 1; i <= row; i++) {
		m.Set(i, oldJ, -1);
	}
	for (int j = 1; j <= col; j++) {
		m.Set(oldI, j, -1);
	}
	m.Seth(Geth());
	m.Set(oldJ, oldI, -1);
	m.simplify();
	for (int i = 0; i < row * 2; i++) {
		m.Set_path(i, Get_path_point(i));
	}
	m.Set_path(Get_num() * 2, oldI);
	m.Set_path(Get_num() * 2 + 1, oldJ);
	m.Set_num(Get_num() + 1);
	return m;
}

int main()
{
	srand(time(NULL));
	const auto start_time = high_resolution_clock::now();
	double T = 5.73622;
	thread t(get_h_now, T);
	t.detach();
	
	//int n1 = 21;
	int n1 = 35;
	matrix m1(n1, n1);
	//cout << m1 << endl;
	m1.InitRnd();
	m1.simplify();
	stack <matrix> que;
	que.push(m1);
	matrix mBest(n1, n1); // save the best path
	mBest.Init0();

	while (!que.empty()) {
		matrix mTemp = que.top();
		que.pop();
		if (mTemp.Geth() > h_upper_bound)
			continue;

		////////////check if this node is pass all its way or not, and change the best path if it's prize is min
		int flagEmpty = 0;
		for (int i = 1; i <= mTemp.rows(); i++) {
			for (int j = 1; j <= mTemp.columns(); j++) {
				if (mTemp.Get(i, j) != -1) {
					flagEmpty = 1;
					break;
				}
				else continue;
			}
		}
		if (flagEmpty == 0 && mTemp.Get_num() != mTemp.rows())
			continue;
		if (flagEmpty == 0) {
			if (mTemp.Geth() < h_upper_bound) {
				mBest.Seth(mTemp.Geth());
				mBest.Set_num(mTemp.Get_num());
				for (int i = 0; i < mBest.rows() * 2; i++) {
					mBest.Set_path(i, mTemp.Get_path_point(i));
				}
				h_upper_bound = mTemp.Geth();
			}
			continue;
		}
		////////////

		////////////find the path need to use or not, use this index i,j to cut matrix into two part
		double tempSum;
		double maxSum = -1;
		int iToCut = 0, jToCut = 0;
		for (int i = 1; i <= mTemp.rows(); i++) {
			for (int j = 1; j <= mTemp.columns(); j++) {
				if (mTemp.Get(i, j) == 0) {
					tempSum = mTemp.find_min_sum(i, j);
					if (tempSum > maxSum) {
						maxSum = tempSum;
						iToCut = i;
						jToCut = j;
					}
				}
			}
		}
		////////////

		////////////check if it is useful to save set which have this path (i,j)
		int flag = 0;
		int flag_save = 1;
		queue<int> qFlag;
		qFlag.push(jToCut);
		int qtemp;
		int step = 0;
		while (!qFlag.empty()) {
			qtemp = qFlag.front();
			qFlag.pop();
			int i = 0;
			flag = 0;
			for (i = 0; i < mTemp.rows(); i++) {
				if (qtemp == mTemp.Get_path_point(i * 2)) {
					flag = 1;
					step++;
					qFlag.push(mTemp.Get_path_point(i * 2 + 1));
					break;
				}
			}
			if (flag == 0) {
				flag_save = 1;
				break;
			}
			else if (iToCut != mTemp.Get_path_point(i * 2 + 1)) {
				continue;
			}
			if (step < mTemp.rows() - 1) {
				flag_save = 0;
				break;
			}
			else {
				flag_save = 1;
				break;
			}
		}
		////////////

		////////////analyis two set, need to push them into stack
		////////////very important to push the min now prize set late, so that we can find the answer quickly
		if (flag_save == 1) {
			matrix mRight = mTemp.cut_matr_right(iToCut, jToCut);
			if (mRight.Geth() <= h_upper_bound) {
				if (maxSum + mTemp.Geth() <= h_upper_bound) {
					matrix mLeft = mTemp.cut_matr_left(iToCut, jToCut);
					if (mLeft.Geth() < mRight.Geth()) {
						que.push(mRight);
						que.push(mLeft);
					}
					else {
						que.push(mLeft);
						que.push(mRight);
					}
				}
				else {
					que.push(mRight);
				}
			}
			else if (maxSum + mTemp.Geth() <= h_upper_bound) {
				matrix mLeft = mTemp.cut_matr_left(iToCut, jToCut);
				que.push(mLeft);
			}
		}
		else if (maxSum + mTemp.Geth() <= h_upper_bound) {
			matrix mLeft = mTemp.cut_matr_left(iToCut, jToCut);
			que.push(mLeft);
		}
		////////////
	} 

	////////////output the best path and min prize
	cout << "the best tour:" << endl;
	int qtemp = mBest.Get_path_point(0);
	int j = 0;
	for (int i = 0; i < mBest.Get_num(); i++) {
		j = 0;
		for (j = 0; j < mBest.rows(); j++) {
			if (qtemp == mBest.Get_path_point(j * 2)) {
				cout << qtemp << "->";
				break;
			}
			else
				continue;
		}
		qtemp = mBest.Get_path_point(j * 2 + 1);
	}
	cout << qtemp << endl;
	cout << "min prize:" << mBest.Geth() << endl;
	////////////

	const auto end_time = high_resolution_clock::now();
	cout << "Time: " << duration_cast<duration<double, milli>>(end_time - start_time).count() << " ms" << endl;
	return 0;
}
