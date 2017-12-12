#include <iostream>
using namespace std;
#define Error 1e-3

void  scanf_data(double upper_lower[], double ** &x_fx, int &m,int &n);		//输入函数
double  solution(double upper_lower[], double ** x_fx, int m, int n);
double **calc_a(double ** x_fx, int  m, int n);
double data_a(double **x_fx, int m, int n);
double * calc_b(double ** x_fx, int  m, int n);
double data_b(double **x_fx, int m, int n);
double * solution_x(double **a, double *b, int n);
void LU_decom(double **&a, double **&l, double **&u, int n);
double * Calc_X(double **&l, double **&u, double *&b, int n);

int main()
{
	double upper_lower[2] = { 0 };				// 上下界
	double **x_fx = NULL;						// x,fx的二维数组
	int n = 1;									// 方程组的次数
	int m = 0;									// 观测数据的个数
	scanf_data(upper_lower, x_fx, m,n);

	while (solution(upper_lower, x_fx, m, n) > Error)
	{
		n++;
		
	}
	cout << "第" << n << "次满足误差" << endl;
	
	getchar();
	getchar();
	return 0;
}

void  scanf_data(double upper_lower[], double ** &x_fx, int &m, int &n)
{

	cout << "********请输入数值的个数:**********\n";
	cin >> m;

	cout << "********请输入方程的次数:**********\n";
	cin >> n;

	x_fx = new double *[2];

	cout << "********请先输入" << m << "个x,然后输入" << m << "个f(x):**********\n";


	for (int i = 0; i < 2; i++)
	{
		x_fx[i] = new double[m];
		for (int j = 0; j < m; j++)
		{
			cin >> x_fx[i][j];
		}
	}
	upper_lower[0] = x_fx[1][0];
	upper_lower[1] = x_fx[1][m - 1];
	//cout << "\n\n********请输入积分上下界:**********\n";

	//for (size_t i = 0; i < 2; i++)
	//{
	//	cin >> upper_lower[i];
	//}
}

double solution(double upper_lower[], double ** x_fx, int m, int n)
{
	double **a = NULL;
	double *b = NULL;
	double *x = NULL;
	cout << "----------"<<n<<"--------\n";
	a = calc_a(x_fx, m, n);
	b = calc_b(x_fx, m, n);
	cout << "\n**a**\n";
	for (int i = 0; i < n + 1; i++)
	{
		for (int j = 0; j < n + 1; j++)
		{
			cout << a[i][j] << " ";
		}
		cout << endl;
	}

	cout << "*****b*****\n";
	for (int i = 0; i < n + 1; i++)
	{
		cout << b[i] << " ";
	}
	cout << endl;

	x = solution_x(a, b, n + 1);


	double error = 0;

	for (int  i = 0; i < m; i++)
	{
		error += pow(x_fx[1][i], 2);
	}

	double temp_error = 0;
	for ( int i = 0; i <n+1 ; i++)
	{
		temp_error += x[i] * b[i];
	}
	cout << "\nerror = "<<error-temp_error << endl;

	return error - temp_error;
	
}

double* calc_b(double ** x_fx, int  m, int n)
{
	double *b = NULL;
	int degree = n + 1;
	b = new double[degree];
	for (int i = 0; i < degree; i++)
	{
		b[i] = data_b(x_fx, m, i);
	}
	return b;
}

double data_b(double **x_fx, int m, int n)
{
	double sum = 0;
	for (int i = 0; i < m; i++)
	{
		sum += x_fx[1][i] * pow(x_fx[0][i], n);
	}
	return sum;
}

double** calc_a(double ** x_fx, int  m, int n)
{
	double **a = NULL;
	int degree = n + 1;
	a = new double *[degree];
	for (int i = 0; i < degree; i++)
	{
		a[i] = new double[degree];
	}

	//cout << "\n*****a*****\n ";
	for (int i = 0; i < degree; i++)
	{
		for (int j = 0; j < degree; j++)
		{
			a[i][j] = data_a(x_fx, m, i + j);
			//cout << a[i][j] << " ";
		}
		//cout << endl;
	}
	return a;
}

double data_a(double **x_fx, int m, int n)
{
	double sum = 0;
	for (int i = 0; i < m; i++)
	{
		sum += pow(x_fx[0][i], n);
	}
	return sum;
}

double * solution_x(double **a, double *b, int n)
{
	double **L = NULL, **U = NULL,*x = NULL;
	//LU分解
	LU_decom(a, L, U, n);
	//根据LU计算x的值
	x = Calc_X(L, U, b, n);
	return x;
}

void LU_decom(double **&a, double **&l, double **&u, int n)
{
	//初始化空间l u
	l = new double *[n];
	u = new double *[n];
	for (int i = 0; i < n; i++)
	{
		l[i] = new double[n];
		u[i] = new double[n];
		for (int j = 0; j < n; j++)
		{
			l[i][j] = 0;
			u[i][j] = 0;
			if (i == j)
				l[i][j] = 1;
		}
	}

	//计算l u矩阵的结果
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			if (j == 0 && i != 0)
				l[i][j] = a[i][j] / u[0][0];
			else if (i == 0)
				u[i][j] = a[i][j];
			else {
				int min_temp = i > j ? j : i;

				if (i > j)
				{
					double sum = 0;
					for (int k = 0; k < min_temp; k++)
					{
						sum += l[i][k] * u[k][j];
					}
					l[i][j] = (a[i][j] - sum) / u[min_temp][min_temp];
				}
				else {
					double sum = 0;
					for (int k = 0; k < min_temp; k++)
					{
						sum += l[i][k] * u[k][j];
					}
					u[i][j] = (a[i][j] - sum);
				}
			}
			//cout << l[i][j] << " ";
		}
		//cout << endl;
	}

}

double * Calc_X(double **&l, double **&u, double *&b, int n)
{
	double *y = NULL, *x = NULL;
	x = new double[n];
	y = new double[n];
	//计算y
	for (int i = 0; i < n; i++)
	{

		double sum_x = 0;
		for (int j = 0; j < i; j++)
		{
			sum_x += l[i][j] * y[j];
		}
		y[i] = b[i] - sum_x;
		//cout << y[i] << " ";
	}

	//计算x

	for (int i = n - 1; i >= 0; i--)
	{

		double sum_y = 0;
		for (int j = n - 1; j > i; j--)
		{
			sum_y += u[i][j] * x[j];
		}
		x[i] = (y[i] - sum_y) / u[i][i];
	}

	/*cout << "x的值为：\n";
	for (int i = 0; i <n; i++)
	{
		cout << "\tx" << "[" << i << "]=" << x[i];
	}*/
	return x;
}