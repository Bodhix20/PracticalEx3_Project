#ifndef _LSOLVE
#define _LSOLVE
#include "matrix.h"


// class encapsulating the solution to a linear systen of equations, uses LU factorisation with pivoting
template<typename T>
class LinearSolve {
private:
	int n; // size of the syatem;
	Matrix<double> a; // matrix
	Vector<T> b;  // right hand side vector
	Vector<T> x;  // solution vector

	// factorisation and solve methods
	void lu_fact(Matrix<double>& l, Matrix<double>& u, const Matrix<double>& p) const;
	void lu_solve(const Matrix<double>& l, const Matrix<double>& u, const Matrix<double>& p);

	// used to create the permutation matrix for pivoting the equations
	Matrix<double> reorder() const;
public:
	// constructors
	LinearSolve();
	LinearSolve(const Matrix<double>& mat, const Vector<T>& rhs, int num);

	// returns the solution vector
	Vector<T> getSolution() const;
};

// constructors
template<typename T>
LinearSolve<T>::LinearSolve() : n(0), a(), b(), x() {}

// code within this constructor generates the solution vector x
template<typename T>
LinearSolve<T>::LinearSolve(const Matrix<double>& mat, const Vector<T>& rhs, int num) : n(num), a(mat), b(rhs), x(n) 
{
	Matrix<double> l(n, n), u(n, n);
	Matrix<double> p = reorder();
	lu_fact(l, u, p);
	lu_solve(l, u, p);
}

// accessor method
template<typename T>
Vector<T> LinearSolve<T>::getSolution() const
{
	return x;
}


/* LU factorisation routine
* Takes in a matrix a of size n and produces the lower (l) and
* upper (u) triangular matrices that factorise a
*/
template<typename T>
void LinearSolve<T>::lu_fact(Matrix<double>& l, Matrix<double>& u, const Matrix<double>& p) const
{
	double mult;

	Matrix<double> temp(p*a);

	// LU decomposition without pivoting
	for (int k = 0; k < n - 1; k++) {
		for (int i = k + 1; i < n; i++) {
			if (fabs(temp[k][k]) < 1.e-07) throw std::out_of_range("zero pivot found");
			mult = temp[i][k]/temp[k][k];
			temp[i][k] = mult;                      // entries of L are saved in temp
			for (int j = k + 1; j < n; j++) { 
				temp[i][j] -= mult*temp[k][j];      // entries of U are saved in temp
				if (fabs(temp[i][i]) < 1.e-07) throw std::out_of_range("zero pivot found");
			}
		}
	}

	// create l and u from temp
	for (int i=0; i<n; i++) l[i][i] = 1.0;
	for (int i=1; i<n; i++) 
	  for (int j=0; j<i; j++) l[i][j] = temp[i][j];

	for (int i=0; i<n; i++)
		for (int j=i; j<n; j++) u[i][j] = temp[i][j];
}


/*
* Solves the equation LUx = b by performing forward and backward
* substitution. Output is the solution vector x
*/
template<typename T>
void LinearSolve<T>::lu_solve(const Matrix<double>& l, const Matrix<double>& u, const Matrix<double>& p)
{
	Vector<T> temp(p*b);

	// forward substitution for L y = b.
	for (int i = 1; i < n; i++) 
		for (int j = 0; j < i; j++) 
			temp[i] = temp[i] - l[i][j]*temp[j];
	
  
	// back substitution for U x = y.  
	for (int i = n - 1; i >= 0; i--) {
		for (int j = i + 1; j < n; j++) temp[i] = temp[i] - u[i][j]*temp[j];
		temp[i] = temp[i]/u[i][i];
	}

	// copy solution into x
	for (int i=0; i<n; i++) x[i] = temp[i];
}

/*
* Computes the permutation matrix P such that the matrix PA can be
* factorised into LU and the system PA = Pb can be solved by forward and
* backward substitution. Output is the permutation matrix P
*/
template<typename T>
Matrix<double> LinearSolve<T>::reorder() const
{

	// Note: pivoting information is stored in temperary vector pvt
	int pvtk, pvti;
	double aet, tmp, mult;
	Matrix<double> p(n, n);
	Vector<int> pvt(n);
	Matrix<double> temp(a);

	for (int k = 0; k < n; k++) pvt[k] = k;

	Vector<double> scale(n);             // find scale vector
	for (int k = 0; k < n; k++) {
		scale[k] = 0;
		for (int j = 0; j < n; j++) 
			if (fabs(scale[k]) < fabs(temp[k][j])) scale[k] = fabs(temp[k][j]);
	} 

	for (int k = 0; k < n - 1; k++) {            // main elimination loop

	// find the pivot in column k in rows pvt[k], pvt[k+1], ..., pvt[n-1]
		int pc = k; 
		aet = fabs(temp[pvt[k]][k]/scale[k]);
		for (int i = k + 1; i < n; i++) {
			tmp = fabs(temp[pvt[i]][k]/scale[pvt[i]]); 
			if (tmp > aet) {
				aet = tmp; 
				pc = i;
			}
		}
		if (fabs(aet) < 1.e-07) throw std::out_of_range("zero pivot found"); 
		if (pc != k) {                      // swap pvt[k] and pvt[pc]
			int ii = pvt[k];
			pvt[k] = pvt[pc];
			pvt[pc] = ii;
		}

		// now eliminate the column entries logically below mx[pvt[k]][k]
		pvtk = pvt[k];                           // pivot row
		for (int i = k + 1; i < n; i++) {
			pvti = pvt[i];
			if (temp[pvti][k] != 0) {
				mult = temp[pvti][k]/temp[pvtk][k]; 
				temp[pvti][k] = mult;
				for (int j = k + 1; j < n; j++) temp[pvti][j] -= mult*temp[pvtk][j];
			}
		} 
	}
	for (int i=0; i<n; i++) p[i][pvt[i]]=1.0;

	return p;
}

#endif