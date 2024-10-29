#include "bsplinecdata.h"
#include "linearsolve.h"
const double tol = 0.000001;



// constructors for the BsplineCdata class
BsplineCdata::BsplineCdata() : cdat(), b() {}

BsplineCdata::BsplineCdata(const Cdata& data) : cdat(data), b() 
{
	computeBsplineInterp();
}

// accessor methods
BsplineCurve BsplineCdata::getBspline() const
{
	return b;
}


// computes the approximating Bspline curve fitting the data points
Matrix<double> BsplineCdata::computeLeastSquaresMatrix(int nseg, int ord) const
{

	Vector<double> param = cdat.getParam();
	int num = cdat.getNum();

	Vector<double> temp(nseg + 1);
	for (int i = 0; i <= nseg; i++) temp[i] = param[0] + i*(param[num - 1] - param[0]) / nseg;

	Vector<double> knots = computeKnots(temp, nseg + 1, ord);
	int dim = nseg - 1 + ord;
	Matrix<double> vmat(num, dim);
	for (int i = 0; i<num; i++) {
		Vector<double> bval(createVectorInterp(dim + ord, ord, knots, param[i]));
		for (int j = 0; j<dim; j++) vmat[i][j] = bval[j];
	}
	double sum = 0.0;
	Matrix<double> nvmat(dim, dim);
	for (int i = 0; i<dim; i++)
	for (int j = 0; j<dim; j++) {
		for (int k = 0; k<num; k++) {
			double temp = vmat[k][i]*vmat[k][j];
			sum = sum + temp;
		}
		nvmat[i][j] = sum;
		sum = 0.0;
	}

	return nvmat;
}


// computes the right hand side vector for the least squares linear equations
Vector<Point3D> BsplineCdata::computeRHSVectorLeastSq(int nseg, int ord) const
{
	Vector<Point3D> pts = cdat.getPoints();
	Vector<double> param = cdat.getParam();
	int num = cdat.getNum();

	Vector<double> temp(nseg + 1);
	for (int i = 0; i <= nseg; i++) temp[i] = param[0] + i*(param[num - 1] - param[0]) / nseg;
	Vector<double> knots = computeKnots(temp, nseg + 1, ord);
	int dim = nseg - 1 + ord;
	Matrix<double> mat(num, dim);
	for (int i = 0; i<num; i++) {
		Vector<double> bval(createVectorInterp(dim + ord, ord, knots, param[i]));
		for (int j = 0; j<dim; j++) mat[i][j] = bval[j];
	}

	Vector<Point3D> b(dim);
	Point3D sum1 = 0.0;
	for (int i = 0; i<dim; i++) {
		for (int j = 0; j<num; j++) {
			Point3D temp = mat[j][i] * pts[j];
			sum1 = sum1 + temp;
		}
		b[i] = sum1;
		sum1 = 0.0;
	}
	return b;
}


// computes the matrix for cubic Bspline interpolation
Matrix<double> BsplineCdata::computeInterpolationMatrix() const
{
	Vector<double> param = cdat.getParam();
	int num = cdat.getNum();
	int dim = num + 2;
	Matrix<double> mat(dim, dim);
	Vector<double> knots = computeKnots(param, num, 4);

	mat[0][0] = -3.0;
	mat[0][1] = 3.0;
	mat[dim - 1][dim - 2] = -3.0;
	mat[dim - 1][dim - 1] = 3.0;

	for (int i = 1; i<dim - 1; i++) {
		Vector<double> temp = createVectorInterp(dim + 4, 4, knots, param[i - 1]);
		for (int j = 0; j<dim; j++) mat[i][j] = temp[j];
	}
	return mat;
}


// computes the right hand side vector for cubic Bspline interpolation
Vector<Point3D> BsplineCdata::computeRHSVectorInterp() const
{
	Vector<double> param = cdat.getParam();
	Vector<Point3D> pts = cdat.getPoints();
	Vector<Point3D> tang = cdat.getTang();
	int num = cdat.getNum();
	int dim = num + 2;
	Vector<Point3D> rhs(dim);

	rhs[0] = (param[1] - param[0])*tang[0];
	rhs[dim - 1] = (param[num - 1] - param[num - 2])*tang[num - 1];

	for (int i = 1; i <= dim - 2; i++) rhs[i] = pts[i - 1];

	return rhs;
}


// computes all the basis function values at the given point val
Vector<double> BsplineCdata::createVectorInterp(int num, int ord, const Vector<double>& knots, double val) const
{
	Matrix<double> v(ord, ord);
	Vector<double> dp(ord), dm(ord);

	int ind = findIndex(ord, num - ord, knots, val);

	v[0][0] = 1.0;
	for (int j = 0; j<ord - 1; j++) {
		dp[j] = knots[ind + j + 1] - val;
		dm[j] = val - knots[ind + 1 - j - 1];
		for (int i = 0; i <= j; i++) {
			double m = v[i][j] / (dp[i] + dm[j - i]);
			v[i][j + 1] = v[i][j + 1] + dp[i] * m;
			v[i + 1][j + 1] = dm[j - i] * m;
		}
	}
	Vector<double> res(num - ord);
	for (int i = 0; i<ord; i++) res[i + ind - ord + 1] = v[i][ord - 1];
	return res;
}

// find the segment number where t0 lies
int BsplineCdata::findIndex(int ord, int dim, const Vector<double>& knots, double t0) const
{
	/* find the correct segment for t0 */
	int i = 0;

	/* check to see that is in the required range
	interval */
	if ((t0 - knots[ord - 1]) <= tol) return ord - 1;
	else if ((knots[dim] - t0) <= tol) return dim - 1;

	else {

		/* traverse the array of knots until knots[i] <= t0 */
		/* ( ==> (t0 - knots[i]) <= tol). */

		do
		{
			i++;
		} while ((knots[i - 1] - t0) <= tol);
		return (i - 2);
	}
}


// computes the knot vector from the parameter values 
Vector<double> BsplineCdata::computeKnots(const Vector<double>& param, int n, int ord) const
{
	int i;

	/* compute knots from the parameter values of the data points */
	Vector<double> knots(n + 2 * (ord - 1));

	/* left end knots */
	for (i = 0; i <= ord - 1; i++) knots[i] = param[0];

	/* interior knots */
	for (i = 1; i <= n - 2; i++) knots[i + ord - 1] = param[i];

	/* right end knots */
	for (i = 0; i <= ord - 1; i++) knots[n - 2 + ord + i] = param[n - 1];

	return knots;
}


// computes the cubic BsplineCurve that interpolates the data set
void BsplineCdata::computeBsplineInterp()
{
	// TO DO
	// 1. declare a Matrix<double> named mat and initiase with the result of calling the method computeInterpolationMatrix()

	// 2. declare a Vector<Point3D> named rhs and initialise with the result of calling the method computeRHSVectorInterp()

	// 3 declare a LinearSOlve object named lsolve and initialise with mat and rhs from 1 and 2

	// 4. declare a Vector<Point3D> named cpoints and initialise with the result of calling getSolution() for lsolve

	// 5. declare a Vector<double> named knots and initialise with the result of calling computeKnots(cdat.getParam(), cdat.getNum(), 4)

	// 6. assign the data member b to the result of calling the BsplineCurve constructor : BsplineCurve(knots[3], knots[knots.size() - 4], 4, knots.size() - 4, cpoints, knots)
}


// computes the least squares BsplineCurve, consisting of nseg segments and order given by ord, that approixmates the data set
void BsplineCdata::computeBsplineLeastSq(int nseg, int ord)
{
	// TO DO
}