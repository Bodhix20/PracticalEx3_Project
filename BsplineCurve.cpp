#include "BsplineCurve.h"
const double tol = 0.000001;


// Bspline constructors
BsplineCurve::BsplineCurve() : AbsCurve<Point3D>(0.0, 0.0), ord(0), dim(0), d(), knots() {}

BsplineCurve::BsplineCurve(double t1, double t2, int k, int num, const Vector<Point3D>& cpts, const Vector<double>& kts) : AbsCurve<Point3D>(t1, t2),
ord(k), dim(num), d(cpts), knots(kts) {}



// accessor methods
Vector<Point3D> BsplineCurve::getCpoints() const
{
	return d;
}

Vector<double> BsplineCurve::getKnots() const
{
	return knots;
}


// overloaded operator for evaluation, uses algorithm in the notes
Point3D BsplineCurve::operator()(double t0) const
{
	Vector<Point3D> temp(ord);
	double lbd;

	/* find the interval for tau0  */
	int ind = findIndex(t0);

	/* compute the point according to de Boor algorithm */
	for (int j = 0; j<ord; j++) temp[j] = d[ind - ord + j + 1];

	for (int i = 0; i<ord - 1; i++)
	for (int j = 0; j<ord - i - 1; j++) {
		lbd = (t0 - knots[ind - ord + i + j + 2]) / (knots[ind + j + 1] - knots[ind - ord + i + j + 2]);
		temp[j] = lbd*temp[j + 1] + (1.0 - lbd)*temp[j];
	}
	return temp[0];
}

// helper method for finding the segment number where the evaluation point t0 lies
int BsplineCurve::findIndex(double t0) const
{
	/* find the correct segment for tau0 */
	int i = 0;

	/* check to see that is in the required range of the B-spline
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

// read/write methods
void BsplineCurve::read(std::istream& is)			// read from the keyboard
{
	std::cout << "input the left and right limits for the curve\n";
	double t1, t2;
	is >> t1 >> t2;
	std::cout << "input the order\n";
	int k;
	is >> k;
	std::cout << "how many control points\n";
	int num;
	is >> num;
	std::cout << "input the points\n";
	Vector<Point3D> cpts(num);
	is >> cpts;
	std::cout << "input the knot values\n";
	Vector<double> t(num + k);
	is >> t;
	(*this) = BsplineCurve(t1, t2, k, num, cpts, t);
}


void BsplineCurve::readfile(std::ifstream& ifs)		// read from a file
{
	double t1, t2;
	ifs >> t1 >> t2;
	int k;
	ifs >> k;
	int num;
	ifs >> num;
	Vector<Point3D> cpts(num);
	ifs >> cpts;
	Vector<double> t(num + k);
	ifs >> t;
	(*this) = BsplineCurve(t1, t2, k, num, cpts, t);
}


void BsplineCurve::write(std::ostream& os)	const	// write to the screen
{
	std::cout << "the left and right limits for the curve are\n";
	std::cout << getLeftLimit() << " " << getRightLimit() << "\n";
	std::cout << "the order of the BsplineCurve is " << ord << "\n";
	std::cout << "there are " << dim << " control points\n";
	std::cout << "the control points are\n";
	std::cout << d;
	std::cout << "the knot values are\n";
	std::cout << knots;
}


void BsplineCurve::writefile(std::ofstream& ofs)	const	// write to a file
{
	ofs << "bspline\n";
	ofs << ord << " " << dim << "\n";
	ofs << "knots\n";
	ofs << knots;
	ofs << "cpoints\n";
	ofs << d;
}

