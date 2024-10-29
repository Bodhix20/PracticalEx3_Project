
#include "point.h"
#include <cmath>

// constructors
Point3D::Point3D() : x(0.0), y(0.0), z(0.0) {}
Point3D::Point3D(double d) : x(d), y(d), z(d) {}
Point3D::Point3D(double X, double Y, double Z) : x(X), y(Y), z(Z) {}
Point3D::Point3D(const Point3D& pt) : x(pt.getX()), y(pt.getY()), z(pt.getZ()) {}

// accessor methods
double Point3D::getX() const { return x; }
double Point3D::getY() const { return y; }
double Point3D::getZ() const { return z; }

// overloaded operators
Point3D& Point3D::operator=(const Point3D& p)		// assignment
{
	x = p.x;
	y = p.y;
	z = p.z;
	return *this;
}


Point3D operator*(double d, const Point3D& p)			// multiplication by a double on the left
{
	return Point3D(d*p.getX(),d*p.getY(),d*p.getZ());
}

Point3D operator*(const Point3D& p, double d)			// multiplication by a double on the right
{
	return Point3D(d*p.getX(),d*p.getY(),d*p.getZ());
}

Point3D operator+(const Point3D& p1, const Point3D& p2)		// addition of two point3D objects
{
	return Point3D(p1.getX()+p2.getX(),p1.getY()+p2.getY(),p1.getZ()+p2.getZ());
}

Point3D operator-(const Point3D& p1, const Point3D& p2)   // subtraction of two point3D objects
{
	return Point3D(p1.getX()-p2.getX(),p1.getY()-p2.getY(),p1.getZ()-p2.getZ());
}

Point3D operator/(const Point3D& p, double d)    // division of a point3D object by a double
{
	return Point3D(p.getX()/d,p.getY()/d,p.getZ()/d);
}


// read/write methods
void Point3D::write(std::ostream& os) const		// write to the screen
{
	os << x << " " << y << " " << z << " " << 1.0;
}



void Point3D::read(std::istream& is)		// read from the keyboard
{
	std::cout << "input x, y, z\n";
	is >> x >> y >> z;
}

void Point3D::writefile(std::ofstream& ofs)	const	// write to a file
{
	ofs << x << " " << y << " " << z << " " << 1.0; 
}

void Point3D::readfile(std::ifstream& ifs)		// read from a file
{
	ifs >> x >> y >> z;
}

