#ifndef _POINT
#define _POINT

#include "TextObject.h"
#include "vector.h"

// class encapsulating a 3D (x,y,z) point
class Point3D :  public TextObject {
private:
	double x;
	double y;
	double z;
public:
	// constructors
	Point3D();
	Point3D(double X, double Y, double Z);
	Point3D(double d);
	Point3D(const Point3D& p);

	// accessor methods
	double getX() const;
	double getY() const;
	double getZ() const;

	// overloaded operators
	Point3D& operator=(const Point3D& p);
	friend Point3D operator*(double d, const Point3D& p);
	friend Point3D operator*(const Point3D& p, double d);
	friend Point3D operator+(const Point3D& p1, const Point3D& p2);
	friend Point3D operator-(const Point3D& p1, const Point3D& p2);
	friend Point3D operator/(const Point3D& p1, double d);

	// read/write methods
	virtual void readfile(std::ifstream& ifs);
	virtual void writefile(std::ofstream& ofs) const;
	virtual void write(std::ostream& os) const;
	virtual void read(std::istream &is);
};

#endif
