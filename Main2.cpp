#include "BsplineCdata.H"

using namespace std;

int main() {


	ifstream ifs("pointset2.dat");
	int num1, num2;
	ofstream ofs("results.dat");

	//reading the first two numbers (size of the matrix)
	ifs >> num1 >> num2;

	//Declaring a matrix of 3d points
	Matrix<Point3D> mat(num1,num2);

	//reading the file and putting the data in the matrix
	for (int i = 0; i < num1; i++) {
		for (int j = 0; j < num2; j++) {
			Point3D pt;
			ifs >> pt;
			mat[i][j] = pt;
		}
	}

	Vector<Point3D> vec(num2);

	//then we put the matrix info in a vector of 3d points
	for (int i = 0; i < num1; i++) {
		cout << "processing row " << i << "of the matrix\n";
		for (int j = 0; j < num2; j++) {
			vec[j] = mat[i][j];
		}		
		Cdata cdat(num2, vec);
		BsplineCdata bsp(cdat);
		ofs << bsp.getBspline();
	}



	return 0;
}