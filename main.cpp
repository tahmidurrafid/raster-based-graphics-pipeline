#include<bits/stdc++.h>

using namespace std;
#define pi (2*acos(0.0))

enum TRANSFORM{T_translate, T_scale, T_rotate};

class Transform;

class Vector{
public:
	double x, y, z;
	Vector(){
	}
	Vector(double px, double py, double pz){
		x = px, y = py, z = pz;
	}
	void copyIt(Vector a){
		x = a.x, y = a.y, z = a.z;
	}

	Vector cross(Vector a){
		return Vector(y*a.z - z*a.y, z*a.x - x*a.z, x*a.y - y*a.x);
	}
	Vector add(Vector a){
		return Vector(x + a.x, y + a.y, z + a.z);
	}

    double dotProduct(Vector a){
        return x*a.x + y*a.y + z*a.z; 
    }

    Vector reflection(Vector a){
        a = a.multiply( - 2 * dotProduct(a) ).add(*this);
        return a;
    }

	Vector multiply(double val){
		return Vector(val*x, val*y, val*z);
	}

    double distance(Vector a){
        return sqrt( (x - a.x)*(x - a.x) + (y-a.y)*(y-a.y) + (z-a.z)*(z-a.z) );
    }

	void shiftInDirection(Vector to, float steps){
		x += (to.x)*steps;
		y += (to.y)*steps;
		z += (to.z)*steps;
	}

	void rotate(Vector per, float angle){
		Vector t = cross(per);
		Vector m = *this;
		m = m.multiply(cos(angle*pi/180.0));
		t = t.multiply(sin(angle*pi/180.0));
		m = t.add(m);
		copyIt(m);
	}

    void unitify(){
        double mul = sqrt(x*x + y*y + z*z );
        x = x/mul;
        y = y/mul;
        z = z/mul;
    }

	Vector transform(Transform *t);

    void input(ifstream &inp){
        inp >> x >> y >> z;
    }

	void print(){
		cout << x << " " << y << " " << z << "\n";
	}

    bool operator<(const Vector &a) const {
        return x < a.x;
    }
};

class Glu{
public:
    Vector eye, look, up;
    double fovY, aspectRatio, near, far;

    void input(ifstream &inp){
        eye.input(inp);
        look.input(inp);
        up.input(inp);
        inp >> fovY >> aspectRatio >> near >> far;
    }
};

class Triangle{
public:
	Vector p1, p2, p3;

    void input(ifstream &inp){
		p1.input(inp);
		p2.input(inp);
		p3.input(inp);
	}

	Triangle transform(Transform* t){
		Triangle x = *this;
		x.p1 = x.p1.transform(t);
		x.p2 = x.p2.transform(t);
		x.p3 = x.p3.transform(t);
		return x;
	}

	void print(){
		p1.print();
		p2.print();
		p3.print();
	}
};


class Transform{
public:
	vector< vector<double> > matrix;

	Transform(){
		matrix = vector<vector<double>>(4, vector<double>(4, 0));
		matrix[3][3] = 1;
	}

	void setTranslate(Vector t){
		matrix[0][3] = t.x;
		matrix[1][3] = t.y;
		matrix[2][3] = t.z;
	}
	
	void setScale(Vector t){
		matrix[0][0] = t.x;
		matrix[1][1] = t.y;
		matrix[2][2] = t.z;
	}

	void addMatrix(vector<vector<double>> x){
		for(int i = 0; i < 4; i++){
			for(int j = 0; j < 4; j++){
				matrix[i][j] += x[i][j];
			}
		}
	}

	vector<vector<double>> multiMatrix(vector<vector<double>> mat, double x){
		for(int i = 0; i < 4; i++){
			for(int j = 0; j < 4; j++){
				mat[i][j] = x*mat[i][j];
			}
		}
		return mat;
	}

	void multiMatrix(Transform t){
		vector<vector<double>> m(4, vector<double>(4, 0));
		for(int i = 0; i < 4; i++){
			for(int j = 0; j < 4; j++){
				for(int k = 0; k < 4; k++){
					m[i][j] += matrix[i][k] * t.matrix[k][j];
				}
			}
		}
		matrix = m;
	}

	void setRotate(double angle, Vector t){
		t.unitify();
		vector<vector<double>> m1(4, vector<double>(4, 0));
		vector<vector<double>> m2(4, vector<double>(4, 0));
		vector<vector<double>> m3(4, vector<double>(4, 0));
		angle = pi*angle/180.0;

		m1[0][0] = m1[1][1] = m1[2][2] = 1;
		
		m2[0][0] = t.x*t.x , m2[0][1] = t.x*t.y, m2[0][2] = t.x*t.z;
		m2[1][0] = t.x*t.y , m2[1][1] = t.y*t.y, m2[1][2] = t.y*t.z;
		m2[2][0] = t.x*t.z , m2[2][1] = t.y*t.z, m2[2][2] = t.z*t.z;

		m3[0][0] = 0 , m3[0][1] = -t.z, m3[0][2] = t.y;
		m3[1][0] = t.z , m3[1][1] = 0, m3[1][2] = -t.x;
		m3[2][0] = -t.y , m3[2][1] = t.x, m3[2][2] = 0;

		addMatrix( multiMatrix(m1, cos(angle) ) );
		addMatrix( multiMatrix(m2, 1 - cos(angle) ) );
		addMatrix( multiMatrix(m3, sin(angle) ) );
	}

	void input(ifstream &inp, TRANSFORM tr){
		if(tr == T_translate){
			Vector t;
			t.input(inp);
			setTranslate(t);
		}else if(tr == T_scale){
			Vector t;
			t.input(inp);
			setScale(t);
		}else if(tr == T_rotate){
			double angle;
			Vector t;
			inp >> angle;
			t.input(inp);
			setRotate(angle, t);
		}
	}
};

Vector Vector::transform(Transform *t){
	vector<double> a(3, 1);
	for(int i = 0; i < 3; i++){
		a[i] = t->matrix[i][0] * x + t->matrix[i][1] * y + t->matrix[i][2] * z; 
	}
	return Vector(a[0], a[1], a[2]);
}	


int main(){
    ifstream scene;
    scene.open ("scene.txt");
	Glu glu;

	glu.input(scene);
	string command;

	stack<Transform> gluTrans;

	Triangle triangle;
	Transform transform, current;

	while(true){
		scene >> command;
		if(command == "triangle"){
			triangle.input(scene);
			triangle.transform(&current).print();
			cout << "-----------------\n";
		}else if(command == "translate"){
			transform.input(scene, T_translate);
			current.multiMatrix(transform);
		}else if(command == "scale"){
			transform.input(scene, T_scale);
			current.multiMatrix(transform);
		}else if(command == "rotate"){
			transform.input(scene, T_rotate);
			current.multiMatrix(transform);			
		}else if(command == "push"){
			gluTrans.push(current);
		}else if(command == "pop"){
			current = gluTrans.top();
			gluTrans.pop();
		}else if(command == "end"){
			break;
		}else{
			cout << "invalid command: " << command << "\n";
		}
	}	
    return 0;
}