#include<bits/stdc++.h>

using namespace std;
#define pi (2*acos(0.0))

enum TRANSFORM{T_translate, T_scale, T_rotate};

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

	void setRotate(double angle, Vector t){

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

int main(){
    ifstream scene;
    scene.open ("scene.txt");
	Glu glu;

	glu.input(scene);
	string command;

	Triangle triangle;

	while(true){
		scene >> command;
		if(command == "triangle"){
			triangle.input(scene);
			triangle.p1.print();
			triangle.p2.print();
			triangle.p3.print();
		}else if(command == "translate"){
		}else if(command == "scale"){
		}else if(command == "rotate"){
		}else if(command == "push"){
		}else if(command == "pop"){
		}else if(command == "end"){
			break;
		}else{
			cout << "invalid command: " << command << "\n";
		}
	}	
    return 0;
}