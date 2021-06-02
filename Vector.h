
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

	Vector transform(Transform *t){
		vector<double> a(3, 1);
		for(int i = 0; i < 3; i++){
			a[i] = t->matrix[i][0] * x + t->matrix[i][1] * y + t->matrix[i][2] * z; 
		}
		return Vector(a[0], a[1], a[2]);
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
