package etomica.space;



/*
 * History
 * Created on Jan 24, 2005 by kofke
 */
public abstract class Vector implements java.io.Serializable, Cloneable { 

        public Object clone() {
    		try {
    			return super.clone();
    		}
    		catch (CloneNotSupportedException ex) {
    			throw new InternalError(ex.toString());
    		}
    	}
        public abstract int D();                              //dimension of the space occupied by vector
        public abstract void assignTo(double[] array);        //converts components to elements of the given array
        public abstract double[] toArray();                   //converts components to array of double
        public abstract boolean equals(Vector v);               //return true if all corresponding elements of this and the given vector are equal
        public abstract void sphericalCoordinates(double[] result); //computes the spherical coordinate representation of the vector; return in the given array to avoid construction of new array with each call
        public double[] toSphericalCoordinateArray() {        //computes spherical coordinates and returns them in a new array
            double[] array = new double[D()];
            sphericalCoordinates(array);
            return array;
        }
        // @deprecated use x 
  //      public double component(int i) {return x(i);}
        // @deprecated use setX 
  //      public void setComponent(int i, double d) {setX(i, d);}
        public abstract double x(int i);              //vector component corresponding to the index i (e.g., i=0, x-component)
        public abstract void setX(int i, double d);   //sets ith component of vector to d
        public abstract void E(Vector u);                     //sets each element of the vector equal to the elements of the vector u
        public abstract void E(double a);                     //sets all components of the vector equal to the constant a
        public abstract void E(double[] a);                   //sets elements of vector to values given in array
        public abstract void E(int[] a);                      //sets elements of vector to values given in array
        public abstract void PE(Vector u);                    //adds (PE is +=) the vector u to this vector
        public abstract void PE(int i, double a);             //adds (+=) a to component i of this vector
        public abstract void PE(double a);                    //adds a constant value to all elements
        public abstract void ME(Vector u);                    //subtracts (-=)
        public abstract void TE(Vector u);                    //multiplies (*=) component-by-component
        public abstract void DE(Vector u);                    //divide (/=) component-by-component
        public abstract void TE(double a);                    //multipies all components by the constant a
        public abstract void TE(int i, double a);             //multiplies "a" times the component i of this vector
        public abstract void DE(double a);                    //divides all components by a
        public abstract void Ea1Tv1(double a, Vector u);      //sets this vector to a*u
        public abstract void Ev1Pa1Tv2(Vector v1, double a1, Vector v2);//sets this vector to v1 + a1*v2
        public abstract void PEa1Tv1(double a, Vector u);     //adds a*u to this vector
        public abstract void Ev1Pv2(Vector u1, Vector u2);    //sets equal to sum of vectors
        public abstract void Ev1Mv2(Vector u1, Vector u2);    //sets equal to difference between vectors
        public abstract Vector P(Vector u);       //adds (+) component-by-component and returns result in another vector
        public abstract Vector M(Vector u);       //subtracts component-by-component and returns result in another vector
        public abstract Vector T(Vector u);       //multiplies component-by-component and returns result in another vector
        public abstract Vector D(Vector u);       //divides component-by-component and returns result in another vector
		public abstract double Mv1Squared(Vector u);    //square of difference between this and given vector
        public abstract void abs();		                      //replaces each component with its absolute value
        public abstract void mod(double a);                   //each component replaced with itself modulo a
        public abstract void mod(Vector u);             //each component replaced with itself modulo component of vector
        public abstract void EModShift(Vector r, Vector u);
        public abstract void EMod2Shift(Vector r, Vector u);
        public abstract double min();                         //returns minimum of all components
        public abstract double max();                         //returns maximum of all components
        public abstract double squared();                     //square-magnitude of vector (e.g., x^2 + y^2)
        public abstract void normalize();                     //scales the vector to unit length
        public abstract boolean isZero();					  //returns true if all elements of vector are zero
        public abstract double dot(Vector u);                 //dot product of this vector with the vector u
        public abstract etomica.space3d.Vector3D cross(etomica.space2d.Vector2D u);       //cross product of this vector with u
        public abstract etomica.space3d.Vector3D cross(etomica.space3d.Vector3D u);       //cross product of this vector with u
        public abstract void XE(etomica.space3d.Vector3D u);            //replaces this vector with its cross product (project result into plane if appropriate)
        public abstract void transform(Tensor A);             //applies the given tensor transformation to this vector
        public abstract void transform(Boundary b, Vector r0, Tensor A);  //applies the transformation to (this - r0)
        public abstract void setRandom(double d);             //
        public abstract void setRandomSphere();               //random point on sphere
        public abstract void setRandomCube();                 //random point in a unit cube
        public abstract void setRandomInSphere();			  //random point in unit sphere
        public abstract void randomStep(double d);            //random step of selected uniformly in cube of edge 2d (adds up to +/- d to present value)
        public final void PEa1Tv1(double[] a, Vector[] u) {   //adds several terms of form a*u to this vector
            for(int i=a.length-1; i>=0; i--) {PEa1Tv1(a[i],u[i]);}
        }
        public etomica.space3d.Vector3D cross(Vector u) {
            if(u instanceof Vector) {return cross((Vector)u);}
            else if(u instanceof Vector) {return cross((Vector)u);}
            else return null;
        }
        public final void setRandomDirection() {setRandomSphere(); normalize();}
        public abstract void randomRotate(double thetaStep);
        public abstract double productOfElements();
    }
