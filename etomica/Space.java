package etomica;

/* History of changes
 * 09/01/02 (DAK) added accelerateTo method to Coordinate
 * 07/10/03 (DAK) added resetV method to CoordinatePair
 * 08/27/03 (DAK) added isZero method to Vector
 * 12/09/03 (DAK) added setRandomInSphere method in Vector
 * 01/22/04 (DAK) added (then removed) CoordinateGroup interface
 */

public abstract class Space implements java.io.Serializable {
    
    public static String VERSION = "01.07.09";

    public final int D;
    private final double rD; // reciprocal of D
    
    public Space(int d) {
        if(d < 1 || d > 3) throw new IllegalArgumentException("Illegal dimension for space");
        D = d;
        rD = 1.0/D;
    }
    
    public int D() {return D;}
    /**
     * Returns the given value raised to the Dth power, where D is the dimension of the space.
     */
    public abstract int powerD(int a);
    /**
     * Returns the given value raised to the Dth power, where D is the dimension of the space.
     */
    public abstract double powerD(double a);
    
    /**
     * Returns the Dth root of the given value, a^(1/D), where D is the dimension of the space.
     */
    public double rootD(double a) {return Math.pow(a, rD);}
    
    public abstract Vector origin();
    
    public abstract Vector makeVector();      //Space.Vector
    public abstract Orientation makeOrientation();
    public abstract Tensor makeTensor();
    public abstract Tensor makeRotationTensor();
    public abstract Coordinate makeCoordinate(Atom a);
    public abstract CoordinatePair makeCoordinatePair();
    public abstract Boundary makeBoundary();  //makes boundary of default type
    public abstract Boundary makeBoundary(Boundary.Type type);
    public abstract Boundary.Type[] boundaryTypes();
    public boolean requiresSpecialBoundary() {return false;}
    /**
     * Returns an array of dimension D, with each element equal to the given value.
     */
    public abstract int[] makeArrayD(int i);
    /**
     * Returns an array of dimension D, with each element equal to the given value.
     */
    public abstract double[] makeArrayD(double d);

    public double sphereVolume(double r) {
        switch(D) {
            case 1: return 2*r;
            case 2: return Math.PI*r*r;
            default:
            case 3: return (Math.PI*4.0*r*r*r/3.0);
        }
    }
    public double sphereArea(double r)  {
        switch(D) {
            case 1: return 2;
            case 2: return 2*Math.PI*r;
            default:
            case 3: return 4*Math.PI*r*r;
        }
    }
    
    /**
     * Returns the square distance between the two vectors, using the given boundary condition.
     */
    public static double r2(Vector u1, Vector u2, Boundary b) { //square distance between two vectors, subject to boundary b
        if(u1.D() != u2.D()) throw new IllegalArgumentException("Space.r2:  Dimension of vectors not equal to each other");
        switch(u1.D()) {
            case 1: return Space1D.r2((Space1D.Vector)u1, (Space1D.Vector)u2, (Space1D.Boundary)b);
            case 2: return Space2D.r2((Space2D.Vector)u1, (Space2D.Vector)u2, (Space2D.Boundary)b);
            case 3: return Space3D.r2((Space3D.Vector)u1, (Space3D.Vector)u2, (Space3D.Boundary)b);
            default: throw new IllegalArgumentException("Space.r2: Unknown vector dimension");
        }
    }
    /**
     * Returns a Vector from the space of the given dimension.
     */
    public static Vector makeVector(int D) {
        switch(D) {
            case 1:  return new Space1D.Vector();
            case 2:  return new Space2D.Vector();
            case 3:  return new Space3D.Vector();
            default: throw new IllegalArgumentException("Space.makeVector: Requested dimension not implemented");
        }
    }
    /**
     * Returns a Vector initialized to the given set of values in the array.
     * Spatial dimension of the Vector is determined by the length of a.
     */
    public static Vector makeVector(double[] a) {
        switch(a.length) {
            case 1:  return new Space1D.Vector(a);
            case 2:  return new Space2D.Vector(a);
            case 3:  return new Space3D.Vector(a);
            default: throw new IllegalArgumentException("Space.makeVector: Requested dimension not implemented");
        }
    }
    
    /**
     * Returns a Vector initialized to the given set of values in the array (cast to double).
     * Spatial dimension of the Vector is determined by the length of a.
     */
    public static Vector makeVector(int[] k) {
        double[] a = new double[k.length];
        for(int i=0; i<k.length; i++) {a[i] = (double)k[i];}
        return makeVector(a);
    }
            
    
//  Vector contains what is needed to describe a point in the space
    public static abstract class Vector implements java.io.Serializable { 
    /* construct planned to replace the displaceBy methods in coordinte
        private Vector saveVector;
        public void save() {
            if(saveVector == null) saveVector = makeVector(D());
            saveVector.E(this);
        }
        public void restore() {this.E(saveVector);}
        */
        public abstract int length();                         //number of components to vector; equal to the dimension of the space
        public abstract int D();                              //dimension of the space occupied by vector
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
        public abstract void PEa1Tv1(double a, Vector u);     //adds a*u to this vector
        public abstract void Ev1Pv2(Vector u1, Vector u2);    //sets equal to sum of vectors
        public abstract void Ev1Mv2(Vector u1, Vector u2);    //sets equal to difference between vectors
        public abstract Space.Vector P(Space.Vector u);       //adds (+) component-by-component and returns result in another vector
        public abstract Space.Vector M(Space.Vector u);       //subtracts component-by-component and returns result in another vector
        public abstract Space.Vector T(Space.Vector u);       //multiplies component-by-component and returns result in another vector
        public abstract Space.Vector D(Space.Vector u);       //divides component-by-component and returns result in another vector
		public abstract double Mv1Squared(Space.Vector u);    //square of difference between this and given vector
        public abstract void abs();		                      //replaces each component with its absolute value
        public abstract void mod(double a);                   //each component replaced with itself modulo a
        public abstract void mod(Space.Vector u);             //each component replaced with itself modulo component of vector
        public abstract double min();                         //returns minimum of all components
        public abstract double max();                         //returns maximum of all components
        public abstract double squared();                     //square-magnitude of vector (e.g., x^2 + y^2)
        public abstract void normalize();                     //scales the vector to unit length
        public abstract boolean isZero();					  //returns true if all elements of vector are zero
        public abstract double dot(Vector u);                 //dot product of this vector with the vector u
        public abstract Space3D.Vector cross(Space2D.Vector u);       //cross product of this vector with u
        public abstract Space3D.Vector cross(Space3D.Vector u);       //cross product of this vector with u
        public abstract void XE(Space3D.Vector u);            //replaces this vector with its cross product (project result into plane if appropriate)
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
        public Space3D.Vector cross(Space.Vector u) {
            if(u instanceof Space3D.Vector) {return cross((Space3D.Vector)u);}
            else if(u instanceof Space2D.Vector) {return cross((Space2D.Vector)u);}
            else return null;
        }
        public final void setRandomDirection() {setRandomSphere(); normalize();}
        public abstract void randomRotate(double thetaStep);
    }
    
//    public static abstract class Tensor implements java.io.Serializable {
//   declare as interface for RotationTensor
    public interface Tensor {
        public abstract int length();
        public abstract double component(int i, int j);
        public abstract void setComponent(int i, int j, double d);
        public abstract void E(Tensor t);
        public abstract void E(Vector u1, Vector u2);
        public abstract void E(double a);
        public abstract void PE(Tensor t);
        public abstract void PE(int i, int j, double a);
        public abstract void PE(Vector u1, Vector u2);
        public abstract double trace();
        public abstract void TE(double a);
    }
    
    public interface RotationTensor extends Tensor {
        public void reset();
        public void invert();
        public void setAxial(int i, double theta);
        public void setAngles(double[] angles);
    }

//  Coordinate collects all vectors needed to describe point in phase space -- position and (maybe) momentum
    public static abstract class Coordinate implements java.io.Serializable {
        public final Atom atom;  //atom that has this as its coordinate        
        protected double mass, rm;    //mass and its reciprocal
        protected Coordinate(Atom a) {atom = a;}          //constructor
        
        public final Atom atom() {return atom;}
        
        public abstract void transform(Vector r0, Tensor A);
        public abstract Vector position();
        public abstract Vector momentum();
        public abstract double position(int i);
        public abstract double momentum(int i);
        public Vector positionCOM() {return position();} //override in CoordinateGroup
        public abstract double kineticEnergy();
        public abstract void freeFlight(double t);
        public abstract void inflate(double scale);
        public abstract void inflate(Space.Vector scale);
        public abstract void translateBy(Space.Vector u);
        public abstract void translateBy(double d, Space.Vector u);
        public abstract void translateTo(Space.Vector u);
        public void translateCOMTo(Space.Vector u) {translateTo(u);}//override in CoordinateGroup
        public abstract void displaceBy(Space.Vector u);
        public abstract void displaceBy(double d, Space.Vector u);
        public abstract void displaceTo(Space.Vector u);
        public abstract void displaceToRandom(etomica.Phase p);
        public abstract void replace();
        public abstract void accelerateBy(Space.Vector u);
        public abstract void accelerateBy(double d, Space.Vector u);
        public abstract void accelerateTo(Space.Vector u);
        public abstract void displaceWithin(double d);
        public abstract void randomizeMomentum(double temperature);
 //       public abstract void save();
 //       public abstract void restore();

        public final void translateToRandom(etomica.Phase p) {translateTo(p.boundary().randomPosition());}

        public final void scaleMomentum(double scale) {momentum().TE(scale);}

        /**
        * @return mass of the atom, in Daltons
        */
        public double mass() {return mass;}

        /**
        * @return reciprocal of the mass of the atom
        */
        public double rm() {return rm;}
        
        public void setMass(double m) {
            mass = m;
            rm = 1.0/m;
        }
        public double getMass() {return mass;}

        /**
        * Sets the atom to be stationary or movable.
        * The atom does not enforce the condition of being stationary, in that it does not
        * override attempts to move it if it is set as stationary.  Rather, this is a flag
        * that can be set and checked by an integrator when deciding how or whether to 
        * move the atom.
        * 
        * @param b If true, the atom becomes stationary, if false it can be moved.
        */
        public void setStationary(boolean b) {
            stationary = b;
            if(stationary) scaleMomentum(0.0);
        }

        /**
        * Flag to indicate of the atom can or cannot be moved
        * 
        * @see setStationary
        */
        public final boolean isStationary() {return stationary;}

        /**
        * Flag indicating whether atom is stationary or mobile.
        * Default is false (atom is mobile)
        */
        private boolean stationary;
        
        public interface Angular {
            public Orientation orientation();
            public Space3D.Vector angularMomentum(); //angular momentum vector in space-fixed frame
            public Space3D.Vector angularVelocity(); //angular velocity vector in space-fixed frame
            public void angularAccelerateBy(Space3D.Vector v);
            public double kineticEnergy();
            public void freeFlight(double t);
        }
    }//end of Space.Coordinate
    
//    public interface CoordinateGroup {
//    	public Space.Vector positionCOM();
//    	public void translateCOMTo(Space.Vector u);
//    }
    
    public static abstract class Orientation {
        public abstract void E(Orientation o); //copies the given orientation to this
        public abstract Vector[] bodyFrame();//body-frame axes in the space-fixed frame
        public abstract double[] angle();//set of angles describing the orientation
        public abstract void rotateBy(double[] t); //rotate all angles by amounts in t array
        public abstract void rotateBy(int i, double dt); //rotate angle i by given amount
        public abstract void randomRotation(double t); //rotate by random amount in solid angle theta on present position
        public abstract void convertToBodyFrame(Vector v); //changes the components of v from space frame to body frame
        public abstract void convertToSpaceFrame(Vector v);//changes the components of v from body frame to space frame
    }
    
    public static abstract class CoordinatePair implements Cloneable, java.io.Serializable {
        public double r2;
        public CoordinatePair() {}  //null constructor
        public abstract void reset();
        public abstract void resetV(); //computes dv
        public abstract void reset(Space.Coordinate c1, Space.Coordinate c2);
        public abstract double v2();
        public abstract double vDotr();
        public abstract double vDot(Space.Vector u);
        public abstract void push(double impulse);  //impart equal and opposite impulse to momenta
        public abstract void setSeparation(double r2New);  //set square-distance between pair to r2New, by moving them along line joining them, keeping center of mass unchanged
        public abstract double r2();
//        public final double r2() {return r2;}
        public abstract Space.Vector dr();   //separation vector
        public abstract double dr(int i);    //component of separation vector
        public abstract double dv(int i);    //component of velocity-difference vector
        public abstract void setBoundary(Space.Boundary b);
        public abstract Space.Boundary getBoundary();
        /**
        * Clones this coordinatePair without cloning the objects it contains
        * The returned coordinatePair refers to the same pair of coordinates as the original
        * Call it "copy" instead of "clone" because fields are not cloned
        */
        public CoordinatePair copy() {
            try {
                return (CoordinatePair)super.clone();
            } catch(CloneNotSupportedException e) {return null;}
        }
    }
    
 /*   public static abstract class Coordinate3 implements java.io.Serializable {
        public Coordinate3() {}
        public abstract void reset();
        public abstract void reset(Space.Coordinate c1, Space.Coordinate c2, Space.Coordinate c3);
        public abstract double r2_12();
        public abstract double r2_13();
        public abstract double r2_23();
        public abstract double r12Dotr13();
    }
*/
    public static abstract class Boundary implements java.io.Serializable {
        protected final float[][] shift0 = new float[0][0];//cannot be static because several phases may be using at once
        protected float[][] shift;
        private Phase phase;
        protected final PhaseEvent inflateEvent = new PhaseEvent(this, PhaseEvent.BOUNDARY_INFLATE);
        public Boundary() {}
        public Boundary(Phase p) {phase = p;}
        public Phase phase() {return phase;}
        public void setPhase(Phase p) {phase = p;}
        public abstract Space.Boundary.Type type();
        public abstract boolean centralImage(Vector r);//returns true if r is changed by applying central image
        public abstract boolean centralImage(Coordinate c);// {return centralImage(c.position());}
        public abstract void nearestImage(Space.Vector dr);
        public abstract double volume();
 //       public void setVolume(double newVolume) {inflate(Math.pow(newVolume/volume(),1.0/D()));}
 //       public double getVolume() {return volume();}
        /**
         * Returns a copy of the dimensions, as a Vector.  Manipulation of this copy
         * will not cause any change to the boundary's dimensions.
         */
        public abstract Vector dimensions();
        public abstract void setDimensions(Vector v);
        public abstract Vector randomPosition();
        public abstract float[][] getOverflowShifts(Vector r, double distance);
        public abstract void inflate(double s);
        public abstract void inflate(Vector s);
       /** Set of vectors describing the displacements needed to translate the central image
        *  to all of the periodic images.  Returns a two dimensional array of doubles.  The
        *  first index specifies each perioidic image, while the second index indicates the
        *  x and y components of the translation vector.
        *
        *  @param nShells the number of shells of images to be computed
        */
        public abstract double[][] imageOrigins(int nShells);
        
        public static abstract class Type extends Constants.TypedConstant {
            protected Type(String label) {super(label);}
        }
        
      /**
       * Interface for a class that can make a boundary
       */
        public interface Maker extends java.io.Serializable {
            /**
             * Returns the boundary made by the object that implements this interface
             */
            public Boundary makeBoundary(Type t);
            /**
             * Returns an array containing a descriptive string for the boundary.
             * This is used by boundary editors to present boundary choices.
             */
            public Type[] boundaryTypes();
            /**
             * Flag indicating if object requires its special boundary to function.
             */
            public boolean requiresSpecialBoundary();
        }
        /**
         * Marker interface for a periodic boundary.
         */
         public interface Periodic {}
        
        /**
         * Placeholder boundary that performs no actions andreturns null or zero 
         * from every method.
         */
        public static final Boundary NULL = new Boundary() {
            public Space.Vector dimensions() {return null;}
            public Space.Boundary.Type type() {return null;}
            public final void nearestImage(Space.Vector dr) {}
            public final boolean centralImage(Coordinate c) {return false;}
            public final boolean centralImage(Space.Vector r) {return false;}
            public double volume() {return 0.0;}
            public void inflate(double s) {}
            public void inflate(Space.Vector s) {}
            public void setDimensions(Space.Vector v) {}
            public double[][] imageOrigins(int nShells) {return null;}
            public float[][] getOverflowShifts(Space.Vector r, double distance) {return null;}
            public Space.Vector randomPosition() {return null;}
        };//end of NULL
        
    }//end of Space.Boundary
}//end of Space    