package etomica;

import java.awt.Graphics;

public abstract class Space implements Space.Boundary.Maker, java.io.Serializable {
    
    public static String VERSION = "01.07.09";

    public Space() {}
    public abstract int D();
    
    public abstract Vector makeVector();      //Space.Vector
    public abstract Orientation makeOrientation();
    public abstract Tensor makeTensor();
    public abstract Coordinate makeCoordinate(Atom a);
    public abstract CoordinatePair makeCoordinatePair(Phase p);
    public abstract Boundary makeBoundary();  //makes boundary of default type
    public abstract Boundary makeBoundary(Boundary.Type type);
    public abstract Boundary.Type[] boundaryTypes();
    public boolean requiresSpecialBoundary() {return false;}
    public abstract double sphereVolume(double r);  //volume of a sphere of radius r
    public abstract double sphereArea(double r);    //surface area of sphere of radius r (used for differential shell volume)
    
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
//            case 3:  return new Space3D.Vector();
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
//            case 3:  return new Space3D.Vector(a);
            default: throw new IllegalArgumentException("Space.makeVector: Requested dimension not implemented");
        }
    }
            
    
//  Vector contains what is needed to describe a point in the space
    public static abstract class Vector implements java.io.Serializable { 
        public abstract int length();                         //number of components to vector; equal to the dimension of the space
        public abstract int D();                              //dimension of the space occupied by vector
        public abstract double[] toArray();                   //converts components to array of double
        public abstract void sphericalCoordinates(double[] result); //computes the spherical coordinate representation of the vector; return in the given array to avoid construction of new array with each call
        public double[] toSphericalCoordinateArray() {        //computes spherical coordinates and returns them in a new array
            double[] array = new double[D()];
            sphericalCoordinates(array);
            return array;
        }
        public abstract double component(int i);              //vector component corresponding to the index i (e.g., i=0, x-component)
        public abstract void setComponent(int i, double d);   //sets ith component of vector to d
        public abstract void E(Vector u);                     //sets each element of the vector equal to the elements of the vector u
//        public abstract void E(int i, double a);              //sets component i of this vector equal to a
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
        public abstract Space.Vector P(Space.Vector u);       //adds (+) component-by-component and returns result in another vector
        public abstract Space.Vector M(Space.Vector u);       //subtracts component-by-component and returns result in another vector
        public abstract Space.Vector T(Space.Vector u);       //multiplies component-by-component and returns result in another vector
        public abstract Space.Vector D(Space.Vector u);       //divides component-by-component and returns result in another vector
        public abstract Space.Vector abs();                   //takes absolute value of each component and returns result in another vector
        public abstract double min();                         //returns minimum of all components
        public abstract double max();                         //returns maximum of all components
        public abstract double squared();                     //square-magnitude of vector (e.g., x^2 + y^2)
        public abstract void normalize();                     //scales the vector to unit length
        public abstract double dot(Vector u);                 //dot product of this vector with the vector u
        public abstract Space3D.Vector cross(Space2D.Vector u);       //cross product of this vector with u
        public abstract Space3D.Vector cross(Space3D.Vector u);       //cross product of this vector with u
        public abstract void XE(Space3D.Vector u);            //replaces this vector with its cross product (project result into plane if appropriate)
        public abstract void setRandom(double d);             //
        public abstract void setRandomSphere();               //random point in unit sphere
        public abstract void setRandomCube();                 //random point in a unit cube
        public abstract void randomStep(double d);            //random step of selected uniformly in cube of edge 2d (adds up to +/- d to present value)
        public final void PEa1Tv1(double[] a, Vector[] u) {   //adds several terms of form a*u to this vector
            for(int i=a.length-1; i>=0; i--) {PEa1Tv1(a[i],u[i]);}
        }
        public Space3D.Vector cross(Space.Vector u) {
            if(u instanceof Space3D.Vector) {return cross((Space3D.Vector)u);}
            else if(u instanceof Space2D.Vector) {return cross((Space2D.Vector)u);}
            else return null;
        }
    }
    
    public static abstract class Tensor implements java.io.Serializable {
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

//  Coordinate collects all vectors needed to describe point in phase space -- position and (maybe) momentum
    public static abstract class Coordinate implements java.io.Serializable {
        protected final Atom atom;  //atom that has this as its coordinate        
        private double mass, rm;    //mass and its reciprocal
        Coordinate(Atom a) {atom = a;}          //constructor
        
        public final Atom atom() {return atom;}

        public abstract Vector position();
        public abstract Vector momentum();
        public abstract double position(int i);
        public abstract double momentum(int i);
        public abstract double kineticEnergy();
        public abstract void freeFlight(double t);
        public abstract void translateBy(Space.Vector u);
        public abstract void translateBy(double d, Space.Vector u);
        public abstract void translateTo(Space.Vector u);   
        public abstract void displaceBy(Space.Vector u);
        public abstract void displaceBy(double d, Space.Vector u);
        public abstract void displaceTo(Space.Vector u);
        public abstract void displaceToRandom(etomica.Phase p);
        public abstract void replace();
        public abstract void accelerateBy(Space.Vector u);
        public abstract void accelerateBy(double d, Space.Vector u);
        public abstract void displaceWithin(double d);
        public abstract void randomizeMomentum(double temperature);

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
            if(!stationary) scaleMomentum(0.0);
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
    
    public static abstract class CoordinateGroup extends Coordinate {
        public CoordinateGroup(AtomGroup a) {super(a);}
        
        public abstract void addCoordinate(Coordinate coord);
    }
    
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
    
    public static abstract class CoordinatePair implements Cloneable, java.io.Serializable, java.util.Observer {
        public double r2;
        public CoordinatePair() {}  //null constructor
        public abstract void reset();
        public abstract void reset(Space.Coordinate c1, Space.Coordinate c2);
        public abstract double v2();
        public abstract double vDotr();
        public abstract double vDot(Space.Vector u);
        public abstract void push(double impulse);  //impart equal and opposite impulse to momenta
        public abstract void setSeparation(double r2New);  //set square-distance between pair to r2New, by moving them along line joining them, keeping center of mass unchanged
        public final double r2() {return r2;}
        public abstract Space.Vector dr();   //separation vector
        public abstract double dr(int i);    //component of separation vector
        public abstract double dv(int i);    //component of velocity-difference vector
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
        
        public abstract void update(java.util.Observable obs, Object arg);
            
    }

    public static abstract class Boundary implements java.io.Serializable {
        private Phase phase;
        public Boundary() {}
        public Boundary(Phase p) {phase = p;}
        public Phase phase() {return phase;}
        public void setPhase(Phase p) {phase = p;}
        public abstract Space.Boundary.Type type();
        public abstract void centralImage(Vector r);
        public void centralImage(Coordinate c) {centralImage(c.position());}
        public abstract void nearestImage(Space.Vector dr);
        public abstract double volume();
 //       public void setVolume(double newVolume) {inflate(Math.pow(newVolume/volume(),1.0/D()));}
 //       public double getVolume() {return volume();}
        public abstract Vector dimensions();
        public abstract Vector randomPosition();
        public abstract double[][] getOverflowShifts(Vector r, double distance);
        public abstract void inflate(double s);
        public abstract void draw(Graphics g, int[] origin, double scale);
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
        
    }//end of Space.Boundary
        
    public void draw(Graphics g, int[] origin, double scale) {}
}//end of Space    