package simulate;

import java.awt.Graphics;

public abstract class Space {
    
    public abstract int D();
    
    public abstract Vector makeVector();      //Space.Vector
    public abstract Coordinate makeCoordinate(Occupant o);
    public abstract CoordinatePair makeCoordinatePair(Boundary b);
    public abstract Boundary makeBoundary(int iBoundary);
    
    /**
     * Something that occupies a Space, and therefore has a coordinate
     * Usually an Atom or a Molecule
     */
    interface Occupant {
        public Coordinate coordinate();
        public Phase parentPhase();
        public double mass();
        public double rm();
    }
    
//  Vector contains what is needed to describe a point in the space
    public static abstract class Vector {    //probably want this to be an abstract class
        public abstract double component(int i);              //vector component corresponding to the index i (e.g., i=0, x-component)
        public abstract void setComponent(int i, double d);   //sets ith component of vector to d
        public abstract void E(Vector u);                     //sets each element of the vector equal to the elements of the vector u
        public abstract void E(int i, double a);              //sets component i of this vector equal to a
        public abstract void E(double a);                     //sets all components of the vector equal to the constant a
        public abstract void PE(Vector u);                    //adds (PE is +=) the vector u to this vector
        public abstract void PE(int i, double a);             //adds (+=) a to component i of this vector
        public abstract void ME(Vector u);                    //subtracts (-=)
        public abstract void TE(Vector u);                    //multiplies (*=) component-by-component
        public abstract void DE(Vector u);                    //divide (/=) component-by-component
        public abstract void TE(double a);                    //multipies all components by the constant a
        public abstract void TE(int i, double a);             //multiplies "a" times the component i of this vector
        public abstract void DE(double a);                    //divides all components by a
        public abstract void Ea1Tv1(double a, Vector u);      //sets this vector to a*u
        public abstract void PEa1Tv1(double a, Vector u);     //adds a*u to this vector
        public abstract double squared();                     //square-magnitude of vector (e.g., x^2 + y^2)
        public abstract double dot(Vector u);                 //dot product of this vector with the vector u
        public abstract void setRandom(double d);             //
        public abstract void setRandomSphere();               //random point in unit sphere
        public abstract void setRandomCube();                 //random point in a unit cube
    }

//  Coordinate collects all vectors needed to describe point in phase space -- position and (maybe) momentum
    public static abstract class Coordinate {
        protected final Space.Occupant parent;        //parent is the "Space-occupant" (e.g, Atom or Molecule) that has this as its coordinate        
        Coordinate(Occupant p) {parent = p;}          //constructor
        public final Space.Occupant parent() {return parent;}
        public final Phase parentPhase() {return parent.parentPhase();}
        public abstract Vector makeVector();
        public abstract Vector position();
        public abstract Vector momentum();
        public abstract double position(int i);
        public abstract double momentum(int i);
        public abstract double kineticEnergy(double mass);
        public void scaleMomentum(double scale) {momentum().TE(scale);}
    }
    
    public static abstract class CoordinatePair implements Cloneable {
        public double r2;
//        public Potential potential;
        public CoordinatePair() {}  //null constructor
        public abstract void reset();
        public abstract void reset(Space.Coordinate c1, Space.Coordinate c2);
        public abstract double v2();
        public abstract double vDotr();
        public abstract void push(double impulse);  //impart equal and opposite impulse to momenta
        public abstract void setSeparation(double r2New);  //set square-distance between pair to r2New, by moving them along line joining them, keeping center of mass unchanged
        public final double r2() {return r2;}
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
    }

    public static abstract class Boundary {
        public static final int NONE = 0;
        public static final int PERIODIC = 1;
        public static final int DEFAULT = PERIODIC; //default PBC is periodic
        public abstract void centralImage(Vector r);
        public abstract double volume();
        public abstract Vector dimensions();
        public abstract Vector randomPosition();
        public abstract double[][] getOverflowShifts(Vector r, double distance);
        public abstract void inflate(double s);
        public abstract void draw(Graphics g, int[] origin, double scale);
    /** Set of vectors describing the displacements needed to translate the central image
        *  to all of the periodic images.  Returns a two dimensional array of doubles.  The
        *  first index specifies each perioidic image, while the second index indicates the
        *  x and y components of the translation vector.
        *  Likely to override in subclasses.
        *
        *  @param nShells the number of shells of images to be computed
        */
        public abstract double[][] imageOrigins(int nShells);
    }
    
    public Potential makePotential() {return new PotentialIdealGas();}  //default  
    
    public void draw(Graphics g, int[] origin, double scale) {}
}    