package etomica.space;

import java.io.Serializable;

import etomica.Constants;
import etomica.NearestImageTransformer;
import etomica.Phase;
import etomica.PhaseEvent;
import etomica.Space;
import etomica.Constants.TypedConstant;
import etomica.Space.Boundary.Type;



/*
 * History
 * Created on Jan 24, 2005 by kofke
 */
public abstract class Boundary implements NearestImageTransformer, java.io.Serializable {
        protected final float[][] shift0 = new float[0][0];//cannot be static because several phases may be using at once
        protected float[][] shift;
        private Phase phase;
        protected final PhaseEvent inflateEvent = new PhaseEvent(this, PhaseEvent.BOUNDARY_INFLATE);
        public Boundary() {}
        public Boundary(Phase p) {phase = p;}
        public Phase phase() {return phase;}
        public void setPhase(Phase p) {phase = p;}
        public abstract Boundary.Type type();
        public abstract boolean centralImage(Vector r);//returns true if r is changed by applying central image
        public abstract boolean centralImage(Coordinate c);// {return centralImage(c.position());}
        public abstract void nearestImage(Vector dr);
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
            public Vector dimensions() {return null;}
            public Boundary.Type type() {return null;}
            public final void nearestImage(Vector dr) {}
            public final boolean centralImage(Coordinate c) {return false;}
            public final boolean centralImage(Vector r) {return false;}
            public double volume() {return 0.0;}
            public void inflate(double s) {}
            public void inflate(Vector s) {}
            public void setDimensions(Vector v) {}
            public double[][] imageOrigins(int nShells) {return null;}
            public float[][] getOverflowShifts(Vector r, double distance) {return null;}
            public Vector randomPosition() {return null;}
        };//end of NULL
        
    }