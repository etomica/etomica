package etomica.space;

import java.io.Serializable;

import etomica.Atom;
import etomica.Phase;



/*
 * History
 * Created on Jan 24, 2005 by kofke
 */
public abstract class Coordinate implements java.io.Serializable {
        //TODO consider possibility of removing atom reference, perhaps needed only in CoordinateGroup
        public final Atom atom;  //atom that has this as its coordinate        
        protected Coordinate(Atom a) {atom = a;}          //constructor
        public abstract void transform(Vector r0, Tensor A);
        public abstract Vector position();
        public abstract Vector momentum();
        public abstract Vector truePosition(double falseTime);
        public abstract double position(int i);
        public abstract double momentum(int i);
        public abstract double truePosition(int i, double falseTime);
        public Vector positionCOM() {return position();} //override in CoordinateGroup
        public abstract double kineticEnergy();
        public abstract void freeFlight(double t);
        public abstract void inflate(double scale);
        public abstract void inflate(Vector scale);
        public abstract void translateBy(Vector u);
        public abstract void translateBy(double d, Vector u);
        public abstract void translateTo(Vector u);
        public void translateCOMTo(Vector u) {translateTo(u);}//override in CoordinateGroup
        public abstract void displaceBy(Vector u);
        public abstract void displaceBy(double d, Vector u);
        public abstract void displaceTo(Vector u);
        public abstract void displaceToRandom(etomica.Phase p);
        public abstract void replace();
        public abstract void accelerateBy(Vector u);
        public abstract void accelerateBy(double d, Vector u);
        public abstract void accelerateTo(Vector u);
        public abstract void trueAccelerateTo(Vector u, double falseTime);
        public abstract void displaceWithin(double d);
        public abstract void randomizeMomentum(double temperature);
 //       public abstract void save();
 //       public abstract void restore();

        public final void translateToRandom(etomica.Phase p) {translateTo(p.boundary().randomPosition());}

        public final void scaleMomentum(double scale) {momentum().TE(scale);}

        /**
        * @return mass of the atom, in Daltons
        */
        public double mass() {return atom.type.getMass();}//mass;}

        /**
        * @return reciprocal of the mass of the atom
        */
        public double rm() {return atom.type.rm();}//rm;}
        
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
            public Vector angularMomentum(); //angular momentum vector in space-fixed frame
            public Vector angularVelocity(); //angular velocity vector in space-fixed frame
            public void angularAccelerateBy(etomica.space3d.Vector v);
            public double kineticEnergy();
            public void freeFlight(double t);
        }
    }
