package etomica.space;

import java.io.Serializable;

import etomica.NearestImageTransformer;




/*
 * History
 * Created on Jan 24, 2005 by kofke
 */
public abstract class CoordinatePair implements Cloneable, java.io.Serializable {
        public CoordinatePair() {}  //null constructor
        public abstract void reset();
        public abstract void resetV(); //computes dv
        public abstract void reset(Coordinate c1, Coordinate c2);
        public abstract void trueReset(Coordinate coord1, Coordinate coord2, double falseTime);
        public abstract void trueReset(double falseTime);
        public abstract double v2();
        public abstract double vDotr();
        public abstract double vDot(Vector u);
        public abstract void push(double impulse);  //impart equal and opposite impulse to momenta
        public abstract void truePush(Vector u, double falseTime);
        public abstract void nudge(double rDelta);  //nudge pair apart by rDelta away from each other
        public abstract void setSeparation(double r2New);  //set square-distance between pair to r2New, by moving them along line joining them, keeping center of mass unchanged
        public abstract double r2();
//        public final double r2() {return r2;}
        public abstract Vector dr();   //separation vector
        public abstract double dr(int i);    //component of separation vector
        public abstract double dv(int i);    //component of velocity-difference vector
        public abstract void setNearestImageTransformer(NearestImageTransformer b);
        public abstract NearestImageTransformer getNearestImageTransformer();
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