package etomica.space2d;

import etomica.Constants;
import etomica.Phase;
import etomica.Constants.TypedConstant;
import etomica.Space2D.Boundary.Type;



/*
 * History
 * Created on Jan 24, 2005 by kofke
 */
public abstract class Boundary extends Boundary {
        public static class Type extends Boundary.Type {
            public Type(String label) {super(label);}
            public Constants.TypedConstant[] choices() {return TYPES;}
        }
        public static final Type NONE = new Type("None");
        public static final Type PERIODIC_SQUARE = new Type("Periodic Square");
        public static final Type HSLIT = new Type("Periodic HSlit");
//        public static final Type SLIDING_BRICK = new Type("Sliding Brick");
        public static final Type[] TYPES = {NONE,PERIODIC_SQUARE, HSLIT/*,SLIDING_BRICK*/};
        public Boundary() {super();}
        public Boundary(Phase p) {super(p);}
        public abstract void nearestImage(Vector dr);
        public abstract boolean centralImage(Vector r);
        public abstract boolean centralImage(Coordinate c);
        public boolean centralImage(Coordinate c) {return centralImage((Coordinate)c);}
    }