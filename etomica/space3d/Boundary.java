package etomica.space3d;

import etomica.Constants;
import etomica.Phase;
import etomica.Constants.TypedConstant;
import etomica.Space3D.Boundary.Type;



/*
 * History
 * Created on Jan 24, 2005 by kofke
 */
public abstract class Boundary extends Boundary {
        public static class Type extends Boundary.Type {
            private Type(String label) {super(label);}
            public Constants.TypedConstant[] choices() {return TYPES;}
        }
        public static final String[] TAGS = {"None", "Periodic Square"/*, "Sliding Brick"*/};
        public static final Type NONE = new Type("None");
        public static final Type PERIODIC_SQUARE = new Type("Periodic Square");
  //      public static final Type SLIDING_BRICK = new Type("Sliding Brick");
        public static final Type[] TYPES = {NONE,PERIODIC_SQUARE/*,SLIDING_BRICK*/};
        public Boundary() {super();}
        public Boundary(Phase p) {super(p);}
        public abstract void nearestImage(Vector dr);
        public abstract void nearestImage(Vector dr, Vector shift);
        public abstract boolean centralImage(Vector r);
        public boolean centralImage(Coordinate c) {return centralImage((Coordinate)c);}
        public abstract boolean centralImage(Coordinate c);
    }