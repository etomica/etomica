package etomica.space2d;

import etomica.Atom;
import etomica.Phase;
import etomica.space.Boundary;
import etomica.statmech.MaxwellBoltzmann;



/*
 * History
 * Created on Jan 24, 2005 by kofke
 */
public class Coordinate extends Coordinate {
    public Coordinate nextCoordinate, previousCoordinate;
    public final Vector r = new Vector();  //Cartesian coordinates
    public final Vector p = new Vector();  //Momentum vector
    public final Vector rLast = new Vector();  //vector for saving position
    protected final Vector work = new Vector();
    public Coordinate(Atom a) {super(a);}
    
    public void setNextAtom(Atom a) {
        if(a == null) nextCoordinate = null;
        else {
            nextCoordinate = (Coordinate)a.coord;
            ((Coordinate)a.coord).previousCoordinate = this;
        }
    }
    public Atom nextAtom() {return nextCoordinate!=null ? nextCoordinate.atom : null;}
    public Atom previousAtom() {return previousCoordinate!=null ? previousCoordinate.atom : null;}
    public void clearPreviousAtom() {previousCoordinate = null;}
            
    public Vector position() {return r;}
    public Vector truePosition(double falseTime) {
        work.E(r);
        work.PEa1Tv1(falseTime*rm(),p);
        return work;
    }
    public Vector momentum() {return p;}
    public double position(int i) {return r.x(i);}
    public double truePosition(int i, double falseTime) {return r.x(i)+falseTime*rm()*p.x(i);}
    public double momentum(int i) {return p.x(i);}
    public double kineticEnergy() {return 0.5*p.squared()*rm();}
    public void freeFlight(double t) {
        double tM = t*rm(); // t/mass
        r.x += p.x*tM;
        r.y += p.y*tM;
    }
    
    public void inflate(double s) {r.x *= s; r.y *= s;}
    public void inflate(Vector s) {Vector u = (Vector)s; r.x *= u.x; r.y *= u.y;}
    
    public void transform(Vector r0, Tensor A) {
        r.transform((Boundary)atom.node.parentPhase().boundary(), (Vector)r0, (Tensor)A);
    }
    /**
    * Moves the atom by some vector distance
    * 
    * @param u
    */
    public void translateBy(Vector u) {
        r.PE((Vector)u);
    }
    /**
    * Moves the atom by some vector distance
    * 
    * @param u
    */
    public void translateBy(double d, Vector u) {
        r.PEa1Tv1(d,(Vector)u);
    }
    /**
    * Moves the atom by some vector distance
    * 
    * @param u
    */
    public void translateTo(Vector u) {
        r.E((Vector)u);
    }      
    public void displaceBy(Vector u) {rLast.E(r); translateBy((Vector)u);}
    public void displaceBy(double d, Vector u) {rLast.E(r); translateBy(d,(Vector)u);}
    public void displaceTo(Vector u) {rLast.E(r); translateTo((Vector)u);}  
    public void displaceWithin(double d) {work.setRandomCube(); displaceBy(d,work);}
    public void displaceToRandom(etomica.Phase p) {rLast.E(r); translateToRandom(p);}
    public void replace() {
        r.E(rLast);
    }
//    public final void inflate(double s) {r.TE(s);}

    public void accelerateBy(Vector u) {p.PE(u);}
    public void accelerateBy(double d, Vector u) {p.PEa1Tv1(d,u);}
    public void accelerateTo(Vector u) {p.E(u);}
    public void trueAccelerateTo(Vector u, double falseTime) {
        double tm = falseTime*rm();
        r.x -= tm * (((Vector)u).x - p.x);
        r.y -= tm * (((Vector)u).y - p.y);
        p.x = ((Vector)u).x;
        p.y = ((Vector)u).y;
    }

    public void randomizeMomentum(double temperature) {
        if(isStationary()) {p.E(0.0); return;}
        p.setX(0,MaxwellBoltzmann.randomMomentumComponent(temperature,mass()));
        p.setX(1,MaxwellBoltzmann.randomMomentumComponent(temperature,mass()));
    }
}