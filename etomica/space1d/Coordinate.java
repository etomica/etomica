package etomica.space1d;

import etomica.Atom;
import etomica.Phase;
import etomica.statmech.MaxwellBoltzmann;



/*
 * History
 * Created on Jan 24, 2005 by kofke
 */
public class Coordinate extends etomica.space.Coordinate {
    public Coordinate nextCoordinate, previousCoordinate;
    public final Vector r = new Vector();  //Cartesian coordinates
    public final Vector p = new Vector();  //Momentum vector
    public final Vector rLast = new Vector();  //vector for saving position
    protected final Vector work = new Vector();
    public Coordinate(Atom atom) {super(atom);}
    
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

    public void transform(etomica.space.Vector r0, etomica.space.Tensor A) {
        r.transform((Boundary)atom.node.parentPhase().boundary(), (Vector)r0, (Tensor)A);
    }
    public etomica.space.Vector position() {return r;}
    public etomica.space.Vector truePosition(double falseTime) {
        work.E(r);
        work.PEa1Tv1(falseTime*rm(),p);
        return work;
    }
    public etomica.space.Vector momentum() {return p;}
    public double position(int i) {return r.x(i);}
    public double truePosition(int i, double falseTime) {return r.x(i)+falseTime*rm()*p.x(i);}
    public double momentum(int i) {return p.x(i);}
    public double kineticEnergy() {return 0.5*p.squared()*rm();}
    public void freeFlight(double t) {r.x += p.x*t*rm();}
    public void inflate(double s) {r.x *= s;}
    public void inflate(etomica.space.Vector s) {r.x *= ((Vector)s).x;}

    /**
    * Moves the atom by some vector distance
    * 
    * @param u
    */
    public void translateBy(etomica.space.Vector u) {
        r.PE((Vector)u);
    }
    /**
    * Moves the atom by some vector distance
    * 
    * @param u
    */
    public void translateBy(double d, etomica.space.Vector u) {
        r.PEa1Tv1(d,u);
    }
    /**
    * Moves the atom by some vector distance
    * 
    * @param u
    */
    public void translateTo(etomica.space.Vector u) {
        r.E((Vector)u);
    }      
    public void replace() {
        r.E(rLast);
    }
    public void displaceBy(etomica.space.Vector u) {rLast.E(r); translateBy(u);}
    public void displaceBy(double d, etomica.space.Vector u) {rLast.E(r); translateBy(d,u);}
    public void displaceTo(etomica.space.Vector u) {rLast.E(r); translateTo(u);}  
    public void displaceWithin(double d) {work.setRandomCube(); displaceBy(d,work);}
    public void displaceToRandom(etomica.Phase phase) {rLast.E(r); translateToRandom(phase);}
//    public final void inflate(double s) {r.TE(s);}

    public void accelerateBy(etomica.space.Vector u) {p.PE(u);}
    public void accelerateBy(double d, etomica.space.Vector u) {p.PEa1Tv1(d,u);}
    public void accelerateTo(etomica.space.Vector u) {p.E(u);}
    public void trueAccelerateTo(etomica.space.Vector u, double falseTime) {
        r.x -= falseTime*rm() * (((Vector)u).x - p.x);
        p.x = ((Vector)u).x;
    }
    
    public void randomizeMomentum(double temperature) {
        if(isStationary()) {p.E(0.0); return;}
        p.setX(0,MaxwellBoltzmann.randomMomentumComponent(temperature,mass()));
    }
}
