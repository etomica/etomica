package etomica.space;
import etomica.*;

public class Coordinate extends Space.Coordinate {
    public final Space.Vector r, p, rLast, work;
//    public final etomica.space.continuum.Vector r, p, rLast, work;
    public Coordinate(Space space, Atom atom) {
        super(atom);
        r = space.makeVector();
        p = space.makeVector();
        rLast = space.makeVector();
        work = space.makeVector();
/*        r = (etomica.space.continuum.Vector)space.makeVector();
        p = (etomica.space.continuum.Vector)space.makeVector();
        rLast = (etomica.space.continuum.Vector)space.makeVector();
        work = (etomica.space.continuum.Vector)space.makeVector();
 */   }
        
    public void transform(Space.Vector r0, Space.Tensor A) {
        r.transform(atom.node.parentPhase().boundary(), r0, A);
    }
    
    public Space.Vector position() {return r;}
    public Space.Vector momentum() {return p;}
    public double position(int i) {return r.component(i);}
    public double momentum(int i) {return p.component(i);}
    public double kineticEnergy() {return 0.5*p.squared()*rm();}
    public void freeFlight(double t) {r.PEa1Tv1(t*rm(),p);}

    /**
    * Moves the atom by some vector distance
    * 
    * @param u
    */
    public void translateBy(Space.Vector u) {r.PE(u);}
    /**
    * Moves the atom by some vector distance
    * 
    * @param u
    */
    public void translateBy(double d, Space.Vector u) {r.PEa1Tv1(d,u);}
    /**
    * Moves the atom by some vector distance
    * 
    * @param u
    */
    public void translateTo(Space.Vector u) {r.E(u);}      
    public void displaceBy(Space.Vector u) {rLast.E(r); translateBy(u);}
    public void displaceBy(double d, Space.Vector u) {rLast.E(r); translateBy(d,u);}
    public void displaceTo(Space.Vector u) {rLast.E(r); translateTo(u);}  
    public void displaceWithin(double d) {work.setRandomCube(); displaceBy(d,work);}
    public void displaceToRandom(etomica.Phase p) {rLast.E(r); translateToRandom(p);}
    public void replace() {r.E(rLast);}
//    public final void inflate(double s) {r.TE(s);}

    public void accelerateBy(Space.Vector u) {p.PE(u);}
    public void accelerateBy(double d, Space.Vector u) {p.PEa1Tv1(d,u);}

    public void randomizeMomentum(double temperature) {  //not very sophisticated; random only in direction, not magnitude
        if(isStationary()) {p.E(0.0); return;}
        double magnitude = Math.sqrt(mass()*temperature*(double)p.D());  //need to divide by sqrt(m) to get velocity
        momentum().setRandomDirection();
        momentum().TE(magnitude);
        //for debugging
    //      momentum().E(position());
    //      momentum().TE(magnitude/30.);
    }
}
