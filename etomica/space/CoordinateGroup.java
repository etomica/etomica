package etomica.space;
import etomica.*;

public class CoordinateGroup extends Coordinate {

    private final AtomIterator childIterator;
    
    public CoordinateGroup(Space space, Atom a) {
        super(space,a);
        childIterator = a.node.parentSimulation().iteratorFactory.makeGroupIteratorSimple();
    }

    /**
        * Applies transformation to COM of group, keeping all internal atoms at same relative
        * positions.
        */
    public void transform(Space.Vector r0, Space.Tensor A) {
        work.E(position()); //work = r
        work.transform(atom.node.parentPhase().boundary(), r0, A);
        work.ME(r);//now work vector contains translation vector for COM
        translateBy(work);
    }
    public Space.Vector position() {
        r.E(0.0); double massSum = 0.0;
        childIterator.reset();
        while(childIterator.hasNext()) {
            Atom a = childIterator.next();
            r.PEa1Tv1(a.coord.mass(), a.coord.position()); 
            massSum += a.coord.mass();
        }
        r.DE(massSum);
        return r;
    }
    public Space.Vector momentum() {
        p.E(0.0);
        childIterator.reset();
        while(childIterator.hasNext()) {
            p.PE(childIterator.next().coord.momentum());
        }
        return p;
    }
    public double position(int i) {
        double sum = 0.0; double massSum = 0.0;
        childIterator.reset();
        while(childIterator.hasNext()) {
            Atom a = childIterator.next();
            sum += a.coord.mass()*a.coord.position(i); 
            massSum += a.coord.mass();
        }
        sum /= massSum;
        return sum;
    }
    public double momentum(int i) {
        double sum = 0.0;
        childIterator.reset();
        while(childIterator.hasNext()) {
            Atom a = childIterator.next();
            sum += a.coord.mass()*a.coord.momentum(i);
        }
        return sum;
    }
    public double kineticEnergy() {
        double sum = 0.0;
        childIterator.reset();
        while(childIterator.hasNext()) {
            sum += childIterator.next().coord.kineticEnergy();
        }
        return sum;
    }
    public void freeFlight(double t) {
        double sum = 0.0;
        childIterator.reset();
        while(childIterator.hasNext()) {
            childIterator.next().coord.freeFlight(t);
        }
    }
    public void translateBy(Space.Vector u) {
        childIterator.reset();
        while(childIterator.hasNext()) {
            childIterator.next().coord.translateBy(u);
        }
    }
    public void translateBy(double d, Space.Vector u) {
        childIterator.reset();
        while(childIterator.hasNext()) {
            childIterator.next().coord.translateBy(d, u);
        }
    }
    public void translateTo(Space.Vector u) {
        work.Ea1Tv1(-1,position()); //position() uses work, so need this first
        work.PE(u);
        translateBy(work);
    }
    public void displaceBy(Space.Vector u) {
        childIterator.reset();
        while(childIterator.hasNext()) {
            childIterator.next().coord.displaceBy(u);
        }
    }
    public void displaceBy(double d, Space.Vector u) {
        childIterator.reset();
        while(childIterator.hasNext()) {
            childIterator.next().coord.displaceBy(d, u);
        }
    }
    public void displaceTo(Space.Vector u) {
        work.Ea1Tv1(-1,position()); //position() uses work, so need this first
        work.PE(u);
        displaceBy(work);
    }
    public void displaceToRandom(etomica.Phase p) {
        displaceTo(p.boundary().randomPosition());
    }
    public void replace() {
        childIterator.reset();
        while(childIterator.hasNext()) {
            childIterator.next().coord.replace();
        }
    }
    public void accelerateBy(Space.Vector u) {
        childIterator.reset();
        while(childIterator.hasNext()) {
            childIterator.next().coord.accelerateBy(u);
        }
    }
    public void accelerateBy(double d, Space.Vector u) {
        childIterator.reset();
        while(childIterator.hasNext()) {
            childIterator.next().coord.accelerateBy(d, u);
        }
    }
    public final void displaceWithin(double d) {work.setRandomCube(); displaceBy(d,work);}
        
    public void randomizeMomentum(double temperature) {
        switch(((AtomGroup)atom).node.childAtomCount()) {
            case 0: return;
            case 1: ((AtomTreeNodeGroup)atom.node).firstChildAtom().coord.randomizeMomentum(temperature);//do not zero COM momentum if only one child atom
                    return;
            default://multi-atom group
                work.E(0.0); double sum=0.0;
                childIterator.reset();
                while(childIterator.hasNext()) {
                    Atom a = childIterator.next();
                    a.coord.randomizeMomentum(temperature);
                    work.PE(a.coord.momentum());
                    sum++;
                }
                work.DE(-sum);
                childIterator.reset();
                while(childIterator.hasNext()) {
                    childIterator.next().coord.accelerateBy(work);
                }
        }//end switch
    }//end randomizeMomentum
    
}//end of CoordinateGroup
