package etomica.space3d;

import etomica.Atom;
import etomica.Phase;
import etomica.atom.AtomTreeNodeGroup;
import etomica.atom.iterator.AtomIteratorListSimple;
import etomica.space.Coordinate;



/*
 * History
 * Created on Jan 24, 2005 by kofke
 */
public class CoordinateGroup extends Coordinate {

    private final AtomIteratorListSimple childIterator = new AtomIteratorListSimple();
    private Atom firstAtom;
    
    public CoordinateGroup(Atom a) {
        super(a);
        childIterator.setList(((AtomTreeNodeGroup)a.node).childList);
    }

    /**
     * Applies transformation to position of group, keeping all internal atoms
     * at same relative positions.
     */
    public void transform(Vector r0, Tensor A) {
        work.E(position()); //work = r
        work.transform(atom.node.parentPhase().boundary(), r0, A);
        work.ME(r);//now work vector contains translation vector for COM
        translateBy(work);
    }
    /**
     * Returns position of group as the position of the first atom in the group.
     * @see etomica.space.Coordinate#position()
     */
    public Vector position() {
		if(firstAtom == null) firstAtom = ((AtomTreeNodeGroup)atom.node).childList.getFirst(); //DAK 
        if(firstAtom == null) {r.E(0.0); return r;}
		return firstAtom.coord.position();
    }
    /**
     * Returns center-of-mass positions of group.  Somewhat computationally
     * expensive, so not appropriate for repeated use.
     * @see etomica.Space.CoordinateGroup#positionCOM()
     */
	public Vector positionCOM() {
        r.E(0.0); double massSum = 0.0;
        childIterator.reset();
        while(childIterator.hasNext()) {
            Atom a = childIterator.nextAtom();
            r.PEa1Tv1(a.coord.mass(), a.coord.positionCOM()); 
            massSum += a.coord.mass();
        }
        if(massSum == 0) r.E(0.0);//added this 08/13/03 (DAK)
        else r.DE(massSum);
        return r;
	}
    public Vector momentum() {
        p.E(0.0);
        childIterator.reset();
        while(childIterator.hasNext()) {
            p.PE(childIterator.nextAtom().coord.momentum());
        }
        return p;
    }
    public double position(int i) {
        double sum = 0.0; double massSum = 0.0;
        childIterator.reset();
        while(childIterator.hasNext()) {
            Atom a = childIterator.nextAtom();
            sum += a.coord.mass()*a.coord.position(i); 
            massSum += a.coord.mass();
        }
        if(massSum == 0) sum = 0.0;
        else sum /= massSum;
        return sum;
    }
    public double momentum(int i) {
        double sum = 0.0;
        childIterator.reset();
        while(childIterator.hasNext()) {
            Atom a = childIterator.nextAtom();
            sum += a.coord.mass()*a.coord.momentum(i);
        }
        return sum;
    }
    public double kineticEnergy() {
        double sum = 0.0;
        childIterator.reset();
        while(childIterator.hasNext()) {
            sum += childIterator.nextAtom().coord.kineticEnergy();
        }
        return sum;
    }
    public void freeFlight(double t) {
        childIterator.reset();
        while(childIterator.hasNext()) {
            childIterator.nextAtom().coord.freeFlight(t);
        }
    }
    public void inflate(double scale) {
        work.E(position());
        work.TE(scale-1.0);
        translateBy(work);
    }
    public void inflate(Vector scale) {
        scale.PE(-1.0);
        work.E(position());
        work.TE(scale);
        translateBy(work);
        scale.PE(1.0);
    }
    
    protected void translateBy(Vector u) {
        childIterator.reset();
        while(childIterator.hasNext()) {
            childIterator.nextAtom().coord.translateBy(u);
        }
    }
    public void translateBy(double d, Vector u) {
        childIterator.reset();
        while(childIterator.hasNext()) {
            childIterator.nextAtom().coord.translateBy(d, u);
        }
    }
    /**
 	 * Translates group so that first atom is at the given position.
     */
    public void translateTo(Vector u) {
        work.Ea1Tv1(-1,position()); //position() uses work, so need this first
        work.PE(u);
        translateBy(work);
    }
    /**
     * Translates center of mass of group to the given position.
     * @param u new position for the COM
     */
	public void translateCOMTo(Vector u) {
		work.Ea1Tv1(-1,positionCOM()); //position() uses work, so need this first
		work.PE(u);
		translateBy(work);
	}
    public void displaceBy(Vector u) {
        childIterator.reset();
        while(childIterator.hasNext()) {
            childIterator.nextAtom().coord.displaceBy(u);
        }
    }
    public void displaceBy(double d, Vector u) {
        childIterator.reset();
        while(childIterator.hasNext()) {
            childIterator.nextAtom().coord.displaceBy(d, u);
        }
    }
    public void displaceTo(Vector u) {
        work.Ea1Tv1(-1,position()); //position() uses work, so need this first
        work.PE(u);
        displaceBy(work);
    }
    public void displaceToRandom(etomica.Phase phase) {
        displaceTo(phase.boundary().randomPosition());
    }
    public void replace() {
        childIterator.reset();
        while(childIterator.hasNext()) {
            childIterator.nextAtom().coord.replace();
        }
    }
    public void accelerateBy(Vector u) {
        childIterator.reset();
        while(childIterator.hasNext()) {
            childIterator.nextAtom().coord.accelerateBy(u);
        }
    }
    public void accelerateBy(double d, Vector u) {
        childIterator.reset();
        while(childIterator.hasNext()) {
            childIterator.nextAtom().coord.accelerateBy(d, u);
        }
    }
    public void accelerateTo(Vector u) {
        work.Ea1Tv1(-1.0/childIterator.size(),momentum());//probably need this first
        work.PE(u);
        accelerateBy(work);
    }
    public final void displaceWithin(double d) {work.setRandomCube(); displaceBy(d,work);}
        
    public void randomizeMomentum(double temperature) {
        switch(((AtomTreeNodeGroup)atom.node).childAtomCount()) {
            case 0: return;
            case 1: ((AtomTreeNodeGroup)atom.node).firstChildAtom().coord.randomizeMomentum(temperature);//do not zero COM momentum if only one child atom
                    return;
            default://multi-atom group
   //             work.E(0.0); double sum=0.0;
                childIterator.reset();
                while(childIterator.hasNext()) {
                    Atom a = childIterator.nextAtom();
                    a.coord.randomizeMomentum(temperature);
   //                 work.PE(a.coord.momentum());
   //                 sum++;
                }
   //             work.DE(-sum);
   //             childIterator.reset();
   //             while(childIterator.hasNext()) {
    //                childIterator.next().coord.accelerateBy(work);
   //             }
        }//end switch
    }//end randomizeMomentum
    
}