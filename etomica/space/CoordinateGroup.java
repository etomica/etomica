package etomica.space;
import etomica.*;

public class CoordinateGroup extends Coordinate {

    public CoordinateGroup(Space space, Atom a) {super(space,a);}

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
        Atom lastChild = atom.node.lastChildAtom();
        for(Atom a=atom.node.firstChildAtom(); a!=null; a=a.seq.nextAtom()) {
            r.PEa1Tv1(a.coord.mass(), a.coord.position()); massSum += a.coord.mass();
            if(a == lastChild) break;
        }
        r.DE(massSum);
        return r;
    }
    public Space.Vector momentum() {
        p.E(0.0);
        Atom lastChild = atom.node.lastChildAtom();
        for(Atom a=atom.node.firstChildAtom(); a!=null; a=a.nextAtom()) {
            p.PE(a.coord.momentum());
            if(a == lastChild) break;
        }
        return p;
    }
    public double position(int i) {
        double sum = 0.0; double massSum = 0.0;
        Atom lastChild = atom.node.lastChildAtom();
        for(Atom a=atom.node.firstChildAtom(); a!=null; a=a.nextAtom()) {
            sum += a.coord.mass()*a.coord.position(i); massSum += a.coord.mass();
            if(a == lastChild) break;
        }
        sum /= massSum;
        return sum;
    }
    public double momentum(int i) {
        double sum = 0.0;
        Atom lastChild = atom.node.lastChildAtom();
        for(Atom a=atom.node.firstChildAtom(); a!=null; a=a.nextAtom()) {
            sum += a.coord.mass()*a.coord.momentum(i);
            if(a == lastChild) break;
        }
        return sum;
    }
    public double kineticEnergy() {
        double sum = 0.0;
        Atom lastChild = atom.node.lastChildAtom();
        for(Atom a=atom.node.firstChildAtom(); a!=null; a=a.nextAtom()) {
            sum += a.coord.kineticEnergy();
            if(a == lastChild) break;
        }
        return sum;
    }
    public void freeFlight(double t) {
        double sum = 0.0;
        Atom lastChild = atom.node.lastChildAtom();
        for(Atom a=atom.node.firstChildAtom(); a!=null; a=a.nextAtom()) {
            a.coord.freeFlight(t);
            if(a == lastChild) break;
        }
    }
    public void translateBy(Space.Vector u) {
        Atom lastChild = atom.node.lastChildAtom();
        for(Atom a=atom.node.firstChildAtom(); a!=null; a=a.nextAtom()) {
            a.coord.translateBy(u);
            if(a == lastChild) break;
        }
    }
    public void translateBy(double d, Space.Vector u) {
        Atom lastChild = atom.node.lastChildAtom();
        for(Atom a=atom.node.firstChildAtom(); a!=null; a=a.nextAtom()) {
            a.coord.translateBy(d, u);
            if(a == lastChild) break;
        }
    }
    public void translateTo(Space.Vector u) {
        work.Ea1Tv1(-1,position()); //position() uses work, so need this first
        work.PE(u);
        translateBy(work);
    }
    public void displaceBy(Space.Vector u) {
        Atom lastChild = atom.node.lastChildAtom();
        for(Atom a=atom.node.firstChildAtom(); a!=null; a=a.nextAtom()) {
            a.coord.displaceBy(u);
            if(a == lastChild) break;
        }
    }
    public void displaceBy(double d, Space.Vector u) {
        Atom lastChild = atom.node.lastChildAtom();
        for(Atom a=atom.node.firstChildAtom(); a!=null; a=a.nextAtom()) {
            a.coord.displaceBy(d, u);
            if(a == lastChild) break;
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
        Atom lastChild = atom.node.lastChildAtom();
        for(Atom a=atom.node.firstChildAtom(); a!=null; a=a.nextAtom()) {
            a.coord.replace();
            if(a == lastChild) break;
        }
    }
    public void accelerateBy(Space.Vector u) {
        Atom lastChild = atom.node.lastChildAtom();
        for(Atom a=atom.node.firstChildAtom(); a!=null; a=a.nextAtom()) {
            a.coord.accelerateBy(u);
            if(a == lastChild) break;
        }
    }
    public void accelerateBy(double d, Space.Vector u) {
        Atom lastChild = atom.node.lastChildAtom();
        for(Atom a=atom.node.firstChildAtom(); a!=null; a=a.nextAtom()) {
            a.coord.accelerateBy(d, u);
            if(a == lastChild) break;
        }
    }
    public final void displaceWithin(double d) {work.setRandomCube(); displaceBy(d,work);}
        
    public void randomizeMomentum(double temperature) {
        switch(((AtomGroup)atom).node.childAtomCount()) {
            case 0: return;
            case 1: atom.node.firstChildAtom().coord.randomizeMomentum(temperature);//do not zero COM momentum if only one child atom
                    return;
            default://multi-atom group
                work.E(0.0); double sum=0.0;
                Atom lastChild = atom.node.lastChildAtom();
                for(Atom a=atom.node.firstChildAtom(); a!=null; a=a.nextAtom()) {
                    a.coord.randomizeMomentum(temperature);
                    work.PE(a.coord.momentum());
                    sum++;
                    if(a == lastChild) break;
                }
                work.DE(-sum);
                for(Atom a=atom.node.firstChildAtom(); a!=null; a=a.nextAtom()) {
                    a.coord.accelerateBy(work);
                    if(a == lastChild) break;
                }
        }//end switch
    }//end randomizeMomentum
    
}//end of CoordinateGroup
