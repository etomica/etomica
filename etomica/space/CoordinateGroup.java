package etomica.space;
import etomica.*;

public class CoordinateGroup extends Coordinate implements Space.CoordinateGroup {
    public Coordinate firstChild, lastChild;
    public CoordinateGroup(Space space, AtomGroup a) {super(space,a);}

    public final Atom firstAtom() {return (firstChild != null) ? firstChild.atom : null;}
    public final void setFirstAtom(Atom a) {firstChild = (a != null) ? (Coordinate)a.coord : null;}
    public final Atom lastAtom() {return (lastChild != null) ? lastChild.atom : null;}
    public final void setLastAtom(Atom a) {lastChild = (a != null) ? (Coordinate)a.coord : null;}
        
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
        for(Coordinate coord=firstChild; coord!=null; coord=coord.nextCoordinate) {
            r.PEa1Tv1(coord.mass(), coord.position()); massSum += coord.mass();
            if(coord == lastChild) break;
        }
        r.DE(massSum);
        return r;
    }
    public Space.Vector momentum() {
        p.E(0.0);
        for(Coordinate coord=firstChild; coord!=null; coord=coord.nextCoordinate) {
            p.PE(coord.momentum());
            if(coord == lastChild) break;
        }
        return p;
    }
    public double position(int i) {
        double sum = 0.0; double massSum = 0.0;
        for(Coordinate coord=firstChild; coord!=null; coord=coord.nextCoordinate) {
            sum += coord.mass()*coord.position(i); massSum += coord.mass();
            if(coord == lastChild) break;
        }
        sum /= massSum;
        return sum;
    }
    public double momentum(int i) {
        double sum = 0.0;
        for(Coordinate coord=firstChild; coord!=null; coord=coord.nextCoordinate) {
            sum += coord.mass()*coord.momentum(i);
            if(coord == lastChild) break;
        }
        return sum;
    }
    public double kineticEnergy() {
        double sum = 0.0;
        for(Coordinate coord=firstChild; coord!=null; coord=coord.nextCoordinate) {
            sum += coord.kineticEnergy();
            if(coord == lastChild) break;
        }
        return sum;
    }
    public void freeFlight(double t) {
        double sum = 0.0;
        for(Coordinate coord=firstChild; coord!=null; coord=coord.nextCoordinate) {
            coord.freeFlight(t);
            if(coord == lastChild) break;
        }
    }
    public void translateBy(Space.Vector u) {
        for(Coordinate coord=firstChild; coord!=null; coord=coord.nextCoordinate) {
            coord.translateBy(u);
            if(coord == lastChild) break;
        }
    }
    public void translateBy(double d, Space.Vector u) {
        for(Coordinate coord=firstChild; coord!=null; coord=coord.nextCoordinate) {
            coord.translateBy(d, u);
            if(coord == lastChild) break;
        }
    }
    public void translateTo(Space.Vector u) {
        work.Ea1Tv1(-1,position()); //position() uses work, so need this first
        work.PE(u);
        translateBy(work);
    }
    public void displaceBy(Space.Vector u) {
        for(Coordinate coord=firstChild; coord!=null; coord=coord.nextCoordinate) {
            coord.displaceBy(u);
            if(coord == lastChild) break;
        }
    }
    public void displaceBy(double d, Space.Vector u) {
        for(Coordinate coord=firstChild; coord!=null; coord=coord.nextCoordinate) {
            coord.displaceBy(d, u);
            if(coord == lastChild) break;
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
        for(Coordinate coord=firstChild; coord!=null; coord=coord.nextCoordinate) {
            coord.replace();
            if(coord == lastChild) break;
        }
    }
    public void accelerateBy(Space.Vector u) {
        for(Coordinate coord=firstChild; coord!=null; coord=coord.nextCoordinate) {
            coord.accelerateBy(u);
            if(coord == lastChild) break;
        }
    }
    public void accelerateBy(double d, Space.Vector u) {
        for(Coordinate coord=firstChild; coord!=null; coord=coord.nextCoordinate) {
            coord.accelerateBy(d, u);
            if(coord == lastChild) break;
        }
    }
    public final void displaceWithin(double d) {work.setRandomCube(); displaceBy(d,work);}
        
    public void randomizeMomentum(double temperature) {
        switch(((AtomGroup)atom).node.childAtomCount()) {
            case 0: return;
            case 1: firstChild.randomizeMomentum(temperature);//do not zero COM momentum if only one child atom
                    return;
            default://multi-atom group
                work.E(0.0); double sum=0.0;
                for(Coordinate coord=firstChild; coord!=null; coord=coord.nextCoordinate) {
                    coord.randomizeMomentum(temperature);
                    work.PE(coord.momentum());
                    sum++;
                    if(coord == lastChild) break;
                }
                work.DE(-sum);
                for(Coordinate coord=firstChild; coord!=null; coord=coord.nextCoordinate) {
                    coord.accelerateBy(work);
                    if(coord == lastChild) break;
                }
        }//end switch
    }//end randomizeMomentum
}//end of CoordinateGroup
