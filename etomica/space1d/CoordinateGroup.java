package etomica.space1d;

import etomica.Atom;
import etomica.Phase;
import etomica.atom.AtomTreeNodeGroup;
import etomica.space.Boundary;
import etomica.space.Coordinate;



/*
 * History
 * Created on Jan 24, 2005 by kofke
 */
public class CoordinateGroup extends Coordinate {
    public Coordinate firstChild, lastChild;//to be removed
    public Atom firstAtom;
    public CoordinateGroup(Atom a) {super(a);}

	//these may be removed at some point, to use childIterator instead (as in other space classes)
    public final Atom firstAtom() {return (firstChild != null) ? firstChild.atom : null;}
    public final void setFirstAtom(Atom a) {firstChild = (a != null) ? (Coordinate)a.coord : null;}
    public final Atom lastAtom() {return (lastChild != null) ? lastChild.atom : null;}
    public final void setLastAtom(Atom a) {lastChild = (a != null) ? (Coordinate)a.coord : null;}
    
    /**
     * Applies transformation to COM of group, keeping all internal atoms at same relative
     * positions.
     */
    public void transform(Vector r0, Tensor A) {
        work.E(position()); //work = r
        work.transform((Boundary)atom.node.parentPhase().boundary(), (Vector)r0, (Tensor)A);
        work.ME(r);//now work vector contains translation vector for COM
        translateBy(work);
    }
    public Vector position() {
		if(firstAtom == null) firstAtom = ((AtomTreeNodeGroup)atom.node).childList.getFirst(); //DAK 
		if(firstAtom == null) {r.E(0.0); return r;}
		return firstAtom.coord.position();
	}
	public Vector positionCOM() {
        r.E(0.0); double massSum = 0.0;
        for(Coordinate coord=firstChild; coord!=null; coord=coord.nextCoordinate) {
            r.PEa1Tv1(coord.mass(), coord.positionCOM()); massSum += coord.mass();
            if(coord == lastChild) break;
        }
        r.DE(massSum);
        return r;
    }
    public Vector momentum() {
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
        for(Coordinate coord=firstChild; coord!=null; coord=coord.nextCoordinate) {
            coord.freeFlight(t);
            if(coord == lastChild) break;
        }
    }
    public void inflate(double scale) {
        work.E(position());
        work.TE(scale-1.0);
        displaceBy(work);
    }
    public void inflate(Vector scale) {
        scale.PE(-1.0);
        work.E(position());
        work.TE(scale);
        displaceBy(work);
        scale.PE(1.0);
    }
    public void translateBy(Vector u) {translateBy((Vector)u);}
    public void translateBy(Vector u0) {
        for(Coordinate coord=firstChild; coord!=null; coord=coord.nextCoordinate) {
            coord.translateBy(u0);
            if(coord == lastChild) break;
        }
    }
    public void translateBy(double d, Vector u) {
        Vector u0 = (Vector)u;
        for(Coordinate coord=firstChild; coord!=null; coord=coord.nextCoordinate) {
            coord.translateBy(d, u0);
            if(coord == lastChild) break;
        }
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
	/**
	 * Translates group so that first atom is at the given position.
	 */
    public void translateTo(Vector u) {
        work.Ea1Tv1(-1,position()); //position() uses work, so need this first
        work.PE((Vector)u);
        translateBy(work);
    }
    public void displaceBy(Vector u) {
        Vector u0 = (Vector)u;
        for(Coordinate coord=firstChild; coord!=null; coord=coord.nextCoordinate) {
            coord.displaceBy(u0);
            if(coord == lastChild) break;
        }
    }
    public void displaceBy(double d, Vector u) {
        Vector u0 = (Vector)u;
        for(Coordinate coord=firstChild; coord!=null; coord=coord.nextCoordinate) {
            coord.displaceBy(d, u0);
            if(coord == lastChild) break;
        }
    }
    public void displaceTo(Vector u) {
        work.Ea1Tv1(-1,position()); //position() uses work, so need this first
        work.PE((Vector)u);
        displaceBy(work);
    }
    public void displaceToRandom(etomica.Phase phase) {
        displaceTo(phase.boundary().randomPosition());
    }
    public void replace() {
        for(Coordinate coord=firstChild; coord!=null; coord=coord.nextCoordinate) {
            coord.replace();
            if(coord == lastChild) break;
        }
    }
    public void accelerateBy(Vector u) {
        Vector u0 = (Vector)u;
        for(Coordinate coord=firstChild; coord!=null; coord=coord.nextCoordinate) {
            coord.accelerateBy(u0);
            if(coord == lastChild) break;
        }
    }
    public void accelerateBy(double d, Vector u) {
        Vector u0 = (Vector)u;
        for(Coordinate coord=firstChild; coord!=null; coord=coord.nextCoordinate) {
            coord.accelerateBy(d, u0);
            if(coord == lastChild) break;
        }
    }
    public void accelerateTo(Vector u) {
//        throw new RuntimeException("Space1D.CoordinateGroup.accelerateTo not implemented");
        //ugly
        int sum = 0;
        for(Coordinate coord=firstChild; coord!=null; coord=coord.nextCoordinate) {
            sum ++;
            if(coord == lastChild) break;
        }

        work.Ea1Tv1(-1.0/(double)sum,momentum());//probably need this first
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
    //            work.E(0.0); double sum=0.0;
                for(Coordinate coord=firstChild; coord!=null; coord=coord.nextCoordinate) {
                    coord.randomizeMomentum(temperature);
    //                work.PE(coord.momentum());
    //                sum++;
    //                if(coord == lastChild) break;
                }
    //            work.DE(-sum);
    //            for(Coordinate coord=firstChild; coord!=null; coord=coord.nextCoordinate) {
    //                coord.accelerateBy(work);
    //                if(coord == lastChild) break;
    //            }
        }//end switch
    }//end randomizeMomentum
}