package simulate.electrostatics;

import simulate.*;

public class LatticeSum extends MeterFunction {
    
    AtomPair.Iterator.All iteratorAll;
    Atom.Iterator atomIterator;
    private double z2Sum;
    private double m2Max;
      
    public void setNPoints(int n) {
        setX(-0.5,(double)n-0.5,n);
    }
    
    public void setPhase(Phase p) {
        super.setPhase(p);
        iteratorAll = p.iterator.makeAtomPairIteratorAll();
        atomIterator = p.iterator.makeAtomIteratorUp();
    }
    
    public double[] currentValue() {
        cubic();
        return y;
    }
    
    /**
     * Performs a sum of the energy over periodic images such that each additional
     * shell yields a cubic shaped region
     * Assumes a 2-D cubic central shell (same length and width)
     */
    private void cubic() { 
        //Compute sum of square charges to add in self-image interactions
        atomIterator.reset();
        z2Sum = 0.0;
        while(atomIterator.hasNext()) {
            Atom a = atomIterator.next();
            double z = ((Monopole)a.type.electroType()).z();
            z2Sum += z*z;
        }
        //Do lattice sum
        Space2D.Vector M = (Space2D.Vector)phase.parentSimulation().space.makeVector();
        double L = phase.boundary().dimensions().component(0);
        m2Max = 1.001*((nPoints-1)*(nPoints-1)*L*L);
        double peSum = energy(M);
        y[0] = peSum;
        for(int m=1; m<nPoints; m++) {
            M.x = m*L;
            M.y = m*L;
            peSum += energy(M);
            M.x = -m*L;
            peSum += energy(M);
            M.y = -m*L;
            peSum += energy(M);
            M.x = +m*L;
            peSum += energy(M);
            for(int i=-m+1; i<m; i++) {
                M.x = +m*L;
                M.y = +i*L;
                peSum += energy(M);
                M.x = -m*L;
                peSum += energy(M);
                M.x = +i*L;
                M.y = +m*L;
                peSum += energy(M);
                M.y = -m*L;
                peSum += energy(M);
            }
            y[m] = peSum;
        }
     }
     
     //Computes all interactions between atoms in central image and those in another
     //image displaced from the origin by the vector M
     private double energy(Space2D.Vector M) {
        double m2 = M.squared();
        if(m2 > m2Max) return 0.0;
        double pe = 0.0;
        iteratorAll.reset();
        while(iteratorAll.hasNext()) {
            AtomPair pair = iteratorAll.next();
            if(m2 > 0) ((Space2D.CoordinatePair)pair.cPair).reset(M);
            double energy = phase.parentSimulation().getPotential(pair).energy(pair);
            if(energy == Double.MAX_VALUE) return Double.MAX_VALUE;
            pe += energy;
        }
        double self = (m2 > 0.0) ? z2Sum/Math.sqrt(m2) : 0.0;  //self-image interaction
        return 2*pe + self;
     }  
}