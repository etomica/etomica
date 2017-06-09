package etomica.potential;

import etomica.api.IPotentialMolecular;
import etomica.box.Box;
import etomica.molecule.IMoleculeList;
import etomica.molecule.IMoleculePositionDefinition;
import etomica.space.Boundary;
import etomica.space.Space;
import etomica.space.Vector;


/**
 * Wraps a soft-spherical potential to apply a truncation to it.  The energy
 * is switched from fully-on to 0 over a short range (i.e. from 0.95*rC to rC).
 */
public class P2MoleculeTruncated extends PotentialMolecular {
    
    public P2MoleculeTruncated(IPotentialMolecular potential, double truncationRadius, Space _space, IMoleculePositionDefinition positionDefinition) {
        super(2, _space);
        this.potential = potential;
        this.positionDefinition = positionDefinition;
        setTruncationRadius(truncationRadius);
        dr = space.makeVector();
    }
    
    /**
     * Returns the wrapped potential.
     */
    public IPotentialMolecular getWrappedPotential() {
        return potential;
    }

    /**
     * Mutator method for the radial cutoff distance.
     */
    public void setTruncationRadius(double rCut) {
        rCutoff = rCut;
        r2Cutoff = rCut*rCut;
        System.out.println("in P2MoleculeTruncated rCut:"+rCutoff);
    }

    /**
     * Accessor method for the radial cutoff distance.
     */
    public double getTruncationRadius() {return rCutoff;}
    
    /**
     * Returns the truncation radius.
     */
    public double getRange() {
        return rCutoff;
    }
    
    public double energy(IMoleculeList atoms) {
    	
        dr.E(positionDefinition.position(atoms.getMolecule(1)));
        dr.ME(positionDefinition.position(atoms.getMolecule(0)));
 //       System.out.println("in P2MoleculeTruncated,before nearest image, distance between two molecules: "+Math.sqrt(dr.squared()));
        boundary.nearestImage(dr);
        
        if (dr.squared() > r2Cutoff) {
        //	System.out.println("out of truncation");
        //	System.out.println("in P2MoleculeTruncated, distance between two molecules: "+Math.sqrt(dr.squared()));
            return 0;
        }
     //   System.out.println("in P2MoleculeTruncated, cutoff is:"+Math.sqrt(r2Cutoff));
      //  System.out.println("in P2MoleculeTruncated, distance between two molecules: "+Math.sqrt(dr.squared()));
        double u = potential.energy(atoms);
        //System.out.println("in P2MoleculeTruncated,potential:"+u);
        return u;
    }
    
    /**
     * Returns the dimension (length) of the radial cutoff distance.
     */
    public etomica.units.Dimension getTruncationRadiusDimension() {return etomica.units.Length.DIMENSION;}
    
    public void setBox(Box newBox) {
        potential.setBox(newBox);
        boundary = newBox.getBoundary();
        if (cutoffRatio > 0){
        	Vector vectorBox = boundary.getBoxSize();
        	double minBoxSize = vectorBox.getX(0);
        	for (int i = 1;i<vectorBox.getD();i++){
        		if (vectorBox.getX(i) < minBoxSize){
        			minBoxSize = vectorBox.getX(i);
        		}
        	}
        	setTruncationRadius(minBoxSize * cutoffRatio);
        }
    }
    
    public void setCutoffRatio(double newCutoffRatio) {
    	cutoffRatio = newCutoffRatio;
    }
    
    private static final long serialVersionUID = 1L;
    protected double rCutoff, r2Cutoff;
    protected final IPotentialMolecular potential;
    protected final Vector dr;
    protected Boundary boundary;
    protected double cutoffRatio;
    protected final IMoleculePositionDefinition positionDefinition;
}
