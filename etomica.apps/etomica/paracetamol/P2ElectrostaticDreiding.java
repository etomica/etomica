package etomica.paracetamol;

import etomica.atom.AtomSet;
import etomica.atom.IAtomPositioned;
import etomica.simulation.ISimulation;
import etomica.space.IVector;
import etomica.space.Space;

/**
 * Empirical Isotropic atom-atom repulsion-dispersion potential
 * Given formula:
 * 				
 * 				U(r) = A*exp(-r/B) - C /r^6	
 * 
 * 
 *	    A is in eV
 *      B is in angstrom
 *      C is in eV.angstrom^6
 *
 * @author Tai Tan
 */

public class P2ElectrostaticDreiding extends etomica.potential.P2Exp6 {
	
	public P2ElectrostaticDreiding(ISimulation sim) {
        this(sim.getSpace(), 1.0, 1.0, 1.0);
       
    }
	
    public P2ElectrostaticDreiding(Space space, double AA, double BB, double CC) {
        super(space, AA, BB, CC);
    }
    
    public double energy(AtomSet atomSet) {
    	
        IAtomPositioned atom0 = (IAtomPositioned)atomSet.getAtom(0);
        IAtomPositioned atom1 = (IAtomPositioned)atomSet.getAtom(1);
    	dr01.Ev1Mv2(atom1.getPosition(), atom0.getPosition());
    	nearestImageTransformer.nearestImage(dr01);
        double r2 = dr01.squared();
        
        int index0 = atom0.getIndex();
        int index1 = atom1.getIndex();
        
        if (true){
        	return 0;
        }
        return constant*SpeciesParacetamol.Echarge[index0]*SpeciesParacetamol.Echarge[index1]/Math.sqrt(r2)
        		+ u(r2);        
        		
    }
    
    public IVector [] gradient(AtomSet atomSet) {
    	
        IAtomPositioned atom0 = (IAtomPositioned)atomSet.getAtom(0);
        IAtomPositioned atom1 = (IAtomPositioned)atomSet.getAtom(1);
    	dr01.Ev1Mv2(atom1.getPosition(), atom0.getPosition());
    	nearestImageTransformer.nearestImage(dr01);
        double r2 = dr01.squared();
        
        int index0 = atom0.getIndex();
        int index1 = atom1.getIndex();
        
        double sumU = du(r2) - constant*SpeciesParacetamol.Echarge[index0]*SpeciesParacetamol.Echarge[index1]/Math.sqrt(r2);
        
        gradient[1].Ea1Tv1(sumU/r2,dr);
        gradient[0].Ea1Tv1(-1,gradient[1]);
        
        return gradient;        		
    }
    
    private double constant = 162.0678e3; //conversion factor
	private static final long serialVersionUID = 1L;
}