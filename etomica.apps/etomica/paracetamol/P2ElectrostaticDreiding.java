package etomica.paracetamol;

import etomica.EtomicaElement;
import etomica.atom.AtomSet;
import etomica.atom.IAtomPositioned;
import etomica.simulation.Simulation;
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

public class P2ElectrostaticDreiding extends etomica.potential.P2IsotropicRepulsionDispersion implements EtomicaElement {
	
	public P2ElectrostaticDreiding(Simulation sim) {
        this(sim.getSpace(), sim.getDefaults().potentialWell, sim.getDefaults().atomSize, sim.getDefaults().potentialWell);
       
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
        
        return constant*AtomParacetamol.Echarge[index0]*AtomParacetamol.Echarge[index1]/Math.sqrt(r2)
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
        
        double sumU = du(r2) - constant*AtomParacetamol.Echarge[index0]*AtomParacetamol.Echarge[index1]/Math.sqrt(r2);
        
        gradient[1].Ea1Tv1(sumU/r2,dr);
        gradient[0].Ea1Tv1(-1,gradient[1]);
        
        return gradient;        		
    }
    
    private double constant = 162.0678e3;
	private static final long serialVersionUID = 1L;
}