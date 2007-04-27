package etomica.paracetamol;

import etomica.EtomicaElement;
import etomica.atom.AtomLeaf;
import etomica.atom.AtomSet;
import etomica.simulation.Simulation;
import etomica.space.IVector;
import etomica.space.Space;

/**
 * Dreiding: Lennard-Jones non-bonding potential.
 * Spherically symmetric potential of the form u(r) = D0*[(r0/r)^12 - 2*(r0/r)^6]
 * where D0 describes the van der Waals well depth [unit Kelvin ], 
 * and r0 is the van der Waals bond length [unit Amstrom ].
 * 
 * r0 = sigma; D0 = epsilon
 *
 * @author Tai Tan
 */

public class P2ElectrostaticDreiding extends etomica.potential.P2LennardJonesDreiding implements EtomicaElement {
	
	public P2ElectrostaticDreiding(Simulation sim) {
        this(sim.getSpace(), sim.getDefaults().atomSize, sim.getDefaults().potentialWell);
       
    }
	
    public P2ElectrostaticDreiding(Space space, double sigma, double epsilon) {
        super(space, sigma, epsilon);
    }
    
    public double energy(AtomSet atomSet) {
    	
    	AtomLeaf atom0 = (AtomLeaf)atomSet.getAtom(0);
    	AtomLeaf atom1 = (AtomLeaf)atomSet.getAtom(1);
    	dr01.Ev1Mv2(atom1.getPosition(), atom0.getPosition());
    	nearestImageTransformer.nearestImage(dr01);
        double r2 = dr01.squared();
        
        int index0 = atom0.getIndex();
        int index1 = atom1.getIndex();
        
        return constant*AtomParacetamol.Echarge[index0]*AtomParacetamol.Echarge[index1]/Math.sqrt(r2)
        		+ u(r2);        
        		
    }
    
    public IVector [] gradient(AtomSet atomSet) {
    	
    	AtomLeaf atom0 = (AtomLeaf)atomSet.getAtom(0);
    	AtomLeaf atom1 = (AtomLeaf)atomSet.getAtom(1);
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
