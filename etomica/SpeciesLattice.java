package etomica;

import etomica.lattice.*;

/**
 * Species in which each molecule is an instance of a Bravais lattice.  
 * Expected that only one molecule would be constructed in a phase (although
 * nothing enforces this).
 *
 * @author David Kofke
 */
 
 /* History of changes
  * 09/04/02 (DAK) new
  */

public class SpeciesLattice extends Species {
    
    public SpeciesLattice() {
        this(Simulation.instance);
    }
    /**
     * Constructs species with default factory that makes a cubic lattice of spherical atoms.
     */
    public SpeciesLattice(Simulation sim) {
        this(sim, 
             new BravaisLattice.Factory(sim, new AtomFactoryMono(sim), sim.space.makeArrayD(3), 
             new PrimitiveCubic(sim,10.)));
    }
    public SpeciesLattice(Simulation sim, BravaisLattice.Factory factory) {
        super(sim, factory);
        setNMolecules(1);
    }
    
    
}//end of SpeciesLattice