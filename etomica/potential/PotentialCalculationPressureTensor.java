package etomica.potential;

import etomica.atom.AtomLeaf;
import etomica.atom.AtomTypeLeaf;
import etomica.atom.IAtom;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.atom.iterator.AtomsetIterator;
import etomica.integrator.IntegratorPhase;
import etomica.phase.Phase;
import etomica.space.ICoordinateKinetic;
import etomica.space.Space;
import etomica.space.Tensor;

/**
 * Calculates the pressure tensor by calculating the force on each atom, along
 * with including the kinetic portion (from the velocities or an Integrator).
 * If the simulation is non-dynamic (MC), the Integrator must be provided.
 */
public class PotentialCalculationPressureTensor extends PotentialCalculation {
        
    private static final long serialVersionUID = 1L;
    protected final Tensor pressureTensor;
    protected final Tensor workTensor;
    protected final Space space;
    protected final AtomIteratorLeafAtoms atomIterator;
    protected IntegratorPhase integrator;
    protected boolean warningPrinted;
    
    public PotentialCalculationPressureTensor(Space space) {
        this.space = space;
        pressureTensor = space.makeTensor();
        workTensor = space.makeTensor();
        atomIterator = new AtomIteratorLeafAtoms();
    }
    
    /**
	 * Adds the pressure tensor contribution based on the forces acting on each
     * pair of atoms produced by the iterator.
	 */
	public void doCalculation(AtomsetIterator iterator, Potential potential) {
		PotentialSoft potentialSoft = (PotentialSoft)potential;

		iterator.reset();
		while(iterator.hasNext()) {
            potentialSoft.gradient(iterator.next(), pressureTensor);
		}
	}
    
    public void setPhase(Phase newPhase) {
        atomIterator.setPhase(newPhase);
    }
    
    public void zeroSum() {
        pressureTensor.E(0);
    }

    /**
     * Sets an integrator to use a source for the temperature to compute the
     * kinetic portion of the pressure.  If running a dynamic simulation
     * (where the Atoms have velocities), this method should not be called.
     */
    public void setIntegrator(IntegratorPhase newIntegrator) {
        integrator = newIntegrator;
    }
    
    /**
     * Returns the pressure tensor based on a previous call to 
     * PotentialMaster.calculate
     */
    public Tensor getPressureTensor() {
        // now handle the kinetic part
        workTensor.E(0);
        atomIterator.reset();

        if (!atomIterator.hasNext()) {
            return pressureTensor;
        }
        
        AtomLeaf firstAtom = (AtomLeaf)atomIterator.nextAtom();
            
        if (firstAtom instanceof ICoordinateKinetic) {
            if (integrator != null) {
                warningPrinted = true;
                System.out.println("Ignoring Integrator's temperature and using actual Atom velocities.  You shouldn't have given me an Integrator.");
            }
        }
        else if (integrator == null) {
            throw new RuntimeException("Need an IntegratorPhase to provide temperature since this is a non-dynamic simulation");
        }
        else {
            int D = integrator.getPhase().getSpace().D();
            for (int i=0; i<D; i++) {
                pressureTensor.PE(atomIterator.size()*integrator.getTemperature());
            }
            return pressureTensor;
        }

        // simulation is dynamic, use the velocities
        atomIterator.reset();

        while (atomIterator.hasNext()) {
            IAtom atom = atomIterator.nextAtom();
            ICoordinateKinetic coord = (ICoordinateKinetic)atom;
            workTensor.Ev1v2(coord.getVelocity(), coord.getVelocity());
            workTensor.TE(((AtomTypeLeaf)atom.getType()).getMass());
            pressureTensor.PE(workTensor);
        }
        
        return pressureTensor;
    }
}
