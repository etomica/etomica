package etomica.data.meter;

import etomica.Atom;
import etomica.AtomTypeLeaf;
import etomica.EtomicaInfo;
import etomica.Phase;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.atom.iterator.AtomIteratorPhaseDependent;
import etomica.data.DataSourceScalar;
import etomica.space.ICoordinateKinetic;
import etomica.units.Dimension;

/**
 * Meter for the total kinetic energy in a phase
 * Computes total KE by summing values of KE returned by every atom in the phase.
 * A different phase-dependent atom integrator may be set to permit calculation
 * over a particular set of atoms in the phase.
 */
 
 /* History of changes
  * 7/03/02  Added non-registering constructor (space argument)
  */
public class MeterKineticEnergy extends DataSourceScalar implements Meter
{
    private AtomIteratorPhaseDependent iterator;
    
    public MeterKineticEnergy() {
        super("Kinetic Energy",Dimension.ENERGY);
        setIterator(new AtomIteratorLeafAtoms());
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Total kinetic energy of molecular motion in a phase");
        return info;
    }

    /**
     * Returns the iterator that defines the atoms summed for their
     * kinetic energy.
     */
	public AtomIteratorPhaseDependent getIterator() {
		return iterator;
	}
	
	/**
	 * Sets the iterator that defines the atoms which are summed for
	 * their total kinetic energy.  Default is a leaf-atom iterator,
	 * giving all leaf atoms in the phase.
	 * @param iterator
	 */
	public void setIterator(AtomIteratorPhaseDependent iterator) {
		this.iterator = iterator;
	}
	
	/**
	 * Returns the total kinetic energy summed over all atoms produced by
	 * the iterator when applied to the given phase.  Does not include contributions
     * from atoms having infinite mass (it assumes they are stationary).
	 */
    public double getDataAsScalar() {
        if (phase == null) throw new IllegalStateException("must call setPhase before using meter");
        double ke = 0.0;
        iterator.setPhase(phase);
        iterator.reset();
        while(iterator.hasNext()) {    //consider doing this with an allAtoms call
            Atom atom = iterator.nextAtom();
            double mass = ((AtomTypeLeaf)atom.type).getMass();
            if(mass == Double.POSITIVE_INFINITY) continue;
            ke += 0.5*mass*((ICoordinateKinetic)atom.coord).velocity().squared();
        }
        return ke;
    }//end of getDataAsScalar
    
    /**
     * @return Returns the phase.
     */
    public Phase getPhase() {
        return phase;
    }
    /**
     * @param phase The phase to set.
     */
    public void setPhase(Phase phase) {
        this.phase = phase;
    }

    private Phase phase;
 }