package etomica.potential;

import etomica.Potential;
import etomica.Space;
import etomica.Species;

/**
 * Zero-body potential implementing the long-range correction 
 * that compensates for truncation of the potential.  Normally, the
 * concrete instance of this class is defined as an inner class to
 * a PotentialTruncation class.  The PotentialTruncation class defines
 * the scheme used to truncate the potential, and its inner Potential0Lrc
 * subclass implements the formulas needed to integrate the potential
 * and its derivatives over the range affected by the truncation.<br>
 *
 * The components and procedures related to the long-range correction 
 * are as follows.  The Potential2 constructor defines a PotentialTruncation
 * that determines whether and how the pair potential may be truncated.
 * If the PotentialTruncation is not passed as a parameter to the constructor,
 * a default is used as indicated by the Default.TRUNCATE_POTENTIALS boolean
 * flag.  By default, all pair potentials are truncated using the 
 * PotentialTruncationSimple scheme; if TRUNCATE_POTENTIALS is set to
 * false, the PotentialTruncation.Null truncation is applied to all new pair potentials.  The
 * PotentialTruncation for a potential cannot be changed after the potential
 * is instantiated.<br>
 *
 * Each PotentialTruncation class defines an inner Potential0Lrc subclass that 
 * provides the appropriate methods for computing the long-range correction 
 * to the energy and its first two derivatives.  This class is instantiated
 * by the PotentialTruncation class in its constructor.  Upon its instantiation,
 * the Potential0Lrc class is added to the group of long-range correction potentials
 * that is kept by a single Potential0GroupLrc instance in the PotentialMaster.<br>
 *
 * Before the calculate method of PotentialMaster is called to compute something,
 * its set(Phase) method must have been previously called, which identifies to
 * all potentials (including Potential0GroupLrc) which phase is subject to the 
 * ensuing calculation.  Potential0Group ignores this notification if the
 * given phase is the same as the one specified in the previous call; otherwise
 * it passes the identified phase to all the set(Phase) methods (inherited from Potential0)
 * of the Potential0Lrc classes it holds.<br>
 *
 * Then when the calculate(IteratorDirective, PotentialCalculation) method of
 * PotentialMaster is invoked, it passes the call on to the calculate methods of
 * its child potentials, including the Potential0GroupLrc instance if present.
 * Potential0GroupLrc checks that a phase has been specified, that its
 * enableLrc flag is <b>true</b> (the default), and that the given iteratorDirective's
 * includeP0Lrc flag is also <b>true</b> (default is <b>false</b>).  If so, it 
 * calls the calculate methods of all child Potential0Lrc classes.
 *
 * The Potential0Lrc class will use the volume from the specified phase and the
 * size method of the iterator of its associated potential to determine the
 * pair density.  The Potential0Lrc methods are called if the PotentialCalculation
 * implements Potential0Calculation.
 *
 * @author David Kofke
 */
public abstract class Potential0Lrc extends Potential0 {
    
    public final Potential potential;
    protected Species[] species;
    
    /**
     * Constructor requires PotentialMaster argument, and calls superclass
     * constructor such that this class is added dirctly to PotentialMaster's
     * lrcMaster instance, which manages the lrc potentials.
     * @param parent  Potential master for the simulation.
     */
    public Potential0Lrc(Space space, Potential potential) {
        super(space);
        this.potential = potential;
        potential.setLrc(this);
    }  
    
	public void setSpecies(Species[] species) {
		this.species = new Species[species.length];
		System.arraycopy(species, 0, this.species, 0, species.length);
	}
  
    
    /**
     * Long-range correction to the energy u.
     */
    public abstract double uCorrection(double pairDensity);
        
    /**
    * Long-range correction to r*du/dr.
    */
    public abstract double duCorrection(double pairDensity);
    
    /**
    * Long-range correction to r^2 d2u/dr2.
    */
    public abstract double d2uCorrection(double pairDensity);
    

}//end of Potential0Lrc
