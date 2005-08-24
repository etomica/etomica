package etomica.data.meter;

import etomica.EtomicaInfo;
import etomica.action.AtomActionTranslateTo;
import etomica.atom.Atom;
import etomica.data.DataSourceScalar;
import etomica.phase.Phase;
import etomica.potential.PotentialMaster;
import etomica.space.Space;
import etomica.species.Species;
import etomica.units.Dimension;
import etomica.util.Default;

/**
 * Meter to measure the chemical potential (as its exponent: exp(-mu/kT)) of a
 * species via the Widom insertion method. Call to getDataAsScalar returns a
 * Widom-insertion average (i.e., sum of exp(-energy/kT) using a configurable
 * number of insertion trials) for a molecule of the desired species in the
 * given phase in its current configuration. It is expected that repeated calls
 * will be averaged to give an overall insertion average for many configurations
 * of the phase.
 * <br>
 * Meter can be configured to measure the residual chemical potential -- above
 * the ideal-gas value -- or the full chemical potential.  If not configured for
 * the residual chemical potential, the insertion average will be multiplied by
 * the species number density in the phase in its current configuration.  Default
 * behavior will give residual chemical potential.
 * 
 * @author David Kofke
 */
public class MeterWidomInsertion extends DataSourceScalar implements Meter {

	public MeterWidomInsertion(Space space, PotentialMaster potentialMaster) {
		super("exp(-\u03BC/kT)", Dimension.NULL);//"\u03BC" is Unicode for greek "mu"
		energyMeter = new MeterPotentialEnergy(potentialMaster);
		nInsert = 100;
		setResidual(true);
		setTemperature(Default.TEMPERATURE);
        atomTranslator = new AtomActionTranslateTo(space); 

	}

	public static EtomicaInfo getEtomicaInfo() {
		EtomicaInfo info = new EtomicaInfo(
				"Chemical potential via Widom's ghost-particle insertion method");
		return info;
	}

	/**
	 * Constructor used if desired to display inserted positions in the given
	 * DisplayConfiguration object
	 */
	/*
	 * public MeterWidomInsertion(Species s, DisplayPhase d) { this(); display =
	 * d; setSpecies(s); setActive(false); }
	 */

	/**
	 * Sets flag specifying if full or residual chemical potential is computed
	 * Default is <code>true</code> (only residual is computed)
	 */
	public void setResidual(boolean b) {
		residual = b;
	}

	/**
	 * Accessor for flag specifying if full or residual chemical potential is
	 * computed
	 */
	public boolean isResidual() {
		return residual;
	}

	/**
	 * Sets the species, takes a prototype molecule, and gets handle to
	 * appropriate species agent in phase
	 */
	public void setSpecies(Species s) {
		species = s;
		testMolecule = s.moleculeFactory().makeAtom();
	}

	/**
	 * Accessor for the species for which chemical potential is evaluated
	 */
	public Species getSpecies() {
		return species;
	}

	/**
	 * Number of Widom insertions attempted with each call to currentValue
	 */
	public void setNInsert(int n) {
		nInsert = n;
	}

	/**
	 * Accessor to number of Widom insertions attempted with each call to
	 * currentValue
	 */
	public int getNInsert() {
		return nInsert;
	}

	/**
	 * @return Returns the temperature.
	 */
	public double getTemperature() {
		return temperature;
	}

	/**
	 * @param temperature
	 *            The temperature to set.
	 */
	public void setTemperature(double temperature) {
		this.temperature = temperature;
	}

	/**
	 * @return Dimension.TEMPERATURE
	 */
	public Dimension getTemperatureDimension() {
		return Dimension.TEMPERATURE;
	}

	/**
	 * Performs a Widom insertion average, doing nInsert insertion attempts
	 * Temperature used to get exp(-uTest/kT) is that of the integrator for the
	 * phase
	 * 
	 * @return the sum of exp(-uTest/kT)/nInsert, multiplied by n <sub>i
	 *         </sub>/V if <code>residual</code> is false
	 */
	public double getDataAsScalar() {
        if (phase == null) throw new IllegalStateException("must call setPhase before using meter");
		double sum = 0.0; //sum for local insertion average
        phase.addMolecule(testMolecule, phase.getAgent(species));
		energyMeter.setTarget(testMolecule);
		for (int i = nInsert; i > 0; i--) { //perform nInsert insertions
            atomTranslator.setDestination(phase.randomPosition());
            atomTranslator.actionPerformed(testMolecule);
			//            if(display != null && i % 10 ==0) display.repaint();
			double u = energyMeter.getDataAsScalar();
			if (u < Double.POSITIVE_INFINITY) //add to test-particle average
				sum += Math.exp(-u / temperature);
		}

        phase.removeMolecule(testMolecule);

		if (!residual)
			sum *= phase.volume() / phase.getAgent(species).moleculeCount(); //multiply
		// by
		// V/N
		return sum / nInsert; //return average
	}
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
        energyMeter.setPhase(phase);
    }

    private Phase phase;

	/**
	 * Number of insertions attempted in each call to currentValue Default is
	 * 100
	 */
	private int nInsert;
	private Species species;
	private Atom testMolecule; //prototype insertion molecule
	private double temperature;
	private boolean residual; //flag to specify if total or residual chemical
								// potential evaluated. Default true
	private AtomActionTranslateTo atomTranslator;
    
	MeterPotentialEnergy energyMeter;
	//  private DisplayPhase display; //used to show location of inserted atoms
	// in the display

	// /*
	//    public static void main(String[] args) {
	//        etomica.simulations.HSMD2D sim = new etomica.simulations.HSMD2D();
	//        MeterWidomInsertion meter = new MeterWidomInsertion();
	//        etomica.graphics.DisplayBox box = new
	// etomica.graphics.DisplayBox((DatumSource)meter);
	//        box.setWhichValue(MeterAbstract.AVERAGE);
	//        sim.elementCoordinator.go();
	//        etomica.graphics.SimulationGraphic.makeAndDisplayFrame(sim);
	//    }
	//    */
}//end of MeterWidomInsertion
