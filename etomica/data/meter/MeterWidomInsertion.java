package etomica.data.meter;

import etomica.EtomicaInfo;
import etomica.action.AtomActionTranslateTo;
import etomica.atom.IAtom;
import etomica.atom.SpeciesAgent;
import etomica.data.DataSourceScalar;
import etomica.integrator.IntegratorPhase;
import etomica.phase.Phase;
import etomica.space.Space;
import etomica.species.Species;
import etomica.units.Null;

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
public class MeterWidomInsertion extends DataSourceScalar {

    public MeterWidomInsertion(Space space) {
        super("exp(-\u03BC/kT)", Null.DIMENSION);//"\u03BC" is Unicode for greek "mu"
        setNInsert(100);
        setResidual(true);
        atomTranslator = new AtomActionTranslateTo(space); 
    }

    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo(
            "Chemical potential via Widom's ghost-particle insertion method");
        return info;
    }

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
     * Performs a Widom insertion average, doing nInsert insertion attempts
     * Temperature used to get exp(-uTest/kT) is that of the integrator for the
     * phase
     * 
     * @return the sum of exp(-uTest/kT)/nInsert, multiplied by V/N if 
     * <code>residual</code> is false
     */
    public double getDataAsScalar() {
        if (integrator == null) throw new IllegalStateException("must call setPhase before using meter");
        Phase phase = integrator.getPhase();
        double sum = 0.0; //sum for local insertion average
        SpeciesAgent agent = phase.getAgent(species);
        agent.addChildAtom(testMolecule);
        energyMeter.setTarget(testMolecule);
        for (int i = nInsert; i > 0; i--) { //perform nInsert insertions
            atomTranslator.setDestination(phase.getBoundary().randomPosition());
            atomTranslator.actionPerformed(testMolecule);
            double u = energyMeter.getDataAsScalar();
            sum += Math.exp(-u / integrator.getTemperature());
        }

        agent.removeChildAtom(testMolecule);

        if (!residual) {
            // multiply by V/N
            sum *= phase.volume() / phase.getAgent(species).getNMolecules();
        }
        return sum / nInsert; //return average
    }

    /**
     * Returns the integrator associated with this class.  The phase, potentialMaster
     * and temperature are taken from the integrator.
     */
    public IntegratorPhase getIntegrator() {
        return integrator;
    }

    /**
     * Sets the integrator associated with this class.  The phase, potentialMaster
     * and temperature are taken from the integrator.
     */
    public void setIntegrator(IntegratorPhase newIntegrator) {
        integrator = newIntegrator;
        energyMeter = new MeterPotentialEnergy(integrator.getPotential());
        energyMeter.setPhase(integrator.getPhase());
    }

    private static final long serialVersionUID = 1L;
    private IntegratorPhase integrator;

    /**
     * Number of insertions attempted in each call to currentValue Default is
     * 100
     */
    private int nInsert;
    private Species species;
    private IAtom testMolecule;// prototype insertion molecule
    private boolean residual; // flag to specify if total or residual chemical
                              // potential evaluated. Default true
    private AtomActionTranslateTo atomTranslator;

    private MeterPotentialEnergy energyMeter;
}
