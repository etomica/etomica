package etomica;
import java.awt.Color;
import etomica.units.Dimension;

//Java2 imports
//import java.util.Iterator;

import etomica.utility.Iterator;


/**
 * Meter to measure the chemical potential of a species via the Widom insertion method
 */
public class MeterWidomInsertion extends Meter implements Molecule.Container, EtomicaElement {

    /**
     * Number of insertions attempted in each call to currentValue
     * Default is 100
     */
    private int nInsert;
    /**
     * Species for which chemical potential is evaluated
     */
    private Species species;
    private Species.Agent speciesAgent;
    private Molecule molecule;  //prototype insertion molecule
    private boolean residual;   //flag to specify if total or residual chemical potential evaluated. Default true
    private DisplayPhase display; //used to show location of inserted atoms in the display
    
    public MeterWidomInsertion() {
        this(Simulation.instance);
    }
    public MeterWidomInsertion(Simulation sim) {
        super(sim);
        setLabel("exp(\u03BC/kT)");  //"\u03BC" is Unicode for greek "mu"
        nInsert = 100;
        setResidual(true);
        
		//add mediator so that by default first species added to simulation is used for sampling
        sim.elementCoordinator.addMediatorPair(new Mediator.MeterSpecies(sim.elementCoordinator) {
            public void add(Species species) {
                if(MeterWidomInsertion.this.wasAdded()) MeterWidomInsertion.this.setSpecies(species);
            }
            public void add(MeterAbstract meter) {
                if(meter != MeterWidomInsertion.this) return;
                for(Iterator ip=mediator.parentSimulation().speciesList.iterator(); ip.hasNext(); ) {
                    Species species = (Species)ip.next();
                    if(species.wasAdded())  {//will make first species the one
                        MeterWidomInsertion.this.setSpecies(species);
                        return;
                    }
                }
            }
        });
    }

    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Chemical potential via Widom's ghost-particle insertion method");
        return info;
    }

    /**
     * Declaration that this meter does not use the boundary object of phase when making its measurements
     */
    public final boolean usesPhaseBoundary() {return false;}
    /**
     * Declaration that this meter does not use the iteratorFactory of phase when making its measurements
     */
    public final boolean usesPhaseIteratorFactory() {return false;}
        
    /**
     * Constructor used if desired to display inserted positions in the given DisplayConfiguration object
     */
    public MeterWidomInsertion(Species s, DisplayPhase d) {
        this();
        display = d;
        setSpecies(s);
        setActive(false);
    }
    
    public Dimension getDimension() {return Dimension.NULL;}  //need to modify to check for isResidual
    
    /**
     * Sets flag specifying if full or residual chemical potential is computed
     * Default is <code>true</code> (only residual is computed)
     */
    public void setResidual(boolean b) {residual = b; reset();}
    /**
     * Accessor for flag specifying if full or residual chemical potential is computed
     */
    public boolean isResidual() {return residual;}
    
    /**
     * Sets the phase and gets handle to appropriate species agent
     */
    public void setPhase(Phase p) {
        super.setPhase(p);
        if(species != null) speciesAgent = species.getAgent(phase);
        setActive(true);
    }
    /**
     * Sets the species, takes a prototype molecule, and gets handle to appropriate species agent in phase
     */
    public void setSpecies(Species s) {
        species = s;
        addMolecule(species.getMolecule());
        if(phase != null) speciesAgent = species.getAgent(phase);
    }
    /**
     * Accessor for the species for which chemical potential is evaluated
     */
    public Species getSpecies(Species s) {return species;}
    
    /**
     * Number of Widom insertions attempted with each call to currentValue
     */
    public void setNInsert(int n) {nInsert = n;}
    /**
     * Accessor to number of Widom insertions attempted with each call to currentValue
     */
    public int getNInsert() {return nInsert;}
    
    /**
     * Performs a Widom insertion average, doing nInsert insertion attempts
     * Temperature used to get exp(-uTest/kT) is that of the integrator for the phase
     * @return the sum of exp(-uTest/kT)/nInsert, multiplied by n<sub>i</sub>/V if <code>residual</code> is false
     */
    public double currentValue()
    {
        double sum = 0.0;                         //sum for local insertion average
        Molecule testMolecule = molecule;         //save a handle to the molecule
        phase.addMolecule(molecule,speciesAgent); //place molecule in phase (removing it from this meter, thus setting molecule to null)
        testMolecule.firstAtom().setColor(Color.red);
        for(int i=nInsert; i>0; i--) {            //perform nInsert insertions
            testMolecule.translateTo(phase.randomPosition());  //select random position
            if(display != null && i % 10 ==0) display.repaint();
            double u = phase.energy.meterPotential().currentValue(testMolecule); //compute energy          
            if(u < Double.MAX_VALUE)              //add to test-particle average
                sum += Math.exp(-u/(phase.integrator().temperature()));
        }
        this.addMolecule(testMolecule);           //remove molecule from phase by placing it back in this meter
        if(!residual) sum *= speciesAgent.moleculeCount()/phase.volume(); //multiply by N/V
        return sum/(double)nInsert;               //return average
    }
    
    //Molecule.Container methods
    
    /**
     * Method to officially add the test molecule to this meter
     */
    public void addMolecule(Molecule m) {
        m.container.removeMolecule(m);
        molecule = m;
        m.setContainer(this);
    }
    /**
     * Method to officially remove the test molecule from this meter
     */
    public void removeMolecule(Molecule m) {
        if(m == molecule) molecule = null;
    }
    
}   