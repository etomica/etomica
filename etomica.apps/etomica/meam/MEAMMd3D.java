package etomica.meam;
import etomica.action.PhaseImposePbc;
import etomica.action.activity.ActivityIntegrate;
import etomica.action.activity.Controller;
import etomica.atom.AtomTypeLeaf;
import etomica.atom.AtomTypeSphere;
import etomica.config.Configuration;
import etomica.config.ConfigurationLattice;
import etomica.data.AccumulatorAverage;
import etomica.data.DataInfo;
import etomica.data.DataPump;
import etomica.data.meter.MeterEnergy;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.DisplayPhase;
import etomica.graphics.DisplayPlot;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.integrator.IntervalActionAdapter;
import etomica.lattice.Crystal;
import etomica.lattice.LatticeCrystal;
import etomica.lattice.crystal.BasisBetaSnA5;
import etomica.lattice.crystal.PrimitiveTetragonal;
import etomica.phase.Phase;
import etomica.simulation.Simulation;
import etomica.space3d.Space3D;
import etomica.space3d.Vector3D;
import etomica.species.Species;
import etomica.species.SpeciesSpheresMono;
import etomica.units.Kelvin;

/**
 * Molecular-Dynamics Simulation Using the Modified Embedded-Atom Method 
 * (MEAM) Potential.  
 * 
 * The MEAM potential is intended for use with metallic and covalently-bonded
 * solid systems.
 * 
 * The MEAM potential for an atom is built using terms describing parts of the
 * relationships between the atom and each of its neighbors, the number of which 
 * is determined by a cutoff and/or screening function.  Each type of pair-
 * wise term is summed over all the neighbors, and then used in expressions 
 * describing the embedding energy and the repulsive energy of the atom.   
 * Effectively, the MEAM potential is a many-body potential.  
 * 
 * Currently, there is no framework in Etomica for looping through each atom 
 * and all of its neigbors at once, so that the potential for an atom might be 
 * calculated in a single potential class.  Happily, the basis of the MEAM 
 * potential upon pair-wise terms makes the creation of such a scheme unnecessary. 
 * Instead, the potential was implemented using a method affectionately referred
 * to as Andrew's Grand Master Plan, although there is some debate about whether or
 * not using both "Andrew" and "Grand Master" in the title is redundant.  
 * 
 * Andrew's Grand Master Plan: 
 * 
 * There is an array with an element for each atom in the system.  Each of these 
 * elements is an "agent"; i.e., it holds descriptions or values for the atom to
 * which it corresponds.  In this case, the agent is a wrapper, creatively called 
 * "agents".  Within each atom's wrapper are two arrays, one called "sums", the 
 * other "gradient Sums".  The array "sums" has an element for each sum required to 
 * calculate the potential of the atom.  The array "gradientSums" has an element 
 * for each sum required to calculate the gradient of the atom's potential.  
 * 
 * Within each time increment, the value of every element in these two arrays is 
 * set to zero, or the zero vector, by the MEAMPInitial class.  This class functions 
 * as a one-body potential.  The MEAMP2 class then cycles through each atom pair 
 * just like a two-body potential, calculating the values of the pair-wise terms 
 * and storing them in the appropriate elements of the arrays. As the pairings 
 * between an atom and each of its neighbors are evaluated, the pair-wise terms 
 * relevant to that atom are summed within that atom's "sums" array or 
 * "gradientSums" array.  I like to think of the elements as bins, in which the 
 * values gradually accumulate into the required sums.  Both the MEAMPInitial class 
 * and the MEAMP2 class have energy and gradient methods, but they must contribute
 * nothing to the value of the potential and its gradient, at least not directly. 
 * In both classes, these methods return either zero or the zero vector.  
 * 
 * The gradient and its vector are finally calculated by the MEAMPMany class, which
 * functions as a one-body potential class.  All of the sums required for the
 * computation of an atom's potential and gradient may be accessed within the
 * relevant atom's bins in its wrapper.  
 * 
 * This class was adapted from LjMd3D.java by K.R. Schadel and A. Schultz in July 
 * 2005.  Intitially, it employed a version of the embedded-atom method potential, 
 * and was later adapted in February 2006 to use the modified embedded-atom method
 * potential.
 */
 
public class MEAMMd3D extends Simulation {
    
    public IntegratorVelocityVerlet integrator;
    public SpeciesSpheresMono species;
    public Phase phase;
    public MEAMPInitial pseudoPotential1;
    public MEAMP2 pseudoPotential2;
    public MEAMPMany potential;
    public Controller controller;
    public DisplayPhase display;
    public DisplayPlot plot;
    public MeterEnergy energy;
    public ActivityIntegrate activityIntegrate;
    public DataInfo info2;

    public static void main(String[] args) {
    	MEAMMd3D sim = new MEAMMd3D();
    	MeterPotentialEnergy energyMeter = new MeterPotentialEnergy(sim.potentialMaster);
    	energyMeter.setPhase(sim.phase);
    	AccumulatorAverage energyAccumulator = new AccumulatorAverage(sim);
        DataPump energyManager = new DataPump(energyMeter,energyAccumulator);
        energyAccumulator.setBlockSize(50);
        IntervalActionAdapter adapter = new IntervalActionAdapter(energyManager, sim.integrator);
        adapter.setActionInterval(5);
        SimulationGraphic simgraphic = new SimulationGraphic(sim);
    	simgraphic.makeAndDisplayFrame();
        ColorSchemeByType colorScheme = ((ColorSchemeByType)((DisplayPhase)simgraphic.displayList().getFirst()).getColorScheme());
    	colorScheme.setColor(sim.species.getMoleculeType(),java.awt.Color.red);
    	//sim.activityIntegrate.setMaxSteps(1000);
    	//sim.getController().run();
    	//DataGroup data = (DataGroup)energyAccumulator.getData(); // kmb change type to Data instead of double[]
        //double PE = ((DataDouble)data.getData(AccumulatorAverage.StatType.AVERAGE.index)).x
                    // /sim.species.getAgent(sim.phase).getNMolecules();  // kmb changed 8/3/05
        //double PE = data[AccumulatorAverage.AVERAGE.index]/sim.species.getAgent(sim.phase).getNMolecules();  // orig line
        //System.out.println("PE(eV)="+ElectronVolt.UNIT.fromSim(PE));
    }
    
    public MEAMMd3D() {
        super(Space3D.getInstance()); //INSTANCE); kmb change 8/3/05
        integrator = new IntegratorVelocityVerlet(this);
        integrator.setTimeStep(0.0005);
        integrator.setTemperature(Kelvin.UNIT.toSim(298));
        activityIntegrate = new ActivityIntegrate(this,integrator);
        activityIntegrate.setSleepPeriod(2);
        getController().addAction(activityIntegrate);
        species = new SpeciesSpheresMono(this);
        species.setNMolecules(216);
        ((AtomTypeLeaf)species.getFactory().getType()).setMass(118.71);
        //The diameter given is really the equilibrium distance between the atoms.
        ((AtomTypeSphere)species.getFactory().getType()).setDiameter(3.44); 
        // This forces the MEAMP2 to request an agentIndex from Atom.
        pseudoPotential2 = new MEAMP2(this, ParameterSetMEAM.Sn);
        phase = new Phase(this);
        phase.setDimensions(new Vector3D(5.8318*3, 5.8318*3, 3.1819*6));
        PrimitiveTetragonal primitive = new PrimitiveTetragonal(space, 5.8318, 3.1819);
        LatticeCrystal crystal = new LatticeCrystal(new Crystal(
        		primitive, new BasisBetaSnA5(primitive)));
        Configuration config = new ConfigurationLattice(crystal);
//        phase.setConfiguration(config);  // kmb remove 8/3/05
        config.initializeCoordinates(phase);  // kmb added 8/3/05
        pseudoPotential1 = new MEAMPInitial(space, pseudoPotential2.getPhaseAgentManager());
        potential = new MEAMPMany(space, ParameterSetMEAM.Sn, pseudoPotential2.getPhaseAgentManager());
        this.potentialMaster.addPotential(potential, new Species[]{species});
        this.potentialMaster.addPotential(pseudoPotential2, new Species[]{species,species});
        this.potentialMaster.addPotential(pseudoPotential1, new Species[]{species});    

        integrator.setPhase(phase);
        PhaseImposePbc imposepbc = new PhaseImposePbc();
        imposepbc.setPhase(phase);
        integrator.addListener(imposepbc);
		
        // IntegratorCoordConfigWriter - Displacement output (3/1/06 - MS)
        //IntegratorCoordConfigWriter coordWriter = new IntegratorCoordConfigWriter(space, "MEAMoutput");
        //coordWriter.setPhase(phase);
        //coordWriter.setIntegrator(integrator);
        //coordWriter.setWriteInterval(100);
        
        // Control simulation lengths
        //activityIntegrate.setMaxSteps(500);

		energy = new MeterEnergy(potentialMaster);
    }
    
}