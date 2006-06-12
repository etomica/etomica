package etomica.meam;
import etomica.action.PhaseImposePbc;
import etomica.action.activity.ActivityIntegrate;
import etomica.action.activity.Controller;
import etomica.atom.AtomTypeLeaf;
import etomica.atom.AtomTypeSphere;
import etomica.config.Configuration;
import etomica.config.ConfigurationLattice;
import etomica.data.AccumulatorHistory;
import etomica.data.DataInfo;
import etomica.data.DataPump;
import etomica.data.meter.MeterEnergy;
import etomica.data.meter.MeterKineticEnergy;
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
import etomica.util.HistoryCollapsingAverage;

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
    public PotentialMEAM potentialN;
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
    	AccumulatorHistory energyAccumulator = new AccumulatorHistory(HistoryCollapsingAverage.FACTORY);
        DisplayPlot plot = new DisplayPlot();
        energyAccumulator.setDataSink(plot.getDataSet().makeDataSink());
    	DataPump energyManager = new DataPump(energyMeter,energyAccumulator);
    	MeterKineticEnergy kineticMeter = new MeterKineticEnergy();
    	kineticMeter.setPhase(sim.phase);
    	AccumulatorHistory kineticAccumulator = new AccumulatorHistory(HistoryCollapsingAverage.FACTORY);
     	DisplayPlot plotKE = new DisplayPlot();
    	kineticAccumulator.setDataSink(plotKE.getDataSet().makeDataSink());
    	DataPump kineticManager = new DataPump(kineticMeter, kineticAccumulator);
    	//energyAccumulator.setBlockSize(50);
        IntervalActionAdapter adapter = new IntervalActionAdapter(energyManager, sim.integrator);
        adapter.setActionInterval(1);
        IntervalActionAdapter kineticAdapter = new IntervalActionAdapter(kineticManager, sim.integrator);
        kineticAdapter.setActionInterval(1);
        SimulationGraphic simgraphic = new SimulationGraphic(sim);
    	simgraphic.makeAndDisplayFrame();
    	simgraphic.panel().add(plot.graphic());
    	simgraphic.panel().add(plotKE.graphic());
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
        integrator.setTimeStep(0.001);
        integrator.setTemperature(Kelvin.UNIT.toSim(298));
        integrator.setThermostatInterval(10);
        integrator.setIsothermal(true);
        activityIntegrate = new ActivityIntegrate(this,integrator);
        activityIntegrate.setSleepPeriod(2);
        getController().addAction(activityIntegrate);
        species = new SpeciesSpheresMono(this);
        species.setNMolecules(216);
        
        //Sn
        
        //The value of the atomic weight of Sn used is from the ASM Handbook. 
        ((AtomTypeLeaf)species.getFactory().getType()).setMass(118.69);
        //The "distance of closest approach" for atoms in beta-tin, as given on 
        //p. 639 of Cullity and Strock's "Elements of X-Ray Diffraction" (2001), 
        //is used as the diameter for Sn atoms in the simulation.  This value was
        //taken from Pearson [G.9, Vol.2].  This is essentially the same as the
        //distance given in the ASM Handbook for a beta-tin atom and its four
        //closest neighbors: 3.02 Angstroms.
        ((AtomTypeSphere)species.getFactory().getType()).setDiameter(3.022); 
        // This forces the MEAMP2 to request an agentIndex from Atom.
        //pseudoPotential2 = new MEAMP2(this, ParameterSetMEAM.Sn);
        phase = new Phase(this);
        //The dimensions of the simulation box must be proportional to those of
        //the unit cell to prevent distortion of the lattice.  The values for the 
        //lattice parameters for tin's beta phase (a = 5.8314 angstroms, c = 3.1815 
        //angstroms) are taken from the ASM Handbook. 
        phase.setDimensions(new Vector3D(5.8314*3, 5.8314*3, 3.1815*6));
        PrimitiveTetragonal primitive = new PrimitiveTetragonal(space, 5.8318, 3.1819);
        //Alternatively, using the parameters calculated in Ravelo & Baskes (1997)
        //phase.setDimensions(new Vector3D(5.92*3, 5.92*3, 3.23*6));
        //PrimitiveTetragonal primitive = new PrimitiveTetragonal(space, 5.92, 3.23);
        LatticeCrystal crystal = new LatticeCrystal(new Crystal(
        		primitive, new BasisBetaSnA5(primitive)));
        

        //pseudoPotential1 = new MEAMPInitial(space, pseudoPotential2.getPhaseAgentManager());
        //potential = new MEAMPMany(space, ParameterSetMEAM.Sn, pseudoPotential2.getPhaseAgentManager());
        //this.potentialMaster.addPotential(potential, new Species[]{species});
        //this.potentialMaster.addPotential(pseudoPotential2, new Species[]{species,species});
        //this.potentialMaster.addPotential(pseudoPotential1, new Species[]{species});    

        //Cu
        /**
	    ((AtomTypeLeaf)species.getFactory().getType()).setMass(63.546);//Cullity & Stock
	    ((AtomTypeSphere)species.getFactory().getType()).setDiameter(2.56);
	    phase = new Phase(this);
	    phase.setDimensions(new Vector3D(3.6148*3, 3.6148*3, 3.6148*6));
	    PrimitiveCubic primitive = new PrimitiveCubic(space, 3.6148);
	    LatticeCrystal crystal = new LatticeCrystal(new Crystal(
		        primitive, new BasisCubicFcc(primitive)));
	    **/
        
		//General   
		Configuration config = new ConfigurationLattice(crystal);
        //phase.setConfiguration(config);  // kmb remove 8/3/05
		config.initializeCoordinates(phase);  // kmb added 8/3/05
        
        
        //N-body potential
		potentialN = new PotentialMEAM (space, ParameterSetMEAM.Sn);//Sn
		//potentialN = new PotentialMEAM (space, ParameterSetMEAM.Cu);//Cu
		
		//System.out.println(ParameterSetMEAM.Sn.Ec);
        this.potentialMaster.addPotential(potentialN, new Species[]{species});    
        
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