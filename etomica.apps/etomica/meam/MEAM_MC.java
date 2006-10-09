package etomica.meam;
import etomica.action.PhaseImposePbc;
import etomica.action.activity.ActivityIntegrate;
import etomica.action.activity.Controller;
import etomica.atom.AtomTypeLeaf;
import etomica.atom.AtomTypeSphere;
import etomica.chem.elements.Copper;
import etomica.chem.elements.Silver;
import etomica.chem.elements.Tin;
import etomica.config.Configuration;
import etomica.config.ConfigurationLattice;
import etomica.data.AccumulatorHistory;
import etomica.data.DataInfo;
import etomica.data.DataPump;
import etomica.data.meter.MeterEnergy;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.DisplayPhase;
import etomica.graphics.DisplayPlot;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorMC;
import etomica.integrator.IntervalActionAdapter;
import etomica.integrator.mcmove.MCMoveAtom;
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
 * @author ub2092
 *
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
public class MEAM_MC extends Simulation {
	
	public IntegratorMC integrator;
	public SpeciesSpheresMono sn;
    public SpeciesSpheresMono ag;
    public SpeciesSpheresMono cu;
	public Phase phase;
	public PotentialMEAM potentialN;
	public Controller controller;
	public DisplayPhase display;
	public DisplayPlot plot;
	public MeterEnergy energy;
	public ActivityIntegrate activityIntegrate;
	public DataInfo info2;

	public static void main(String[] args) {
	    MEAM_MC sim = new MEAM_MC();
	    MeterPotentialEnergy energyMeter = new MeterPotentialEnergy(sim.potentialMaster);
	    energyMeter.setPhase(sim.phase);
	    AccumulatorHistory energyAccumulator = new AccumulatorHistory(HistoryCollapsingAverage.FACTORY);
	    DisplayPlot plot = new DisplayPlot();
	    energyAccumulator.setDataSink(plot.getDataSet().makeDataSink());
	    DataPump energyManager = new DataPump(energyMeter,energyAccumulator);
	    //MeterKineticEnergy kineticMeter = new MeterKineticEnergy();
	    //kineticMeter.setPhase(sim.phase);
	    //AccumulatorHistory kineticAccumulator = new AccumulatorHistory(HistoryCollapsingAverage.FACTORY);
	    //DisplayPlot plotKE = new DisplayPlot();
	    //kineticAccumulator.setDataSink(plotKE.getDataSet().makeDataSink());
	    //DataPump kineticManager = new DataPump(kineticMeter, kineticAccumulator);
	    //energyAccumulator.setBlockSize(50);
	    IntervalActionAdapter adapter = new IntervalActionAdapter(energyManager, sim.integrator);
	    adapter.setActionInterval(1);
	    //IntervalActionAdapter kineticAdapter = new IntervalActionAdapter(kineticManager, sim.integrator);
	    //kineticAdapter.setActionInterval(1);
	    SimulationGraphic simgraphic = new SimulationGraphic(sim);
	    simgraphic.makeAndDisplayFrame();
	    simgraphic.panel().add(plot.graphic());
	    //simgraphic.panel().add(plotKE.graphic());
	    ColorSchemeByType colorScheme = ((ColorSchemeByType)((DisplayPhase)simgraphic.displayList().getFirst()).getColorScheme());
	    colorScheme.setColor(sim.sn.getMoleculeType(),java.awt.Color.blue);
    	colorScheme.setColor(sim.ag.getMoleculeType(),java.awt.Color.gray);
    	colorScheme.setColor(sim.cu.getMoleculeType(),java.awt.Color.orange);
	    //sim.activityIntegrate.setMaxSteps(1000);
	    //sim.getController().run();
	    //DataGroup data = (DataGroup)energyAccumulator.getData(); // kmb change type to Data instead of double[]
	    //double PE = ((DataDouble)data.getData(AccumulatorAverage.StatType.AVERAGE.index)).x
	    // /sim.species.getAgent(sim.phase).getNMolecules();  // kmb changed 8/3/05
	    //double PE = data[AccumulatorAverage.AVERAGE.index]/sim.species.getAgent(sim.phase).getNMolecules();  // orig line
	    //System.out.println("PE(eV)="+ElectronVolt.UNIT.fromSim(PE));
	    }
	    
	public MEAM_MC() {
	    super(Space3D.getInstance()); //INSTANCE); kmb change 8/3/05
	    integrator = new IntegratorMC(this);
	    integrator.getMoveManager().addMCMove(new MCMoveAtom(potentialMaster, 0.1, 0.2, false));
	    integrator.setTemperature(Kelvin.UNIT.toSim(298));
	    //integrator.setThermostatInterval(10);
	    integrator.setIsothermal(true);
	    activityIntegrate = new ActivityIntegrate(this,integrator);
	    activityIntegrate.setSleepPeriod(2);
	    getController().addAction(activityIntegrate);
	    sn = new SpeciesSpheresMono(this, Tin.INSTANCE);
        ag = new SpeciesSpheresMono(this, Silver.INSTANCE);
        cu = new SpeciesSpheresMono(this, Copper.INSTANCE);
	    
        /** The following values come from either the ASM Handbook or Cullity & Stock's 
         * "Elements of X-Ray Diffraction" (2001)
         */
        ((AtomTypeSphere)sn.getFactory().getType()).setDiameter(3.022); 
        
        ((AtomTypeSphere)ag.getFactory().getType()).setDiameter(2.8895); 
        
        ((AtomTypeSphere)cu.getFactory().getType()).setDiameter(2.5561); 
        
	    phase = new Phase(this);
        phase.getAgent(sn).setNMolecules(216);
        phase.getAgent(ag).setNMolecules(0);
        phase.getAgent(cu).setNMolecules(0);
	    
	    //beta-Sn phase
	    
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
	    
	    //FCC Cu
        /**
	    phase.setDimensions(new Vector3D(3.6148*3, 3.6148*3, 3.6148*6));
	    PrimitiveCubic primitive = new PrimitiveCubic(space, 3.6148);
	    LatticeCrystal crystal = new LatticeCrystal(new Crystal(
		        primitive, new BasisCubicFcc(primitive)));
	    **/
        
        //FCC Ag
        /**
	    phase.setDimensions(new Vector3D(4.0863*3, 4.0863*3, 4.0863*6));
	    PrimitiveCubic primitive = new PrimitiveCubic(space, 4.0863);
	    LatticeCrystal crystal = new LatticeCrystal(new Crystal(
		        primitive, new BasisCubicFcc(primitive)));
	    **/
	    
		Configuration config = new ConfigurationLattice(crystal);
		config.initializeCoordinates(phase);  
		
		potentialN = new PotentialMEAM(space);
		potentialN.setParameters(sn, ParameterSetMEAM.Sn);
		potentialN.setParameters(ag, ParameterSetMEAM.Ag);
		potentialN.setParameters(cu, ParameterSetMEAM.Cu);
		potentialN.setParametersIMC(cu, ParameterSetMEAM.Cu3Sn);
		potentialN.setParametersIMC(ag, ParameterSetMEAM.Ag3Sn);
        this.potentialMaster.addPotential(potentialN, new Species[]{sn, ag, cu}); 
	        
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
