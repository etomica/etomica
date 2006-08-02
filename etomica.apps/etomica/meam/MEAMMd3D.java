package etomica.meam;
import etomica.action.PhaseImposePbc;
import etomica.action.activity.ActivityIntegrate;
import etomica.action.activity.Controller;
import etomica.atom.AtomTypeLeaf;
import etomica.atom.AtomTypeSphere;
import etomica.config.Configuration;
import etomica.config.ConfigurationLattice;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorHistory;
import etomica.data.Data;
import etomica.data.DataInfo;
import etomica.data.DataProcessor;
import etomica.data.DataPump;
import etomica.data.DataTag;
import etomica.data.AccumulatorAverage.StatType;
import etomica.data.meter.MeterEnergy;
import etomica.data.meter.MeterKineticEnergy;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.types.DataDouble;
import etomica.data.types.DataDouble.DataInfoDouble;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.DisplayBox;
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
import etomica.units.CompoundDimension;
import etomica.units.CompoundUnit;
import etomica.units.Dimension;
import etomica.units.Energy;
import etomica.units.Joule;
import etomica.units.Kelvin;
import etomica.units.Mole;
import etomica.units.Quantity;
import etomica.units.Temperature;
import etomica.units.Unit;
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
    	MEAMMd3D sim = new MEAMMd3D();
    	
    	MeterPotentialEnergy energyMeter = new MeterPotentialEnergy(sim.potentialMaster);
    	MeterKineticEnergy kineticMeter = new MeterKineticEnergy();
    	
    	energyMeter.setPhase(sim.phase);
    	kineticMeter.setPhase(sim.phase);
        
        AccumulatorHistory energyAccumulator = new AccumulatorHistory(HistoryCollapsingAverage.FACTORY);
        AccumulatorHistory kineticAccumulator = new AccumulatorHistory(HistoryCollapsingAverage.FACTORY);
        
        AccumulatorAverage accumulatorAveragePE = new AccumulatorAverage(50);
    	AccumulatorAverage accumulatorAverageKE = new AccumulatorAverage(50);
    	
    	DataPump energyManager = new DataPump(energyMeter,accumulatorAveragePE);   	
    	DataPump kineticManager = new DataPump(kineticMeter, accumulatorAverageKE);
    	
    	accumulatorAveragePE.addDataSink(energyAccumulator, new StatType[]{StatType.MOST_RECENT});
    	accumulatorAverageKE.addDataSink(kineticAccumulator, new StatType[]{StatType.MOST_RECENT});
    	
    	DisplayPlot plot = new DisplayPlot();
        //DisplayPlot plotKE = new DisplayPlot();
    	
        energyAccumulator.setDataSink(plot.getDataSet().makeDataSink());
        kineticAccumulator.setDataSink(plot.getDataSet().makeDataSink());
        //energyAccumulator.setBlockSize(50);
        
    	accumulatorAveragePE.setPushInterval(1);
    	accumulatorAverageKE.setPushInterval(1);
    	
    	sim.register(energyMeter, energyManager);
    	sim.register(kineticMeter, kineticManager);
    	
    	//Heat Capacity (PE)
    	DataProcessor dataProcessorPE = new DataProcessor(){
    		MEAMMd3D sim = new MEAMMd3D();
    		DataDouble data;
    		DataInfoDouble dataInfo;
    		DataTag tag = new DataTag();
    		
    		public Data processData(Data inputData){
    			data.x = ((DataDouble)inputData).x;
    			double systemTemp = sim.integrator.getTemperature();
    			data.x /= systemTemp;
    			data.x *= data.x/sim.phase.moleculeCount();
    			return data;
    		}
    		
    		public DataProcessor getDataCaster(DataInfo inputDataInfo){
    			return null;
    		}
    		
    		public DataInfo processDataInfo(DataInfo inputDataInfo){
    			Dimension cvDimension = new CompoundDimension(new Dimension[]{Energy.DIMENSION, Temperature.DIMENSION, Quantity.DIMENSION}, new double[]{1,-1,-1});
    			data = new DataDouble();
    			dataInfo = new DataInfoDouble("Heat Capacity", cvDimension);
    			return dataInfo;
    		}
    		
    		public DataTag getTag(){
    			return tag;
    		}
    	};
    	
//    	Heat Capacity (KE)
    	DataProcessor dataProcessorKE = new DataProcessor(){
    		MEAMMd3D sim = new MEAMMd3D();
    		DataDouble data;
    		DataInfoDouble dataInfo;
    		DataTag tag = new DataTag();
    		
    		public Data processData(Data inputData){
    			data.x = ((DataDouble)inputData).x;
    			double systemTemp = sim.integrator.getTemperature();
    			data.x /= systemTemp;
    			data.x *= data.x/sim.phase.moleculeCount();
    			return data;
    		}
    		
    		public DataProcessor getDataCaster(DataInfo inputDataInfo){
    			return null;
    		}
    		
    		public DataInfo processDataInfo(DataInfo inputDataInfo){
    			Dimension cvDimension = new CompoundDimension(new Dimension[]{Energy.DIMENSION, Temperature.DIMENSION, Quantity.DIMENSION}, new double[]{1,-1,-1});
    			data = new DataDouble();
    			dataInfo = new DataInfoDouble("Heat Capacity", cvDimension);
    			return dataInfo;
    		}
    		
    		public DataTag getTag(){
    			return tag;
    		}
    	};
    	
    	accumulatorAveragePE.addDataSink(dataProcessorPE, new StatType[]{StatType.STANDARD_DEVIATION});
    	accumulatorAverageKE.addDataSink(dataProcessorKE, new StatType[]{StatType.STANDARD_DEVIATION});
    	  	
        IntervalActionAdapter adapter = new IntervalActionAdapter(energyManager, sim.integrator);
        adapter.setActionInterval(1);
        IntervalActionAdapter kineticAdapter = new IntervalActionAdapter(kineticManager, sim.integrator);
        kineticAdapter.setActionInterval(1);    

        SimulationGraphic simgraphic = new SimulationGraphic(sim);
        
    	DisplayBox cvBoxPE = new DisplayBox();
    	dataProcessorPE.setDataSink(cvBoxPE);
    	cvBoxPE.setUnit(new CompoundUnit(new Unit[]{Joule.UNIT, Kelvin.UNIT, Mole.UNIT}, new double []{1,-1,-1}));
    	cvBoxPE.setLabel("PE Cv contrib.");
    	DisplayBox cvBoxKE = new DisplayBox();
    	dataProcessorKE.setDataSink(cvBoxKE);
    	cvBoxKE.setUnit(new CompoundUnit(new Unit[]{Joule.UNIT, Kelvin.UNIT, Mole.UNIT}, new double []{1,-1,-1}));
    	cvBoxKE.setLabel("KE Cv contrib.");
    	
    	simgraphic.makeAndDisplayFrame();
    	simgraphic.panel().add(plot.graphic());
    	//simgraphic.panel().add(plotKE.graphic());
    	simgraphic.panel().add(cvBoxPE.graphic());
    	simgraphic.panel().add(cvBoxKE.graphic());
    	
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
    
    public MEAMMd3D() {
        super(Space3D.getInstance()); //INSTANCE); kmb change 8/3/05
        integrator = new IntegratorVelocityVerlet(this);
        integrator.setTimeStep(0.001);
        integrator.setTemperature(Kelvin.UNIT.toSim(298));
        integrator.setThermostatInterval(100);
        integrator.setIsothermal(true);
        activityIntegrate = new ActivityIntegrate(this,integrator);
        activityIntegrate.setSleepPeriod(2);
        getController().addAction(activityIntegrate);
        sn = new SpeciesSpheresMono(this);
        ag = new SpeciesSpheresMono(this);
        cu = new SpeciesSpheresMono(this);
        sn.setNMolecules(216);
        ag.setNMolecules(0);
        cu.setNMolecules(0);
        
        /** The following values come from either the ASM Handbook or Cullity & Stock's 
         * "Elements of X-Ray Diffraction" (2001)
         */
        ((AtomTypeLeaf)sn.getFactory().getType()).setMass(118.69);
        ((AtomTypeSphere)sn.getFactory().getType()).setDiameter(3.022); 
        
        ((AtomTypeLeaf)ag.getFactory().getType()).setMass(107.868);
        ((AtomTypeSphere)ag.getFactory().getType()).setDiameter(2.8895); 
        
        ((AtomTypeLeaf)cu.getFactory().getType()).setMass(63.546);
        ((AtomTypeSphere)cu.getFactory().getType()).setDiameter(2.5561); 
        
        
        phase = new Phase(this);
        
        // beta-Sn phase
        
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
	    */
        
        //FCC Ag
        /**
	    phase.setDimensions(new Vector3D(4.0863*3, 4.0863*3, 4.0863*6));
	    PrimitiveCubic primitive = new PrimitiveCubic(space, 4.0863);
	    LatticeCrystal crystal = new LatticeCrystal(new Crystal(
		        primitive, new BasisCubicFcc(primitive)));
	    */
           
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