package etomica.meam;
import etomica.action.activity.ActivityIntegrate;
import etomica.action.activity.Controller;
import etomica.atom.AtomTypeSphere;
import etomica.chem.elements.Copper;
import etomica.chem.elements.Silver;
import etomica.chem.elements.Tin;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorHistory;
import etomica.data.DataInfo;
import etomica.data.DataPump;
import etomica.data.AccumulatorAverage.StatType;
import etomica.data.meter.MeterEnergy;
import etomica.data.meter.MeterKineticEnergy;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.DisplayBox;
import etomica.graphics.DisplayPhase;
import etomica.graphics.DisplayPlot;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.integrator.IntervalActionAdapter;
import etomica.lattice.BravaisLatticeCrystal;
import etomica.lattice.crystal.BasisBetaSnA5;
import etomica.lattice.crystal.BasisCubicFcc;
import etomica.lattice.crystal.PrimitiveCubic;
import etomica.lattice.crystal.PrimitiveTetragonal;
import etomica.nbr.CriterionSimple;
import etomica.nbr.list.PotentialMasterList;
import etomica.phase.Phase;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularSlit;
import etomica.space3d.Space3D;
import etomica.space3d.Vector3D;
import etomica.species.Species;
import etomica.species.SpeciesSpheresMono;
import etomica.units.CompoundUnit;
import etomica.units.Joule;
import etomica.units.Kelvin;
import etomica.units.Mole;
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
 
public class MEAM_3DMDwithSnAgGB extends Simulation {
    
    private static final long serialVersionUID = 1L;
    public IntegratorVelocityVerlet integrator;
    public SpeciesSpheresMono snFixedA;
    public SpeciesSpheresMono snA;
    public SpeciesSpheresMono agA;
    public SpeciesSpheresMono cuA;
    public SpeciesSpheresMono agFixedB;
    public SpeciesSpheresMono snB;
    public SpeciesSpheresMono agB;
    public SpeciesSpheresMono cuB;
    public Phase phase;
    public PotentialMEAM potentialN;
    public Controller controller;
    public DisplayPhase display;
    public DisplayPlot plot;
    public MeterEnergy energy;
    public ActivityIntegrate activityIntegrate;
    public DataInfo info2;

    public static void main(String[] args) {
    	MEAM_3DMDwithSnAgGB sim = new MEAM_3DMDwithSnAgGB();
    	
    	MeterPotentialEnergy energyMeter = new MeterPotentialEnergy(sim.getPotentialMaster());
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
    	DataProcessorCvMD dataProcessorPE = new DataProcessorCvMD();
    	dataProcessorPE.setIntegrator(sim.integrator);
    	
    	//Heat Capacity (KE)
    	DataProcessorCvMD dataProcessorKE = new DataProcessorCvMD();
    	dataProcessorKE.setIntegrator(sim.integrator);
    	
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
    	colorScheme.setColor(sim.snFixedA.getMoleculeType(),java.awt.Color.white);
    	colorScheme.setColor(sim.snA.getMoleculeType(),java.awt.Color.white);
    	colorScheme.setColor(sim.agA.getMoleculeType(),java.awt.Color.gray);
    	colorScheme.setColor(sim.cuA.getMoleculeType(),java.awt.Color.orange);
    	colorScheme.setColor(sim.agFixedB.getMoleculeType(),java.awt.Color.gray);
    	colorScheme.setColor(sim.snB.getMoleculeType(),java.awt.Color.white);
    	colorScheme.setColor(sim.agB.getMoleculeType(),java.awt.Color.gray);
    	colorScheme.setColor(sim.cuB.getMoleculeType(),java.awt.Color.orange);
    	//sim.activityIntegrate.setMaxSteps(1000);
    	//sim.getController().run();
    	//DataGroup data = (DataGroup)energyAccumulator.getData(); // kmb change type to Data instead of double[]
        //double PE = ((DataDouble)data.getData(AccumulatorAverage.StatType.AVERAGE.index)).x
                    // /sim.species.getAgent(sim.phase).getNMolecules();  // kmb changed 8/3/05
        //double PE = data[AccumulatorAverage.AVERAGE.index]/sim.species.getAgent(sim.phase).getNMolecules();  // orig line
        //System.out.println("PE(eV)="+ElectronVolt.UNIT.fromSim(PE));
    }
    
    public MEAM_3DMDwithSnAgGB() {
        super(Space3D.getInstance(), true, new PotentialMasterList(Space3D.getInstance())); //INSTANCE); kmb change 8/3/05
        integrator = new IntegratorVelocityVerlet(this);
        integrator.setTimeStep(0.001);
        integrator.setTemperature(Kelvin.UNIT.toSim(295));
        integrator.setThermostatInterval(100);
        integrator.setIsothermal(true);
        activityIntegrate = new ActivityIntegrate(this,integrator);
        activityIntegrate.setSleepPeriod(2);
        getController().addAction(activityIntegrate);
        Tin SnF = new Tin("SnF", Double.POSITIVE_INFINITY);
        snFixedA = new SpeciesSpheresMono(this, SnF);
        snA = new SpeciesSpheresMono(this, Tin.INSTANCE);
        agA = new SpeciesSpheresMono(this, Silver.INSTANCE);
        cuA = new SpeciesSpheresMono(this, Copper.INSTANCE);
        agFixedB = new SpeciesSpheresMono(this, SnF);
        snB = new SpeciesSpheresMono(this, Tin.INSTANCE);
        agB = new SpeciesSpheresMono(this, Silver.INSTANCE);
        cuB = new SpeciesSpheresMono(this, Copper.INSTANCE);
        
        
        double aA, bA, cA, aB, bB, cB;
        
        int nCellsAx, nCellsAy, nCellsAz, nAMobile, nAFixed, basisA,
		nCellsBx, nCellsBy, nCellsBz, nBMobile, nBFixed, basisB,
		nA, nB, nAImpurity, nBImpurity, nAVacancy, nBVacancy;
        
        nCellsAx = 5; nCellsAy = 5; nCellsAz = 4;
        nCellsBx = 7; nCellsBy = 7; nCellsBz = 4;
        nAImpurity = 0; nAVacancy = 0;
        nBImpurity = 0; nBVacancy = 0;
        
        phase = new Phase(this);
        phase.setBoundary(new BoundaryRectangularSlit(this, 2));
        
        // beta-Sn phase
        
        //The dimensions of the simulation box must be proportional to those of
        //the unit cell to prevent distortion of the lattice.  The values for the 
        //lattice parameters for tin's beta phase (a = 5.8314 angstroms, c = 3.1815 
        //angstroms) are taken from the ASM Handbook. 
        aA = bA = 5.796; cA = 3.1815; basisA = 4;
        PrimitiveTetragonal primitiveA = new PrimitiveTetragonal(space, aA, cA);
        //Alternatively, using the parameters calculated in Ravelo & Baskes (1997)
        //phase.setDimensions(new Vector3D(5.92*3, 5.92*3, 3.23*6));
        //PrimitiveTetragonal primitive = new PrimitiveTetragonal(space, 5.92, 3.23);
        BravaisLatticeCrystal latticeA = new BravaisLatticeCrystal(primitiveA, new BasisBetaSnA5());
        
        aB = bB = cB = 4.14; basisB = 4;
	    PrimitiveCubic primitiveB = new PrimitiveCubic(space, aB);
	    BravaisLatticeCrystal latticeB = new BravaisLatticeCrystal(primitiveB, new BasisCubicFcc());
        
        phase.setDimensions(new Vector3D(aA*nCellsAx, aA*nCellsAy, (cA*nCellsAz)+(cB*nCellsBz)));
		
	    nA = (nCellsAx * nCellsAy * nCellsAz) * basisA;
	    nAFixed = (nCellsAx * nCellsAy * 2) * basisA;
	    nAMobile = nA - nAFixed - nAImpurity - nAVacancy;
	    nB = (nCellsBx * nCellsBy * nCellsBz) * basisB;
	    nBFixed = (nCellsBx * nCellsBy * 2) * basisB;
	    nBMobile = nB - nBFixed - nBImpurity - nBVacancy;
	        
	    snFixedA.getAgent(phase).setNMolecules(nAFixed);
	    snA.getAgent(phase).setNMolecules(nAMobile);
	    agA.getAgent(phase).setNMolecules(0);
	    cuA.getAgent(phase).setNMolecules(0);
	    agFixedB.getAgent(phase).setNMolecules(nBFixed);
	    snB.getAgent(phase).setNMolecules(0);
	    agB.getAgent(phase).setNMolecules(nBMobile);
	    cuB.getAgent(phase).setNMolecules(0);
	    
	        
	    /** The following values come from either the ASM Handbook or Cullity & Stock's 
	     * "Elements of X-Ray Diffraction" (2001)
	     */
	    
	    ((AtomTypeSphere)snFixedA.getFactory().getType()).setDiameter(3.022); 
	    
	    ((AtomTypeSphere)snA.getFactory().getType()).setDiameter(3.022); 
	        
	    ((AtomTypeSphere)agA.getFactory().getType()).setDiameter(2.8895); 
	        
	    ((AtomTypeSphere)cuA.getFactory().getType()).setDiameter(2.5561); 
	    
	    ((AtomTypeSphere)agFixedB.getFactory().getType()).setDiameter(2.8895); 
	        
	    ((AtomTypeSphere)snB.getFactory().getType()).setDiameter(3.022); 
	        
	    ((AtomTypeSphere)agB.getFactory().getType()).setDiameter(2.8895); 
	        
	    ((AtomTypeSphere)cuB.getFactory().getType()).setDiameter(2.5561); 
	     
	    GrainBoundaryConfiguration config = new GrainBoundaryConfiguration(latticeA, latticeB);
	    config.setDimensions(nCellsAx, nCellsAy, nCellsAz, nCellsBx, nCellsBy, 
	    		nCellsBz, aA, bA, cA, aB, bB, cB);
	    config.initializeCoordinates(phase);
        
	    System.out.println("In simulation class  " + snFixedA.getAgent(phase).getNode().firstLeafAtom().getCoord().getPosition());
		potentialN = new PotentialMEAM(space);
		potentialN.setParameters(snFixedA, ParameterSetMEAM.Sn);
		potentialN.setParameters(snA, ParameterSetMEAM.Sn);
		potentialN.setParameters(agA, ParameterSetMEAM.Ag);
		potentialN.setParameters(cuA, ParameterSetMEAM.Cu);
		potentialN.setParameters(agFixedB, ParameterSetMEAM.Ag);
		potentialN.setParameters(snB, ParameterSetMEAM.Sn);
		potentialN.setParameters(agB, ParameterSetMEAM.Ag);
		potentialN.setParameters(cuB, ParameterSetMEAM.Cu);
		potentialN.setParametersIMC(cuA, ParameterSetMEAM.Cu3Sn);
		potentialN.setParametersIMC(agA, ParameterSetMEAM.Ag3Sn);
		potentialN.setParametersIMC(cuB, ParameterSetMEAM.Cu3Sn);
		potentialN.setParametersIMC(agB, ParameterSetMEAM.Ag3Sn);
		potentialN.setParametersIMC(agFixedB, ParameterSetMEAM.Ag3Sn);
        this.potentialMaster.addPotential(potentialN, new Species[]{snFixedA, snA, agA, cuA, agFixedB, snB, agB, cuB});    
        ((PotentialMasterList)potentialMaster).setRange(potentialN.getRange()*1.1);
        ((PotentialMasterList)potentialMaster).setCriterion(potentialN, new CriterionSimple(this, potentialN.getRange(), potentialN.getRange()*1.1));
        integrator.addListener(((PotentialMasterList)potentialMaster).getNeighborManager());
        
        integrator.setPhase(phase);
		
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