package etomica.meam;

import etomica.action.activity.ActivityIntegrate;
import etomica.action.activity.Controller;
import etomica.atom.AtomTypeSphere;
import etomica.chem.elements.Tin;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorHistory;
import etomica.data.DataPump;
import etomica.data.IDataInfo;
import etomica.data.AccumulatorAverage.StatType;
import etomica.data.meter.MeterEnergy;
import etomica.data.meter.MeterKineticEnergy;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.DisplayBox;
import etomica.graphics.DisplayPhase;
import etomica.graphics.DisplayPlot;
import etomica.graphics.SimulationGraphic;
import etomica.graphics.SimulationPanel;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.integrator.IntervalActionAdapter;
import etomica.lattice.BravaisLatticeCrystal;
import etomica.lattice.crystal.BasisBetaSnA5;
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
 
public class MEAM_3DMDwithGB extends Simulation {
    
    private static final long serialVersionUID = 1L;
    private static final String APP_NAME = "MEAM 3DMD w/GB";
    public PotentialMasterList potentialMaster;
    public IntegratorVelocityVerlet integrator;
    public SpeciesSpheresMono snFixedA;
    public SpeciesSpheresMono snA;
//    public SpeciesSpheresMono agA;
//    public SpeciesSpheresMono cuA;
    public SpeciesSpheresMono snFixedB;
    public SpeciesSpheresMono snB;
//    public SpeciesSpheresMono agB;
//    public SpeciesSpheresMono cuB;
    public Phase phase;
    public PotentialMEAM potentialN;
    public Controller controller;
    public DisplayPlot plot;
    public MeterEnergy energy;
    public ActivityIntegrate activityIntegrate;
    public IDataInfo info2;

    public static void main(String[] args) {
    	MEAM_3DMDwithGB sim = new MEAM_3DMDwithGB();
    	
    	MeterPotentialEnergy energyMeter = new MeterPotentialEnergy(sim.potentialMaster);
    	MeterKineticEnergy kineticMeter = new MeterKineticEnergy();
    	
    	energyMeter.setPhase(sim.phase);
    	kineticMeter.setPhase(sim.phase);
        
        AccumulatorHistory energyAccumulator = new AccumulatorHistory(new HistoryCollapsingAverage());
        AccumulatorHistory kineticAccumulator = new AccumulatorHistory(new HistoryCollapsingAverage());
        
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

        SimulationGraphic simgraphic = new SimulationGraphic(sim, SimulationGraphic.TABBED_PANE, APP_NAME);
        
    	DisplayBox cvBoxPE = new DisplayBox();
    	dataProcessorPE.setDataSink(cvBoxPE);
    	cvBoxPE.setUnit(new CompoundUnit(new Unit[]{Joule.UNIT, Kelvin.UNIT, Mole.UNIT}, new double []{1,-1,-1}));
    	cvBoxPE.setLabel("PE Cv contrib.");
    	DisplayBox cvBoxKE = new DisplayBox();
    	dataProcessorKE.setDataSink(cvBoxKE);
    	cvBoxKE.setUnit(new CompoundUnit(new Unit[]{Joule.UNIT, Kelvin.UNIT, Mole.UNIT}, new double []{1,-1,-1}));
    	cvBoxKE.setLabel("KE Cv contrib.");

    	simgraphic.getPanel().plotPanel.add(plot.graphic(), SimulationPanel.getVertGBC());
    	//simgraphic.panel().add(plotKE.graphic());
    	simgraphic.getPanel().plotPanel.add(cvBoxPE.graphic(), SimulationPanel.getVertGBC());
    	simgraphic.getPanel().plotPanel.add(cvBoxKE.graphic(), SimulationPanel.getVertGBC());

    	simgraphic.getController().getReinitButton().setPostAction(simgraphic.getDisplayPhasePaintAction(sim.phase));

    	ColorSchemeByType colorScheme = ((ColorSchemeByType)((DisplayPhase)simgraphic.displayList().getFirst()).getColorScheme());
    	colorScheme.setColor(sim.snFixedA.getMoleculeType(),java.awt.Color.blue);
    	colorScheme.setColor(sim.snA.getMoleculeType(),java.awt.Color.red);
//    	colorScheme.setColor(sim.agA.getMoleculeType(),java.awt.Color.gray);
//    	colorScheme.setColor(sim.cuA.getMoleculeType(),java.awt.Color.orange);
    	colorScheme.setColor(sim.snFixedB.getMoleculeType(),java.awt.Color.yellow);
    	colorScheme.setColor(sim.snB.getMoleculeType(),java.awt.Color.green);
//    	colorScheme.setColor(sim.agB.getMoleculeType(),java.awt.Color.gray);
//    	colorScheme.setColor(sim.cuB.getMoleculeType(),java.awt.Color.orange);

    	simgraphic.makeAndDisplayFrame(APP_NAME);

    	//sim.activityIntegrate.setMaxSteps(1000);
    	//sim.getController().run();
    	//DataGroup data = (DataGroup)energyAccumulator.getData(); // kmb change type to Data instead of double[]
        //double PE = ((DataDouble)data.getData(AccumulatorAverage.StatType.AVERAGE.index)).x
                    // /sim.species.getAgent(sim.phase).getNMolecules();  // kmb changed 8/3/05
        //double PE = data[AccumulatorAverage.AVERAGE.index]/sim.species.getAgent(sim.phase).getNMolecules();  // orig line
        //System.out.println("PE(eV)="+ElectronVolt.UNIT.fromSim(PE));
    }
    
    public MEAM_3DMDwithGB() {
        super(Space3D.getInstance(), true); //INSTANCE); kmb change 8/3/05
        potentialMaster = new PotentialMasterList(this);
        integrator = new IntegratorVelocityVerlet(this, potentialMaster);
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
//        agA = new SpeciesSpheresMono(this, Silver.INSTANCE);
//        cuA = new SpeciesSpheresMono(this, Copper.INSTANCE);
        snFixedB = new SpeciesSpheresMono(this, SnF);
        snB = new SpeciesSpheresMono(this, Tin.INSTANCE);
//        agB = new SpeciesSpheresMono(this, Silver.INSTANCE);
//        cuB = new SpeciesSpheresMono(this, Copper.INSTANCE);
        
        
        double aA, bA, cA, aB, bB, cB;
        
        int nCellsAx, nCellsAy, nCellsAz, nAMobile, nAFixed,
		nCellsBx, nCellsBy, nCellsBz, nBMobile, nBFixed,
		nA, nB, nAImpurity, nBImpurity, nAVacancy, nBVacancy;
        
        nCellsAx = 6; nCellsAy = 6; nCellsAz = 4;
        nCellsBx = 6; nCellsBy = 6; nCellsBz = 4;
        nAImpurity = 0; nAVacancy = 0;
        nBImpurity = 0; nBVacancy = 0;
        
        phase = new Phase(new BoundaryRectangularSlit(this, 2));
        addPhase(phase);
        
        // beta-Sn phase
        
        //The dimensions of the simulation box must be proportional to those of
        //the unit cell to prevent distortion of the lattice.  The values for the 
        //lattice parameters for tin's beta phase (a = 5.8314 angstroms, c = 3.1815 
        //angstroms) are taken from the ASM Handbook. 
        aA = bA = 5.8314; cA = 3.1815; int basisA = 4;
        PrimitiveTetragonal primitiveA = new PrimitiveTetragonal(space, aA, cA);
        //Alternatively, using the parameters calculated in Ravelo & Baskes (1997)
        //phase.setDimensions(new Vector3D(5.92*3, 5.92*3, 3.23*6));
        //PrimitiveTetragonal primitive = new PrimitiveTetragonal(space, 5.92, 3.23);
        BravaisLatticeCrystal latticeA = new BravaisLatticeCrystal(primitiveA, new BasisBetaSnA5());
        
        aB = bB = 5.8314; cB = 3.1815; int basisB = 4;
        PrimitiveTetragonal primitiveB = new PrimitiveTetragonal(space, aB, cB);
        BravaisLatticeCrystal latticeB = new BravaisLatticeCrystal(primitiveB, new BasisBetaSnA5());
        
        phase.setDimensions(new Vector3D(aA*nCellsAx, aA*nCellsAy, (cA*nCellsAz)+(cB*nCellsBz)));
        
        
        //FCC Cu
        /**
        aA = bA = cA = 3.6148;
	    phaseA.setDimensions(new Vector3D(aA*4, aA*4, aA*4));
	    PrimitiveCubic primitiveA = new PrimitiveCubic(space, aA);
	    LatticeCrystal latticeA = new LatticeCrystal(new Crystal(
		        primitiveA, new BasisCubicFcc(primitiveA)));
		        
		aB = bB = cB = 3.6148;
	    phaseB.setDimensions(new Vector3D(aB*4, aB*4, aB*4));
	    PrimitiveCubic primitiveB = new PrimitiveCubic(space, aB);
	    LatticeCrystal latticeB = new LatticeCrystal(new Crystal(
		        primitiveB, new BasisCubicFcc(primitiveB)));        
		
	    */
        
        //FCC Ag
        /**
        aA = bA = cA = 4.0863;
	    phaseA.setDimensions(new Vector3D(aA*4, aA*4, aA*4));
	    PrimitiveCubic primitiveA = new PrimitiveCubic(space, aA);
	    LatticeCrystal latticeA = new LatticeCrystal(new Crystal(
		        primitiveA, new BasisCubicFcc(primitiveA)));
		        
		aB = bB = cB = 4.0863;
	    phaseB.setDimensions(new Vector3D(aB*4, aB*4, aB*4));
	    PrimitiveCubic primitiveB = new PrimitiveCubic(space, aB);
	    LatticeCrystal latticeB = new LatticeCrystal(new Crystal(
		        primitiveB, new BasisCubicFcc(primitiveB)));
	    */
		
	    nA = (nCellsAx * nCellsAy * nCellsAz) * basisA;
	    nAFixed = (nCellsAx * nCellsAy * 2) * basisA;
	    nAMobile = nA - nAFixed - nAImpurity - nAVacancy;
	    nB = (nCellsBx * nCellsBy * nCellsBz) * basisB;
	    nBFixed = (nCellsBx * nCellsBy * 2) * basisB;
	    nBMobile = nB - nBFixed - nBImpurity - nBVacancy;

        getSpeciesManager().addSpecies(snFixedA);
        getSpeciesManager().addSpecies(snA);
//        getSpeciesRoot().addSpecies(agA);
//        getSpeciesRoot().addSpecies(cuA);
        getSpeciesManager().addSpecies(snFixedB);
        getSpeciesManager().addSpecies(snB);
//        getSpeciesRoot().addSpecies(agB);
//        getSpeciesRoot().addSpecies(cuB);
        
	    snFixedA.getAgent(phase).setNMolecules(nAFixed);
	    snA.getAgent(phase).setNMolecules(nAMobile);
//	    agA.getAgent(phase).setNMolecules(0);
//	    cuA.getAgent(phase).setNMolecules(0);
	    snFixedB.getAgent(phase).setNMolecules(nBFixed);
	    snB.getAgent(phase).setNMolecules(nBMobile);
//	    agB.getAgent(phase).setNMolecules(0);
//	    cuB.getAgent(phase).setNMolecules(0);
	    
	        
	    /** The following values come from either the ASM Handbook or Cullity & Stock's 
	     * "Elements of X-Ray Diffraction" (2001)
	     */
	    
	    ((AtomTypeSphere)snFixedA.getFactory().getType()).setDiameter(3.022); 
	    
	    ((AtomTypeSphere)snA.getFactory().getType()).setDiameter(3.022); 
	        
//	    ((AtomTypeSphere)agA.getFactory().getType()).setDiameter(2.8895); 
//	        
//	    ((AtomTypeSphere)cuA.getFactory().getType()).setDiameter(2.5561); 
//	    
	    ((AtomTypeSphere)snFixedB.getFactory().getType()).setDiameter(3.022); 
	    
	    ((AtomTypeSphere)snB.getFactory().getType()).setDiameter(3.022); 
	        
//	    ((AtomTypeSphere)agB.getFactory().getType()).setDiameter(2.8895); 
//	        
//	    ((AtomTypeSphere)cuB.getFactory().getType()).setDiameter(2.5561); 
//	     
	    GrainBoundaryConfiguration config = new GrainBoundaryConfiguration(latticeA, latticeB);
	    config.setDimensions(nCellsAx, nCellsAy, nCellsAz, nCellsBx, nCellsBy, 
	    		nCellsBz, aA, bA, cA, aB, bB, cB);
	    config.initializeCoordinates(phase);
        
		potentialN = new PotentialMEAM(space);
		potentialN.setParameters(snFixedA, ParameterSetMEAM.Sn);
		potentialN.setParameters(snA, ParameterSetMEAM.Sn);
//		potentialN.setParameters(agA, ParameterSetMEAM.Ag);
//		potentialN.setParameters(cuA, ParameterSetMEAM.Cu);
		potentialN.setParameters(snFixedB, ParameterSetMEAM.Sn);
		potentialN.setParameters(snB, ParameterSetMEAM.Sn);
//		potentialN.setParameters(agB, ParameterSetMEAM.Ag);
//		potentialN.setParameters(cuB, ParameterSetMEAM.Cu);
//		potentialN.setParametersIMC(cuA, ParameterSetMEAM.Cu3Sn);
//		potentialN.setParametersIMC(agA, ParameterSetMEAM.Ag3Sn);
//		potentialN.setParametersIMC(cuB, ParameterSetMEAM.Cu3Sn);
//		potentialN.setParametersIMC(agB, ParameterSetMEAM.Ag3Sn);
//        this.potentialMaster.addPotential(potentialN, new Species[]{snFixedA, snA, agA, cuA, snFixedB, snB, agB, cuB});    
        this.potentialMaster.addPotential(potentialN, new Species[]{snFixedA, snA, snFixedB, snB});    
        potentialMaster.setRange(potentialN.getRange()*1.1);
        potentialMaster.setCriterion(potentialN, new CriterionSimple(this, potentialN.getRange(), potentialN.getRange()*1.1));
        integrator.addListener(potentialMaster.getNeighborManager(phase));
        
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