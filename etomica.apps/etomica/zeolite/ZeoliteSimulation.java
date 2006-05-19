package etomica.zeolite;

import etomica.action.PhaseImposePbc;
import etomica.action.activity.ActivityIntegrate;
import etomica.atom.Atom;
import etomica.atom.AtomArrayList;
import etomica.atom.AtomLeaf;
import etomica.atom.AtomTreeNodeGroup;
import etomica.atom.AtomType;
import etomica.atom.AtomTypeLeaf;
import etomica.atom.iterator.AtomIteratorTree;
import etomica.config.ConfigurationLattice;
import etomica.config.ConfigurationSequential;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.DeviceNSelector;
import etomica.graphics.DisplayBox;
import etomica.graphics.DisplayBoxesCAE;
import etomica.graphics.DisplayPhase;
import etomica.graphics.SimulationGraphic;
import etomica.graphics.DisplayPlot;
import etomica.graphics.DisplayTable;
import etomica.integrator.IntegratorHard;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.integrator.IntervalActionAdapter;
import etomica.lattice.LatticeCubicFcc;
import etomica.nbr.CriterionSimple;
import etomica.nbr.CriterionSpecies;
import etomica.nbr.NeighborCriterion;
import etomica.nbr.list.NeighborListManager;
import etomica.nbr.list.PotentialMasterList;
import etomica.phase.Phase;
import etomica.potential.P2HardSphere;
import etomica.potential.P1HardBoundary;
import etomica.potential.P2SoftSphericalTruncated;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.Species;
import etomica.species.SpeciesSpheresMono;
import etomica.util.Default;
import etomica.util.HistoryCollapsing;
import etomica.config.ConfigurationFile;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorHistory;
import etomica.data.DataFork;
import etomica.data.DataPump;
import etomica.data.DataSink;
import etomica.data.meter.MeterEnergy;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.meter.MeterKineticEnergy;
import etomica.data.meter.MeterTemperature;
import etomica.potential.P2LennardJones;
import etomica.potential.P2WCA;
import etomica.units.Kelvin;
import etomica.units.Pixel;
import etomica.units.Energy;
import etomica.units.Bar;
import etomica.space.BoundaryRectangularSlit;
import etomica.space.BoundaryRectangularPeriodic;

import java.awt.ActiveEvent;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;


/**
 * 
 * Three-dimensional hard-sphere molecular dynamics simulation, using
 * neighbor listing.  
 * <p>
 * Developed as a prototype and example for the construction of a basic simulation.
 *
 * @author David Kofke and Andrew Schultz
 *
 */
public class ZeoliteSimulation extends Simulation {

    //the following fields are made accessible for convenience to permit simple
    //mutation of the default behavior

    /**
     * The Phase holding the atoms. 
     */
    public final Phase phase;
    /**
     * The Integrator performing the dynamics.
     */
    //public final IntegratorHard integrator;
    public final IntegratorVelocityVerlet integrator;
    /**
     * The single hard-sphere species.
     */
    public final SpeciesSpheresMono[] species;
    /**
     * The hard-sphere potential governing the interactions.
     */
    //public final P2HardSphere potential;
    public final P2LennardJones potentialMM;
    private int nAtomsMeth;
    private AccumulatorAverage temperatureAvg;
    private AccumulatorAverage keAvg;
    private int interval;
    public ActivityIntegrate activityIntegrate;
    public DisplayPlot ePlot;
    /**
     * Sole public constructor, makes a simulation using a 3D space.
     */
    public ZeoliteSimulation() {
        this(new Default());
    }
    
    public ZeoliteSimulation(Default defaults) {
        this(Space3D.getInstance(), defaults);
    }
    
    //we use a second, private constructor to permit the space to
    //appear twice in the call to the superclass constructor; alternatively
    //we could have passed Space3D.getInstance() twice
    private ZeoliteSimulation(Space space, Default defaults) {

        // invoke the superclass constructor
        // "true" is indicating to the superclass that this is a dynamic simulation
        // the PotentialMaster is selected such as to implement neighbor listing
        super(space, true, new PotentialMasterList(space, 1.6), Default.BIT_LENGTH, defaults);

        //Additions for Zeolite Calculations
        //Start by reading the first line, which is number of Atoms
        String fileName = "2unitcell";
        //String fileName = "pbu2";
        ConfigurationFileXYZ config = new ConfigurationFileXYZ(this.space,fileName);
        int[] numAtoms = config.getNumAtoms();
        
        nAtomsMeth = numAtoms[numAtoms.length - 1];
        double neighborRangeFac = 1.2;
        //defaults.makeLJDefaults();
        defaults.timeStep = 0.01;
        defaults.atomMass = 16;
        defaults.pixelUnit = new Pixel(10);
        defaults.temperature = Kelvin.UNIT.toSim(298.0);
        defaults.pressure = Bar.UNIT.toSim(1.01325);
        //Setting sizes of molecules
        double[] atomicSize = new double[3];
        atomicSize[0] = 0.73;
        atomicSize[1] = 1.18;
        atomicSize[2] = 2.088;
        
        double range = 8.035;
        ((PotentialMasterList)potentialMaster).setRange(3.214*neighborRangeFac*2.5);
        
        
        integrator = new IntegratorVelocityVerlet(this);
        integrator.setIsothermal(true);
        integrator.setThermostatInterval(10);
        integrator.setTimeStep(0.00611);
        this.register(integrator);
        
        
        
        
        NeighborListManager nbrManager = ((PotentialMasterList)potentialMaster).getNeighborManager();
        //nbrManager.setRange(range*1.6);
        nbrManager.getPbcEnforcer().setApplyToMolecules(false);
        integrator.addListener(nbrManager);

        activityIntegrate = new ActivityIntegrate(this,integrator);
        activityIntegrate.setDoSleep(true);
        activityIntegrate.setSleepPeriod(2);
        activityIntegrate.setMaxSteps(500);
        getController().addAction(activityIntegrate);
        
        species = new SpeciesSpheresMono[numAtoms.length];
        for(int i=0;i<numAtoms.length;i++){
        	species[i] = new SpeciesSpheresMono(this);
        	species[i].setNMolecules(numAtoms[i]);
        	((etomica.atom.AtomTypeSphere)species[i].getFactory().getType()).setDiameter(atomicSize[i]);
        	if (i!=(numAtoms.length-1)){
        	((etomica.atom.AtomTypeLeaf)species[i].getFactory().getType()).setMass(Double.POSITIVE_INFINITY);
        	}
        }
        //Setting up potential for Methane-Methane interactions
        potentialMM = new P2LennardJones(this);
        //Setting sigma and epsilon (from J. Chem. Soc Faraday Trans 1991
        potentialMM.setEpsilon(Kelvin.UNIT.toSim(147.95));
        //potentialMM.setEpsilon(1);
        potentialMM.setSigma(3.0);
        //Setting up potential for Methane-Oxygen interactions
        P2LennardJones potentialMO = new P2LennardJones(this);
        //Setting sigma and epsilon 
        potentialMO.setEpsilon(Kelvin.UNIT.toSim(133.3));
        //potentialMO.setEpsilon(100);
        potentialMO.setSigma(3.214);
        
        //Setting up Methane - Silicon interactions
        //P2LennardJones potentialMS = potentialMO;
        P2WCA potentialMS = new P2WCA(space,atomicSize[1],potentialMO.getEpsilon());
        
        //Wrap LJ potentials to truncate
        P2SoftSphericalTruncated MM = new P2SoftSphericalTruncated(potentialMM,2.5*potentialMM.getSigma());
        P2SoftSphericalTruncated MO = new P2SoftSphericalTruncated(potentialMO,2.5*potentialMO.getSigma());
        //P2SoftSphericalTruncated MS = new P2SoftSphericalTruncated(potentialMS,2.5*potentialMS.getSigma());
        
        
        NeighborCriterion criterionMM = new CriterionSpecies(new CriterionSimple(this,MM.getRange(), MM.getRange()*neighborRangeFac), species[2], species[2]);
        NeighborCriterion criterionMO = new CriterionSpecies(new CriterionSimple(this,MO.getRange(), MO.getRange()*neighborRangeFac), species[0], species[2]);
        NeighborCriterion criterionMS = new CriterionSpecies(new CriterionSimple(this,potentialMS.getRange(), potentialMS.getRange()*neighborRangeFac*2), species[1], species[2]);
        
        MM.setCriterion(criterionMM);
        MO.setCriterion(criterionMO);
        potentialMS.setCriterion(criterionMS);
        
        
        potentialMaster.addPotential(MM,new Species[]{species[2],species[2]});
        potentialMaster.addPotential(MO,new Species[]{species[0],species[2]});
        potentialMaster.addPotential(potentialMS,new Species[]{species[1],species[2]});
        
        phase = new Phase(this);
        //Initializes the coordinates and positions
        config.initializeCoordinates(phase);
        phase.getBoundary().setDimensions(config.getUpdatedDimensions());
        integrator.setPhase(phase);
        //integrator.addListener(new PhaseImposePbc(phase));
          
        //PARAMETERS For Simulation Run
        //activityIntegrate.setMaxSteps(5000000);
        activityIntegrate.setMaxSteps(1000000);
        double ts = 0.00611;
        integrator.setTimeStep(ts);
        interval = 2000;
        integrator.setThermostatInterval(interval/1000);
        
        //      Adding coordinate writer by Mike Sellars
     
        filename = (numAtoms[2]+"_"+activityIntegrate.getMaxSteps()+"_"+ts+"_"+interval+"_WCA");
        sp = new Species[1];
        sp[0] = (Species)species[2];
        /*
        MSDCoordWriter coordWriter = new MSDCoordWriter(this.space, filename,sp);
        coordWriter.setPhase(this.phase);
        coordWriter.setNatoms(numAtoms[2]);
        coordWriter.setIntegrator(integrator);
        coordWriter.setWriteInterval(interval);
        */
    } //end of constructor
    int getInterval(){
    	return interval;
    }
    int getMethane(){
    	return nAtomsMeth;
    }
    String getFileName(){
    	return filename;
    }
    Species[] getSpeciesRMS(){
    	return sp;
    }
    private Species[] sp;
    private String filename;
    public static void Converter(String inputFile,int meth) {
		// TODO Auto-generated method stub
		String outputFile = inputFile+"__Result.txt";
		MSDProcessor proc = new MSDProcessor(Space3D.getInstance(),inputFile,outputFile);
		
		//proc.setDeltaTmax(1);
		proc.setMethane(meth);
		proc.fillArrays();
		System.out.println("Converter done");
	}
    /**
     * Demonstrates how this class is implemented.
     */
    public static void main(String[] args) {
    	Default defaults = new Default();
        defaults.doSleep = true;
        defaults.ignoreOverlap = true;
        //defaults.temperature = Kelvin.UNIT.toSim(298.0);
        ZeoliteSimulation sim = new ZeoliteSimulation(defaults);
        zeoliteSimGraphic simGraphic = new zeoliteSimGraphic(sim);
        int num = sim.species.length;
        DeviceNSelector nSelector = new DeviceNSelector(sim,sim.phase.getAgent(sim.species[num-1]));
        simGraphic.add(nSelector);
        
        //Energy
        int history = sim.getInterval()*10;
        //Settings
        
        MeterEnergy eMeter = new MeterEnergy(sim.potentialMaster);
        eMeter.setPhase(sim.phase);
        AccumulatorHistory energyHistory = new AccumulatorHistory();
        energyHistory.setHistoryLength(history);
        AccumulatorAverage enAcc = new AccumulatorAverage(sim);
        enAcc.setPushInterval(20);
        DataFork enFork = new DataFork(new DataSink[]{energyHistory, enAcc});
        DataPump energyPump = new DataPump(eMeter, enFork);
        IntervalActionAdapter energyAdapter = new IntervalActionAdapter(energyPump);
        energyAdapter.setActionInterval(10);
        energyHistory.setPushInterval(10);
        sim.integrator.addListener(energyAdapter);
        sim.register(eMeter,energyPump);
		
        /*
        MeterPotentialEnergy peMeter = new MeterPotentialEnergy(sim.potentialMaster);
        peMeter.setPhase(sim.phase);
        AccumulatorHistory peHistory = new AccumulatorHistory();
        peHistory.setHistoryLength(history);
        AccumulatorAverage peAccumulator = new AccumulatorAverage(sim);
        DataFork peFork = new DataFork(new DataSink[]{peHistory, peAccumulator});
        DataPump pePump = new DataPump(peMeter, peFork);
        IntervalActionAdapter peAdapter = new IntervalActionAdapter(pePump);
        peAdapter.setActionInterval(10);
        peHistory.setPushInterval(10);
        sim.register(peMeter,pePump);
        sim.integrator.addListener(peAdapter);
		
		MeterKineticEnergy keMeter = new MeterKineticEnergy();
        keMeter.setPhase(sim.phase);
        AccumulatorHistory keHistory = new AccumulatorHistory();
        keHistory.setHistoryLength(history);
        DataPump kePump = new DataPump(keMeter, keHistory);
        IntervalActionAdapter keAdapter = new IntervalActionAdapter(kePump);
        keAdapter.setActionInterval(10);
        keHistory.setPushInterval(10);
        sim.register(keMeter,kePump);
        sim.integrator.addListener(keAdapter);
        
        
        MeterTemperature temp = new MeterTemperature();
		temp.setPhase(sim.phase);
		AccumulatorHistory tHistory = new AccumulatorHistory();
        tHistory.setHistoryLength(history);
        AccumulatorAverage tAccumulator = new AccumulatorAverage(sim);
        DataFork tFork = new DataFork(new DataSink[]{tHistory, tAccumulator});
        DataPump tPump = new DataPump(temp, tFork);
        IntervalActionAdapter tAdapter = new IntervalActionAdapter(tPump);
        tAdapter.setActionInterval(10);
        tHistory.setPushInterval(10);
        sim.register(temp,tPump);
        sim.integrator.addListener(tAdapter);
        */
        
        DisplayPlot ePlot = sim.ePlot;
        ePlot = new DisplayPlot();
        ePlot.setLabel("Energy");
        energyHistory.setDataSink(ePlot.getDataSet().makeDataSink());
        //peHistory.setDataSink(ePlot.getDataSet().makeDataSink());
        //keHistory.setDataSink(ePlot.getDataSet().makeDataSink());
		//tHistory.setDataSink(ePlot.getDataSet().makeDataSink());
		ePlot.setDoLegend(true);
		simGraphic.add(ePlot);
		
		
		
		/*
		DisplayBoxesCAE energy = new DisplayBoxesCAE();
		energy.setAccumulator(enAcc);
		simGraphic.add(energy);
		*/
		
		
        simGraphic.makeAndDisplayFrame();
        ColorSchemeByType colorScheme = ((ColorSchemeByType)((DisplayPhase)simGraphic.displayList().getFirst()).getColorScheme());
        for(int i=0;i<sim.species.length;i++){
        	switch(i){
        		case 0:
        			colorScheme.setColor(sim.species[i].getMoleculeType(), java.awt.Color.red);
        			break;
        		case 1:
        			colorScheme.setColor(sim.species[i].getMoleculeType(), java.awt.Color.blue);
        			break;
        		default:
        			colorScheme.setColor(sim.species[i].getMoleculeType(), java.awt.Color.white);
        	}
        	
        }
        simGraphic.panel().setBackground(java.awt.Color.yellow);
    }//end of main

}//end of class
