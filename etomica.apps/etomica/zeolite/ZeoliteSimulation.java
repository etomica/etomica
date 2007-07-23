package etomica.zeolite;

import etomica.action.SimulationRestart;
import etomica.action.activity.ActivityIntegrate;
import etomica.box.Box;
import etomica.chem.elements.ElementSimple;
import etomica.data.AccumulatorAverage;
import etomica.data.AccumulatorAverageCollapsing;
import etomica.data.AccumulatorHistory;
import etomica.data.DataFork;
import etomica.data.DataPump;
import etomica.data.DataSink;
import etomica.data.meter.MeterEnergy;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.DeviceNSelector;
import etomica.graphics.DisplayBox;
import etomica.graphics.DisplayPlot;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.nbr.list.NeighborListManager;
import etomica.nbr.list.PotentialMasterList;
import etomica.potential.P2LennardJones;
import etomica.potential.P2SoftSphericalTruncated;
import etomica.potential.P2WCA;
import etomica.simulation.Simulation;
import etomica.space3d.Space3D;
import etomica.species.Species;
import etomica.species.SpeciesSpheresMono;
import etomica.units.Kelvin;


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

    private static final long serialVersionUID = 1L;

    private static final String APP_NAME = "Zeolite Simulation";

    /**
     * The Box holding the atoms. 
     */
    public final Box box;
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
    private int interval;
    public ActivityIntegrate activityIntegrate;
    public DisplayPlot ePlot;
    /**
     * Sole public constructor, makes a simulation using a 3D space.
     */
    //we use a second, private constructor to permit the space to
    //appear twice in the call to the superclass constructor; alternatively
    //we could have passed Space3D.getInstance() twice
    private ZeoliteSimulation() {

        // invoke the superclass constructor
        // "true" is indicating to the superclass that this is a dynamic simulation
        // the PotentialMaster is selected such as to implement neighbor listing
        super(Space3D.getInstance(), true);

        PotentialMasterList potentialMaster = new PotentialMasterList(this, 1.6);
        //Additions for Zeolite Calculations
        //Start by reading the first line, which is number of Atoms
        String fileName = "2unitcell";
        //String fileName = "pbu2";
        ConfigurationFileXYZ config = new ConfigurationFileXYZ(fileName);
        int[] numAtoms = config.getNumAtoms();
        
        nAtomsMeth = numAtoms[numAtoms.length - 1];
        double neighborRangeFac = 1.2;
        //defaults.makeLJDefaults();
        //Setting sizes of molecules
        double[] atomicSize = new double[3];
        atomicSize[0] = 0.73;
        atomicSize[1] = 1.18;
        atomicSize[2] = 2.088;
        
        double range = 8.035;
        potentialMaster.setRange(3.214*neighborRangeFac*2.5);
        
        
        integrator = new IntegratorVelocityVerlet(this, potentialMaster);
        integrator.setIsothermal(true);
        integrator.setThermostatInterval(10);
        integrator.setTimeStep(0.00611);
        integrator.setTemperature(Kelvin.UNIT.toSim(298.0));
        

        activityIntegrate = new ActivityIntegrate(integrator, 2, true);
        activityIntegrate.setMaxSteps(500);
        getController().addAction(activityIntegrate);
        
        box = new Box(this);
        addBox(box);
        NeighborListManager nbrManager = potentialMaster.getNeighborManager(box);
        integrator.addNonintervalListener(nbrManager);
        integrator.addIntervalAction(nbrManager);
        species = new SpeciesSpheresMono[numAtoms.length];
        for(int i=0;i<numAtoms.length;i++){
        	species[i] = new SpeciesSpheresMono(this);
            getSpeciesManager().addSpecies(species[i]);
        	box.setNMolecules(species[i], numAtoms[i]);
        	((etomica.atom.AtomTypeSphere)species[i].getMoleculeType()).setDiameter(atomicSize[i]);
        	if (i!=(numAtoms.length-1)){
                // all elements except the last (methane) are fixed
        	    ((ElementSimple)((etomica.atom.AtomTypeLeaf)species[i].getMoleculeType()).getElement()).setMass(Double.POSITIVE_INFINITY);
        	}
            else {
                ((ElementSimple)((etomica.atom.AtomTypeLeaf)species[i].getMoleculeType()).getElement()).setMass(16);
            }
        }
        //Setting up potential for Methane-Methane interactions
        potentialMM = new P2LennardJones(space);
        //Setting sigma and epsilon (from J. Chem. Soc Faraday Trans 1991
        potentialMM.setEpsilon(Kelvin.UNIT.toSim(147.95));
        //potentialMM.setEpsilon(1);
        potentialMM.setSigma(3.0);
        //Setting up potential for Methane-Oxygen interactions
        P2LennardJones potentialMO = new P2LennardJones(space);
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
        
        
        potentialMaster.addPotential(MM,new Species[]{species[2],species[2]});
        potentialMaster.addPotential(MO,new Species[]{species[0],species[2]});
        potentialMaster.addPotential(potentialMS,new Species[]{species[1],species[2]});
        
        //Initializes the coordinates and positions
        config.initializeCoordinates(box);
        box.getBoundary().setDimensions(config.getUpdatedDimensions());
        integrator.setBox(box);
        //integrator.addListener(new BoxImposePbc(box));
          
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
        sp[0] = species[2];
        /*
        MSDCoordWriter coordWriter = new MSDCoordWriter(this.space, filename,sp);
        coordWriter.setBox(this.box);
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
    public static void Converter(String inputFile) {
		// TODO Auto-generated method stub
		String outputFile = inputFile+"__Result.txt";
		MSDProcessor proc = new MSDProcessor(Space3D.getInstance(),inputFile,outputFile);
		
		//proc.setDeltaTmax(1);
		proc.fillArrays();
		System.out.println("Converter done");
	}
    /**
     * Demonstrates how this class is implemented.
     */
    public static void main(String[] args) {
        //defaults.temperature = Kelvin.UNIT.toSim(298.0);
        ZeoliteSimulation sim = new ZeoliteSimulation();
        zeoliteSimGraphic simGraphic = new zeoliteSimGraphic(sim, APP_NAME);
        int num = sim.species.length;
        DeviceNSelector nSelector = new DeviceNSelector(sim.getController());
        nSelector.setResetAction(new SimulationRestart(sim));
        nSelector.setSpecies(sim.species[num-1]);
        nSelector.setBox(sim.box);
        simGraphic.add(nSelector);
        
        //Energy
        int history = sim.getInterval()*10;
        //Settings
        
        MeterEnergy eMeter = new MeterEnergy(sim.integrator.getPotential());
        eMeter.setBox(sim.box);
        AccumulatorHistory energyHistory = new AccumulatorHistory();
        energyHistory.getHistory().setHistoryLength(history);
        AccumulatorAverage enAcc = new AccumulatorAverageCollapsing();
        enAcc.setPushInterval(20);
        DataFork enFork = new DataFork(new DataSink[]{energyHistory, enAcc});
        DataPump energyPump = new DataPump(eMeter, enFork);
        sim.integrator.addIntervalAction(energyPump);
        sim.integrator.setActionInterval(energyPump, 10);
        energyHistory.setPushInterval(10);
        simGraphic.getController().getDataStreamPumps().add(energyPump);
		
        /*
        MeterPotentialEnergy peMeter = new MeterPotentialEnergy(sim.potentialMaster);
        peMeter.setBox(sim.box);
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
        keMeter.setBox(sim.box);
        AccumulatorHistory keHistory = new AccumulatorHistory();
        keHistory.setHistoryLength(history);
        DataPump kePump = new DataPump(keMeter, keHistory);
        IntervalActionAdapter keAdapter = new IntervalActionAdapter(kePump);
        keAdapter.setActionInterval(10);
        keHistory.setPushInterval(10);
        sim.register(keMeter,kePump);
        sim.integrator.addListener(keAdapter);
        
        
        MeterTemperature temp = new MeterTemperature();
		temp.setBox(sim.box);
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


        simGraphic.makeAndDisplayFrame(APP_NAME);

        ColorSchemeByType colorScheme = ((ColorSchemeByType)((DisplayBox)simGraphic.displayList().getFirst()).getColorScheme());
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

    }//end of main

}//end of class
