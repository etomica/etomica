package etomica.paracetamol;

import etomica.action.PDBWriter;
import etomica.action.WriteConfiguration;
import etomica.action.WriteConfigurationDLPOLY;
import etomica.action.activity.ActivityIntegrate;
import etomica.action.activity.Controller;
import etomica.atom.AtomTypeGroup;
import etomica.box.Box;
import etomica.data.DataLogger;
import etomica.data.DataPump;
import etomica.data.DataTableWriter;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveMolecule;
import etomica.integrator.mcmove.MCMoveRotateMolecule3D;
import etomica.integrator.mcmove.MCMoveStepTracker;
import etomica.lattice.BravaisLatticeCrystal;
import etomica.lattice.crystal.PrimitiveOrthorhombic;
import etomica.normalmode.MCMoveHarmonicStep;
import etomica.normalmode.NormalModesFromFile;
import etomica.normalmode.WaveVectorFactory;
import etomica.normalmode.CoordinateDefinition.BasisCell;
import etomica.potential.PotentialDLPOLY;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space.Space;
import etomica.species.Species;
import etomica.units.Kelvin;

/**
 * 
 * Three-dimensional soft-sphere MC simulation for paracetamol molecule
 *  using the potential model from DL_MULTI package
 *  with inclusion of multipole ewald-summation
 * 
 * Orthorhombic Crystal
 * 
 *  * The constructor takes in:
 * 		1. Space
 * 		2. Number of molecules
 * 		3. Temperature
 * 		4. Type of simulation 
 * 
 * Selection of type of simulation:
 * 		0. MC Simulation (MCMoveMolecule & MCMoveRotateMolecule3D)
 * 		1. Equilibration (MCMoveHarmonicStep)
 * 		2. SimCalcS      (MCMoveRotateMolecule3D & MCMoveMoleculeCoupledDLPOLY)
 * 
 * @author Tai Tan
 *
 */
public class MCParacetamolOrthorhombicDLMULTI extends Simulation {

	private static final long serialVersionUID = 1L;
//	private final static String APP_NAME = "MC Paracetamol Orthorhombic";
    public Box box;
    public IntegratorMC integrator;
    public MCMoveMolecule mcMoveMolecule;
    public MCMoveRotateMolecule3D mcMoveRotateMolecule;
    public MCMoveHarmonicStep mcMoveHarmonicStep;
    public MCMoveMoleculeCoupledDLPOLY mcMoveMoleculeCoupledDLPOLY;
    public SpeciesParacetamol species;
    public Controller controller;
    public static final int MCSIMULATION = 0;
    public static final int EQUILIBRATION = 1;
    public static final int SIMCALCS = 2;
  
    public MCParacetamolOrthorhombicDLMULTI(Space space, int numMolecules, double temperature, int simType) {
    	super(space, true);
        potentialMaster = new PotentialMaster(space);
        this.simType = simType;
    	
    	/*
    	 * Orthorhombic Crystal
    	 */
    	
        primitive = new PrimitiveOrthorhombic(space, 17.248, 12.086, 7.382);
        // 17.248, 12.086, 7.382
        BasisOrthorhombicParacetamol basis = new BasisOrthorhombicParacetamol();
        lattice = new BravaisLatticeCrystal(primitive, basis); 
        
        integrator = new IntegratorMC(this, potentialMaster);
        integrator.setIsothermal(false);
        //integrator.setThermostatInterval(1);
        integrator.setTemperature(Kelvin.UNIT.toSim(temperature));

        actionIntegrate = new ActivityIntegrate(integrator, 0, false);
        getController().addAction(actionIntegrate);
        
        ConformationParacetamolOrthorhombic conformation = new ConformationParacetamolOrthorhombic(space);
        species = new SpeciesParacetamol(this);
        ((AtomTypeGroup)species.getMoleculeType()).setConformation(conformation);
        getSpeciesManager().addSpecies(species);
        
        box = new Box(this);
        addBox(box);
        box.setDimensions(Space.makeVector(new double[] {25,25,25}));
        box.setNMolecules(species, numMolecules);        

        bdry =  new BoundaryRectangularPeriodic(space, getRandom(), 1); //unit cell
        bdry.setDimensions(Space.makeVector(new double []{2*17.248, 3*12.086, 4*7.382}));
        box.setBoundary(bdry);

        coordDef = new CoordinateDefinitionParacetamol(box, primitive, basis);
        coordDef.setBasisOrthorhombic();
        
   
        if (simType == 0){
        	coordDef.initializeCoordinates(new int []{2, 3, 4});
            mcMoveMolecule = new MCMoveMolecule(this, potentialMaster);
            mcMoveMolecule.setStepSize(0.2077);  //Step size to input
            ((MCMoveStepTracker)mcMoveMolecule.getTracker()).setTunable(true);
            ((MCMoveStepTracker)mcMoveMolecule.getTracker()).setNoisyAdjustment(true);
            
            mcMoveRotateMolecule = new MCMoveRotateMolecule3D(potentialMaster, random);
            mcMoveRotateMolecule.setStepSize(0.0922);
            ((MCMoveStepTracker)mcMoveRotateMolecule.getTracker()).setNoisyAdjustment(true);
         
            integrator.getMoveManager().setEquilibrating(true);
            integrator.getMoveManager().addMCMove(mcMoveMolecule);
            integrator.getMoveManager().addMCMove(mcMoveRotateMolecule);
        } else 
        	if (simType == 1){
        		coordDef.initializeCoordinates(new int []{2, 3, 4});
        		mcMoveHarmonicStep = new MCMoveHarmonicStep(potentialMaster, getRandom());
                mcMoveHarmonicStep.setCoordinateDefinition(coordDef);
                mcMoveHarmonicStep.setBox(box);
                
                int[] modes = new int[45];
                for (int i=0; i<45; i++){
                	modes[i] = i+3;
                }
                
                NormalModesFromFile normalModes = new NormalModesFromFile("Normal_Modes_Paracetamol_FormII_10.0K",3);
                WaveVectorFactory waveVectorFactory = normalModes.getWaveVectorFactory();
                waveVectorFactory.makeWaveVectors(box);
                
                mcMoveHarmonicStep.setModes(modes);
                mcMoveHarmonicStep.setEigenVectors(normalModes.getEigenvectors(box)[0]);
                mcMoveHarmonicStep.setStepSize(0.006);
                ((MCMoveStepTracker)mcMoveHarmonicStep.getTracker()).setAdjustInterval(5);
                ((MCMoveStepTracker)mcMoveHarmonicStep.getTracker()).setNoisyAdjustment(true);
               
                integrator.getMoveManager().addMCMove(mcMoveHarmonicStep);
                integrator.getMoveManager().setEquilibrating(true);
                
        	} else 
        		if (simType == 2){
        			
        	        mcMoveRotateMolecule = new MCMoveRotateMolecule3D(potentialMaster, random);
        	        mcMoveRotateMolecule.setStepSize(0.0922);
        	        ((MCMoveStepTracker)mcMoveRotateMolecule.getTracker()).setNoisyAdjustment(true);
        	      
        	        mcMoveMoleculeCoupledDLPOLY = new MCMoveMoleculeCoupledDLPOLY(potentialMaster, getRandom());
        	        
        	        //integrator.getMoveManager().addMCMove(mcMoveMolecule);
        	        integrator.getMoveManager().addMCMove(mcMoveMoleculeCoupledDLPOLY);
        	        integrator.getMoveManager().addMCMove(mcMoveRotateMolecule);
        	        integrator.getMoveManager().setEquilibrating(true);
        	        
        	        coordDef.initializeCoordinates(new int []{1, 2, 2});
        	        double[] u = new double[]{-8.519546134513288, -6.357161053952702, -7.449621419035304, -0.003932863936184616, -0.0016776660890108416, -0.013522565626184084, -8.726086341054327, -5.722820924623122, -7.451563204983988, -0.005551030882875035, -0.008781239888660329, 0.005654659111204717, -8.737722525835544, -6.358584368530293, -7.312474753722753, -0.004632407803175042, 0.0073903160411125865, -0.0022154770629279168, -8.513251161366203, -5.718672928462283, -7.308964948201192, -0.008253140936357385, -0.0016033969986892181, -0.007589294260920792, -8.726675083680286, -5.711819873904577, -7.345829560494344, 0.0012238914448621102, 0.009773882446509393, 0.0029289771411575016, -8.4909577926205, -6.372620681124473, -7.328828834258954, 0.008259322829339406, 7.612944355468755E-4, 0.015416013389648568, -8.536258682509295, -5.721912267935009, -7.4175804780220504, 0.0048799658666174375, 0.01669560616620193, -0.004313176827043905, -8.741502278415917, -6.38040790144197, -7.441136801296486, 0.0043516618097694075, 0.010742961480577531, 0.0048663615734551415}; 
        	        BasisCell[] cell = coordDef.getBasisCells() ;
        	        for (int i=0; i<cell.length; i++){
        	        	coordDef.setToU(cell[i].molecules, u);
        	        }
        	        
        		} else {
        			
        			throw new RuntimeException ("Need to input type of simulation!!!!");
        		}
        
        WriteConfigurationDLPOLY configDLPOLY = new WriteConfigurationDLPOLY();
        configDLPOLY.setConfName("CONFIG");
        configDLPOLY.setBox(box);
        configDLPOLY.getElementHash().put(HydrogenP.INSTANCE, "HP");
  
        PotentialDLPOLY potentialDLPOLY = new PotentialDLPOLY(space);
        potentialDLPOLY.setConfigDLPOLY(configDLPOLY);
        potentialMaster.addPotential(potentialDLPOLY, new Species[0]);

        integrator.setBox(box);
        
    } //end of constructor
    
   
    public static void main(String[] args) {
    	
    	int numMolecules = 192;
    	double temperature = Kelvin.UNIT.toSim(10);
    	long simSteps = 100;
    	int simType = 0;
    	
        String filename = "Paracetamol_Orthorhombic_"+ Kelvin.UNIT.fromSim(temperature)+"K";
        if (args.length > 0) {
            filename = args[0];
        }
        if (args.length > 1) {
            simSteps = Long.parseLong(args[1]);
        }
        if (args.length > 2) {
            temperature = Kelvin.UNIT.toSim(Double.parseDouble(args[2]));
        }
        if (args.length > 3) {
            simType = Integer.parseInt(args[3]);
        }

        System.out.println("Running "+ "Orthorhombic Paracetamol simulation");
        System.out.println("Type of simulation: "+simType);
        System.out.println(numMolecules + " molecules " +" and temperature "+ Kelvin.UNIT.fromSim(temperature) +"K");
        System.out.println(simSteps+ " steps");
        System.out.println("output data to " + filename);
        
    	MCParacetamolOrthorhombicDLMULTI sim = 
        	new MCParacetamolOrthorhombicDLMULTI(Space.getInstance(3), numMolecules, temperature, simType);
       
       sim.actionIntegrate.setMaxSteps(simSteps);
       MeterPotentialEnergy meterPE = new MeterPotentialEnergy(sim.potentialMaster);
       meterPE.setBox(sim.box);
       
       DataLogger dataLoggerPE = new DataLogger();
       dataLoggerPE.setWriteInterval(5);
       dataLoggerPE.setFileName("Paracetamol_Orthorhombic_"+Kelvin.UNIT.fromSim(temperature)+"K");
       dataLoggerPE.setDataSink(new DataTableWriter());
       dataLoggerPE.setAppending(true);
       dataLoggerPE.setCloseFileEachTime(true);
       
       DataPump PEpump = new DataPump(meterPE, dataLoggerPE);
       sim.getController().getEventManager().addListener(dataLoggerPE);
       
       sim.integrator.addIntervalAction(PEpump);
       sim.integrator.setActionInterval(PEpump, 1);
       
       
       sim.getController().actionPerformed();
       
       WriteConfiguration writeConfig = new WriteConfiguration();
       writeConfig.setConfName("FinalCoord_Paracetamol_Orthorhombic_"+Kelvin.UNIT.fromSim(temperature)+"K");
       writeConfig.setBox(sim.box);
       writeConfig.setDoApplyPBC(false);
       writeConfig.actionPerformed();
        
       PDBWriter pdbWriter = new PDBWriter(sim.box);
       pdbWriter.setFileName("GraphicView_Orthorhombic_"+ Kelvin.UNIT.fromSim(temperature)+"K.pdb");
       pdbWriter.actionPerformed();
        
    }//end of main

    public PotentialMaster potentialMaster;
    public BravaisLatticeCrystal lattice;
    public BoundaryRectangularPeriodic bdry;
    public ActivityIntegrate actionIntegrate;
    public PrimitiveOrthorhombic primitive;
    public CoordinateDefinitionParacetamol coordDef;
    public int simType;
}//end of class