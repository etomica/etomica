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
 *  To get the lowest minimum cofiguration with MCMoveHarmonicStep
 * 
 * Orthorhombic Crystal
 * 
 * @author Tai Tan
 *
 */
public class MCParacetamolOrthorhombicDLMULTIEquilibration extends Simulation {

	private static final long serialVersionUID = 1L;
//	private final static String APP_NAME = "MC Paracetamol Orthorhombic";
    public Box box;
    public IntegratorMC integrator;
    public MCMoveMolecule mcMoveMolecule;
    public MCMoveRotateMolecule3D mcMoveRotateMolecule;
    public SpeciesParacetamol species;
    public Controller controller;

  
    public MCParacetamolOrthorhombicDLMULTIEquilibration(Space space, int numMolecules, double temperature) {
    	super(space, true);
        potentialMaster = new PotentialMaster(space);
    	
    	/*
    	 * Orthorhombic Crystal
    	 */
    	
        PrimitiveOrthorhombic primitive = new PrimitiveOrthorhombic(space, 17.248, 12.086, 7.382);
        // 17.248, 12.086, 7.382
        BasisOrthorhombicParacetamol basis = new BasisOrthorhombicParacetamol();
        lattice = new BravaisLatticeCrystal(primitive, basis); 
        
        integrator = new IntegratorMC(this, potentialMaster);
        integrator.setIsothermal(false);
        //integrator.setThermostatInterval(1);
        integrator.setTemperature(Kelvin.UNIT.toSim(temperature));

        /*
        mcMoveMolecule = new MCMoveMolecule(this, potentialMaster);
        mcMoveMolecule.setStepSize(0.2077);  //Step size to input
        ((MCMoveStepTracker)mcMoveMolecule.getTracker()).setTunable(true);
        ((MCMoveStepTracker)mcMoveMolecule.getTracker()).setNoisyAdjustment(true);
        
        mcMoveRotateMolecule = new MCMoveRotateMolecule3D(potentialMaster, random);
        mcMoveRotateMolecule.setStepSize(0.0922);
        ((MCMoveStepTracker)mcMoveRotateMolecule.getTracker()).setNoisyAdjustment(true);
      	*/
      
        //integrator.getMoveManager().addMCMove(mcMoveMolecule);
        //integrator.getMoveManager().addMCMove(mcMoveRotateMolecule);

        
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

        CoordinateDefinitionParacetamol coordDef = new CoordinateDefinitionParacetamol(box, primitive, basis);
        coordDef.setBasisOrthorhombic();
        coordDef.initializeCoordinates(new int []{2, 3, 4});
   
        /*
         * For Equilibration Purposes
         */
        MCMoveHarmonicStep mcMoveHarmonicStep = new MCMoveHarmonicStep(potentialMaster, getRandom());
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
        ///////////////////////////////////////////////////////////////////////////////////////
        
        WriteConfigurationDLPOLY configDLPOLY = new WriteConfigurationDLPOLY();
        configDLPOLY.setConfName("CONFIG");
        configDLPOLY.setBox(box);
        configDLPOLY.getElementHash().put(HydrogenP.INSTANCE, "HP");
//        configDLPOLY.actionPerformed();
//        System.exit(1); 
        
        PotentialDLPOLY potentialDLPOLY = new PotentialDLPOLY(space);
        potentialDLPOLY.setConfigDLPOLY(configDLPOLY);
        potentialMaster.addPotential(potentialDLPOLY, new Species[0]);
        
//        XYZWriter pdbWriter = new XYZWriter(box);
//        pdbWriter.setFileName("TestPoly_CONFIG_Paracetamol");
//        pdbWriter.actionPerformed();
        
        
        integrator.setBox(box);
        //BoxImposePbc pbc = new BoxImposePbc(box);
        //pbc.actionPerformed();
        //pbc.setApplyToMolecules(true);
        //integrator.addListener(new IntervalActionAdapter(pbc));

        
    } //end of constructor
    
   
    public static void main(String[] args) {
    	
    	int numMolecules = 192;
    	double temperature = Kelvin.UNIT.toSim(10);
    	long simSteps = 100;
    	
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

        System.out.println("Running "+ "Orthorhombic Paracetamol Equilibration simulation");
        System.out.println(numMolecules + " molecules " +" and temperature "+ Kelvin.UNIT.fromSim(temperature) +"K");
        System.out.println(simSteps+ " steps");
        System.out.println("output data to " + filename);
        
    	etomica.paracetamol.MCParacetamolOrthorhombicDLMULTIEquilibration sim = 
        	new etomica.paracetamol.MCParacetamolOrthorhombicDLMULTIEquilibration(Space.getInstance(3), numMolecules, temperature);
       
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
   	private String CoordFile;
}//end of class