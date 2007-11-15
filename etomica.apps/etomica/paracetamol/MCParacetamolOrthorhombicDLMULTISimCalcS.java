package etomica.paracetamol;

import etomica.action.WriteConfiguration;
import etomica.action.WriteConfigurationDLPOLY;
import etomica.action.activity.ActivityIntegrate;
import etomica.action.activity.Controller;
import etomica.atom.AtomTypeGroup;
import etomica.box.Box;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveRotateMolecule3D;
import etomica.integrator.mcmove.MCMoveStepTracker;
import etomica.lattice.BravaisLatticeCrystal;
import etomica.lattice.crystal.PrimitiveOrthorhombic;
import etomica.normalmode.MeterNormalMode;
import etomica.normalmode.WaveVectorFactorySimple;
import etomica.normalmode.WriteS;
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
 * MC simulation of Paracetamol molecules in Form II (Orthorhombic) with 
 *  tabulation of the collective-coordinate S-matrix.
 * No graphic display
 * 
 * Orthorhombic Crystal
 * 
 * @author Tai Tan
 *
 */
public class MCParacetamolOrthorhombicDLMULTISimCalcS extends Simulation {

	private static final long serialVersionUID = 1L;
    public Box box;
    public IntegratorMC integrator;
    //public MCMoveMolecule mcMoveMolecule;
    public MCMoveMoleculeCoupledDLPOLY mcMoveMoleculeCoupledDLPOLY;
    public MCMoveRotateMolecule3D mcMoveRotateMolecule;
    public SpeciesParacetamol species;
    public Controller controller;

  
    public MCParacetamolOrthorhombicDLMULTISimCalcS(Space space, int numMolecules, double temperature) {
    	super(space, true);
        potentialMaster = new PotentialMaster(space);
    	
    	/*
    	 * Orthorhombic Crystal
    	 */
    	
        primitive = new PrimitiveOrthorhombic(space, 17.248, 12.086, 7.382);
        // 17.248, 12.086, 7.382
        BasisOrthorhombicParacetamol basis = new BasisOrthorhombicParacetamol();
        lattice = new BravaisLatticeCrystal(primitive, basis);
        
        integrator = new IntegratorMC(this, potentialMaster);
        integrator.setIsothermal(false);
        integrator.setTemperature(Kelvin.UNIT.toSim(temperature));
        
        mcMoveRotateMolecule = new MCMoveRotateMolecule3D(potentialMaster, random);
        mcMoveRotateMolecule.setStepSize(0.0922);
        ((MCMoveStepTracker)mcMoveRotateMolecule.getTracker()).setNoisyAdjustment(true);
      
        mcMoveMoleculeCoupledDLPOLY = new MCMoveMoleculeCoupledDLPOLY(potentialMaster, getRandom());
        
        //integrator.getMoveManager().addMCMove(mcMoveMolecule);
        integrator.getMoveManager().addMCMove(mcMoveMoleculeCoupledDLPOLY);
        integrator.getMoveManager().addMCMove(mcMoveRotateMolecule);
        integrator.getMoveManager().setEquilibrating(true);
        
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
        bdry.setDimensions(Space.makeVector(new double []{1*17.248, 2*12.086, 2*7.382}));
        box.setBoundary(bdry);

    	//ConfigurationFile configFile = new ConfigurationFile("FinalCoord_Min_Paracetamol_Orthorhombic");
        coordinateDefinition = new CoordinateDefinitionParacetamol(box, primitive, basis);
        coordinateDefinition.setBasisOrthorhombic();
        //coordinateDefinition.setConfiguration(configFile);
        coordinateDefinition.initializeCoordinates(new int []{1, 2, 2});
        double[] u = new double[]{-8.519546134513288, -6.357161053952702, -7.449621419035304, -0.003932863936184616, -0.0016776660890108416, -0.013522565626184084, -8.726086341054327, -5.722820924623122, -7.451563204983988, -0.005551030882875035, -0.008781239888660329, 0.005654659111204717, -8.737722525835544, -6.358584368530293, -7.312474753722753, -0.004632407803175042, 0.0073903160411125865, -0.0022154770629279168, -8.513251161366203, -5.718672928462283, -7.308964948201192, -0.008253140936357385, -0.0016033969986892181, -0.007589294260920792, -8.726675083680286, -5.711819873904577, -7.345829560494344, 0.0012238914448621102, 0.009773882446509393, 0.0029289771411575016, -8.4909577926205, -6.372620681124473, -7.328828834258954, 0.008259322829339406, 7.612944355468755E-4, 0.015416013389648568, -8.536258682509295, -5.721912267935009, -7.4175804780220504, 0.0048799658666174375, 0.01669560616620193, -0.004313176827043905, -8.741502278415917, -6.38040790144197, -7.441136801296486, 0.0043516618097694075, 0.010742961480577531, 0.0048663615734551415}; 
        BasisCell[] cell = coordinateDefinition.getBasisCells() ;
        for (int i=0; i<cell.length; i++){
        	coordinateDefinition.setToU(cell[i].molecules, u);
        }
   
        WriteConfigurationDLPOLY configDLPOLY = new WriteConfigurationDLPOLY();
        configDLPOLY.setConfName("CONFIG");
        configDLPOLY.setBox(box);
        configDLPOLY.getElementHash().put(HydrogenP.INSTANCE, "HP");
        
        PotentialDLPOLY potentialDLPOLY = new PotentialDLPOLY(space);
        potentialDLPOLY.setConfigDLPOLY(configDLPOLY);
        potentialMaster.addPotential(potentialDLPOLY, new Species[0]);
        //mcMoveMoleculeCoupled.setPotential(potentialMaster.getPotential(new AtomType[]{new S}));
        
        integrator.setBox(box);

        
    } //end of constructor
    
   
    public static void main(String[] args) {
    	
    	int numMolecules = 32;
    	double temperature = Kelvin.UNIT.toSim(10);
    	long simSteps = 10000;
    	
        String filename = "SimCalcS_Paracetamol_Orthorhombic_"+ Kelvin.UNIT.fromSim(temperature)+"K";
        if (args.length > 0) {
            filename = args[0];
        }
        if (args.length > 1) {
            simSteps = Long.parseLong(args[1]);
        }
        if (args.length > 2) {
            temperature = Kelvin.UNIT.toSim(Double.parseDouble(args[2]));
        }

        System.out.println("Running "+ " Sim Calculate-S Orthorhombic Paracetamol simulation");
        System.out.println(numMolecules + " molecules " +" and temperature "+ Kelvin.UNIT.fromSim(temperature) +"K");
        System.out.println(simSteps+ " steps");
        System.out.println("output data to " + filename);
        
    	MCParacetamolOrthorhombicDLMULTISimCalcS sim = 
        	new MCParacetamolOrthorhombicDLMULTISimCalcS(Space.getInstance(3), numMolecules, temperature);
        
        sim.actionIntegrate.setMaxSteps(simSteps);
        PrimitiveOrthorhombic primitive = sim.primitive;
       
        //Set up Normal-Mode Meter
        MeterNormalMode meterNormalMode = new MeterNormalMode();
        meterNormalMode.setCoordinateDefinition(sim.coordinateDefinition);
        WaveVectorFactorySimple waveVectorFactory = new WaveVectorFactorySimple(primitive);
       
        meterNormalMode.setWaveVectorFactory(waveVectorFactory);
        meterNormalMode.setBox(sim.box);
       
        sim.integrator.addIntervalAction(meterNormalMode);
        sim.integrator.setActionInterval(meterNormalMode, 300);
    
        //Write S-Vectors
        WriteS sWriter = new WriteS();
        sWriter.setFilename(filename);
        sWriter.setMeter(meterNormalMode);
        sWriter.setWaveVectorFactory(waveVectorFactory);
        sWriter.setTemperature(temperature);
        sWriter.setOverwrite(true);
       
        sim.integrator.addIntervalAction(sWriter);
        sim.integrator.setActionInterval(sWriter, 1000);
       
        sim.getController().actionPerformed();
        
        WriteConfiguration writeConfig = new WriteConfiguration();
        writeConfig.setConfName("FinalCoord_SimCalcS_Paracetamol_Orthorhombic_"+Kelvin.UNIT.fromSim(temperature)+"K");
        writeConfig.setBox(sim.box);
        writeConfig.setDoApplyPBC(false);
        writeConfig.actionPerformed();
        
        
        
    }//end of main

    public PotentialMaster potentialMaster;
    public BravaisLatticeCrystal lattice;
    public BoundaryRectangularPeriodic bdry;
    public CoordinateDefinitionParacetamol coordinateDefinition;
    public ActivityIntegrate actionIntegrate;
    public PrimitiveOrthorhombic primitive;
}//end of class