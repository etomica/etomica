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
import etomica.potential.P2Exp6;
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
    public SpeciesParacetamol species;
    public P2Exp6 potentialCC , potentialCHy , potentialHyHy;
    public P2Exp6 potentialCN , potentialNO  , potentialNN  ;
    public P2Exp6 potentialHyN, potentialHyO , potentialOO  ;
    public P2Exp6 potentialCO , potentialHpHp, potentialCHp ;
    public P2Exp6 potentialHpN, potentialOHp , potentialHyHp;
    public Controller controller;

  
    public MCParacetamolOrthorhombicDLMULTI(Space space, int numMolecules, double temperature) {
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

        System.out.println("Running "+ "Orthorhombic Paracetamol simulation");
        System.out.println(numMolecules + " molecules " +" and temperature "+ Kelvin.UNIT.fromSim(temperature) +"K");
        System.out.println(simSteps+ " steps");
        System.out.println("output data to " + filename);
        
    	MCParacetamolOrthorhombicDLMULTI sim = 
        	new MCParacetamolOrthorhombicDLMULTI(Space.getInstance(3), numMolecules, temperature);
       
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
}//end of class