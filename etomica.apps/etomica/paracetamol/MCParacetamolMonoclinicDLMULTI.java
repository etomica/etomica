
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
import etomica.lattice.crystal.PrimitiveMonoclinic;
import etomica.potential.P2Exp6;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.BoundaryDeformableLattice;
import etomica.space.Space;
import etomica.units.Kelvin;

/**
 * 
 * Three-dimensional soft-sphere MC simulation for paracetamol molecule
 *  using the potential model from DL_MULTI package
 *  with inclusion of multipole ewald-summation
 *  
 * To get the lowest minimum cofiguration with MCMoveHarmonicStep
 * 
 * Monoclinic Crystal
 * 
 * @author Tai Tan
 *
 */
public class MCParacetamolMonoclinicDLMULTI extends Simulation{
	private static final long serialVersionUID = 1L;
//	private final static String APP_NAME = "MC Paracetamol Monoclinic";
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

  
    public MCParacetamolMonoclinicDLMULTI(Space space, int numMolecules, double temperature) {
    	super(space, true);
    	potentialMaster = new PotentialMaster(space);
    	
    	/*
    	 * Monoclinic Crystal
    	 */
        
    	PrimitiveMonoclinic primitive = new PrimitiveMonoclinic(space,  12.119, 8.944, 7.278,  1.744806);
    	//8.944, 12.119, 7.277, 1.74533
    	BasisMonoclinicParacetamol basis = new BasisMonoclinicParacetamol();
    	lattice = new BravaisLatticeCrystal (primitive, basis);
        
        integrator = new IntegratorMC(this, potentialMaster);
        integrator.setIsothermal(false);
        //integrator.setThermostatInterval(1);
        integrator.setTemperature(Kelvin.UNIT.toSim(20));
        

        integrator.getMoveManager().setEquilibrating(true);

        mcMoveMolecule = new MCMoveMolecule(this, potentialMaster);
        mcMoveMolecule.setStepSize(0.1747);  //Step size to input
        ((MCMoveStepTracker)mcMoveMolecule.getTracker()).setTunable(true);
        ((MCMoveStepTracker)mcMoveMolecule.getTracker()).setNoisyAdjustment(true);
        
        mcMoveRotateMolecule = new MCMoveRotateMolecule3D(potentialMaster, random);
        mcMoveRotateMolecule.setStepSize(0.068);
        ((MCMoveStepTracker)mcMoveRotateMolecule.getTracker()).setNoisyAdjustment(true);
        
        
        integrator.getMoveManager().setEquilibrating(true);
        integrator.getMoveManager().addMCMove(mcMoveMolecule);
        integrator.getMoveManager().addMCMove(mcMoveRotateMolecule);
        
        actionIntegrate = new ActivityIntegrate(integrator, 0, false);
        getController().addAction(actionIntegrate);
        
        ConformationParacetamolMonoclinic conformation = new ConformationParacetamolMonoclinic(space);
        species = new SpeciesParacetamol(this);
        ((AtomTypeGroup)species.getMoleculeType()).setConformation(conformation);
        getSpeciesManager().addSpecies(species);
        
        box = new Box(this);
        addBox(box);
        box.setDimensions(Space.makeVector(new double[] {25,25,25}));
        box.setNMolecules(species, numMolecules);        
        
        bdry =  new BoundaryDeformableLattice( primitive, getRandom(), new int []{2, 3, 4});
        bdry.setDimensions(Space.makeVector(new double []{2*12.119, 3*8.944, 4*7.278}));
        box.setBoundary(bdry);
        
        CoordinateDefinitionParacetamol coordDef = new CoordinateDefinitionParacetamol(box, primitive, basis);
        coordDef.setBasisMonoclinic();
        coordDef.initializeCoordinates(new int []{2, 3, 4});
        
        WriteConfigurationDLPOLY configDLPOLY = new WriteConfigurationDLPOLY();
        configDLPOLY.setConfName("CONFIG");
        configDLPOLY.setBox(box);
        configDLPOLY.getElementHash().put(HydrogenP.INSTANCE, "HP");
       // configDLPOLY.actionPerformed();
       // System.exit(1);
   
       //PotentialDLPOLY potentialDLPOLY = new PotentialDLPOLY(space);
       //potentialDLPOLY.setConfigDLPOLY(configDLPOLY);
       //potentialMaster.addPotential(potentialDLPOLY, new Species[0]);

        integrator.setBox(box);
    } //end of constructor
    
    public static void main(String[] args) {
    	
    	int numMolecules = 96;
    	double temperature = Kelvin.UNIT.toSim(123);
    	long simSteps = 100;
      
        String filename = "Paracetamol_Monoclinic_"+ Kelvin.UNIT.fromSim(temperature)+"K";
        if (args.length > 0) {
            filename = args[0];
        }
        if (args.length > 1) {
            simSteps = Long.parseLong(args[1]);
        }
        if (args.length > 2) {
            temperature = Kelvin.UNIT.toSim(Double.parseDouble(args[2]));
        }
        
        System.out.println("Running "+ "Monoclinic Paracetamol simulation");
        System.out.println(numMolecules + " molecules " +" and temperature "+ Kelvin.UNIT.fromSim(temperature) +"K");
        System.out.println(simSteps+ " steps");
        System.out.println("output data to " + filename);
        
    	MCParacetamolMonoclinicDLMULTI sim = 
        	new MCParacetamolMonoclinicDLMULTI(Space.getInstance(3), numMolecules, temperature);

        sim.actionIntegrate.setMaxSteps(simSteps);
        MeterPotentialEnergy meterPE = new MeterPotentialEnergy(sim.potentialMaster);
        meterPE.setBox(sim.box);
        
        DataLogger dataLoggerPE = new DataLogger();
        dataLoggerPE.setWriteInterval(5);
        dataLoggerPE.setFileName("Paracetamol_Monoclinic_"+Kelvin.UNIT.fromSim(temperature)+"K");
        dataLoggerPE.setDataSink(new DataTableWriter());
        dataLoggerPE.setAppending(true);
        dataLoggerPE.setCloseFileEachTime(true);
        
        DataPump PEpump = new DataPump(meterPE, dataLoggerPE);
        sim.getController().getEventManager().addListener(dataLoggerPE);
        
        sim.integrator.addIntervalAction(PEpump);
        sim.integrator.setActionInterval(PEpump, 1);
        
        
        sim.getController().actionPerformed();
        
        WriteConfiguration writeConfig = new WriteConfiguration();
        writeConfig.setConfName("FinalCoord_Paracetamol_Monoclinic_"+Kelvin.UNIT.fromSim(temperature)+"K");
        writeConfig.setBox(sim.box);
        writeConfig.setDoApplyPBC(false);
        writeConfig.actionPerformed();
         
        PDBWriter pdbWriter = new PDBWriter(sim.box);
        pdbWriter.setFileName("GraphicView_Monoclinic_"+ Kelvin.UNIT.fromSim(temperature)+"K.pdb");
        pdbWriter.actionPerformed();
        
    }//end of main
    
    public PotentialMaster potentialMaster;
    public BravaisLatticeCrystal lattice;
    public BoundaryDeformableLattice bdry;
    public ActivityIntegrate actionIntegrate;

}//end of class