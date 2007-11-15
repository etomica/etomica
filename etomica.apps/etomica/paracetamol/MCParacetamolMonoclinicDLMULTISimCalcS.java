
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
import etomica.lattice.crystal.PrimitiveMonoclinic;
import etomica.normalmode.MeterNormalMode;
import etomica.normalmode.WaveVectorFactorySimple;
import etomica.normalmode.WriteS;
import etomica.normalmode.CoordinateDefinition.BasisCell;
import etomica.potential.PotentialDLPOLY;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.BoundaryDeformableLattice;
import etomica.space.Space;
import etomica.species.Species;
import etomica.units.Kelvin;

/**
 * 
 * Three-dimensional soft-sphere MC simulation for paracetamol molecule
 *  using the potential model from DL_MULTI package
 *  with inclusion of multipole ewald-summation
 * 
 * MC simulation of Paracetamol molecules in Form I (Monoclinic) with 
 *  tabulation of the collective-coordinate S-matrix.
 * No graphic display
 *  
 * Monoclinic Crystal
 * 
 * @author Tai Tan
 *
 */
public class MCParacetamolMonoclinicDLMULTISimCalcS extends Simulation{
	private static final long serialVersionUID = 1L;
	public Box box;
    public IntegratorMC integrator;
    public MCMoveMoleculeCoupledDLPOLY mcMoveMoleculeCoupledDLPOLY;
    public MCMoveRotateMolecule3D mcMoveRotateMolecule;
    public SpeciesParacetamol species;
    public Controller controller;

  
    public MCParacetamolMonoclinicDLMULTISimCalcS(Space space, int numMolecules, double temperature) {
    	super(space, true);
    	potentialMaster = new PotentialMaster(space);
    	
    	/*
    	 * Monoclinic Crystal
    	 */
        
    	primitive = new PrimitiveMonoclinic(space,  12.119, 8.944, 7.278,  1.744806);
    	//8.944, 12.119, 7.277, 1.74533
    	BasisMonoclinicParacetamol basis = new BasisMonoclinicParacetamol();
    	lattice = new BravaisLatticeCrystal (primitive, basis);
        
        integrator = new IntegratorMC(this, potentialMaster);
        integrator.setIsothermal(false);
        //integrator.setThermostatInterval(1);
        integrator.setTemperature(Kelvin.UNIT.toSim(123));
        
        mcMoveRotateMolecule = new MCMoveRotateMolecule3D(potentialMaster, random);
        mcMoveRotateMolecule.setStepSize(0.068);
        ((MCMoveStepTracker)mcMoveRotateMolecule.getTracker()).setNoisyAdjustment(true);
      
        mcMoveMoleculeCoupledDLPOLY = new MCMoveMoleculeCoupledDLPOLY(potentialMaster, getRandom());
        
        integrator.getMoveManager().addMCMove(mcMoveMoleculeCoupledDLPOLY);
        integrator.getMoveManager().addMCMove(mcMoveRotateMolecule);
        integrator.getMoveManager().setEquilibrating(true);
        
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
        
        bdry =  new BoundaryDeformableLattice( primitive, getRandom(), new int []{2, 2, 2});
        bdry.setDimensions(Space.makeVector(new double []{2*12.119, 2*8.944, 2*7.278}));
        box.setBoundary(bdry);
        
        //ConfigurationFile configFile = new ConfigurationFile("FinalCoord_Min_Paracetamol_Monoclinic");
        coordinateDefinition = new CoordinateDefinitionParacetamol(box, primitive, basis);
        coordinateDefinition.setBasisMonoclinic();
        //coordinateDefinition.setConfiguration(configFile);
        coordinateDefinition.initializeCoordinates(new int []{2, 2, 2});
        double[] u = new double[]{0.19923028993600184, -0.10831241063851138, -0.3374206766423242, -0.056318915244514295, -0.08373011094015517, -0.20967425215989952, -0.15036406963107662, -0.06903390627642114, 0.2911023015207981, -0.062482609873595024, -0.0836259970912052, -0.17490727322325597, -0.19043886713958358, 0.09585326949021145, 0.29982577023219825, -0.052978731871725596, 0.0794450448846585, 0.20078353728718995, 0.1791946947288432, 0.07473598884040555, -0.33875779999181965, -0.06542404448301108, 0.0856334706980024, 0.20954733660556393};
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
        
        
        integrator.setBox(box);
    } //end of constructor
    
    public static void main(String[] args) {
    	
    	int numMolecules = 32;
    	double temperature = Kelvin.UNIT.toSim(123);
    	long simSteps = 10000;
      
        String filename = "SimCalcS_Paracetamol_Monoclinic_"+ Kelvin.UNIT.fromSim(temperature)+"K";
        if (args.length > 0) {
            filename = args[0];
        }
        if (args.length > 1) {
            simSteps = Long.parseLong(args[1]);
        }
        if (args.length > 2) {
            temperature = Kelvin.UNIT.toSim(Double.parseDouble(args[2]));
        }
        
        System.out.println("Running "+ "Sim Calculate-S Monoclinic Paracetamol simulation");
        System.out.println(numMolecules + " molecules " +" and temperature "+ Kelvin.UNIT.fromSim(temperature) +"K");
        System.out.println(simSteps+ " steps");
        System.out.println("output data to " + filename);
        
    	MCParacetamolMonoclinicDLMULTISimCalcS sim = 
        	new MCParacetamolMonoclinicDLMULTISimCalcS(Space.getInstance(3), numMolecules, temperature);

        sim.actionIntegrate.setMaxSteps(simSteps);
        PrimitiveMonoclinic primitive = sim.primitive;
        
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
        writeConfig.setConfName("FinalCoord_SimCalcS_Paracetamol_Monoclinic_"+Kelvin.UNIT.fromSim(temperature)+"K");
        writeConfig.setBox(sim.box);
        writeConfig.setDoApplyPBC(false);
        writeConfig.actionPerformed();
        
    }//end of main
    
    public PotentialMaster potentialMaster;
    public BravaisLatticeCrystal lattice;
    public BoundaryDeformableLattice bdry;
    public CoordinateDefinitionParacetamol coordinateDefinition;
    public ActivityIntegrate actionIntegrate;
    public PrimitiveMonoclinic primitive;

}//end of class