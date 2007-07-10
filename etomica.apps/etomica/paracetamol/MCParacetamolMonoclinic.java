
package etomica.paracetamol;

import java.util.ArrayList;

import etomica.action.activity.ActivityIntegrate;
import etomica.action.activity.Controller;
import etomica.atom.AtomType;
import etomica.atom.AtomTypeGroup;
import etomica.data.DataPump;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.DisplayTextBox;
import etomica.graphics.DisplayBox;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveMolecule;
import etomica.integrator.mcmove.MCMoveRotateMolecule3D;
import etomica.integrator.mcmove.MCMoveStepTracker;
import etomica.lattice.BravaisLattice;
import etomica.lattice.BravaisLatticeCrystal;
import etomica.lattice.crystal.Basis;
import etomica.lattice.crystal.Primitive;
import etomica.lattice.crystal.PrimitiveMonoclinic;
import etomica.box.Box;
import etomica.potential.P2SoftSphericalTruncated;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.BoundaryDeformableLattice;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.units.ElectronVolt;
import etomica.units.Kelvin;
import etomica.units.Pixel;

/**
 * 
 * Three-dimensional soft-sphere MC simulation for paracetamol molecule
 * 
 * Monoclinic Crystal
 * 
 * @author Tai Tan
 *
 */
public class MCParacetamolMonoclinic extends Simulation {

	private static final long serialVersionUID = 1L;
	private final static String APP_NAME = "MC Paracetamol Monoclinic";
    public Box box;
    public IntegratorMC integrator;
    public MCMoveMolecule mcMoveMolecule;
    public MCMoveRotateMolecule3D mcMoveRotateMolecule;
    public SpeciesParacetamol species;
    public P2ElectrostaticDreiding potentialCC , potentialCHy , potentialHyHy;
    public P2ElectrostaticDreiding potentialCN , potentialNO  , potentialNN  ;
    public P2ElectrostaticDreiding potentialHyN, potentialHyO , potentialOO  ;
    public P2ElectrostaticDreiding potentialCO , potentialHpHp, potentialCHp ;
    public P2ElectrostaticDreiding potentialHpN, potentialOHp , potentialHyHp;
    public Controller controller;

  
    public MCParacetamolMonoclinic() {
        this(96);
    }
    
    private MCParacetamolMonoclinic(int numMolecules) {

    	super(Space3D.getInstance(), false);
    	potentialMaster = new PotentialMaster(space);
    	
    	/*
    	 * Monoclinic Crystal
    	 */
        
    	primitive = new PrimitiveMonoclinic(space,  12.119, 8.944, 7.278,  1.744806);
    	//8.944, 12.119, 7.277, 1.74533
    	basis = new BasisMonoclinicParacetamol();
    	lattice = new BravaisLatticeCrystal (primitive, basis);
        
        integrator = new IntegratorMC(this, potentialMaster);
        integrator.setIsothermal(false);
        //integrator.setThermostatInterval(1);
        integrator.setTemperature(Kelvin.UNIT.toSim(20));
        
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
        
        actionIntegrate = new ActivityIntegrate(integrator, false, false);
        //actionIntegrate.setMaxSteps(1000000);
        getController().addAction(actionIntegrate);
        
        
        ConformationParacetamolMonoclinic conformation = new ConformationParacetamolMonoclinic(space);
        species = new SpeciesParacetamol(this);
        ((AtomTypeGroup)species.getMoleculeType()).setConformation(conformation);
        getSpeciesManager().addSpecies(species);
        
        box = new Box(this);
        addBox(box);
        box.setDimensions(Space.makeVector(new double[] {25,25,25}));
        box.setNMolecules(species, numMolecules);        
        
        
        
        /*
         * Intermolecular Potential
         */
    
        double truncationRadiusCC   = 3.0* 3.395524116;
        double truncationRadiusCHy  = 3.0* 2.670105986;
        double truncationRadiusHyHy = 3.0* 2.099665865;
        double truncationRadiusCN   = 3.0* 3.237739512;
        double truncationRadiusNO   = 3.0* 3.035146951;
        double truncationRadiusNN   = 3.0* 3.087286897;
        double truncationRadiusHyN  = 3.0* 2.546030404;
        double truncationRadiusHyO  = 3.0* 2.503031484;
        double truncationRadiusOO   = 3.0* 2.983887553;
        double truncationRadiusCO   = 3.0* 3.183058614;
        double truncationRadiusHpHp = 3.0* 1.543178334;
        double truncationRadiusCHp  = 3.0* 2.289082817;
        double truncationRadiusHpN  = 3.0* 2.182712943;
        double truncationRadiusOHp  = 3.0* 2.145849932;
        double truncationRadiusHyHp = 3.0* 1.800044389;
        
        potentialCC   = new P2ElectrostaticDreiding(space, ElectronVolt.UNIT.toSim(3832.14700),
        		0.277778, ElectronVolt.UNIT.toSim(25.286949));
        potentialCHy  = new P2ElectrostaticDreiding(space, ElectronVolt.UNIT.toSim( 689.53672),
        		0.272480, ElectronVolt.UNIT.toSim( 5.978972));
        potentialHyHy = new P2ElectrostaticDreiding(space, ElectronVolt.UNIT.toSim( 124.07167), 
        		0.267380, ElectronVolt.UNIT.toSim( 1.413698));
        potentialCN   = new P2ElectrostaticDreiding(space, ElectronVolt.UNIT.toSim(3179.51460),
        		0.271003, ElectronVolt.UNIT.toSim(19.006710));
        potentialNO   = new P2ElectrostaticDreiding(space, ElectronVolt.UNIT.toSim(2508.04480),
        		0.258398, ElectronVolt.UNIT.toSim(12.898341));
        potentialNN   = new P2ElectrostaticDreiding(space, ElectronVolt.UNIT.toSim(2638.02850),
        		0.264550, ElectronVolt.UNIT.toSim(14.286224));
        potentialHyN  = new P2ElectrostaticDreiding(space, ElectronVolt.UNIT.toSim( 572.10541),
        		0.265957, ElectronVolt.UNIT.toSim( 4.494041));
        potentialHyO  = new P2ElectrostaticDreiding(space, ElectronVolt.UNIT.toSim( 543.91604),
        		0.259740, ElectronVolt.UNIT.toSim( 4.057452));
        potentialOO   = new P2ElectrostaticDreiding(space, ElectronVolt.UNIT.toSim(2384.46580),
        		0.252525, ElectronVolt.UNIT.toSim(11.645288));
        potentialCO   = new P2ElectrostaticDreiding(space, ElectronVolt.UNIT.toSim(3022.85020), 
        		0.264550, ElectronVolt.UNIT.toSim(17.160239));
        potentialHpHp = new P2ElectrostaticDreiding(space, ElectronVolt.UNIT.toSim(  52.12899), 
        		0.214592, ElectronVolt.UNIT.toSim( 0.222819));
        potentialCHp  = new P2ElectrostaticDreiding(space, ElectronVolt.UNIT.toSim( 446.95185), 
        		0.242131, ElectronVolt.UNIT.toSim( 2.373693));
        potentialHpN  = new P2ElectrostaticDreiding(space, ElectronVolt.UNIT.toSim( 370.83387), 
        		0.236967, ElectronVolt.UNIT.toSim( 1.784166));
        potentialOHp  = new P2ElectrostaticDreiding(space, ElectronVolt.UNIT.toSim( 352.56176), 
        		0.232019, ElectronVolt.UNIT.toSim( 1.610837));
        potentialHyHp = new P2ElectrostaticDreiding(space, ElectronVolt.UNIT.toSim(  80.42221), 
        		0.238095, ElectronVolt.UNIT.toSim( 0.561248));
        
        // CA-CA
        if(truncationRadiusCC > 0.5*box.getBoundary().getDimensions().x(0)) {
            throw new RuntimeException("Truncation radius too large.  " +
            		"Max allowed is"+0.5*box.getBoundary().getDimensions().x(0));
            }
        P2SoftSphericalTruncated interpotentialCC = new P2SoftSphericalTruncated (potentialCC, truncationRadiusCC); 
        potentialMaster.addPotential(interpotentialCC, new AtomType[]{(
        		(AtomFactoryParacetamol)species.getMoleculeFactory()).cType, ((AtomFactoryParacetamol)species.getMoleculeFactory()).cType} );
        
        // CA-HY
        if(truncationRadiusCHy > 0.5*box.getBoundary().getDimensions().x(0)) {
            throw new RuntimeException("Truncation radius too large.  " +
            		"Max allowed is"+0.5*box.getBoundary().getDimensions().x(0));
            }
        P2SoftSphericalTruncated interpotentialCHy = new P2SoftSphericalTruncated (potentialCHy, truncationRadiusCHy); 
        potentialMaster.addPotential(interpotentialCHy, new AtomType[]{(
        		(AtomFactoryParacetamol)species.getMoleculeFactory()).cType, ((AtomFactoryParacetamol)species.getMoleculeFactory()).hyType} );
        
        // HY-HY
        if(truncationRadiusHyHy > 0.5*box.getBoundary().getDimensions().x(0)) {
            throw new RuntimeException("Truncation radius too large.  " +
            		"Max allowed is"+0.5*box.getBoundary().getDimensions().x(0));
            }
        P2SoftSphericalTruncated interpotentialHyHy = new P2SoftSphericalTruncated (potentialHyHy, truncationRadiusHyHy); 
        potentialMaster.addPotential(interpotentialHyHy, new AtomType[]{(
        		(AtomFactoryParacetamol)species.getMoleculeFactory()).hyType, ((AtomFactoryParacetamol)species.getMoleculeFactory()).hyType} );
               
        // CA-NI
        if(truncationRadiusCN > 0.5*box.getBoundary().getDimensions().x(0)) {
            throw new RuntimeException("Truncation radius too large.  " +
            		"Max allowed is"+0.5*box.getBoundary().getDimensions().x(0));
            }
        P2SoftSphericalTruncated interpotentialCN = new P2SoftSphericalTruncated (potentialCN, truncationRadiusCN); 
        potentialMaster.addPotential(interpotentialCN, new AtomType[]{(
        		(AtomFactoryParacetamol)species.getMoleculeFactory()).cType, ((AtomFactoryParacetamol)species.getMoleculeFactory()).nType} );
        
        // NI-OX
        if(truncationRadiusNO > 0.5*box.getBoundary().getDimensions().x(0)) {
            throw new RuntimeException("Truncation radius too large.  " +
            		"Max allowed is"+0.5*box.getBoundary().getDimensions().x(0));
            }
        P2SoftSphericalTruncated interpotentialNO = new P2SoftSphericalTruncated (potentialNO, truncationRadiusNO); 
        potentialMaster.addPotential(interpotentialNO, new AtomType[]{(
        		(AtomFactoryParacetamol)species.getMoleculeFactory()).nType, ((AtomFactoryParacetamol)species.getMoleculeFactory()).oType} );
        
        //NI-NI
        if(truncationRadiusNN > 0.5*box.getBoundary().getDimensions().x(0)) {
            throw new RuntimeException("Truncation radius too large.  " +
            		"Max allowed is"+0.5*box.getBoundary().getDimensions().x(0));
            }
        P2SoftSphericalTruncated interpotentialNN = new P2SoftSphericalTruncated (potentialNN, truncationRadiusNN); 
        potentialMaster.addPotential(interpotentialNN, new AtomType[]{(
        		(AtomFactoryParacetamol)species.getMoleculeFactory()).nType, ((AtomFactoryParacetamol)species.getMoleculeFactory()).nType} );
        
        // HY-NI
        if(truncationRadiusHyN > 0.5*box.getBoundary().getDimensions().x(0)) {
            throw new RuntimeException("Truncation radius too large.  " +
            		"Max allowed is"+0.5*box.getBoundary().getDimensions().x(0));
            }
        P2SoftSphericalTruncated interpotentialHyN = new P2SoftSphericalTruncated (potentialHyN, truncationRadiusHyN); 
        potentialMaster.addPotential(interpotentialHyN, new AtomType[]{(
        		(AtomFactoryParacetamol)species.getMoleculeFactory()).hyType, ((AtomFactoryParacetamol)species.getMoleculeFactory()).nType} );
        
        // HY-OX
        if(truncationRadiusHyO > 0.5*box.getBoundary().getDimensions().x(0)) {
            throw new RuntimeException("Truncation radius too large. " +
            		" Max allowed is"+0.5*box.getBoundary().getDimensions().x(0));
            }
        P2SoftSphericalTruncated interpotentialHyO = new P2SoftSphericalTruncated (potentialHyO, truncationRadiusHyO); 
        potentialMaster.addPotential(interpotentialHyO, new AtomType[]{(
        		(AtomFactoryParacetamol)species.getMoleculeFactory()).hyType, ((AtomFactoryParacetamol)species.getMoleculeFactory()).oType} );
             
        // OX-OX
        if(truncationRadiusOO > 0.5*box.getBoundary().getDimensions().x(0)) {
            throw new RuntimeException("Truncation radius too large. " +
            		" Max allowed is"+0.5*box.getBoundary().getDimensions().x(0));
            }
        P2SoftSphericalTruncated interpotentialOO = new P2SoftSphericalTruncated (potentialOO, truncationRadiusOO); 
        potentialMaster.addPotential(interpotentialOO, new AtomType[]{(
        		(AtomFactoryParacetamol)species.getMoleculeFactory()).oType, ((AtomFactoryParacetamol)species.getMoleculeFactory()).oType} );
        
        // CA-OX
        if(truncationRadiusCO > 0.5*box.getBoundary().getDimensions().x(0)) {
            throw new RuntimeException("Truncation radius too large. " +
            		" Max allowed is"+0.5*box.getBoundary().getDimensions().x(0));
            }
        P2SoftSphericalTruncated interpotentialCO = new P2SoftSphericalTruncated (potentialCO, truncationRadiusCO); 
        potentialMaster.addPotential(interpotentialCO, new AtomType[]{(
        		(AtomFactoryParacetamol)species.getMoleculeFactory()).cType, ((AtomFactoryParacetamol)species.getMoleculeFactory()).oType} );
        
        // HP-HP
        if(truncationRadiusHpHp > 0.5*box.getBoundary().getDimensions().x(0)) {
            throw new RuntimeException("Truncation radius too large. " +
            		" Max allowed is"+0.5*box.getBoundary().getDimensions().x(0));
            }
        P2SoftSphericalTruncated interpotentialHpHp = new P2SoftSphericalTruncated (potentialHpHp, truncationRadiusHpHp); 
        potentialMaster.addPotential(interpotentialHpHp, new AtomType[]{(
        		(AtomFactoryParacetamol)species.getMoleculeFactory()).hpType, ((AtomFactoryParacetamol)species.getMoleculeFactory()).hpType} );
        
        // CA-HP
        if(truncationRadiusCHp > 0.5*box.getBoundary().getDimensions().x(0)) {
            throw new RuntimeException("Truncation radius too large. " +
            		" Max allowed is"+0.5*box.getBoundary().getDimensions().x(0));
            }
        P2SoftSphericalTruncated interpotentialCHp = new P2SoftSphericalTruncated (potentialCHp, truncationRadiusCHp); 
        potentialMaster.addPotential(interpotentialCHp, new AtomType[]{(
        		(AtomFactoryParacetamol)species.getMoleculeFactory()).cType, ((AtomFactoryParacetamol)species.getMoleculeFactory()).hpType} );
               
        // HP-NI
        if(truncationRadiusHpN > 0.5*box.getBoundary().getDimensions().x(0)) {
            throw new RuntimeException("Truncation radius too large.  " +
            		"Max allowed is"+0.5*box.getBoundary().getDimensions().x(0));
            }
        P2SoftSphericalTruncated interpotentialHpN = new P2SoftSphericalTruncated (potentialHpN, truncationRadiusHpN); 
        potentialMaster.addPotential(interpotentialHpN, new AtomType[]{(
        		(AtomFactoryParacetamol)species.getMoleculeFactory()).hpType, ((AtomFactoryParacetamol)species.getMoleculeFactory()).nType} );
        
        // OX-HP
        if(truncationRadiusOHp > 0.5*box.getBoundary().getDimensions().x(0)) {
            throw new RuntimeException("Truncation radius too large. " +
            		" Max allowed is"+0.5*box.getBoundary().getDimensions().x(0));
            }
        P2SoftSphericalTruncated interpotentialOHp = new P2SoftSphericalTruncated (potentialOHp, truncationRadiusOHp); 
        potentialMaster.addPotential(interpotentialOHp, new AtomType[]{(
        		(AtomFactoryParacetamol)species.getMoleculeFactory()).oType, ((AtomFactoryParacetamol)species.getMoleculeFactory()).hpType} );
        
        // HY-HP
        if(truncationRadiusHyHp > 0.5*box.getBoundary().getDimensions().x(0)) {
            throw new RuntimeException("Truncation radius too large.  " +
            		"Max allowed is"+0.5*box.getBoundary().getDimensions().x(0));
            }
        P2SoftSphericalTruncated interpotentialHyHp = new P2SoftSphericalTruncated (potentialHyHp, truncationRadiusHyHp); 
        potentialMaster.addPotential(interpotentialHyHp, new AtomType[]{(
        		(AtomFactoryParacetamol)species.getMoleculeFactory()).hyType, ((AtomFactoryParacetamol)species.getMoleculeFactory()).hpType} );
  
        potentialMaster.lrcMaster().setEnabled(false);
       /*
        *
        */
        
        bdry =  new BoundaryDeformableLattice( primitive, getRandom(), new int []{2, 3, 4});
        //bdry.setDimensions(Space.makeVector(new double []{2*12.119, 3*8.944, 4*7.278}));
        box.setBoundary(bdry);
        CoordinateDefinitionParacetamol coordDef = new CoordinateDefinitionParacetamol(box, primitive, basis);
        coordDef.setBasisMonoclinic();
        coordDef.initializeCoordinates(new int []{2, 3, 4});
        
        integrator.setBox(box);
        
    } //end of constructor
    
   
    public static void main(String[] args) {
    	int numMolecules = 96;
        etomica.paracetamol.MCParacetamolMonoclinic sim = new etomica.paracetamol.MCParacetamolMonoclinic(numMolecules);
        SimulationGraphic simGraphic = new SimulationGraphic(sim, APP_NAME);
        //sim.getController().actionPerformed();
        
   /*****************************************************************************/    
        
        MeterPotentialEnergy meterPE = new MeterPotentialEnergy(sim.potentialMaster);
        meterPE.setBox(sim.box);
        DisplayTextBox PEbox = new DisplayTextBox();
        DataPump PEpump = new DataPump(meterPE, PEbox);
        ArrayList dataStreamPumps = simGraphic.getController().getDataStreamPumps();
        dataStreamPumps.add(PEpump);
        
        sim.integrator.addIntervalAction(PEpump);
        sim.integrator.setActionInterval(PEpump, 500);
        
 /**********************************************************************/   
        simGraphic.add(PEbox);
        
        
        simGraphic.makeAndDisplayFrame(APP_NAME);
        simGraphic.getDisplayBox(sim.box).setPixelUnit(new Pixel(10));
        ColorSchemeByType colorScheme = ((ColorSchemeByType)((DisplayBox)simGraphic.displayList().getFirst()).getColorScheme());
        AtomTypeGroup atomType = (AtomTypeGroup)sim.species.getMoleculeType();
        colorScheme.setColor(atomType.getChildTypes()[0], java.awt.Color.red);
        colorScheme.setColor(atomType.getChildTypes()[1], java.awt.Color.gray);
        colorScheme.setColor(atomType.getChildTypes()[2], java.awt.Color.blue);
        colorScheme.setColor(atomType.getChildTypes()[3], java.awt.Color.white);
        colorScheme.setColor(atomType.getChildTypes()[4], java.awt.Color.white);
        
        simGraphic.getDisplayBox(sim.box).repaint();
        
    }//end of main

    public PotentialMaster potentialMaster;
    public BravaisLattice lattice;
    public BoundaryDeformableLattice bdry;
    public CoordinateDefinitionParacetamol coordinateDefinition;
    public ActivityIntegrate actionIntegrate;
    public Primitive primitive;
    public Basis basis;
}//end of class