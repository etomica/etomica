
package etomica.paracetamol;

import etomica.action.activity.ActivityIntegrate;
import etomica.atom.AtomType;
import etomica.atom.AtomTypeGroup;
import etomica.data.DataPump;
import etomica.data.meter.MeterEnergy;
import etomica.data.meter.MeterKineticEnergy;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.meter.MeterTemperature;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.DisplayBox;
import etomica.graphics.DisplayPhase;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.integrator.IntervalActionAdapter;
import etomica.lattice.BravaisLattice;
import etomica.lattice.BravaisLatticeCrystal;
import etomica.lattice.crystal.PrimitiveMonoclinic;
import etomica.nbr.CriterionInterMolecular;
import etomica.nbr.CriterionNone;
import etomica.nbr.list.NeighborListManager;
import etomica.nbr.list.PotentialMasterList;
import etomica.phase.Phase;
import etomica.potential.P2Dreiding;
import etomica.potential.P2Exp6;
import etomica.potential.P2SoftSphericalTruncated;
import etomica.potential.P3BondAngleDreiding;
import etomica.potential.P4TorsionDreiding;
import etomica.potential.PotentialGroup;
import etomica.simulation.Simulation;
import etomica.space.BoundaryDeformableLattice;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.Species;
import etomica.units.Kelvin;
import etomica.units.Pixel;
import etomica.util.Default;

/**
 * 
 * Three-dimensional soft-sphere molecular dynamics simulation for paracetamol molecule, using
 * neighbor listing.  
 * 
 * Monoclinic Crystal
 * 
 * @author Tai Tan
 *
 */
public class MDParacetamolMonoclinic extends Simulation {

	private static final String APP_NAME = "MD Paracetamol Monoclinic";
	private static final int PIXEL_SIZE = 15;

    /**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	//the following fields are made accessible for convenience to permit simple
    //mutation of the default behavior

    public final PotentialMasterList potentialMaster;
    
    /**
     * The Phase holding the atoms. 
     */
    public final Phase phase;
    /**
     * The Integrator performing the dynamics.
     */
    public final IntegratorVelocityVerlet integrator;
    /**
     * The single soft-sphere species.
     */
    public final SpeciesParacetamol species;

    
    /**
     * Sole public constructor, makes a simulation using a 3D space.
     */
    public MDParacetamolMonoclinic() {
        this(new Default());
    }
    
    public MDParacetamolMonoclinic(Default defaults) {
        this(Space3D.getInstance(), defaults);
    }
    
    /* 
     * we use a second, private constructor to permit the space to
     * appear twice in the call to the superclass constructor; alternatively
     * we could have passed Space3D.getInstance() twice
     *
     */
    private MDParacetamolMonoclinic(Space space, Default defaults) {

        /*
         * invoke the superclass constructor
         *	"true" is indicating to the superclass that this is a dynamic simulation
         * the PotentialMaster is selected such as to implement neighbor listing
         */
    	super(space, true, Default.BIT_LENGTH, defaults);
        
        potentialMaster = new PotentialMasterList(this, 1.6);

    	/*
    	 * Monoclinic Crystal
    	 */
        
    	PrimitiveMonoclinic primitive = new PrimitiveMonoclinic(space,  12.119, 8.944, 7.278,  1.744806);
    	//8.944, 12.119, 7.277, 1.74533
    	BasisMonoclinicParacetamol basis = new BasisMonoclinicParacetamol();
    	lattice = new BravaisLatticeCrystal (primitive, basis);
    	configMonoLattice = new ConfigurationMonoclinicLattice(lattice);
    	
        double neighborRangeFac = 1.6;
        potentialMaster.setRange(neighborRangeFac*defaults.atomSize);

        integrator = new IntegratorVelocityVerlet(this, potentialMaster);
        integrator.setIsothermal(false);
        //integrator.setThermostatInterval(1);
        integrator.setTimeStep(0.001); //1 = pico sec
        this.register(integrator);
        integrator.setTemperature(Kelvin.UNIT.toSim(20));

        ActivityIntegrate activityIntegrate = new ActivityIntegrate(this,integrator);
        //activityIntegrate.setMaxSteps(100000);
        activityIntegrate.setDoSleep(true);
        activityIntegrate.setSleepPeriod(1);
        getController().addAction(activityIntegrate);
        
        phase = new Phase(this);
        addPhase(phase);
        phase.setDimensions(Space.makeVector(new double[] {25,25,25}));
        species = new SpeciesParacetamol(this);
        getSpeciesManager().addSpecies(species);
        species.getAgent(phase).setNMolecules(96);
        
        NeighborListManager nbrManager = potentialMaster.getNeighborManager(phase);
        nbrManager.setRange(defaults.atomSize*1.6);
        integrator.addListener(nbrManager);
               
        PotentialGroup intramolecularpotential = potentialMaster.makePotentialGroup(1);
        potentialMaster.addPotential(intramolecularpotential, new Species[]{species});
        
        /*
         *  Bond Stretch Potential             
         * 
      	 *	Equalibrium Radius [unit Amstrom]; Pre-factor [unit Kelvin]
         */
 
        P2Dreiding potentialC4O1 = new P2Dreiding(space, 1.352 , .17612581e6); 
        intramolecularpotential.addPotential(potentialC4O1, 
        		new ApiIndexList(new int[][] {{AtomParacetamol.indexC4,AtomParacetamol.indexO1}}));
        
        P2Dreiding potentialring1 = new P2Dreiding(space, 1.395, .264188715e6);
        intramolecularpotential.addPotential(potentialring1, 
        		new ApiIndexList(new int[][]{{AtomParacetamol.indexC1,AtomParacetamol.indexC6},
        									 {AtomParacetamol.indexC1,AtomParacetamol.indexC2},
        									 {AtomParacetamol.indexC3,AtomParacetamol.indexC4},
        									 {AtomParacetamol.indexC4,AtomParacetamol.indexC5},
        		}));
      
        P2Dreiding potentialring2 = new P2Dreiding(space, 1.385, .264188715e6); 
        intramolecularpotential.addPotential(potentialring2, 
        		new ApiIndexList(new int[][]{{AtomParacetamol.indexC2,AtomParacetamol.indexC3},
        									 {AtomParacetamol.indexC5,AtomParacetamol.indexC6}
        		}));
        
        P2Dreiding potentialC7C8 = new P2Dreiding(space, 1.503, .17612581e6); 
        intramolecularpotential.addPotential(potentialC7C8, 
        		new ApiIndexList(new int[][]{{AtomParacetamol.indexC7,AtomParacetamol.indexC8}}));
        
        P2Dreiding potentialC1N1 = new P2Dreiding(space, 1.394, .17612581e6); 
        intramolecularpotential.addPotential(potentialC1N1, 
        		new ApiIndexList(new int[][]{{AtomParacetamol.indexC1,AtomParacetamol.indexN}}));
 
        P2Dreiding potentialC7N1 = new P2Dreiding(space, 1.366, .17612581e6); 
        intramolecularpotential.addPotential(potentialC7N1, 
        		new ApiIndexList(new int[][]{{AtomParacetamol.indexC7,AtomParacetamol.indexN}}));

        P2Dreiding potentialC7O2 = new P2Dreiding(space, 1.226 , .35225162e6);
        intramolecularpotential.addPotential(potentialC7O2, 
        		new ApiIndexList(new int[][]{{AtomParacetamol.indexC7,AtomParacetamol.indexO2}}));

      /*
       * Bond Angle Potential
       * 
       * Equilibrium Angle [unit radians]; Pre-factor [unit Kelvin]
       */

        P3BondAngleDreiding potential120 = new P3BondAngleDreiding(space, 2.094395102, .3354777333e5);
        intramolecularpotential.addPotential(potential120, 
        		new AtomsetIteratorIndexList(new int[][]{{AtomParacetamol.indexC3,AtomParacetamol.indexC4,AtomParacetamol.indexC5},
        											     {AtomParacetamol.indexC4,AtomParacetamol.indexC5,AtomParacetamol.indexC6},
        											     {AtomParacetamol.indexC5,AtomParacetamol.indexC6,AtomParacetamol.indexC1},
        											     {AtomParacetamol.indexC6,AtomParacetamol.indexC1,AtomParacetamol.indexC2},
        											     {AtomParacetamol.indexC1,AtomParacetamol.indexC2,AtomParacetamol.indexC3},
        											     {AtomParacetamol.indexC2,AtomParacetamol.indexC3,AtomParacetamol.indexC4}									     				     
        		}));
        
        P3BondAngleDreiding potentialO1C4C5 = new P3BondAngleDreiding(space, 2.14675498, .3354777333e5); 
        intramolecularpotential.addPotential(potentialO1C4C5, 
        		new AtomsetIteratorIndexList(new int[][]{{AtomParacetamol.indexO1,AtomParacetamol.indexC4,AtomParacetamol.indexC5}}));
       
        P3BondAngleDreiding potentialO1C4C3 = new P3BondAngleDreiding(space, 2.042035225, .3354777333e5); 
        intramolecularpotential.addPotential(potentialO1C4C3, 
        		new AtomsetIteratorIndexList(new int[][]{{AtomParacetamol.indexO1,AtomParacetamol.indexC4,AtomParacetamol.indexC3}}));
        
        P3BondAngleDreiding potentialC6C1N1 = new P3BondAngleDreiding(space, 2.059488517, .3354777333e5); 
        intramolecularpotential.addPotential(potentialC6C1N1, 
        		new AtomsetIteratorIndexList(new int[][]{{AtomParacetamol.indexC6,AtomParacetamol.indexC1,AtomParacetamol.indexN }}));
  
        P3BondAngleDreiding potentialC2C1N1 = new P3BondAngleDreiding(space, 2.129301687, .3354777333e5); 
        intramolecularpotential.addPotential(potentialC2C1N1, 
        		new AtomsetIteratorIndexList(new int[][]{{AtomParacetamol.indexC2,AtomParacetamol.indexC1,AtomParacetamol.indexN }}));
  
        P3BondAngleDreiding potentialN1C7C8 = new P3BondAngleDreiding(space, 1.989675347, .3354777333e5); 
        intramolecularpotential.addPotential(potentialN1C7C8, 
        		new AtomsetIteratorIndexList(new int[][]{{AtomParacetamol.indexN ,AtomParacetamol.indexC7,AtomParacetamol.indexC8}}));
  
        P3BondAngleDreiding potentialN1C7O2 = new P3BondAngleDreiding(space, 2.181661565, .3354777333e5); 
        intramolecularpotential.addPotential(potentialN1C7O2, 
        		new AtomsetIteratorIndexList(new int[][]{{AtomParacetamol.indexN ,AtomParacetamol.indexC7,AtomParacetamol.indexO2}}));
       
        P3BondAngleDreiding potentialO2C7C8 = new P3BondAngleDreiding(space, 2.111848395, .3354777333e5); 
        intramolecularpotential.addPotential(potentialO2C7C8, 
        		new AtomsetIteratorIndexList(new int[][]{{AtomParacetamol.indexO2,AtomParacetamol.indexC7,AtomParacetamol.indexC8}}));
      
        P3BondAngleDreiding potentialC1N1C7 = new P3BondAngleDreiding(space, 2.251474735, .2742552176e5); 
        intramolecularpotential.addPotential(potentialC1N1C7, 
        		new AtomsetIteratorIndexList(new int[][]{{AtomParacetamol.indexC1,AtomParacetamol.indexN ,AtomParacetamol.indexC7}}));
  
        /*
         * Torsional Angle Potential
         * 
         * Equilibrium Phi [unit radians]; Pre-factor [unit Kelvin]; Periodicity [integer]
         */
   
        P4TorsionDreiding potentialCRCR = new P4TorsionDreiding(space, 3.141592654, 6290.2075, 2);
        intramolecularpotential.addPotential(potentialCRCR, 
        		new Atomset4IteratorIndexList(new int[][]{{AtomParacetamol.indexO1, AtomParacetamol.indexC4,
        												   AtomParacetamol.indexC5, AtomParacetamol.indexC6 },
        												  {AtomParacetamol.indexO1, AtomParacetamol.indexC4,
        												   AtomParacetamol.indexC3, AtomParacetamol.indexC2 },
        												  {AtomParacetamol.indexC4, AtomParacetamol.indexC5,
            											   AtomParacetamol.indexC6, AtomParacetamol.indexC1 },
            											  {AtomParacetamol.indexC4, AtomParacetamol.indexC3,
            											   AtomParacetamol.indexC2, AtomParacetamol.indexC1 },
            											  {AtomParacetamol.indexC5, AtomParacetamol.indexC6,
            											   AtomParacetamol.indexC1, AtomParacetamol.indexN  },
            											  {AtomParacetamol.indexC3, AtomParacetamol.indexC2,
            											   AtomParacetamol.indexC1, AtomParacetamol.indexN  }
        		
        		}));

        P4TorsionDreiding potentialCRN3 = new P4TorsionDreiding(space, 0.0, 251.6083, 6);
        intramolecularpotential.addPotential(potentialCRN3, 
        		new Atomset4IteratorIndexList(new int[][]{{AtomParacetamol.indexC6, AtomParacetamol.indexC1,
        												   AtomParacetamol.indexN , AtomParacetamol.indexC7 },
        												  {AtomParacetamol.indexC2, AtomParacetamol.indexC1,
            											   AtomParacetamol.indexN , AtomParacetamol.indexC7 }
        		}));

        
        P4TorsionDreiding potentialN3C3 = new P4TorsionDreiding(space, 3.141592654, 503.2166, 3);
        intramolecularpotential.addPotential(potentialN3C3, 
        		new Atomset4IteratorIndexList(new int[][]{{AtomParacetamol.indexC1, AtomParacetamol.indexN ,
        												   AtomParacetamol.indexC7, AtomParacetamol.indexO2  },
        												  {AtomParacetamol.indexC1, AtomParacetamol.indexN ,
        												   AtomParacetamol.indexC7, AtomParacetamol.indexC8  }
        		}));
       
        
        /*
         * Intra-nonbonded Potential
         * 
         * Equilibrium Radius [unit Amstrom]; Pre-factor [unit Kelvin]
         */
   
        P2Exp6 potentialCC = new P2Exp6(space, 3832.147000*11604.45728, 0.277778, 25.286949*11604.45728);
        intramolecularpotential.addPotential(potentialCC, 
        		new ApiIndexList(new int[][]{{AtomParacetamol.indexC4,AtomParacetamol.indexC7},
        									 {AtomParacetamol.indexC4,AtomParacetamol.indexC8},
        									 {AtomParacetamol.indexC5,AtomParacetamol.indexC7},
        									 {AtomParacetamol.indexC5,AtomParacetamol.indexC8},
        									 {AtomParacetamol.indexC3,AtomParacetamol.indexC7},
        									 {AtomParacetamol.indexC3,AtomParacetamol.indexC8},
        									 {AtomParacetamol.indexC6,AtomParacetamol.indexC8},
        									 {AtomParacetamol.indexC2,AtomParacetamol.indexC8}
        		}));

        P2Exp6 potentialCO = new P2Exp6(space, 3022.850200*11604.45728, 0.264550, 17.160239*11604.45728);
        intramolecularpotential.addPotential(potentialCO, 
        		new ApiIndexList(new int[][]{{AtomParacetamol.indexC1,AtomParacetamol.indexO1},
        									 {AtomParacetamol.indexC7,AtomParacetamol.indexO1},
        									 {AtomParacetamol.indexC8,AtomParacetamol.indexO1},
        									 {AtomParacetamol.indexC4,AtomParacetamol.indexO2},
        									 {AtomParacetamol.indexC5,AtomParacetamol.indexO2},
        									 {AtomParacetamol.indexC3,AtomParacetamol.indexO2},
        									 {AtomParacetamol.indexC6,AtomParacetamol.indexO2},
        									 {AtomParacetamol.indexC2,AtomParacetamol.indexO2}
        		}));
        
        P2Exp6 potentialON = new P2Exp6(space, 2508.044800*11604.45728, 0.258398, 12.898341*11604.45728);
        intramolecularpotential.addPotential(potentialON, 
        		new ApiIndexList(new int[][]{{AtomParacetamol.indexO1,AtomParacetamol.indexN }}));
        
        P2Exp6 potentialCN = new P2Exp6(space, 3179.514600*11604.45728, 0.271003, 19.006710*11604.45728);
        intramolecularpotential.addPotential(potentialCN, 
        		new ApiIndexList(new int[][]{{AtomParacetamol.indexC4,AtomParacetamol.indexN }}));
        
        P2Exp6 potentialO1O2 = new P2Exp6(space, 2384.465800*11604.45728, 0.252525, 11.645288*11604.45728);
        intramolecularpotential.addPotential(potentialO1O2, 
        		new ApiIndexList(new int[][]{{AtomParacetamol.indexO1,AtomParacetamol.indexO2}}));
       
      
        
        /*
         * Intermolecular Potential
         */
    
        double truncationRadiusCC = 3.0* 3.472990473;
        double truncationRadiusCO = 3.0* 3.253072125;
        double truncationRadiusCN = 3.0* 3.369296217;
        double truncationRadiusON = 3.0* 3.149377868;
        double truncationRadiusOO = 3.0* 3.033153776;
        double truncationRadiusNN = 3.0* 3.262560196;
        
        if(truncationRadiusCC > 0.5*phase.getBoundary().getDimensions().x(0)) {
            throw new RuntimeException("Truncation radius too large.  Max allowed is"+0.5*phase.getBoundary().getDimensions().x(0));
            }
        potentialMaster.setCellRange(2);
        potentialMaster.setRange(1.2*truncationRadiusCC);
        P2SoftSphericalTruncated interpotentialCC = new P2SoftSphericalTruncated (new P2ElectrostaticDreiding(space, 3832.14700*11604.45728, 0.277778, 25.286949*11604.45728), truncationRadiusCC); 
        potentialMaster.addPotential(interpotentialCC, new AtomType[]{((AtomFactoryParacetamol)species.getFactory()).cType, ((AtomFactoryParacetamol)species.getFactory()).cType} );
        
        if(truncationRadiusCO > 0.5*phase.getBoundary().getDimensions().x(0)) {
            throw new RuntimeException("Truncation radius too large.  Max allowed is"+0.5*phase.getBoundary().getDimensions().x(0));
            }
        potentialMaster.setCellRange(2);
        potentialMaster.setRange(1.2*truncationRadiusCO);
        P2SoftSphericalTruncated interpotentialCO = new P2SoftSphericalTruncated (new P2ElectrostaticDreiding(space, 3022.850200*11604.45728, 0.264550, 17.160239*11604.45728), truncationRadiusCO); 
        potentialMaster.addPotential(interpotentialCO, new AtomType[]{((AtomFactoryParacetamol)species.getFactory()).cType, ((AtomFactoryParacetamol)species.getFactory()).oType} );
        
        if(truncationRadiusCN > 0.5*phase.getBoundary().getDimensions().x(0)) {
            throw new RuntimeException("Truncation radius too large.  Max allowed is"+0.5*phase.getBoundary().getDimensions().x(0));
            }
        potentialMaster.setCellRange(2);
        potentialMaster.setRange(1.2*truncationRadiusCN);
        P2SoftSphericalTruncated interpotentialCN = new P2SoftSphericalTruncated (new P2ElectrostaticDreiding(space, 3179.514600*11604.45728, 0.271003, 19.006710*11604.45728), truncationRadiusCN); 
        potentialMaster.addPotential(interpotentialCN, new AtomType[]{((AtomFactoryParacetamol)species.getFactory()).cType, ((AtomFactoryParacetamol)species.getFactory()).nType} );
        
        if(truncationRadiusON > 0.5*phase.getBoundary().getDimensions().x(0)) {
            throw new RuntimeException("Truncation radius too large.  Max allowed is"+0.5*phase.getBoundary().getDimensions().x(0));
            }
        potentialMaster.setCellRange(2);
        potentialMaster.setRange(1.2*truncationRadiusON);
        P2SoftSphericalTruncated interpotentialON = new P2SoftSphericalTruncated (new P2ElectrostaticDreiding(space, 2508.044800*11604.45728, 0.258398, 12.898341*11604.45728), truncationRadiusON); 
        potentialMaster.addPotential(interpotentialON, new AtomType[]{((AtomFactoryParacetamol)species.getFactory()).oType, ((AtomFactoryParacetamol)species.getFactory()).nType} );
        
        if(truncationRadiusOO > 0.5*phase.getBoundary().getDimensions().x(0)) {
            throw new RuntimeException("Truncation radius too large.  Max allowed is"+0.5*phase.getBoundary().getDimensions().x(0));
            }
        potentialMaster.setCellRange(2);
        potentialMaster.setRange(1.2*truncationRadiusOO);
        P2SoftSphericalTruncated interpotentialOO = new P2SoftSphericalTruncated (new P2ElectrostaticDreiding(space, 2384.465800*11604.45728, 0.252525, 11.645288*11604.45728), truncationRadiusOO); 
        potentialMaster.addPotential(interpotentialOO, new AtomType[]{((AtomFactoryParacetamol)species.getFactory()).oType, ((AtomFactoryParacetamol)species.getFactory()).oType} );
        
        if(truncationRadiusNN > 0.5*phase.getBoundary().getDimensions().x(0)) {
            throw new RuntimeException("Truncation radius too large.  Max allowed is"+0.5*phase.getBoundary().getDimensions().x(0));
            }
        potentialMaster.setCellRange(2);
        potentialMaster.setRange(1.2*truncationRadiusNN);
        P2SoftSphericalTruncated interpotentialNN = new P2SoftSphericalTruncated (new P2ElectrostaticDreiding(space, 2638.028500*11604.45728, 0.264550, 14.286224*11604.45728), truncationRadiusNN); 
        potentialMaster.addPotential(interpotentialNN, new AtomType[]{((AtomFactoryParacetamol)species.getFactory()).nType, ((AtomFactoryParacetamol)species.getFactory()).nType} );
        
        ((CriterionInterMolecular)potentialMaster.getCriterion(interpotentialCC)).setIntraMolecularCriterion(new CriterionNone());
        ((CriterionInterMolecular)potentialMaster.getCriterion(interpotentialCO)).setIntraMolecularCriterion(new CriterionNone());
        ((CriterionInterMolecular)potentialMaster.getCriterion(interpotentialCN)).setIntraMolecularCriterion(new CriterionNone());
        ((CriterionInterMolecular)potentialMaster.getCriterion(interpotentialON)).setIntraMolecularCriterion(new CriterionNone());
        ((CriterionInterMolecular)potentialMaster.getCriterion(interpotentialOO)).setIntraMolecularCriterion(new CriterionNone());
        ((CriterionInterMolecular)potentialMaster.getCriterion(interpotentialNN)).setIntraMolecularCriterion(new CriterionNone());
        
        
        bdry =  new BoundaryDeformableLattice(primitive, getRandom(), new int []{2, 3, 4});
        // bdry.setDimensions(Space.makeVector(new double []{3*12.119, 4*8.944, 4*7.278}));
        phase.setBoundary(bdry);
        configMonoLattice.initializeCoordinates(phase);
       	
        CoordinateDefinitionParacetamol coordDef = new CoordinateDefinitionParacetamol(phase, primitive, basis);

        integrator.setPhase(phase);
        
    } //end of constructor
    
   //////////////////////////////////////////////////////////////////////////////////////////////////////////

    /**
     * Demonstrates how this class is implemented.
     */
    public static void main(String[] args) {
        Default defaults = new Default();
        defaults.doSleep = false;
        defaults.ignoreOverlap = true;
        etomica.paracetamol.MDParacetamolMonoclinic sim = new etomica.paracetamol.MDParacetamolMonoclinic(defaults);
        SimulationGraphic simGraphic = new SimulationGraphic(sim, APP_NAME);
        
   /*****************************************************************************/    
        MeterKineticEnergy meterKE = new MeterKineticEnergy();
        meterKE.setPhase(sim.phase);
        DisplayBox KEbox = new DisplayBox();
        DataPump KEpump = new DataPump(meterKE, KEbox);
        new IntervalActionAdapter(KEpump, sim.integrator);
        
        MeterPotentialEnergy meterPE = new MeterPotentialEnergy(sim.potentialMaster);
        meterPE.setPhase(sim.phase);
        DisplayBox PEbox = new DisplayBox();
        DataPump PEpump = new DataPump(meterPE, PEbox);
        new IntervalActionAdapter(PEpump, sim.integrator);
        
        MeterEnergy meterTotal = new MeterEnergy(sim.potentialMaster);
        meterTotal.setPhase(sim.phase);   
        DisplayBox meterTotalbox = new DisplayBox();
        DataPump meterTotalpump = new DataPump(meterTotal, meterTotalbox);
        new IntervalActionAdapter(meterTotalpump, sim.integrator);
           
        MeterTemperature meterTemp = new MeterTemperature();
        meterTemp.setPhase(sim.phase);
        DisplayBox tempBox = new DisplayBox();
        tempBox.setUnit(Kelvin.UNIT);
        DataPump tempPump = new DataPump(meterTemp, tempBox);
        new IntervalActionAdapter(tempPump, sim.integrator);
  
        
 /**********************************************************************/   
        simGraphic.add(KEbox);
        simGraphic.add(PEbox);
        simGraphic.add(meterTotalbox);
        simGraphic.add(tempBox);
        
        
        simGraphic.getDisplayPhase(sim.phase).setPixelUnit(new Pixel(PIXEL_SIZE));
        ColorSchemeByType colorScheme = ((ColorSchemeByType)((DisplayPhase)simGraphic.displayList().getFirst()).getColorScheme());
        AtomTypeGroup atomType = (AtomTypeGroup)sim.species.getMoleculeType();
        colorScheme.setColor(atomType.getChildTypes()[0], java.awt.Color.red);
        colorScheme.setColor(atomType.getChildTypes()[1], java.awt.Color.gray);
        colorScheme.setColor(atomType.getChildTypes()[2], java.awt.Color.blue);
        colorScheme.setColor(atomType.getChildTypes()[3], java.awt.Color.white);
        colorScheme.setColor(atomType.getChildTypes()[4], java.awt.Color.white);

        simGraphic.makeAndDisplayFrame(APP_NAME);

        simGraphic.getDisplayPhase(sim.phase).repaint();

    }//end of main
    
    public BravaisLattice lattice;
    public BoundaryDeformableLattice bdry;
    public CoordinateDefinitionParacetamol coordinateDefinition;
    public ConfigurationMonoclinicLattice configMonoLattice;
}//end of class