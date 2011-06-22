package etomica.modules.adsorption;
import etomica.action.activity.ActivityIntegrate;
import etomica.api.IAtomType;
import etomica.api.IBox;
import etomica.api.IVectorMutable;
import etomica.box.Box;
import etomica.chem.elements.Carbon;
import etomica.config.ConfigurationLattice;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorHard;
import etomica.integrator.IntegratorMC;
import etomica.lattice.LatticeCubicFcc;
import etomica.nbr.PotentialMasterHybrid;
import etomica.potential.P2SquareWell;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularSlit;
import etomica.space.ISpace;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;
import etomica.util.RandomNumberGenerator;

/**
 * Simulation for Adsorption module.
 * Design by Lev Gelb
 *
 * @author Andrew Schultz
 */
public class Adsorption extends Simulation {
    
    private static final long serialVersionUID = 1L;
    public final SpeciesSpheresMono species;
    public final IBox box;
    public final IntegratorHard integratorMD;
    public final IntegratorMC integratorMC;
    public final IntegratorHybrid integratorHybrid;
    public final ActivityIntegrate activityIntegrate;
    public final P2SquareWell p2;
    public final P1Wall p1Wall;
    public final MyMCMove mcMoveID;
    
    public Adsorption(ISpace _space) {
        super(_space);
        setRandom(new RandomNumberGenerator(1));
        PotentialMasterHybrid potentialMaster = new PotentialMasterHybrid(this, 2, space); //List(this, 2.0);
        
        //controller and integrator
	    integratorMD = new IntegratorHard(this, potentialMaster, space);
	    integratorMD.setTimeStep(0.002);
	    integratorMD.setIsothermal(false);
	    
	    integratorMC = new IntegratorMC(potentialMaster, random, 2);
	    integratorMC.setTemperature(1);
	    
	    integratorHybrid = new IntegratorHybrid(potentialMaster, integratorMD, integratorMC, 2);
	    
        activityIntegrate = new ActivityIntegrate(integratorHybrid);
        getController().addAction(activityIntegrate);

        double sigma = 1;
        double lambda = 1.5;
        double epsilon = 1.0;
        double epsilonWF = 5.0;
        
	    //species and potentials
        species = new SpeciesSpheresMono(space, Carbon.INSTANCE);
        species.setIsDynamic(true);
        addSpecies(species);

        //construct box
        box = new Box(new BoundaryRectangularSlit(1, 20.0, space), space);
        addBox(box);

        mcMoveID = new MyMCMove(integratorMC, random, space, 0.1, sigma, 1);
        mcMoveID.setMu(-12);
        integratorMC.getMoveManager().addMCMove(mcMoveID);
        integratorHybrid.setMCMoveInsertDelete(mcMoveID);
        mcMoveID.setSpecies(species);
        mcMoveID.setBox(box);

        p2 = new P2SquareWell(space, sigma, lambda, epsilon, false);
        potentialMaster.addPotential(p2,new IAtomType[]{species.getLeafType(), species.getLeafType()});


        p1Wall = new P1Wall(space);
        p1Wall.setSigma(sigma);
        p1Wall.setRange(sigma/2);
        p1Wall.setEpsilon(epsilonWF);
        p1Wall.setThermalize(integratorMC, 0.1, random);

        potentialMaster.addPotential(p1Wall, new IAtomType[]{species.getLeafType()});
        
        integratorMD.getEventManager().addListener(potentialMaster.getNeighborManager(box));

        IVectorMutable dim = space.makeVector();
        dim.E(8*sigma);
        dim.setX(1, 12*sigma);
        box.getBoundary().setBoxSize(dim);
        
        box.setNMolecules(species, 40);

        ConfigurationLattice config = new ConfigurationLattice(new LatticeCubicFcc(space), space);
        config.initializeCoordinates(box);
        
        integratorHybrid.setBox(box);
    }
    
    public static void main(String[] args) {
        ISpace space = Space3D.getInstance();
        
        Adsorption sim = new Adsorption(space);
        SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.GRAPHIC_ONLY, "Catalysis", 1, sim.space, sim.getController());
        simGraphic.makeAndDisplayFrame();
    }
}
