package etomica.modules.catalysis;
import etomica.action.activity.ActivityIntegrate;
import etomica.api.IAtomType;
import etomica.api.IBox;
import etomica.api.IVectorMutable;
import etomica.atom.AtomTypeSphere;
import etomica.box.Box;
import etomica.chem.elements.Carbon;
import etomica.chem.elements.ElementSimple;
import etomica.chem.elements.Oxygen;
import etomica.graphics.SimulationGraphic;
import etomica.integrator.IntegratorHard;
import etomica.integrator.IntegratorMD.ThermostatType;
import etomica.nbr.list.PotentialMasterList;
import etomica.potential.P1HardBoundary;
import etomica.potential.P2SquareWell;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularSlit;
import etomica.space.Space;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresMono;
import etomica.units.Calorie;
import etomica.units.ElectronVolt;
import etomica.units.Kelvin;
import etomica.units.Mole;

/**
 * Simulation for Catalysis module.
 * Design by Ken Benjamin
 *
 * @author Andrew Schultz
 */
public class Catalysis extends Simulation {
    
    private static final long serialVersionUID = 1L;
    public final SpeciesSpheresMono speciesO, speciesC, speciesSurface;
    public final IBox box;
    public final IntegratorHard integrator;
    public final ActivityIntegrate activityIntegrate;
    public final P2SquareWellBonding potentialOO;
    public final P2SquareWellBondingCO potentialCO;
    public final P2SquareWell potentialCC;
    public final P2SquareWellSurface potentialCS, potentialOS;
    public final ConfigurationCatalysis config;
    public final InteractionTracker interactionTracker;
    
    public Catalysis(Space _space) {
        super(_space);
        PotentialMasterList potentialMaster = new PotentialMasterList(this, 9, space); //List(this, 2.0);
        
        //controller and integrator
	    integrator = new IntegratorHard(this, potentialMaster, space);
	    integrator.setTimeStep(0.001);
	    integrator.setTemperature(Kelvin.UNIT.toSim(600));
	    integrator.setIsothermal(false);
        integrator.setThermostat(ThermostatType.ANDERSEN_SINGLE);
        integrator.setThermostatInterval(1);
        activityIntegrate = new ActivityIntegrate(integrator);
        getController().addAction(activityIntegrate);

        double sigmaO = 3.6;
        double sigmaC = 3.8;
        double sigmaS = 3.7;
        double epsilonO = 1000*Mole.UNIT.fromSim(Calorie.UNIT.toSim(0.15));
        double epsilonC = 1000*Mole.UNIT.fromSim(Calorie.UNIT.toSim(0.08));
        double epsilonOS = 10*epsilonO;
        double epsilonCS = 10*epsilonC;
        
	    //species and potentials
	    speciesO = new SpeciesSpheresMono(this, space, Oxygen.INSTANCE);
	    ((AtomTypeSphere)speciesO.getLeafType()).setDiameter(sigmaO);
        getSpeciesManager().addSpecies(speciesO);
        speciesC = new SpeciesSpheresMono(this, space, Carbon.INSTANCE);
        ((AtomTypeSphere)speciesC.getLeafType()).setDiameter(sigmaC);
        getSpeciesManager().addSpecies(speciesC);
        speciesSurface = new SpeciesSpheresMono(this, space, new ElementSimple("Surface", Double.POSITIVE_INFINITY));
        ((AtomTypeSphere)speciesSurface.getLeafType()).setDiameter(sigmaS);
        getSpeciesManager().addSpecies(speciesSurface);

        //construct box
        box = new Box(new BoundaryRectangularSlit(1, 20.0, space), space);
        addBox(box);
        interactionTracker = new InteractionTracker(box, speciesSurface);

	    potentialOO = new P2SquareWellBonding(space, interactionTracker.getAgentManager(), sigmaO, 1.3, epsilonO, 3, 1, 1, 7.4);
        potentialMaster.addPotential(potentialOO,new IAtomType[]{speciesO.getLeafType(), speciesO.getLeafType()});

        potentialCO = new P2SquareWellBondingCO(space, interactionTracker.getAgentManager(), 0.5*(sigmaO+sigmaC), 1.1, Math.sqrt(epsilonC*epsilonO), 20, 1, 3, 7.4);
        potentialMaster.addPotential(potentialCO,new IAtomType[]{speciesO.getLeafType(), speciesC.getLeafType()});

        potentialCC = new P2SquareWell(space, sigmaC, 1.3, epsilonC, false);
        potentialMaster.addPotential(potentialCC,new IAtomType[]{speciesC.getLeafType(), speciesC.getLeafType()});

        potentialOS = new P2SquareWellSurface(space, interactionTracker.getAgentManager(), 0.5*(sigmaO+sigmaS), 1.3, epsilonOS, 3);
        potentialMaster.addPotential(potentialOS,new IAtomType[]{speciesO.getLeafType(), speciesSurface.getLeafType()});
        
        potentialCS = new P2SquareWellSurface(space, interactionTracker.getAgentManager(), 0.5*(sigmaC+sigmaS), 1.3, epsilonCS, 3);
        potentialMaster.addPotential(potentialCS,new IAtomType[]{speciesC.getLeafType(), speciesSurface.getLeafType()});
        
        P1HardBoundary p1HardWallO = new P1HardBoundary(space, true);
        p1HardWallO.setActive(0, true, false);
        p1HardWallO.setActive(0, false, false);
        p1HardWallO.setActive(1, true, false);
        p1HardWallO.setActive(1, false, true);
        p1HardWallO.setActive(2, true, false);
        p1HardWallO.setActive(2, false, false);
        potentialMaster.addPotential(p1HardWallO, new IAtomType[]{speciesO.getLeafType()});
        P1HardBoundary p1HardWallC = new P1HardBoundary(space, true);
        p1HardWallC.setActive(0, true, false);
        p1HardWallC.setActive(0, false, false);
        p1HardWallC.setActive(1, true, false);
        p1HardWallC.setActive(1, false, true);
        p1HardWallC.setActive(2, true, false);
        p1HardWallC.setActive(2, false, false);
        potentialMaster.addPotential(p1HardWallC, new IAtomType[]{speciesC.getLeafType()});
        
        integrator.addCollisionListener(interactionTracker);
        ReactionManagerCO reactionManager = new ReactionManagerCO(this);
        integrator.getEventManager().addListener(reactionManager);
        
        integrator.getEventManager().addListener(potentialMaster.getNeighborManager(box));

        IVectorMutable dim = space.makeVector();
        dim.E(150);
        box.getBoundary().setBoxSize(dim);
        int nCO = 40, nO2 = 40;
        
        config = new ConfigurationCatalysis(this, space, speciesSurface, speciesC, speciesO, interactionTracker.getAgentManager());
        config.setNCellsX(40);
        config.setNCellsZ(25);
        config.setCellSizeX(sigmaS);
        config.setCellSizeZ(sigmaS*Math.sqrt(3));
        config.setNumCO(nCO);
        config.setNumO2(nO2);
        config.initializeCoordinates(box);
        
        integrator.setBox(box);
    }
    
    public static void main(String[] args) {
        System.out.println(1000*Mole.UNIT.fromSim(Calorie.UNIT.toSim(0.15)));
//        System.out.println(Kelvin.UNIT.fromSim(ElectronVolt.UNIT.toSim(1)));
//        System.out.println();
//        System.exit(1);
        Space space = Space3D.getInstance();
        if(args.length != 0) {
            try {
                int D = Integer.parseInt(args[0]);
                if (D == 3) {
                    space = Space3D.getInstance();
                }
            } catch(NumberFormatException e) {}
        }
            
        Catalysis sim = new Catalysis(space);
        SimulationGraphic simGraphic = new SimulationGraphic(sim, SimulationGraphic.GRAPHIC_ONLY, "Catalysis", 1, sim.space, sim.getController());
        simGraphic.makeAndDisplayFrame();
    }
}
