package etomica.modules.catalysis;
import etomica.action.activity.ActivityIntegrate;
import etomica.api.IAtomType;
import etomica.api.IBox;
import etomica.api.IVectorMutable;
import etomica.box.Box;
import etomica.chem.elements.ElementSimple;
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

/**
 * Simulation for Catalysis module.
 * Design by Ken Benjamin
 *
 * @author Andrew Schultz
 */
public class Catalysis extends Simulation {
    
    private static final long serialVersionUID = 1L;
    public final SpeciesSpheresMono speciesO2, speciesCO, speciesO, speciesCO2, speciesSurface;
    public final IBox box;
    public final IntegratorHard integrator;
    public final ActivityIntegrate activityIntegrate;
    public final P2SquareWell potentialO2O2, potentialO2CO, potentialCOCO;
    public final P2SquareWell potentialSurfaceO2, potentialSurfaceCO;
    public final ConfigurationCatalysis config;
    
    public Catalysis(Space _space) {
        super(_space);
        PotentialMasterList potentialMaster = new PotentialMasterList(this, 2.5, space); //List(this, 2.0);
        
        int N = 256;  //number of atoms
        
        //controller and integrator
	    integrator = new IntegratorHard(this, potentialMaster, space);
//	    integrator.setTimeStep(0.01);
	    integrator.setTemperature(1);
	    integrator.setIsothermal(false);
        integrator.setThermostat(ThermostatType.ANDERSEN_SINGLE);
        integrator.setThermostatInterval(1);
        activityIntegrate = new ActivityIntegrate(integrator);
        getController().addAction(activityIntegrate);

	    //species and potentials
	    speciesO2 = new SpeciesSpheresMono(this, space, new ElementSimple("O2", 15.9994*2));
        getSpeciesManager().addSpecies(speciesO2);
        speciesCO = new SpeciesSpheresMono(this, space, new ElementSimple("CO", 12.011+15.9994));
        getSpeciesManager().addSpecies(speciesCO);
        speciesO = new SpeciesSpheresMono(this, space, new ElementSimple("O", 15.9994));
        getSpeciesManager().addSpecies(speciesO);
        speciesCO2 = new SpeciesSpheresMono(this, space, new ElementSimple("CO2", 15.9994*2+12.011));
        getSpeciesManager().addSpecies(speciesCO2);
        speciesSurface = new SpeciesSpheresMono(this, space, new ElementSimple("Surface", Double.POSITIVE_INFINITY));
        getSpeciesManager().addSpecies(speciesSurface);
        
        //instantiate several potentials for selection in combo-box
	    potentialO2O2 = new P2SquareWell(space, 1, 1.5, 1, true);
        potentialMaster.addPotential(potentialO2O2,new IAtomType[]{speciesO2.getLeafType(), speciesO2.getLeafType()});

        potentialO2CO = new P2SquareWell(space, 1, 1.5, 1, true);
        potentialMaster.addPotential(potentialO2CO,new IAtomType[]{speciesO2.getLeafType(), speciesCO.getLeafType()});

        potentialCOCO = new P2SquareWell(space, 1, 1.5, 1, true);
        potentialMaster.addPotential(potentialCOCO,new IAtomType[]{speciesCO.getLeafType(), speciesCO.getLeafType()});

        potentialSurfaceO2 = new P2SquareWell(space, 1, 1.3, 9, true);
        potentialMaster.addPotential(potentialSurfaceO2,new IAtomType[]{speciesO2.getLeafType(), speciesSurface.getLeafType()});
        
        potentialSurfaceCO = new P2SquareWell(space, 1, 1.3, 6, true);
        potentialMaster.addPotential(potentialSurfaceCO,new IAtomType[]{speciesCO.getLeafType(), speciesSurface.getLeafType()});
        
        P1HardBoundary p1HardWallO2 = new P1HardBoundary(space, true);
        p1HardWallO2.setActive(0, true, false);
        p1HardWallO2.setActive(0, false, false);
        p1HardWallO2.setActive(1, true, false);
        p1HardWallO2.setActive(1, false, true);
        p1HardWallO2.setActive(2, true, false);
        p1HardWallO2.setActive(2, false, false);
        potentialMaster.addPotential(p1HardWallO2, new IAtomType[]{speciesO2.getLeafType()});
        P1HardBoundary p1HardWallCO = new P1HardBoundary(space, true);
        p1HardWallO2.setActive(0, true, false);
        p1HardWallO2.setActive(0, false, false);
        p1HardWallO2.setActive(1, true, false);
        p1HardWallO2.setActive(1, false, true);
        p1HardWallO2.setActive(2, true, false);
        p1HardWallO2.setActive(2, false, false);
        potentialMaster.addPotential(p1HardWallCO, new IAtomType[]{speciesCO.getLeafType()});
	    
        //construct box
	    box = new Box(new BoundaryRectangularSlit(1, 20.0, space), space);
        addBox(box);
        integrator.getEventManager().addListener(potentialMaster.getNeighborManager(box));
        IVectorMutable dim = space.makeVector();
        dim.E(30);
        box.getBoundary().setBoxSize(dim);
        box.setNMolecules(speciesO2, N/2);
        box.setNMolecules(speciesCO, N/2);
        
        config = new ConfigurationCatalysis(this, space, speciesSurface);
        config.setNCellsX(30);
        config.setNCellsZ(20);
        config.setCellSizeX(1);
        config.setCellSizeZ(Math.sqrt(3));
        config.initializeCoordinates(box);
        
        integrator.setBox(box);
    }
    
    public static void main(String[] args) {
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
