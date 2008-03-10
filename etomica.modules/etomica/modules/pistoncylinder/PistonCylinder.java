package etomica.modules.pistoncylinder;

import etomica.action.activity.ActivityIntegrate;
import etomica.api.IBox;
import etomica.api.IPotentialMaster;
import etomica.api.ISpecies;
import etomica.api.IVector;
import etomica.atom.AtomTypeSphere;
import etomica.box.Box;
import etomica.chem.elements.ElementSimple;
import etomica.config.ConfigurationLattice;
import etomica.integrator.IntegratorMD.ThermostatType;
import etomica.lattice.LatticeCubicFcc;
import etomica.lattice.LatticeOrthorhombicHexagonal;
import etomica.potential.P1HardBoundary;
import etomica.potential.P1HardMoleculeMonatomic;
import etomica.potential.P1HardMovingBoundary;
import etomica.potential.P2HardMoleculeMonatomic;
import etomica.potential.P2SquareWell;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space2d.Vector2D;
import etomica.space3d.Vector3D;
import etomica.species.SpeciesSpheresMono;
import etomica.units.Bar;

/**
 * Simple hard-sphere MD in piston-cylinder apparatus
 */
public class PistonCylinder extends Simulation {
    
    private static final long serialVersionUID = 1L;
    private final int INIT_NUM_MOLECULES = 100;
    public IntegratorHardPiston integrator;
    public SpeciesSpheresMono species;
    public IBox box;
    public P2HardMoleculeMonatomic potentialWrapper;
    public P1HardBoundary wallPotential;
    public P1HardMoleculeMonatomic pistonPotentialWrapper;
    public ActivityIntegrate ai;
    public double lambda;
    public ConfigurationLattice config;

    public PistonCylinder(int D) {
        super(Space.getInstance(D), true);
        IPotentialMaster potentialMaster = new PotentialMaster(space);
        lambda = 2.0;
        double sigma = 4.0;
        species = new SpeciesSpheresMono(this);
        ((ElementSimple)species.getLeafType().getElement()).setMass(16);
        ((AtomTypeSphere)species.getLeafType()).setDiameter(sigma);
        getSpeciesManager().addSpecies(species);
        box = new Box(new BoundaryPistonCylinder(space, getRandom()), space);
        addBox(box);
        box.setNMolecules(species, INIT_NUM_MOLECULES);
        IVector newDim;
        if (space.D() == 2) {
            config = new ConfigurationLattice(new LatticeOrthorhombicHexagonal(), space);
            newDim = new Vector2D(80,150);
        }
        else {
            config = new ConfigurationLattice(new LatticeCubicFcc(), space);
            newDim = new Vector3D(30,80,30);
        }
        config.setBoundaryPadding(sigma);
        box.setDimensions(newDim);
        config.initializeCoordinates(box);
        
        P2SquareWell potentialSW = new P2SquareWell(space,sigma,lambda,31.875,true);
        potentialWrapper = new P2HardMoleculeMonatomic(space,potentialSW);
        potentialMaster.addPotential(potentialWrapper,new ISpecies[]{species,species});
        
        wallPotential = new P1HardBoundary(space, true);
        wallPotential.setCollisionRadius(sigma*0.5); //potential.getCoreDiameter()*0.5);
        potentialMaster.addPotential(new P1HardMoleculeMonatomic(space, wallPotential),new ISpecies[]{species});
        wallPotential.setActive(0,true,true);  // left wall
        wallPotential.setActive(0,false,true); // right wall
        if (D==3) {
            wallPotential.setActive(1,true,true); // top wall
            wallPotential.setActive(1,false,false); // bottom wall
            wallPotential.setActive(2,true,true);  // front wall
            wallPotential.setActive(2,false,true); // back wall
        }
        else {
            wallPotential.setActive(1,true,false); // top wall
            wallPotential.setActive(1,false,true); // bottom wall
        }

        P1HardMovingBoundary pistonPotential = new P1HardMovingBoundary(space,box.getBoundary(),1,400, true);
        pistonPotential.setCollisionRadius(sigma*0.5);
        if (D == 3) {
            pistonPotential.setWallPosition(box.getBoundary().getDimensions().x(1)*0.5);
            pistonPotential.setWallVelocity(-0.5);
            pistonPotential.setPressure(-Bar.UNIT.toSim(1.0));
        }
        else {
            pistonPotential.setWallPosition(-box.getBoundary().getDimensions().x(1)*0.5);
            pistonPotential.setWallVelocity(0.5);
            pistonPotential.setPressure(Bar.UNIT.toSim(100.0));
        }
        pistonPotential.setThickness(1.0);
        pistonPotentialWrapper = new P1HardMoleculeMonatomic(space, pistonPotential);
        potentialMaster.addPotential(pistonPotentialWrapper, new ISpecies[]{species});
        ((BoundaryPistonCylinder)box.getBoundary()).setPistonPotential(pistonPotential);
        
        integrator = new IntegratorHardPiston(this,potentialMaster,pistonPotentialWrapper, space);
        integrator.setBox(box);
        integrator.setIsothermal(true);
        integrator.setThermostatInterval(1);
        integrator.setThermostat(ThermostatType.ANDERSEN_SINGLE);
        integrator.setTimeStep(1.0);
        ai = new ActivityIntegrate(integrator,0,true);
        getController().addAction(ai);
        
    }
}
