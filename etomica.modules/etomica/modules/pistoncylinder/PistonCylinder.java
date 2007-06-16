package etomica.modules.pistoncylinder;

import etomica.action.activity.ActivityIntegrate;
import etomica.action.activity.Controller;
import etomica.atom.AtomTypeLeaf;
import etomica.atom.AtomTypeSphere;
import etomica.chem.elements.ElementSimple;
import etomica.config.Configuration;
import etomica.config.ConfigurationLattice;
import etomica.integrator.IntegratorMD.ThermostatType;
import etomica.lattice.LatticeCubicFcc;
import etomica.lattice.LatticeOrthorhombicHexagonal;
import etomica.phase.Phase;
import etomica.potential.P1HardBoundary;
import etomica.potential.P1HardMovingBoundary;
import etomica.potential.P2SquareWell;
import etomica.potential.Potential2HardSphericalWrapper;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.IVector;
import etomica.space.Space;
import etomica.space2d.Vector2D;
import etomica.space3d.Vector3D;
import etomica.species.Species;
import etomica.species.SpeciesSpheresMono;
import etomica.units.Bar;

/**
 * Simple hard-sphere MD in piston-cylinder apparatus
 */
public class PistonCylinder extends Simulation {
    
    private static final long serialVersionUID = 1L;
    public IntegratorHardPiston integrator;
    public SpeciesSpheresMono species;
    public Phase phase;
    public Controller controller;
    public Potential2HardSphericalWrapper potentialWrapper;
    public P1HardBoundary wallPotential;
    public P1HardMovingBoundary pistonPotential;
    public ActivityIntegrate ai;
    public double lambda;

    public PistonCylinder(int D) {
        super(Space.getInstance(D), true);
        PotentialMaster potentialMaster = new PotentialMaster(space);
        lambda = 1.5;
        double sigma = 3.0;
        controller = getController();
        species = new SpeciesSpheresMono(this);
        ((ElementSimple)((AtomTypeLeaf)species.getMoleculeType()).getElement()).setMass(16);
        ((AtomTypeSphere)species.getMoleculeType()).setDiameter(sigma);
        getSpeciesManager().addSpecies(species);
        phase = new Phase(new BoundaryPistonCylinder(space, getRandom()));
        addPhase(phase);
        phase.getAgent(species).setNMolecules(112);
        IVector newDim;
        Configuration config;
        if (space.D() == 2) {
            config = new ConfigurationLattice(new LatticeOrthorhombicHexagonal());
            newDim = new Vector2D(80,150);
        }
        else {
            config = new ConfigurationLattice(new LatticeCubicFcc());
            newDim = new Vector3D(80,80,80);
        }
        phase.setDimensions(newDim);
        config.initializeCoordinates(phase);
        
        P2SquareWell potentialSW = new P2SquareWell(space,sigma,lambda,31.875,true);
        potentialWrapper = new Potential2HardSphericalWrapper(space,potentialSW);
//        potential = new P2HardSphere(space,Default.atomSize);
        potentialMaster.addPotential(potentialWrapper,new Species[]{species,species});
        
        wallPotential = new P1HardBoundary(space, true);
        wallPotential.setCollisionRadius(sigma*0.5); //potential.getCoreDiameter()*0.5);
        potentialMaster.addPotential(wallPotential,new Species[]{species});
        wallPotential.setActive(0,true,true);  // left wall
        wallPotential.setActive(0,false,true); // right wall
        wallPotential.setActive(1,true,false); // top wall
        wallPotential.setActive(1,false,true); // bottom wall
        if (D==3) {
            wallPotential.setActive(2,true,true);  // front wall
            wallPotential.setActive(2,false,true); // back wall
        }

        pistonPotential = new P1HardMovingBoundary(space,phase.getBoundary(),1,400, true);
        pistonPotential.setCollisionRadius(sigma*0.5);
        pistonPotential.setWallPosition(-phase.getBoundary().getDimensions().x(1)*0.5);
        pistonPotential.setWallVelocity(0.5);
        if (D == 3) {
            pistonPotential.setPressure(Bar.UNIT.toSim(1.0));
        }
        else {
            pistonPotential.setPressure(Bar.UNIT.toSim(100.0));
        }
        pistonPotential.setThickness(1.0);
        potentialMaster.addPotential(pistonPotential,new Species[]{species});
        ((BoundaryPistonCylinder)phase.getBoundary()).setPistonPotential(pistonPotential);
        
        integrator = new IntegratorHardPiston(this,potentialMaster,pistonPotential);
        integrator.setPhase(phase);
        integrator.setIsothermal(true);
        integrator.setThermostatInterval(1);
        integrator.setThermostat(ThermostatType.ANDERSEN_SINGLE);
        integrator.setTimeStep(1.0);
        ai = new ActivityIntegrate(integrator);
        getController().addAction(ai);
        
    }
}
