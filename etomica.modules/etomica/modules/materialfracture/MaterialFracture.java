package etomica.modules.materialfracture;

import etomica.action.BoxImposePbc;
import etomica.action.activity.ActivityIntegrate;
import etomica.api.IAtomTypeLeaf;
import etomica.api.IAtomTypeSphere;
import etomica.api.IBox;
import etomica.api.IVectorMutable;
import etomica.box.Box;
import etomica.chem.elements.ElementSimple;
import etomica.config.ConfigurationLattice;
import etomica.integrator.IntegratorMD;
import etomica.integrator.IntegratorVelocityVerlet;
import etomica.lattice.BravaisLatticeCrystal;
import etomica.lattice.crystal.BasisOrthorhombicHexagonal;
import etomica.lattice.crystal.PrimitiveGeneral;
import etomica.potential.P2LennardJones;
import etomica.potential.P2SoftSphericalTruncatedForceShifted;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.BoundaryRectangularSlit;
import etomica.space2d.Space2D;
import etomica.species.SpeciesSpheresMono;

public class MaterialFracture extends Simulation {

    public final Box box;
    public final ConfigurationLattice config;
    public final SpeciesSpheresMono species;
    public final IntegratorVelocityVerlet integrator;
    public final PotentialCalculationForceStress pc;
    public final P1Tension p1Tension;

    public MaterialFracture() {
        super(Space2D.getInstance());
        PotentialMaster potentialMaster = new PotentialMaster();
        box = new Box(space);
        box.setBoundary(new BoundaryRectangularSlit(0, space));
        box.getBoundary().setDimensions(space.makeVector(new double[]{90,30}));
        addBox(box);
        integrator = new IntegratorVelocityVerlet(this, potentialMaster, space);
        integrator.setIsothermal(true);
        integrator.setTemperature(300.0);
        integrator.setTimeStep(0.007);
        integrator.setThermostatInterval(400);
        integrator.setThermostat(IntegratorMD.ThermostatType.ANDERSEN_NODRIFT);
        integrator.setBox(box);
        pc = new PotentialCalculationForceStress(space);
        integrator.setForceSum(pc);
        getController().addAction(new ActivityIntegrate(integrator));
        P2LennardJones p2LJ = new P2LennardJones(space, 3, 2000);
        P2SoftSphericalTruncatedForceShifted pt = new P2SoftSphericalTruncatedForceShifted(space, p2LJ, 7.5);

        p1Tension = new P1Tension(space); 
        species = new SpeciesSpheresMono(this, space);
        ((ElementSimple)species.getLeafType().getElement()).setMass(40);
        ((IAtomTypeSphere)species.getLeafType()).setDiameter(3.0);
        getSpeciesManager().addSpecies(species);
        box.setNMolecules(species, 198);

        potentialMaster.addPotential(pt, new IAtomTypeLeaf[]{species.getLeafType(), species.getLeafType()});
        potentialMaster.addPotential(p1Tension, new IAtomTypeLeaf[]{species.getLeafType()});

        PrimitiveGeneral primitive = new PrimitiveGeneral(space, new IVectorMutable[]{space.makeVector(new double[]{Math.sqrt(3),0}), space.makeVector(new double[]{0,1})});
        config = new ConfigurationLattice(new BravaisLatticeCrystal(primitive, new BasisOrthorhombicHexagonal()), space) {
            public void initializeCoordinates(IBox aBox) {
                IVectorMutable d = space.makeVector();
                d.E(aBox.getBoundary().getDimensions());
                d.setX(0, 64.7);
                aBox.getBoundary().setDimensions(d);
                super.initializeCoordinates(aBox);
                d.setX(0, 90);
                aBox.getBoundary().setDimensions(d);
            }
        };
        config.initializeCoordinates(box);

        integrator.addIntervalAction(new BoxImposePbc(box, space));
    }
}
