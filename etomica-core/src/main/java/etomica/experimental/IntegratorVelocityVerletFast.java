package etomica.experimental;

import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.atom.IAtomKinetic;
import etomica.box.Box;
import etomica.data.DataSourceScalar;
import etomica.integrator.IntegratorMD;
import etomica.potential.*;
import etomica.space.BoundaryRectangularPeriodic;
import etomica.space3d.Space3D;
import etomica.space3d.Vector3D;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Energy;
import etomica.util.random.IRandom;
import jdk.incubator.vector.DoubleVector;
import jdk.incubator.vector.Shapes;

import java.util.stream.IntStream;

public class IntegratorVelocityVerletFast extends IntegratorMD {
    private static final DoubleVector.DoubleSpecies<Shapes.S256Bit> SPECIES = VectorSystem3D.SPECIES;
    private final VectorSystem3D positions;
    private final VectorSystem3D forces;
    private final VectorSystem3D velocities;
    private final int[] atomTypes;
    private final AtomType[] types;
    private final Potential2Soft[][] potentials;

    private final BoundaryRectangularPeriodic boundary;

    public IntegratorVelocityVerletFast(PotentialMaster potentialMaster, IRandom random, double timeStep, double temperature, Box box) {
        super(potentialMaster, random, timeStep, temperature, box);
        this.atomTypes = new int[box.getLeafList().size()];
        this.boundary = (BoundaryRectangularPeriodic) box.getBoundary();

        this.positions = new VectorSystem3D(this.box);
        this.forces = new VectorSystem3D(this.box.getLeafList().size());
        this.velocities = new VectorSystem3D(this.box.getLeafList().size());

        for (int i = 0; i < this.box.getLeafList().size(); i++) {
            this.velocities.setVector(i, ((IAtomKinetic) box.getLeafList().get(i)).getVelocity());
            this.atomTypes[i] = box.getLeafList().get(i).getType().getIndex();
        }

        types = new AtomType[]{this.box.getLeafList().get(0).getType()};
        potentials = new Potential2Soft[][]{{new P2SoftSphericalTruncated(new P2LennardJones(Space3D.getInstance()), 3)}};
        computeForces();
    }

    private void computeForces() {
        this.forces.setAll(0);
        for (int i = 0; i < this.positions.getRows(); i++) {
            for (int j = i + 1; j < this.positions.getRows(); j++) {
                Vector3D dr = this.positions.diff(i, j);
                this.box.getBoundary().nearestImage(dr);
                double r2 = dr.squared();
                Potential2Soft potential = potentials[atomTypes[i]][atomTypes[j]];
                double du = potential.du(r2);
                if (du == 0) {
                    continue;
                }

                dr.TE(du / r2);
                this.forces.sub(i, dr);
                this.forces.add(j, dr);
//                System.out.printf("%d, %f, %f, %s%n", j, du, r2, dr);
            }
//            System.out.println(forces.get(i));
//            System.exit(1);
        }
    }

    private static DoubleVector<Shapes.S256Bit> nearestImage(DoubleVector<Shapes.S256Bit> dr, DoubleVector<Shapes.S256Bit> dim, DoubleVector<Shapes.S256Bit> dimHalf) {
        // which elements of dr are greater than dimHalf? Makes a mask like [false, true, true, false] (last lane is always false)
        var mask = dr.greaterThan(dimHalf);
        var image = dr; // does this work?

        while (mask.anyTrue()) {
            // subtract dimHalf for only the elements that need it
            image = image.sub(dimHalf, mask);
            mask = image.greaterThan(dimHalf);
        }

        var negDimHalf = dimHalf.neg();
        mask = image.lessThan(negDimHalf);
        while (mask.anyTrue()) {
            image = image.add(dimHalf, mask);
            mask = image.lessThan(negDimHalf);
        }

        return image;
    }

    private void computeForcesSimd() {
        var arr = new double[4];
        boundary.getBoxSize().assignTo(arr);
        var dim = SPECIES.fromArray(arr, 0);
        var dimHalf = dim.div(2);

        this.forces.setAll(0);
        for (int i = 0; i < this.positions.getRows(); i++) {
            for (int j = i + 1; j < this.positions.getRows(); j++) {
                var dr = this.positions.diffV(i, j);

//                this.box.getBoundary().nearestImage(dr);
//                dr = nearestImage(dr, dim, dimHalf);

                // *** nearestImage (it doesn't inline for some reason?)
                var mask = dr.greaterThan(dimHalf);

                while (mask.anyTrue()) {
                    // subtract dimHalf for only the elements that need it
                    dr = dr.sub(dimHalf, mask);
                    mask = dr.greaterThan(dimHalf);
                }

                var negDimHalf = dimHalf.neg();
                mask = dr.lessThan(negDimHalf);
                while (mask.anyTrue()) {
                    dr = dr.add(dimHalf, mask);
                    mask = dr.lessThan(negDimHalf);
                }
                // ***

//                double r2 = dr.squared();
                double r2 = dr.mul(dr).addAll();

                Potential2Soft potential = potentials[atomTypes[i]][atomTypes[j]];
                double du = potential.du(r2);
                if (du == 0) {
                    continue;
                }

//                dr.TE(du / r2);
                dr = dr.mul(du / r2);

                this.forces.subV(i, dr);
                this.forces.addV(j, dr);
            }
        }
    }

    private void computeForcesParallel() {
        IntStream.range(0, this.positions.getRows()).parallel().forEach(i -> {
            for (int j = i + 1; j < this.positions.getRows(); j++) {
                Vector3D dr = this.positions.diff(i, j);
                this.box.getBoundary().nearestImage(dr);
                double r2 = dr.squared();
                Potential2Soft potential = potentials[atomTypes[i]][atomTypes[j]];
                double du = potential.du(r2);
                if (du == 0) {
                    continue;
                }

                dr.TE(du / r2);
                this.forces.subAtomic(i, dr);
                this.forces.addAtomic(j, dr);
            }
        });
    }

    private double getEnergy() {
        double energy = 0;
        for (int i = 0; i < this.positions.getRows(); i++) {
            for (int j = i + 1; j < this.positions.getRows(); j++) {
                Vector3D dr = this.positions.diff(i, j);
                this.box.getBoundary().nearestImage(dr);
                double r2 = dr.squared();
                Potential2Soft potential = potentials[atomTypes[i]][atomTypes[j]];
                energy += potential.u(r2);
            }
        }
        return energy;

//        return IntStream.range(0, this.positions.getRows()).parallel().mapToDouble(i -> {
//            double energy = 0;
//            for (int j = i + 1; j < this.positions.getRows(); j++) {
//                Vector3D dr = this.positions.diff(i, j);
//                this.box.getBoundary().nearestImage(dr);
//                double r2 = dr.squared();
//                Potential2Soft potential = potentials[atomTypes[i]][atomTypes[j]];
//                energy += potential.u(r2);
//            }
//            return energy;
//        }).sum();
    }

    protected void doStepInternal() {
        super.doStepInternal();

        for (int i = 0; i < this.positions.getRows(); i++) {
            this.velocities.addScaled(i, i, 0.5 * timeStep * this.types[this.atomTypes[i]].rm(), forces);
            this.positions.addScaled(i, i, timeStep, velocities);
        }

        this.computeForcesSimd();

        for (int i = 0; i < this.positions.getRows(); i++) {
            this.velocities.addScaled(i, i, 0.5 * timeStep * this.types[this.atomTypes[i]].rm(), forces);
        }

    }

    public MeterPotentialEnergy getMeter() {
        return new MeterPotentialEnergy();
    }

    public class MeterPotentialEnergy extends DataSourceScalar {

        public MeterPotentialEnergy() {
            super("Potential Energy", Energy.DIMENSION);
        }

        @Override
        public double getDataAsScalar() {
            return getEnergy();
        }
    }
}
