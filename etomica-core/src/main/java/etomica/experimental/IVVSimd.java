package etomica.experimental;

import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.atom.IAtomKinetic;
import etomica.box.Box;
import etomica.data.DataSourceScalar;
import etomica.integrator.IntegratorMD;
import etomica.potential.*;
import etomica.space3d.Space3D;
import etomica.space3d.Vector3D;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Energy;
import etomica.util.random.IRandom;
import jdk.incubator.vector.DoubleVector;
import jdk.incubator.vector.Vector;

import java.util.Arrays;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class IVVSimd extends IntegratorMD implements EnergyMeter {
    private final VecSys3DAlt positions;
    private final VecSys3DAlt forces;
    private final VecSys3DAlt velocities;
    private final int[] atomTypes;
    private final AtomType[] types;
    private final Potential2Soft[][] potentials;

    private static final DoubleVector.DoubleSpecies S = DoubleVector.preferredSpecies();

    public IVVSimd(PotentialMaster potentialMaster, IRandom random, double timeStep, double temperature, Box box) {
        super(potentialMaster, random, timeStep, temperature, box);
        this.atomTypes = new int[box.getLeafList().size()];


        this.randomizeMomenta();
        shiftMomenta();
        scaleMomenta();
        this.positions = new VecSys3DAlt(this.box);
        this.forces = new VecSys3DAlt(this.box.getLeafList().size());
        this.velocities = new VecSys3DAlt(this.box.getLeafList().size());

        for (int i = 0; i < this.box.getLeafList().size(); i++) {
            this.velocities.setVector(i, ((IAtomKinetic) box.getLeafList().get(i)).getVelocity());
            this.atomTypes[i] = box.getLeafList().get(i).getType().getIndex();
        }

        types = new AtomType[]{this.box.getLeafList().get(0).getType()};
        potentials = new Potential2Soft[][]{{new P2SoftSphericalTruncated(new P2LennardJones(Space3D.getInstance()), 3)}};
        computeForces();
    }

    private void computeForces() {
        var dim = this.box.getBoundary().getBoxSize();
        this.forces.setAll(0);
        for (int i = 0; i < this.positions.size(); i++) {
            int j = i + 1;
            int nbrCount = 0;
            int nbrSize = this.positions.size() - j;
            int vectorizedSize = nbrSize & (-S.length()); // rounds down to nearest multiple of vec len
//            vectorizedSize = 0;

            var ix = S.broadcast(this.positions.xs[i]);
            var iy = S.broadcast(this.positions.ys[i]);
            var iz = S.broadcast(this.positions.zs[i]);

            for (; nbrCount < vectorizedSize; j+=S.length(), nbrCount += S.length()) {
                var dx = S.fromArray(this.positions.xs, j).sub(ix);
                var dy = S.fromArray(this.positions.ys, j).sub(iy);
                var dz = S.fromArray(this.positions.zs, j).sub(iz);


                dx = nearestImage(dx, dim.getX(0));
                dy = nearestImage(dy, dim.getX(1));
                dz = nearestImage(dz, dim.getX(2));

//                var r2s = dx.mul(dx).add(dy.mul(dy)).add(dz.mul(dz));
                var r2s = dx.fma(dx, dy.fma(dy, dz.mul(dz)));
                var cutoffMask = r2s.greaterThanEq(3*3);
                if (cutoffMask.allTrue()) {
                    continue;
                }
                var dus = duVec(r2s);
//                if (dus.equal(0).allTrue()) {
//                    continue;
//                }

                var div = dus.div(r2s);
                dx = dx.mul(div);
                dy = dy.mul(div);
                dz = dz.mul(div);

//                if (j >= 26 && j <= 29) {
//                    System.out.println(i + "," + "29" + " "  + dz.toArray()[29 -  j] + " sum: " + this.forces.zs[29]);
//                }
//
//                if (i == 29) {
//                    System.out.println(i + "," + j + " " + Arrays.toString(dz.toArray()) + " sum: " + this.forces.zs[29]);
//                }

                this.forces.xs[i] += dx.addAll();
                this.forces.ys[i] += dy.addAll();
                this.forces.zs[i] += dz.addAll();

                S.fromArray(this.forces.xs, j).sub(dx).intoArray(this.forces.xs, j);
                S.fromArray(this.forces.ys, j).sub(dy).intoArray(this.forces.ys, j);
                S.fromArray(this.forces.zs, j).sub(dz).intoArray(this.forces.zs, j);

            }

            for (; nbrCount < nbrSize; j++, nbrCount++) {
                Vector3D dr = this.positions.diffSingle(j, i);
                this.box.getBoundary().nearestImage(dr);
                double r2 = dr.squared();
                var potential = potentials[0][0];
                double du = potential.du(r2);
                if (du == 0) {
                    continue;
                }

                dr.TE(du / r2);
//                if (j == 29 || i == 29) {
//                    System.out.println(i + "," + j + " "  + dr);
//                }
//                System.out.println(j + " " + dr.getX(0));
                this.forces.addSingle(i, dr);
                this.forces.subSingle(j, dr);
            }

//            if (i == 29) {
//                System.out.println(Arrays.toString(new double[] {this.forces.xs[i], this.forces.ys[i], this.forces.zs[i]}));
//            }
//            System.exit(1);
        }
    }

    private static DoubleVector duVec(DoubleVector r2s) {
        var s2 = S.broadcast(1.0).div(r2s);
        var s6 = s2.mul(s2).mul(s2);
        var du = s6.mul(-48.0).mul(s6.sub(0.5));
        return du.blend(0.0, r2s.greaterThanEq(3*3));

    }

    private static DoubleVector nearestImage(DoubleVector component, double dimComponent) {
        double half = dimComponent / 2;
        var comp = component;
        var mask = comp.greaterThan(half);
//        while (mask.anyTrue()) {
            comp = comp.sub(dimComponent, mask);
//            mask = comp.greaterThan(half);
//        }

        mask = comp.lessThan(-half);
//        while (mask.anyTrue()) {
            comp = comp.add(dimComponent, mask);
//            mask = comp.lessThan(-half);
//        }

        return comp;
    }

    private double getEnergy() {
        double energy = 0;
        for (int i = 0; i < this.positions.size(); i++) {
            for (int j = i + 1; j < this.positions.size(); j++) {
                Vector3D dr = this.positions.diffSingle(i, j);
                this.box.getBoundary().nearestImage(dr);
                double r2 = dr.squared();
                Potential2Soft potential = potentials[atomTypes[i]][atomTypes[j]];
                energy += potential.u(r2);
            }
        }
        return energy;
    }

    protected void doStepInternal() {
        super.doStepInternal();

//        var oldPos = this.positions.toVectors();
//        System.out.println(this.positions.toVectors().stream().map(v -> v.toString()).collect(Collectors.joining("\n")));
//        System.out.println(this.velocities.toVectors().stream().map(v -> "vel: "+ v.toString()).collect(Collectors.joining("\n")));
//        System.out.println("---");
        this.velocities.addScaled(0.5 * timeStep * this.types[0].rm(), forces);
        this.positions.addScaled(timeStep, velocities);
//        var newPos = this.positions.toVectors();
////        for (int i = 0; i < oldPos.size(); i++) {
////            if (oldPos.get(i).Mv1Squared(newPos.get(i)) > 1e-10) {
////                System.out.println(i + ": " + oldPos.get(i) + " " + newPos.get(i));
////            }
////        }
//        System.exit(1);
//        System.out.println("pos "+this.positions.toVectors().get(0));

        this.computeForces();

        this.velocities.addScaled(0.5 * timeStep * this.types[0].rm(), forces);

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
