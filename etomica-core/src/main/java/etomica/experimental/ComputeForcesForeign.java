package etomica.experimental;

import etomica.atom.IAtomKinetic;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.data.DataSourceScalar;
import etomica.integrator.IntegratorMD;
import etomica.potential.P2LennardJones;
import etomica.potential.P2SoftSphericalTruncated;
import etomica.potential.Potential2Soft;
import etomica.potential.PotentialMaster;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.space3d.Vector3D;
import etomica.units.dimensions.Energy;
import etomica.util.random.IRandom;

import java.foreign.Libraries;
import java.foreign.NativeTypes;
import java.foreign.Scope;
import java.foreign.layout.Layout;
import java.foreign.memory.Array;
import java.foreign.memory.Pointer;
import java.io.IOException;
import java.lang.invoke.MethodHandles;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.stream.Collectors;
import java.util.stream.IntStream;


public class ComputeForcesForeign extends IntegratorMD implements EnergyMeter {
    private final Scope sc;
    private final Array<Double> xs;
    private final Array<Double> ys;
    private final Array<Double> zs;

    private final Array<Double> fxs;
    private final Array<Double> fys;
    private final Array<Double> fzs;

    private final Array<Double> vxs;
    private final Array<Double> vys;
    private final Array<Double> vzs;

    private final Potential2Soft[][] potentials;

    static etomica.foreign.computeForces cf = Libraries.bind(etomica.foreign.computeForces.class, Libraries.load(MethodHandles.lookup(), "/Users/alex/workspace/etomica/computeForcesForeign/computeForces.dylib"));

    /**
     * Constructs integrator with a default for non-isothermal sampling.
     *
     * @param potentialMaster PotentialMaster instance used to compute the energy and forces
     * @param random          random number generator used for initial velocities and some thermostats
     * @param timeStep        time step for integration
     * @param temperature     used by thermostat and/or to initialize velocities
     * @param box
     */
    public ComputeForcesForeign(PotentialMaster potentialMaster, IRandom random, double timeStep, double temperature, Box box, Scope sc) {
        super(potentialMaster, random, timeStep, temperature, box);
        this.sc = sc;
        potentials = new Potential2Soft[][]{{new P2SoftSphericalTruncated(new P2LennardJones(Space3D.getInstance()), 3)}};

        // uncomment if using own velocities
        this.randomizeMomenta();
        this.shiftMomenta();
        this.scaleMomenta();

        var atoms = box.getLeafList();
        xs = sc.allocateArray(NativeTypes.DOUBLE, atoms.size());
        ys = sc.allocateArray(NativeTypes.DOUBLE, atoms.size());
        zs = sc.allocateArray(NativeTypes.DOUBLE, atoms.size());

        fxs = sc.allocateArray(NativeTypes.DOUBLE, atoms.size());
        fys = sc.allocateArray(NativeTypes.DOUBLE, atoms.size());
        fzs = sc.allocateArray(NativeTypes.DOUBLE, atoms.size());

        vxs = sc.allocateArray(NativeTypes.DOUBLE, atoms.size());
        vys = sc.allocateArray(NativeTypes.DOUBLE, atoms.size());
        vzs = sc.allocateArray(NativeTypes.DOUBLE, atoms.size());

        for (int i = 0; i < atoms.size(); i++) {
            xs.set(i, atoms.get(i).getPosition().getX(0));
            ys.set(i, atoms.get(i).getPosition().getX(1));
            zs.set(i, atoms.get(i).getPosition().getX(2));

            var v = ((IAtomKinetic) atoms.get(i)).getVelocity();
            vxs.set(i, v.getX(0));
            vys.set(i, v.getX(1));
            vzs.set(i, v.getX(2));

            fxs.set(i, 0.0);
            fys.set(i, 0.0);
            fzs.set(i, 0.0);
        }
        computeForcesJava();
    }

    private void computeForcesJava() {
        cf.computeForces(
                xs.elementPointer(), ys.elementPointer(), zs.elementPointer(),
                fxs.elementPointer(), fys.elementPointer(), fzs.elementPointer(),
                this.box.getLeafList().size(), this.box.getBoundary().getBoxSize().getX(0));

    }

    @Override
    protected void doStepInternal() {
        super.doStepInternal();

//        var atoms = box.getLeafList();
//
//        for (int i = 0; i < atoms.size(); i++) {
//            var a = ((IAtomKinetic) atoms.get(i));
//            var vel = a.getVelocity();
//            if (i == 0) {
//                System.out.println("Before: " + vel);
//            }
//            var fVec = new Vector3D(fxs.get(i), fys.get(i), fzs.get(i));
//            vel.PEa1Tv1(0.5 * timeStep * a.getType().rm(), fVec);
//            xs.set(i, xs.get(i) + a.getVelocity().getX(0) * timeStep);
//            ys.set(i, ys.get(i) + a.getVelocity().getX(1) * timeStep);
//            zs.set(i, zs.get(i) + a.getVelocity().getX(2) * timeStep);
//            if (i == 0) {
//                System.out.println("Force: " + fVec);
//                System.out.println("After: " + vel);
//                System.out.printf("New pos: (%s %s %s)%n", xs.get(i), ys.get(i), zs.get(i));
//            }
//        }
//
//        try {
//            Files.write(Paths.get("stuff_native_" + stepCount + ".txt"),
//                    IntStream.range(0, atoms.size()).<String>mapToObj(i -> {
//                        return String.format("force (%s %s %s) pos (%s %s %s) vel %s", fxs.get(i), fys.get(i), fzs.get(i), xs.get(i), ys.get(i), zs.get(i), ((IAtomKinetic) atoms.get(i)).getVelocity());
//                    }).collect(Collectors.toList())
//            );
//        } catch (IOException e) {
//            e.printStackTrace();
//        }
        cf.updatePreForces(
                xs.elementPointer(), ys.elementPointer(), zs.elementPointer(),
                vxs.elementPointer(), vys.elementPointer(), vzs.elementPointer(),
                fxs.elementPointer(), fys.elementPointer(), fzs.elementPointer(),
                box.getLeafList().size(), box.getLeafList().get(0).getType().rm(), timeStep);

        computeForcesJava();

        cf.updatePostForces(
                vxs.elementPointer(), vys.elementPointer(), vzs.elementPointer(),
                fxs.elementPointer(), fys.elementPointer(), fzs.elementPointer(),
                box.getLeafList().size(), box.getLeafList().get(0).getType().rm(), timeStep);

//        try {
//            Files.write(Paths.get("forces_native.txt"),
//                    IntStream.range(0, atoms.size()).<String>mapToObj(i -> String.format("(%f %f %f)", fxs.get(i), fys.get(i), fzs.get(i))).collect(Collectors.toList())
//            );
//        } catch (IOException e) {
//            e.printStackTrace();
//        }

//        for (int i = 0; i < atoms.size(); i++) {
//            var a = ((IAtomKinetic) atoms.get(i));
//            var fVec = new Vector3D(fxs.get(i), fys.get(i), fzs.get(i));
//            a.getVelocity().PEa1Tv1(0.5 * timeStep * a.getType().rm(), fVec);
//            if (i == 0) {
//                System.out.println("Force: " + fVec + " vel: " + a.getVelocity());
//            }
//        }
//        System.out.println(getEnergy());
//        System.exit(1);
    }

    private double getEnergy() {

        double energy = 0;
        IAtomList atoms = box.getLeafList();
        Vector3D dr = new Vector3D();
        for (int i = 0; i < atoms.size(); i++) {
            for (int j = i + 1; j < atoms.size(); j++) {
                dr.Ev1Mv2(
                        Vector.of(this.xs.get(i), this.ys.get(i), this.zs.get(i)),
                        Vector.of(this.xs.get(j), this.ys.get(j), this.zs.get(j))
                        );
                this.box.getBoundary().nearestImage(dr);
                double r2 = dr.squared();
                Potential2Soft potential = potentials[0][0];
                energy += potential.u(r2);
            }
        }
        return energy;
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
