package etomica.experimental;

import etomica.atom.AtomType;
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

public class IntegratorVelocityVerletLessFast extends IntegratorMD {
    private final Vector[] forces;
    private final int[] atomTypes;
    private final AtomType[] types;
    private final Potential2Soft[][] potentials;

    public IntegratorVelocityVerletLessFast(PotentialMaster potentialMaster, IRandom random, double timeStep, double temperature, Box box) {
        super(potentialMaster, random, timeStep, temperature, box);
        this.atomTypes = new int[box.getLeafList().size()];


        this.forces = new Vector[box.getLeafList().size()];
        for (int i = 0; i < this.forces.length; i++) {
            forces[i] = new Vector3D();
        }

        for (int i = 0; i < this.box.getLeafList().size(); i++) {
            this.atomTypes[i] = box.getLeafList().get(i).getType().getIndex();
        }

        types = new AtomType[]{this.box.getLeafList().get(0).getType()};
        potentials = new Potential2Soft[][]{{new P2SoftSphericalTruncated(new P2LennardJones(Space3D.getInstance()), 3)}};
        computeForces();
    }

    private void computeForces() {
        for (Vector force : this.forces) {
            force.E(0);
        }
        IAtomList atoms = box.getLeafList();
        Vector3D dr = new Vector3D();
        for (int i = 0; i < atoms.size(); i++) {
            for (int j = i + 1; j < atoms.size(); j++) {
                dr.Ev1Mv2(atoms.get(i).getPosition(), atoms.get(j).getPosition());
                this.box.getBoundary().nearestImage(dr);
                double r2 = dr.squared();
                Potential2Soft potential = potentials[atomTypes[i]][atomTypes[j]];
                double du = potential.du(r2);
                if (du == 0) {
                    continue;
                }

                dr.TE(du / r2);
                this.forces[i].ME(dr);
                this.forces[j].PE(dr);
//                System.out.printf("%d, %f, %f, %s%n", j, du, r2, dr);
            }
//            System.out.println(forces[i]);
//            System.exit(1);
        }
    }

    private double getEnergy() {
        double energy = 0;
        IAtomList atoms = box.getLeafList();
        Vector3D dr = new Vector3D();
        for (int i = 0; i < atoms.size(); i++) {
            for (int j = i + 1; j < atoms.size(); j++) {
                dr.Ev1Mv2(atoms.get(i).getPosition(), atoms.get(j).getPosition());
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

        IAtomList atoms = box.getLeafList();
        for (int i = 0; i < atoms.size(); i++) {

            IAtomKinetic a = ((IAtomKinetic) atoms.get(i));
            a.getVelocity().PEa1Tv1(0.5 * timeStep * this.types[this.atomTypes[i]].rm(), forces[i]);
            a.getPosition().PEa1Tv1(timeStep, a.getVelocity());
        }

        this.computeForces();

        for (int i = 0; i < atoms.size(); i++) {

            IAtomKinetic a = ((IAtomKinetic) atoms.get(i));
            a.getVelocity().PEa1Tv1(0.5 * timeStep * this.types[this.atomTypes[i]].rm(), forces[i]);
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
