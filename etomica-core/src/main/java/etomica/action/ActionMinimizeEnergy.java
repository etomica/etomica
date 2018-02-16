package etomica.action;

import etomica.atom.AtomLeafAgentManager;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.atom.IAtomOriented;
import etomica.box.Box;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.integrator.Integrator;
import etomica.potential.IteratorDirective;
import etomica.potential.PotentialCalculationTorqueSum;
import etomica.potential.PotentialMaster;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Orientation3D;

public class ActionMinimizeEnergy implements IAction, AtomLeafAgentManager.AgentSource<ActionMinimizeEnergy.MyForceTorque> {

    protected final int maxSteps;
    protected final Box box;
    protected final PotentialMaster potentialMaster;
    protected final PotentialCalculationTorqueSum pc;
    protected final AtomLeafAgentManager<MyForceTorque> agentManager;
    protected final Space space;
    protected final IteratorDirective id;
    protected double tolerance;

    public ActionMinimizeEnergy(Space space, int maxSteps, PotentialMaster potentialMaster, Box box) {
        this.space = space;
        this.maxSteps = maxSteps;
        this.potentialMaster = potentialMaster;
        this.box = box;
        id = new IteratorDirective();
        pc = new PotentialCalculationTorqueSum();
        agentManager = new AtomLeafAgentManager<>(this, box);
        pc.setAgentManager(agentManager);
    }

    public double getTolerance() {
        return tolerance;
    }

    public void setTolerance(double newTolerance) {
        tolerance = newTolerance;
    }

    protected void takeStep(Vector[] direction, Vector[] rotationDirection, double stepSize) {
        IAtomList atoms = box.getLeafList();
        int numAtoms = atoms.getAtomCount();

        for (int i = 0; i < numAtoms; i++) {
            IAtom atom = atoms.getAtom(i);
            atom.getPosition().PEa1Tv1(stepSize, direction[i]);
            if (atom instanceof IAtomOriented) {
                if (space.getD() == 2) {
                    throw new RuntimeException("hahahahahahah.  you lose");
                }
                ((Orientation3D) ((IAtomOriented) atom).getOrientation()).rotateBy(stepSize, rotationDirection[i]);
            }
        }
    }

    @Override
    public void actionPerformed() {
        final MeterPotentialEnergy meterPotentialEnergy = new MeterPotentialEnergy(potentialMaster, box);
        IAtomList atoms = box.getLeafList();
        int numAtoms = atoms.getAtomCount();
        Vector[] direction = space.makeVectorArray(numAtoms);
        Vector[] rotationDirection = space.makeVectorArray(numAtoms);
        double u0 = meterPotentialEnergy.getDataAsScalar();
        double uLast = u0;
        double step0 = 1e-5;
        for (int iter = 0; iter < maxSteps; iter++) {
            pc.reset();
            potentialMaster.calculate(box, id, pc);
            double dSum = 0;
            for (int i = 0; i < numAtoms; i++) {
                IAtom atom = atoms.getAtom(i);
                direction[i].E(agentManager.getAgent(atom).force);
                dSum += direction[i].squared();
                if (atom instanceof IAtomOriented) {
                    rotationDirection[i].E(agentManager.getAgent(atom).torque);
                    dSum += rotationDirection[i].squared();
                }
            }
            dSum = Math.sqrt(dSum);
            for (int i = 0; i < numAtoms; i++) {
                IAtom atom = atoms.getAtom(i);
                direction[i].TE(1 / dSum);
                if (atom instanceof IAtomOriented) {
                    rotationDirection[i].TE(1 / dSum);
                }
            }
            takeStep(direction, rotationDirection, step0);

            pc.reset();
            potentialMaster.calculate(box, id, pc);
            double dSum2 = 0;
            for (int i = 0; i < numAtoms; i++) {
                IAtom atom = atoms.getAtom(i);
                dSum2 += direction[i].dot(agentManager.getAgent(atom).force);
                if (atom instanceof IAtomOriented) {
                    dSum2 += rotationDirection[i].dot(agentManager.getAgent(atom).torque);
                }
            }
            dSum2 = Math.sqrt(dSum2);
            double step2 = -dSum2 * step0 / (dSum2 - dSum);
            takeStep(direction, rotationDirection, step2);
            double u2 = meterPotentialEnergy.getDataAsScalar();
            if (u2 > uLast + tolerance) {
                break;
            }
        }
    }

    public Box getBox() {
        return box;
    }

    public void releaseAgent(MyForceTorque agent, IAtom atom, Box box) {
    }

    public MyForceTorque makeAgent(IAtom atom, Box box) {
        return new MyForceTorque(space);
    }

    public static class MyForceTorque implements Integrator.Torquable, Integrator.Forcible {

        protected final Vector force, torque;

        public MyForceTorque(Space space) {
            force = space.makeVector();
            torque = space.makeVector();
        }

        @Override
        public Vector torque() {
            return torque;
        }

        @Override
        public Vector force() {
            return force;
        }
    }
}
