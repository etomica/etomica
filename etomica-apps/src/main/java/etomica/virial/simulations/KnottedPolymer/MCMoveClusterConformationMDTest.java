package etomica.virial.simulations.KnottedPolymer;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.data.meter.MeterRadiusGyration;
import etomica.molecule.IMolecule;
import etomica.molecule.MoleculePositionCOM;
import etomica.space.RotationTensor;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.RotationTensor3D;
import etomica.space3d.Vector3D;
import etomica.virial.BoxCluster;
import etomica.virial.MCMoveClusterMolecule;
import etomica.virial.simulations.SimulationVirialOverlap2;

import java.util.ArrayList;
import java.util.List;

public class MCMoveClusterConformationMDTest extends MCMoveClusterMolecule {
    private int nBead;
    protected transient Vector r0;
    protected transient RotationTensor rotationTensor;
    private MoleculePositionCOM molCOM;
    private List<List<Vector3D>> oldPos = new ArrayList<>(2);
    public ArrayList<List<Vector>> CoordComboList = new ArrayList<>();
    private StarPolymerMD simMD;


    public MCMoveClusterConformationMDTest(SimulationVirialOverlap2 sim, Space _space, double temperature, int f, int l,
                                           boolean useNbrs, boolean doCollecting) {
        super(sim, _space);
        rotationTensor = _space.makeRotationTensor();
        r0 = _space.makeVector();
        molCOM = new MoleculePositionCOM(_space);
        nBead = f * l + 1;
        if (doCollecting) {
            double t1 = System.currentTimeMillis();
            equilibrating(f, l, temperature, 0.005, useNbrs);
            collecting(2000);
            double t2 = System.currentTimeMillis();
            System.out.println("time spent on MD: " + (t2 - t1) / 1000.0);
        }
    }

    public void equilibrating(int f, int l, double temperature, double tStep, boolean useNbrs) {
        this.simMD = new StarPolymerMD(f, l, temperature, tStep, useNbrs);
        long step = 1000000;
        simMD.ai.setMaxSteps(step);
        System.out.println("Equilibrating MD for " + step + " steps ....");
        simMD.ai.actionPerformed();
        System.out.print("Finished it!\n");
    }

    public void collecting(int sizeConformations) {
        int steps = 1000;
        System.out.println("Running MD production for " + sizeConformations + " Conformations ...");

        Box boxMD = simMD.getBox(0);

        MeterRadiusGyration meter = new MeterRadiusGyration(simMD.getSpace());
        meter.setBox(boxMD);
        double sumRG = 0;
        int count = 0;

        for (int i = 0; i < sizeConformations; i++) {
            simMD.ai.setMaxSteps(steps);
            simMD.ai.actionPerformed();
            IAtomList atomList1 = boxMD.getMoleculeList().get(0).getChildList();
            ArrayList<Vector> pos1 = new ArrayList<>();
            for (IAtom atom : atomList1) {
                pos1.add(atom.getPosition());
            }
            CoordComboList.add(pos1);

            sumRG += meter.getDataAsScalar();
            count++;
        }
        System.out.println("Finished collecting conformations");
        System.out.println("Averaged squared radius of gyration: " + sumRG / count);
    }

    public ArrayList<List<Vector>> getConformationsList() {
        if (CoordComboList.size() != 0) {
            return CoordComboList;
        } else {
            throw new RuntimeException("Collecting needed");
        }
    }

    public void setConformationsList(ArrayList<List<Vector>> confList) {
        this.CoordComboList = confList;
    }

    @Override
    public boolean doTrial() {
        if (CoordComboList.size() == 0) throw new RuntimeException("Collecting needed before doing trials");

        int nBody = box.getMoleculeList().size();
        uOld = ((BoxCluster) box).getSampleCluster().value((BoxCluster) box);
//        System.out.println("uOld: " + uOld);
        for (int i = 0; i < nBody; i++) {
            IMolecule mol = box.getMoleculeList().get(i);
            List<IAtom> atoms = mol.getChildList().getAtoms();
            List<Vector3D> old = new ArrayList<>();
            Vector oldCOM = space.makeVector();
            oldCOM.E(molCOM.position(mol));

            int rNo = random.nextInt(CoordComboList.size());
            List<Vector> newPos = CoordComboList.get(rNo);
//            CoordComboList.remove(rNo);

            for (int j = 0; j < nBead; j++) {
                //Storing old coordinates
                Vector3D pos = new Vector3D();
                pos.E(atoms.get(j).getPosition());
                old.add(pos);
//                System.out.println("old pos for atom: " + j + pos);

                //Assigning new coordinates
                atoms.get(j).getPosition().E(newPos.get(j));
//                System.out.println("new pos for atom: " + j + atoms.get(j).getPosition());
            }

            Vector newCOM = molCOM.position(mol);

//            System.out.println("old com" + oldCOM);
//            System.out.println("new com" + newCOM);

            groupTranslationVector.E(oldCOM);
            groupTranslationVector.ME(newCOM);
            moveMoleculeAction.actionPerformed(mol);
//            System.out.println("com" + molCOM.position(mol));

            oldPos.add(old);

            molecule = box.getMoleculeList().get(i);
            r0.E(oldCOM);


            // Rotation begins here

//            double dTheta = (2 * random.nextDouble() - 1.0) * stepSize;
//            rotationTensor.setAxial(random.nextInt(3), dTheta);


            double[] quat = new double[4];
            double u1 = random.nextDouble();
            double u2 = 2 * Math.PI * random.nextDouble();
            double u3 = 2 * Math.PI * random.nextDouble();
            double s1 = Math.sqrt(u1);
            double s2 = Math.sqrt(1 - u1);
            quat[0] = s1 * Math.sin(u2);
            quat[1] = s1 * Math.cos(u2);
            quat[2] = s2 * Math.sin(u3);
            quat[3] = s2 * Math.cos(u3);

            rotationTensor = new RotationTensor3D();
            ((RotationTensor3D) rotationTensor).setQuaternions(quat);

            doTransform();

//            if (trialCount-- == 0) {
//                relaxAction.actionPerformed(molecule);
//                trialCount = relaxInterval;
//            }

        }

        ((BoxCluster) box).trialNotify();
        uNew = ((BoxCluster) box).getSampleCluster().value((BoxCluster) box);

//        System.out.println("uNew is "+uNew);
//        if (uNew == 0 ) throw new RuntimeException();
        return true;
    }

    protected void doTransform() {
        IAtomList childList = molecule.getChildList();
        for (int iChild = 0; iChild < childList.size(); iChild++) {
            IAtom a = childList.get(iChild);
            Vector r = a.getPosition();
            r.ME(r0);
            box.getBoundary().nearestImage(r);
            rotationTensor.transform(r);
            r.PE(r0);
        }
    }

    @Override
    public void acceptNotify() {
        ((BoxCluster) box).acceptNotify();
//        System.out.println("accepted");
//        System.out.println(this.getTracker().acceptanceProbability());
//        System.out.println(this.getTracker().acceptanceRatio());
        oldPos.clear();
    }

    @Override
    public void rejectNotify() {
        super.rejectNotify();

        for (int i = 0; i < box.getMoleculeList().size(); i++) {
            IMolecule mol = box.getMoleculeList().get(i);
            List<IAtom> atoms = mol.getChildList().getAtoms();
            List<Vector3D> old = oldPos.get(i);

            for (int j = 0; j < nBead; j++) {
                atoms.get(j).getPosition().E(old.get(j));
            }
        }

        oldPos.clear();
    }
}
