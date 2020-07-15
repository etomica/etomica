package etomica.starpolymer;

import etomica.action.MoleculeAction;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.molecule.IMolecule;
import etomica.molecule.MoleculePositionCOM;
import etomica.simulation.Simulation;
import etomica.space.RotationTensor;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.RotationTensor3D;
import etomica.space3d.Vector3D;
import etomica.virial.BoxCluster;
import etomica.virial.MCMoveClusterMolecule;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class MCMoveClusterUpdateConformation extends MCMoveClusterMolecule {
    final String fileName;
    final int nBead = 201;
    final int conformationSize = 10000;
    protected transient Vector r0;
    protected transient RotationTensor rotationTensor;
    protected MoleculePositionCOM molCOM;
    private List<List<Vector3D>> oldPos = new ArrayList<List<Vector3D>>(2);
    private List<Long> bytePositionList;
    private List<List<Vector3D>> randomCoordList;

    protected int trialCount, relaxInterval = 100;
    protected MoleculeAction relaxAction;

    public MCMoveClusterUpdateConformation(Simulation sim, Space _space, String fileName) {
        super(sim, _space);
        this.fileName = fileName;
        rotationTensor = _space.makeRotationTensor();
        r0 = _space.makeVector();
        molCOM = new MoleculePositionCOM(_space);
        randomCoordList = getRandomCoordListBuffer(conformationSize);
    }

    @Override
    public boolean doTrial() {
        int nBody = box.getMoleculeList().size();

//        if (randomCoordList.isEmpty()){
//            randomCoordList = getRandomCoordList(1000);
//        }
        uOld = ((BoxCluster) box).getSampleCluster().value((BoxCluster) box);
//        System.out.println("uOld: " + uOld);
        for (int i = 0; i < nBody; i++) {
            int rNo = random.nextInt(conformationSize);
//            randomCoordList.remove(randomCoordList.size()-1);
            List<Vector3D> newPos = randomCoordList.get(rNo);
            IMolecule mol = box.getMoleculeList().get(i);
            List<IAtom> atoms = mol.getChildList().getAtoms();
            List<Vector3D> old = new ArrayList<>();
            Vector oldCOM = space.makeVector();
            oldCOM.E(molCOM.position(mol));

            for (int j = 0; j < nBead; j++) {
                //Storing old coordinates
                Vector3D pos = new Vector3D();
                pos.E(atoms.get(j).getPosition());
                old.add(pos);
//                System.out.println("old pos for atom: " + j + pos);

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

    public List<List<Vector3D>> getRandomCoordListBuffer(int size) {
        List<List<Vector3D>> randomCoordList = new ArrayList<>();
        FileReader fileReader;
        try {
            fileReader = new FileReader(fileName);
        } catch (IOException e) {
            throw new RuntimeException("Cannot open " + fileName + ", caught IOException: " + e.getMessage());
        }

        try {
            BufferedReader bufReader = new BufferedReader(fileReader);


            for (int i = 0; i < size; i++) {
                //Start reading coordinates now
                List<Vector3D> coordinates = new ArrayList<>();
                bufReader.readLine();
                bufReader.readLine();
                for (int iBead = 0; iBead < nBead; iBead++) {
                    String[] coordStr = bufReader.readLine().split("[ ]+");
                    double x = Double.parseDouble(coordStr[1]);
                    double y = Double.parseDouble(coordStr[2]);
                    double z = Double.parseDouble(coordStr[3]);
                    Vector3D coord = new Vector3D();
                    coord.E(x, y, z);
                    coordinates.add(coord);
                }

                if (coordinates.size() != nBead) throw new RuntimeException("Error on reading coords!");
                randomCoordList.add(coordinates);
            }

            if (randomCoordList.size() != size) throw new RuntimeException("Didn't read all the coords!");
            fileReader.close();
        } catch (IOException e1) {
            e1.printStackTrace();
        }

        return randomCoordList;
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
