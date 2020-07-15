package etomica.starpolymer;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.integrator.mcmove.MCMoveMolecule;
import etomica.molecule.IMolecule;
import etomica.molecule.MoleculePositionCOM;
import etomica.simulation.Simulation;
import etomica.space.RotationTensor;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.RotationTensor3D;
import etomica.space3d.Vector3D;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class MCMoveUpdateConformation extends MCMoveMolecule {
    final String fileName;
    final int nBead = 201;
    final int conformationSize = 10000;
    protected transient Vector r0;
    protected transient RotationTensor rotationTensor;
    protected MoleculePositionCOM molCOM;
    private List<List<Vector3D>> randomCoordList;
    private int rNo = 0;

    public MCMoveUpdateConformation(Simulation sim, Space _space, String fileName) {
        super(sim, null, _space);
        this.fileName = fileName;
        rotationTensor = _space.makeRotationTensor();
        r0 = _space.makeVector();
        molCOM = new MoleculePositionCOM(_space);
        randomCoordList = getRandomCoordListBuffer(conformationSize);
    }

    @Override
    public boolean doTrial() {
        int nBody = box.getMoleculeList().size();

//        System.out.println("uOld: " + uOld);
        for (int i = 0; i < nBody; i++) {
//            if (rNo < conformationSize - 1){
//                rNo++;
//            }else{
//                throw new RuntimeException("OOps");
//            }
            int rNo = random.nextInt(conformationSize);
//            System.out.println("conformation "+rNo);
            List<Vector3D> newPos = randomCoordList.get(rNo);
            IMolecule mol = box.getMoleculeList().get(i);
            List<IAtom> atoms = mol.getChildList().getAtoms();
            List<Vector3D> old = new ArrayList<>();
            Vector oldCOM = space.makeVector();
            oldCOM.E(molCOM.position(mol));

            for (int j = 0; j < nBead; j++) {

                //Assigning new coordinates
                atoms.get(j).getPosition().E(newPos.get((j + nBead - 1) % nBead));
//                System.out.println("new pos for atom: " + j + atoms.get(j).getPosition());
            }

            Vector newCOM = molCOM.position(mol);

//            System.out.println("old com" + oldCOM);
//            System.out.println("new com" + newCOM);

            groupTranslationVector.E(oldCOM);
            groupTranslationVector.ME(newCOM);
//            moveMoleculeAction.actionPerformed(mol);
////            System.out.println("com" + molCOM.position(mol));


            molecule = box.getMoleculeList().get(i);
            r0.E(oldCOM);


            // Rotation begins here

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
        }

        return true;
    }

    public double getChi(double temperature) {
        return 1;
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
    }

    @Override
    public void rejectNotify() {
        throw new RuntimeException("oops");
    }
}
