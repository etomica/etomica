package etomica.virial;

import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.util.random.IRandom;

/**
 *  Overrides MCMoveAtom to fix the position of 0th atom such that
 *  center of mass (COM) of the system is preserved.
 *
 *  @author Arpit Bansal
 */

public class MCMoveClusterFixCOM extends MCMoveAtom {

    public MCMoveClusterFixCOM(IRandom random, Space _space)
    {
        super(random, null, _space, 0.1, 0.5, false);
    }

    public void setBox(Box box) {
        super.setBox(box);
        if (translationVectors == null) {
            translationVectors = new Vector[box.getLeafList().size() - startAtom];
            for (int i=0; i<translationVectors.length; i++) {
                translationVectors[i] = space.makeVector();
            }
        }
    }

    public void setStartAtom(int newStartAtom) {
        startAtom = newStartAtom;
        if (translationVectors != null && translationVectors.length != box.getLeafList().size() - startAtom) {
            translationVectors = new Vector[box.getLeafList().size() - startAtom];
            for (int i = 0; i < translationVectors.length; i++) {
                translationVectors[i] = space.makeVector();
            }
        }
    }

    public int getStartAtom() {
        return startAtom;
    }

    public boolean doTrial() {
//        System.out.println("Do Trial");
        uOld = ((BoxCluster)box).getSampleCluster().value((BoxCluster)box);
        IAtomList leafAtoms = box.getLeafList();
//        System.out.println(uOld+" Positions before new Move:");
//        for (int i = 0; i < leafAtoms.size(); i++) {
//            System.out.println(leafAtoms.get(i).getPosition());
//        }
        Vector a = leafAtoms.get(0).getPosition();
        for(int i = startAtom; i<leafAtoms.size(); i++) {
            translationVectors[i - startAtom].setRandomCube(random);
            translationVectors[i-startAtom].TE(stepSize);
            a.ME(translationVectors[i - startAtom]);
            Vector r = leafAtoms.get(i).getPosition();
            r.PE(translationVectors[i - startAtom]);
            if (imposePBC) r.PE(box.getBoundary().centralImage(r));
        }
        if (imposePBC) a.PE(box.getBoundary().centralImage(a));
//        System.out.println("Positions proposed by new Move:");
//        for (int i = 0; i < leafAtoms.size(); i++) {
//            System.out.println(leafAtoms.get(i).getPosition());
//        }
        ((BoxCluster)box).trialNotify();
        uNew = ((BoxCluster)box).getSampleCluster().value((BoxCluster)box);
        return true;
    }

    public double getChi(double temperature) {
        return (uOld==0.0) ? Double.POSITIVE_INFINITY : uNew/uOld;
    }

    public void setDoImposePBC(boolean doImposePBC) {
        imposePBC = doImposePBC;
    }

    public void rejectNotify() {
        IAtomList leafAtoms = box.getLeafList();
        Vector a = leafAtoms.get(0).getPosition();
        for(int i = startAtom; i<leafAtoms.size(); i++) {
            a.PE(translationVectors[i - startAtom]);
            Vector r = leafAtoms.get(i).getPosition();
            r.ME(translationVectors[i - startAtom]);
            if (imposePBC) r.PE(box.getBoundary().centralImage(r));
        }
        if (imposePBC) a.PE(box.getBoundary().centralImage(a));
//        System.out.println("Positions after new Move is rejected:");
//        for (int i = 0; i < box.getLeafList().size(); i++) {
//            System.out.println(box.getLeafList().get(i).getPosition());
//        }
        ((BoxCluster)box).rejectNotify();
    }

    public void acceptNotify() {
        super.acceptNotify();
        ((BoxCluster)box).acceptNotify();
//        System.out.println("Positions after new Move is accepted:");
//        for (int i = 0; i < box.getLeafList().size(); i++) {
//            System.out.println(box.getLeafList().get(i).getPosition());
//        }
    }

    protected Vector[] translationVectors;
    protected int startAtom = 1;
    protected boolean imposePBC = false;
}
