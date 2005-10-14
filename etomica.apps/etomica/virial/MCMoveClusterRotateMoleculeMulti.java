/*
 * Created on Oct 4, 2005
 *
 * TODO To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package etomica.virial;

import etomica.action.AtomAction;
import etomica.action.AtomTransform;
import etomica.atom.Atom;
import etomica.atom.AtomList;
import etomica.atom.AtomTreeNodeGroup;
import etomica.integrator.mcmove.MCMoveRotateMolecule3D;
import etomica.phase.Phase;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space.Vector;

/**
 * @author andrew
 *
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
public class MCMoveClusterRotateMoleculeMulti extends MCMoveRotateMolecule3D
        implements MCMoveCluster {

    /**
     * @param potentialMaster
     * @param space
     */
    public MCMoveClusterRotateMoleculeMulti(PotentialMaster potentialMaster,
            Space space, int numMolecules) {
        super(potentialMaster, space);
        weightMeter = new MeterClusterWeight(potential);
        setName("MCMoveClusterMolecule");
        nMolecules = numMolecules;
        selectedMolecules = new Atom[nMolecules];
        oldPositions = new Vector[nMolecules][];
    }
    
    public void setPhase(Phase[] p) {
        super.setPhase(p);
        weightMeter.setPhase(p[0]);
        selectMolecules();
        for (int i=0; i<nMolecules; i++) {
            molecule = selectedMolecules[i];
            oldPositions[i] = new Vector[((AtomTreeNodeGroup)molecule.node).childList.size()-1];
            for (int j=0; j<oldPositions[i].length; j++) {
                oldPositions[i][j] = p[0].space().makeVector();
            }
        }
    }

    public boolean doTrial() {
        uOld = weightMeter.getDataAsScalar();
        boolean doRelax = false;
        if (trialCount-- == 0) {
            doRelax = true;
            trialCount = relaxInterval;
        }
        for (int i=0; i<selectedMolecules.length; i++) {
            molecule = selectedMolecules[i];
            leafAtomIterator.setRoot(molecule);
            leafAtomIterator.reset();
            Atom O = leafAtomIterator.nextAtom();
        
            double dTheta = (2*Simulation.random.nextDouble() - 1.0)*stepSize;
            rotationTensor.setAxial(Simulation.random.nextInt(3),dTheta);
            
            int j=0;
            while (leafAtomIterator.hasNext()) {
                oldPositions[i][j++].E(leafAtomIterator.nextAtom().coord.position());
            }
            leafAtomIterator.reset();
            r0.E(O.coord.position());
//            System.out.println(molecule+" starting at "+molecule.node.lastLeafAtom().coord.position());
            AtomTransform.doTransform(leafAtomIterator, r0, rotationTensor);
//            System.out.println(molecule+" moved to "+molecule.node.lastLeafAtom().coord.position());
            
            
            if (doRelax && relaxAction != null) {
                relaxAction.setAtom(molecule);
                relaxAction.actionPerformed();
            }
        }
            
        uNew = Double.NaN;
        ((PhaseCluster)phases[0]).trialNotify(null);
        return true;
    }
    
    public void selectMolecules() {
        AtomList atomList = ((AtomTreeNodeGroup)phases[0].getSpeciesMaster().firstSpecies().node).childList;
        int total=atomList.size();
        for(int i=total-1; i>total-nMolecules-1; i--) {
            selectedMolecules[total-1-i] = atomList.get(i);
        }
    }

    public double trialRatio() {return 1.0;}
    
    public double lnProbabilityRatio() {
        return Math.log(probabilityRatio());
    }
    
    public double probabilityRatio() {
        uNew = weightMeter.getDataAsScalar();
        return (uOld==0.0) ? Double.POSITIVE_INFINITY : uNew/uOld;
    }
    
    public void acceptNotify() {
        super.acceptNotify();
        ((PhaseCluster)phases[0]).acceptNotify();
    }
    
    public void rejectNotify() {
//        super.rejectNotify();
        for (int i=0; i<selectedMolecules.length; i++) {
            molecule = selectedMolecules[i];
            leafAtomIterator.setRoot(molecule);
            leafAtomIterator.reset();
            // skip first Atom
            leafAtomIterator.nextAtom();
            int j=0;
            while (leafAtomIterator.hasNext()) {
                leafAtomIterator.nextAtom().coord.position().E(oldPositions[i][j++]);
            }
        }
        ((PhaseCluster)phases[0]).rejectNotify();
    }
    
    public void setRelaxAction(AtomAction action) {
        relaxAction = action;
    }
    
    private final MeterClusterWeight weightMeter;
    private final Atom[] selectedMolecules;
    private final Vector[][] oldPositions;
    private final int nMolecules;
    private int trialCount, relaxInterval = 100;
    private AtomAction relaxAction;
}
