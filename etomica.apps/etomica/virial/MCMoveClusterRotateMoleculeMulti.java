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
import etomica.atom.AtomArrayList;
import etomica.atom.AtomLeaf;
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
public class MCMoveClusterRotateMoleculeMulti extends MCMoveRotateMolecule3D {

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
    
    public void setPhase(Phase p) {
        super.setPhase(p);
        weightMeter.setPhase(p);
        selectMolecules();
        for (int i=0; i<nMolecules; i++) {
            molecule = selectedMolecules[i];
            oldPositions[i] = new Vector[((AtomTreeNodeGroup)molecule.getNode()).getChildList().size()];
            for (int j=0; j<oldPositions[i].length; j++) {
                oldPositions[i][j] = p.space().makeVector();
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
            r0.E(molecule.getType().getPositionDefinition().position(molecule));
//            System.out.println(molecule+" before position "+r0);
        
            double dTheta = (2*Simulation.random.nextDouble() - 1.0)*stepSize;
            rotationTensor.setAxial(Simulation.random.nextInt(3),dTheta);
            
            int j=0;
            while (leafAtomIterator.hasNext()) {
                oldPositions[i][j++].E(((AtomLeaf)leafAtomIterator.nextAtom()).coord.position());
            }
            leafAtomIterator.reset();
//            System.out.println(molecule+" starting at "+molecule.node.lastLeafAtom().coord.position());
            AtomTransform.doTransform(leafAtomIterator, r0, rotationTensor);
//            System.out.println(molecule+" moved to "+molecule.node.lastLeafAtom().coord.position());
            
            
            if (doRelax && relaxAction != null) {
                relaxAction.setAtom(molecule);
                relaxAction.actionPerformed();
            }
        }

        ((PhaseCluster)phase).trialNotify();
        uNew = weightMeter.getDataAsScalar();
        return true;
    }
    
    public void selectMolecules() {
        AtomArrayList atomList = ((AtomTreeNodeGroup)((AtomTreeNodeGroup)phase.getSpeciesMaster().getNode()).getChildList().get(0).getNode()).getChildList();
        System.arraycopy(atomList.toArray(),1,selectedMolecules,0,atomList.size()-1);
    }

    public double getB() {
        return 0.0;
    }
    
    public double getA() {
        return (uOld==0.0) ? Double.POSITIVE_INFINITY : uNew/uOld;
    }
    
    public void acceptNotify() {
        super.acceptNotify();
        ((PhaseCluster)phase).acceptNotify();
    }
    
    public void rejectNotify() {
//        super.rejectNotify();
        for (int i=0; i<selectedMolecules.length; i++) {
            molecule = selectedMolecules[i];
            leafAtomIterator.setRoot(molecule);
            leafAtomIterator.reset();
            int j=0;
            while (leafAtomIterator.hasNext()) {
                ((AtomLeaf)leafAtomIterator.nextAtom()).coord.position().E(oldPositions[i][j++]);
            }
        }
        ((PhaseCluster)phase).rejectNotify();
    }
    
    public void setRelaxAction(AtomAction action) {
        relaxAction = action;
    }
    
    private static final long serialVersionUID = 1L;
    private final MeterClusterWeight weightMeter;
    private final Atom[] selectedMolecules;
    private final Vector[][] oldPositions;
    private final int nMolecules;
    private int trialCount, relaxInterval = 100;
    private AtomAction relaxAction;
}
