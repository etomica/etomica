/*
 * Created on Oct 4, 2005
 *
 * TODO To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package etomica.virial;

import etomica.action.AtomTransform;
import etomica.atom.Atom;
import etomica.atom.AtomList;
import etomica.atom.AtomTreeNodeGroup;
import etomica.integrator.mcmove.MCMoveRotateMolecule3D;
import etomica.models.water.AtomTreeNodeWater;
import etomica.phase.Phase;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.RotationTensor;
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
            Space space, int numMolecules, double bondDistance, double bondAngle) {
        super(potentialMaster, space);
        weightMeter = new MeterClusterWeight(potential);
        setName("MCMoveClusterMolecule");
        nMolecules = numMolecules;
        selectedMolecules = new Atom[nMolecules];
        oldPositions = new Vector[nMolecules][2];
        for (int i=0; i<nMolecules; i++) {
            for (int j=0; j<2; j++) {
                oldPositions[i][j] = space.makeVector();
            }
        }
        distance = bondDistance;
        cosAngle = Math.cos(bondAngle);
        sinAngle = Math.sin(bondAngle);
        work = space.makeVector();
    }
    
    public void setPhase(Phase[] p) {
        super.setPhase(p);
        weightMeter.setPhase(p[0]);
        selectMolecules();
    }

    public boolean doTrial() {
        uOld = weightMeter.getDataAsScalar();
        for (int i=0; i<selectedMolecules.length; i++) {
            molecule = selectedMolecules[i];
            AtomTreeNodeWater waterNode= (AtomTreeNodeWater)molecule.node;
            Atom O = waterNode.O;
            Atom H1 = waterNode.H1;
            Atom H2 = waterNode.H2;
//            if (foo-- == 0) {
//                foo = 10000;
//                Vector v1 = (Vector)H1.coord.position().clone();
//                v1.ME(O.coord.position());
//                Vector v2 = (Vector)H2.coord.position().clone();
//                v2.ME(O.coord.position());
//                double d = v1.dot(v2);
//                System.out.println(i+"d = "+(d+0.3338068592338)+" "+(Math.sqrt(v1.squared())-1)+" "+(Math.sqrt(v2.squared())-1));
//            }
        
            double dTheta = (2*Simulation.random.nextDouble() - 1.0)*stepSize;
            rotationTensor.setAxial(Simulation.random.nextInt(3),dTheta);
            oldPositions[i][0].E(H1.coord.position());
            oldPositions[i][1].E(H2.coord.position());
            leafAtomIterator.setRoot(molecule);
            leafAtomIterator.reset();
            r0.E(O.coord.position());
//            System.out.println(molecule+" starting at "+molecule.node.lastLeafAtom().coord.position());
            AtomTransform.doTransform(leafAtomIterator, r0, rotationTensor);
//            System.out.println(molecule+" moved to "+molecule.node.lastLeafAtom().coord.position());
            
            if (foo2-- == 0) {
                foo2 = 100;
                // normalize OH1
                Vector p1 = H1.coord.position();
                p1.ME(O.coord.position());
                p1.TE(1/Math.sqrt(p1.squared()));
                Vector p2 = H2.coord.position();
                p2.ME(O.coord.position());
                p2.TE(1/Math.sqrt(p2.squared()));
                // move H2 to fix bond angle
                double d = p1.dot(p2);
                work.Ev1Pa1Tv2(p2,-d,p1);
                work.TE(1/Math.sqrt(work.squared()));
                p2.Ea1Tv1(sinAngle,work);
                p2.PEa1Tv1(cosAngle,p1);
                p2.TE(distance/Math.sqrt(p2.squared()));
                p1.TE(distance);
                p1.PE(O.coord.position());
                p2.PE(O.coord.position());
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
            AtomTreeNodeWater waterNode = (AtomTreeNodeWater)molecule.node;
            waterNode.H1.coord.position().E(oldPositions[i][0]);
            waterNode.H2.coord.position().E(oldPositions[i][1]);
//            System.out.println(molecule+" back to "+molecule.node.lastLeafAtom().coord.position());
        }
        ((PhaseCluster)phases[0]).rejectNotify();
    }
    
    private final MeterClusterWeight weightMeter;
    private final Atom[] selectedMolecules;
    private final Vector[][] oldPositions;
    private final int nMolecules;
    private int foo = 10000, foo2 = 100;
    private final double distance;
    private final double cosAngle, sinAngle;
    private final Vector work;
}
