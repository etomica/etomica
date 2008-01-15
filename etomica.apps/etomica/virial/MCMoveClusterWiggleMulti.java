package etomica.virial;

import etomica.atom.AtomSet;
import etomica.atom.IAtomPositioned;
import etomica.atom.IMolecule;
import etomica.box.Box;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.integrator.mcmove.MCMoveMolecule;
import etomica.potential.PotentialMaster;
import etomica.simulation.ISimulation;
import etomica.space.IVector;
import etomica.space3d.IVector3D;
import etomica.space3d.Vector3D;
import etomica.util.Debug;
import etomica.util.IRandom;

/**
 * An MC Move for cluster simulations that "wiggles" a chain molecule.  If the 
 * first or last atom in the chain is chosen, it is moved to a new position 
 * with the same bond length as before, but perturbed by some angle from its
 * original position.  If an Atom in the middle of the chain is chosen, a 
 * crankshaft move is performed that maintains its distances with its 
 * neighbors.  If a middle Atom has a bond angle too close to 180 degrees
 * (such that rotation does nothing) the Atom is not moved at all.
 * In each doTrial, wiggle moves are attempted on all molecules in the box. 
 * 
 * @author Andrew Schultz
 */
public class MCMoveClusterWiggleMulti extends MCMoveMolecule {

    private static final long serialVersionUID = 1L;
    private final MeterClusterWeight weightMeter;
    private final MeterPotentialEnergy energyMeter;

    public MCMoveClusterWiggleMulti(ISimulation sim, PotentialMaster potentialMaster, int nAtoms) {
    	this(potentialMaster,sim.getRandom(), 1.0, nAtoms);
        setBondLength(1.0);
    }
    
    /**
     * Constructor for MCMoveAtomMulti.
     * @param parentIntegrator
     * @param nAtoms number of atoms to move in a trial.  Number of atoms in
     * box should be at least one greater than this value (greater
     * because first atom is never moved)
     */
    public MCMoveClusterWiggleMulti(PotentialMaster potentialMaster, 
            IRandom random, double stepSize, int nAtoms) {
        super(potentialMaster,random,stepSize,Double.POSITIVE_INFINITY,false);
        setStepSizeMax(Math.PI);
        weightMeter = new MeterClusterWeight(potential);
        energyMeter = new MeterPotentialEnergy(potential);
        work1 = new Vector3D();
        work2 = new Vector3D();
        work3 = new Vector3D();
    }

    public void setBox(Box p) {
        super.setBox(p);
        selectedAtoms = new IAtomPositioned[box.getMoleculeList().getAtomCount()];
        translationVectors = new Vector3D[box.getMoleculeList().getAtomCount()];
        for (int i=0; i<translationVectors.length; i++) {
            translationVectors[i] = (IVector3D)p.getSpace().makeVector();
        }
        weightMeter.setBox(p);
        energyMeter.setBox(p);
    }
    
    //note that total energy is calculated
    public boolean doTrial() {
        uOld = energyMeter.getDataAsScalar();
        wOld = weightMeter.getDataAsScalar();

        AtomSet moleculeList = box.getMoleculeList();
        for(int i=0; i<moleculeList.getAtomCount(); i++) {
            AtomSet childList = ((IMolecule)moleculeList.getAtom(i)).getChildList();
            int numChildren = childList.getAtomCount();

            int j = random.nextInt(numChildren);
            selectedAtoms[i] = (IAtomPositioned)childList.getAtom(j);
//            System.out.println(selectedAtoms[i]+" "+j+" before "+selectedAtoms[i].coord.position());
            IVector position = selectedAtoms[i].getPosition();
            translationVectors[i].Ea1Tv1(-1,position);
            if (j == 0 || j == numChildren-1) {
//                System.out.println("end"+j+" move");

                //work1 is the current vector from the bonded atom to atom j
                work1.E(position);
                if (j == 0) {
                    work1.ME(((IAtomPositioned)childList.getAtom(j+1)).getPosition());
                    position.E(((IAtomPositioned)childList.getAtom(j+1)).getPosition());
                }
                else {
                    work1.ME(((IAtomPositioned)childList.getAtom(j-1)).getPosition());
                    position.E(((IAtomPositioned)childList.getAtom(j-1)).getPosition());
                }
                //work2 is a vector perpendicular to work1.  it can be any 
                //perpendicular vector, but that just makes it harder!
                if (work1.x(0)*work1.x(0) < 0.5*bondLength*bondLength) {
                    // if work1 doesn't point in the X direction (mostly) then
                    // find a vector in the plane containing the X axis and work1
                    double a = -work1.x(0)/bondLength;
                    work2.Ea1Tv1(a,work1);
                    work2.setX(0,work2.x(0)+bondLength);
                }
                else {
                    // work1 does point in the X direction (mostly) so
                    // find a vector in the plane containing the Y axis and work1
                    double a = -work1.x(1)/bondLength;
                    work2.Ea1Tv1(a,work1);
                    work2.setX(1,work2.x(1)+bondLength);
                }
                //normalize
                work2.TE(bondLength/Math.sqrt(work2.squared()));
                //work3 is a vector normal to both work1 and work2
                work3.E(work1);
                work3.XE(work2);
                work3.TE(1.0/bondLength);
                
                double phi = (random.nextDouble()-0.5)*Math.PI;
                work2.TE(Math.cos(phi));
                work2.PEa1Tv1(Math.sin(phi),work3);
            }
            else {
//                System.out.println("middle move "+j);
                IVector position0 = ((IAtomPositioned)childList.getAtom(j-1)).getPosition();
                IVector position2 = ((IAtomPositioned)childList.getAtom(j+1)).getPosition();
                work2.Ev1Pv2(position0, position2);
                work2.TE(0.5);
                //work1 is vector between the 0-2 midpoint and 1
                work1.Ev1Mv2(position,work2);
                if (work1.squared() < bondLength*0.05) {
                    translationVectors[i].E(0);
                    continue;
                }
                position.E(work2);
                work2.ME(position0);
                work2.TE(-1);
                work2.XE(work1);
                //work2 is vector between the 0-2 midpoint and 0, normalized to
                //to be the same length as work1
                work2.TE(Math.sqrt(work1.squared()/work2.squared()));
            }
            
            double theta = (random.nextDouble()-0.5)*stepSize;
            position.PEa1Tv1(Math.cos(theta),work1);
            position.PEa1Tv1(Math.sin(theta),work2);

            translationVectors[i].PE(position);
            if (Debug.ON && Debug.DEBUG_NOW) {
                for (int k=0; k<numChildren; k++) {
//                    System.out.println(i+" after "+k+" "+((AtomLeaf)childList.get(k)).coord.position());
                    if (k > 0) {
                        work2.E(((IAtomPositioned)childList.getAtom(k)).getPosition());
                        work2.ME(((IAtomPositioned)childList.getAtom(k-1)).getPosition());
                        double d = Math.sqrt(work2.squared());
//                        System.out.println("distance "+d);
                        if (Math.abs(d - bondLength)/bondLength > 0.000001) {
                            throw new IllegalStateException("wiggle "+i+" "+k+" bond length should be close to "+bondLength+" ("+d+")");
                        }
                    }
                }
            }
        }
        ((BoxCluster)box).trialNotify();
        wNew = weightMeter.getDataAsScalar();
        uNew = energyMeter.getDataAsScalar();
        return true;
    }
    
    public void setBondLength(double b) {
        bondLength = b;
    }
	
    public void rejectNotify() {
        AtomSet moleculeList = box.getMoleculeList();
        for(int i=0; i<selectedAtoms.length; i++) {
            AtomSet childList = ((IMolecule)moleculeList.getAtom(i)).getChildList();
            work1.E(translationVectors[i]);
            work1.TE(1.0/childList.getAtomCount());
            for (int k=0; k<childList.getAtomCount(); k++) {
                ((IAtomPositioned)childList.getAtom(k)).getPosition().PE(work1);
            }
            selectedAtoms[i].getPosition().ME(translationVectors[i]);
        }
        ((BoxCluster)box).rejectNotify();
    }

    public double getB() {
        return -(uNew - uOld);
    }
    
    public double getA() {
        return (wOld==0.0) ? Double.POSITIVE_INFINITY : wNew/wOld;
    }
	
    private IAtomPositioned[] selectedAtoms;
    private double bondLength;
    private final IVector3D work1, work2, work3;
    private IVector3D[] translationVectors;
    private double wOld, wNew;
}
