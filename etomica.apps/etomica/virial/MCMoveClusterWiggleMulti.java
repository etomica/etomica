package etomica.virial;

import etomica.atom.Atom;
import etomica.atom.AtomArrayList;
import etomica.atom.AtomLeaf;
import etomica.atom.AtomTreeNodeGroup;
import etomica.atom.iterator.AtomIteratorAllMolecules;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.integrator.mcmove.MCMoveMolecule;
import etomica.phase.Phase;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.Vector;
import etomica.space3d.Vector3D;
import etomica.util.Debug;

/**
 * An MC Move for cluster simulations that "wiggles" a chain molecule.  If the 
 * first or last atom in the chain is chosen, it is moved to a new position 
 * with the same bond length as before, but perturbed by some angle from its
 * original position.  If an Atom in the middle of the chain is chosen, a 
 * crankshaft move is performed that maintains its distances with its 
 * neightobors.  If a middle Atom has a bond angle too close to 180 degrees
 * (such that rotation does nothing) the Atom is not moved at all.
 * In each doTrail, wiggle moves are attempted on all molecules in the phase. 
 * 
 * @author Andrew Schultz
 */
public class MCMoveClusterWiggleMulti extends MCMoveMolecule {

    private static final long serialVersionUID = 1L;
    private final MeterClusterWeight weightMeter;
    private final MeterPotentialEnergy energyMeter;

    public MCMoveClusterWiggleMulti(Simulation sim, int nAtoms) {
    	this(sim.getPotentialMaster(),sim.getDefaults().atomSize, nAtoms);
        setBondLength(1.0);
    }
    
    /**
     * Constructor for MCMoveAtomMulti.
     * @param parentIntegrator
     * @param nAtoms number of atoms to move in a trial.  Number of atoms in
     * phase should be at least one greater than this value (greater
     * because first atom is never moved)
     */
    public MCMoveClusterWiggleMulti(PotentialMaster potentialMaster, 
            double stepSize, int nAtoms) {
        super(potentialMaster,stepSize,Double.POSITIVE_INFINITY,false);
        this.nAtoms = nAtoms;
        setStepSizeMax(Math.PI);
        selectedMolecules = new Atom[nAtoms];
        selectedAtoms = new AtomLeaf[nAtoms];
        translationVectors = new Vector3D[nAtoms];
        for (int i=0; i<nAtoms; i++) {
            translationVectors[i] = new Vector3D();
        }
        weightMeter = new MeterClusterWeight(potential);
        energyMeter = new MeterPotentialEnergy(potential);
        setName("MCMoveClusterMolecule");
        work1 = new Vector3D();
        work2 = new Vector3D();
        work3 = new Vector3D();
    }

    public void setPhase(Phase p) {
        super.setPhase(p);
        weightMeter.setPhase(p);
        energyMeter.setPhase(p);
    }
    
    //note that total energy is calculated
    public boolean doTrial() {
        if (selectedMolecules[0] == null) selectMolecules();
        uOld = energyMeter.getDataAsScalar();
        wOld = weightMeter.getDataAsScalar();

        for(int i=0; i<selectedMolecules.length; i++) {
            AtomArrayList childList = ((AtomTreeNodeGroup)selectedMolecules[i].getNode()).getChildList();
            int numChildren = childList.size();

            int j = Simulation.random.nextInt(numChildren);
            selectedAtoms[i] = (AtomLeaf)childList.get(j);
//            System.out.println(selectedAtoms[i]+" "+j+" before "+selectedAtoms[i].coord.position());
            Vector position = selectedAtoms[i].getCoord().getPosition();
            translationVectors[i].Ea1Tv1(-1,position);
            if (j == 0 || j == numChildren-1) {
//                System.out.println("end"+j+" move");

                //work1 is the current vector from the bonded atom to atom j
                work1.E(position);
                if (j == 0) {
                    work1.ME(((AtomLeaf)childList.get(j+1)).getCoord().getPosition());
                    position.E(((AtomLeaf)childList.get(j+1)).getCoord().getPosition());
                }
                else {
                    work1.ME(((AtomLeaf)childList.get(j-1)).getCoord().getPosition());
                    position.E(((AtomLeaf)childList.get(j-1)).getCoord().getPosition());
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
                
                double phi = (Simulation.random.nextDouble()-0.5)*Math.PI;
                work2.TE(Math.cos(phi));
                work2.PEa1Tv1(Math.sin(phi),work3);
            }
            else {
//                System.out.println("middle move "+j);
                Vector position0 = ((AtomLeaf)childList.get(j-1)).getCoord().getPosition();
                Vector position2 = ((AtomLeaf)childList.get(j+1)).getCoord().getPosition();
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
            
            double theta = (Simulation.random.nextDouble()-0.5)*stepSize;
            position.PEa1Tv1(Math.cos(theta),work1);
            position.PEa1Tv1(Math.sin(theta),work2);

            translationVectors[i].PE(position);
            if (Debug.ON && Debug.DEBUG_NOW) {
                for (int k=0; k<numChildren; k++) {
//                    System.out.println(i+" after "+k+" "+((AtomLeaf)childList.get(k)).coord.position());
                    if (k > 0) {
                        work2.E(((AtomLeaf)childList.get(k)).getCoord().getPosition());
                        work2.ME(((AtomLeaf)childList.get(k-1)).getCoord().getPosition());
                        double d = Math.sqrt(work2.squared());
//                        System.out.println("distance "+d);
                        if (Math.abs(d - bondLength)/bondLength > 0.000001) {
                            throw new IllegalStateException("wiggle "+i+" "+k+" bond length should be close to "+bondLength+" ("+d+")");
                        }
                    }
                }
            }
        }
        ((PhaseCluster)phase).trialNotify();
        wNew = weightMeter.getDataAsScalar();
        uNew = energyMeter.getDataAsScalar();
//        System.out.println(uOld+" => "+uNew+"   "+wOld+" => "+wNew+" "+stepSize);
        return true;
    }
    
    public void setBondLength(double b) {
        bondLength = b;
    }
	
    protected Atom[] selectMolecules() {
        AtomIteratorAllMolecules iterator = new AtomIteratorAllMolecules(phase);
        if (iterator.size() != nAtoms) throw new IllegalStateException("move should work on number of molecules in phase");
        iterator.reset();
        int i=0;
        while (iterator.hasNext()) {
            selectedMolecules[i++] = iterator.nextAtom();
        }
        return selectedMolecules;
    }
	
    public void rejectNotify() {
        for(int i=0; i<selectedMolecules.length; i++) {
            selectedAtoms[i].getCoord().getPosition().ME(translationVectors[i]);
//            System.out.println(selectedAtoms[i]+" rejected => "+selectedAtoms[i].coord.position());
        }
        ((PhaseCluster)phase).rejectNotify();
    }

    public void acceptNotify() {
        ((PhaseCluster)phase).acceptNotify();
        for(int i=0; i<selectedMolecules.length; i++) {
//            System.out.println(selectedAtoms[i]+" accepted => "+selectedAtoms[i].coord.position());
        }
    }
    
    public double getB() {
        return -(uNew - uOld);
    }
    
    public double getA() {
        return (wOld==0.0) ? Double.POSITIVE_INFINITY : wNew/wOld;
    }
	
    private final int nAtoms;
    private final Atom[] selectedMolecules;
    private final AtomLeaf[] selectedAtoms;
    private double bondLength;
    private final Vector3D work1, work2, work3;
    private final Vector3D[] translationVectors;
    private double wOld, wNew;
}
