package etomica.data.meter;

import etomica.EtomicaInfo;
import etomica.api.IBox;
import etomica.api.IVector;
import etomica.atom.AtomSet;
import etomica.atom.AtomTypeLeaf;
import etomica.atom.AtomTypeMoleculeOriented;
import etomica.atom.IAtomKinetic;
import etomica.atom.IAtomOrientedKinetic;
import etomica.atom.IMolecule;
import etomica.atom.MoleculeOrientedDynamic;
import etomica.data.DataSourceScalar;
import etomica.simulation.ISimulation;
import etomica.space.Space;
import etomica.space3d.IOrientationFull3D;
import etomica.space3d.RotationTensor3D;
import etomica.species.ISpecies;
import etomica.units.Energy;

/**
 * Meter for the total kinetic energy in a box
 * Computes total KE by summing values of KE returned by every atom in the box.
 * Kinetic energy for rigid molecules is calculated based on the its
 * translational and rotational energy.
 * 
 * @author Andrew Schultz
 */
public class MeterKineticEnergyRigid extends DataSourceScalar {

    public MeterKineticEnergyRigid(Space space, ISimulation sim) {
        this(space, sim, null);
    }
    
    public MeterKineticEnergyRigid(Space space, ISimulation sim, IBox box) {
        super("Kinetic Energy",Energy.DIMENSION);
        angularVelocity = space.makeVector();
        rotationTensor = (RotationTensor3D)space.makeRotationTensor();
        this.box = box;
        this.sim = sim;
    }
    
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Total kinetic energy of molecular motion in a box");
        return info;
    }

	/**
	 * Returns the total kinetic energy summed over all atoms produced by
	 * the iterator when applied to the given box.  Does not include contributions
     * from atoms having infinite mass (it assumes they are stationary).
	 */
    public double getDataAsScalar() {
        if (box == null) throw new IllegalStateException("must call setBox before using meter");
        double ke = 0.0;
        ISpecies[] species = sim.getSpeciesManager().getSpecies();
        for (int i=0; i<species.length; i++) {
            AtomSet moleculeList = box.getMoleculeList(species[i]);
            if (moleculeList.getAtomCount() == 0) {
                continue;
            }
            IMolecule molecule0 = (IMolecule)moleculeList.getAtom(0);
            if (molecule0 instanceof MoleculeOrientedDynamic) {
                for (int j=0; j<moleculeList.getAtomCount(); j++) {
                    IAtomOrientedKinetic moleculeOrientedKinetic = (IAtomOrientedKinetic)moleculeList.getAtom(j);
                    double mass = ((AtomTypeMoleculeOriented)molecule0.getType()).getMass();
                    if (Double.isInfinite(mass)) {
                        continue;
                    }
                    ke += 0.5*mass*((IAtomKinetic)molecule0).getVelocity().squared();
                    IVector moment = ((AtomTypeMoleculeOriented)molecule0.getType()).getMomentOfInertia();
        
                    angularVelocity.E(moleculeOrientedKinetic.getAngularVelocity());
                    rotationTensor.setOrientation((IOrientationFull3D)moleculeOrientedKinetic.getOrientation());
                    rotationTensor.transform(angularVelocity);
                    angularVelocity.DE(moment);
        
                    angularVelocity.TE(angularVelocity);
                    angularVelocity.TE(moment);
                    ke += angularVelocity.x(0) + angularVelocity.x(1)+ angularVelocity.x(2);
                }
            }
            else {
                for (int j=0; j<moleculeList.getAtomCount(); j++) {
                    IMolecule molecule = (IMolecule)moleculeList.getAtom(j);
                    AtomSet children = molecule.getChildList();
                    for (int iLeaf=0; iLeaf<children.getAtomCount(); iLeaf++) {
                        IAtomKinetic a = (IAtomKinetic)children.getAtom(iLeaf);
                        double mass = ((AtomTypeLeaf)a.getType()).getMass();
                        if(mass == Double.POSITIVE_INFINITY) continue;
        //                    System.out.println("force: "+((MyAgent)a.ia).force.toString());
                        IVector velocity = a.getVelocity();
                        ke += velocity.squared() * mass;
                    }
                }
            }
        }
        ke *= 0.5;
        return ke;
    }

    /**
     * Sets the box to the given box.
     */
    public void setBox(IBox newBox) {
        box = newBox;
    }
    
    /**
     * @return Returns the box.
     */
    public IBox getBox() {
        return box;
    }

    private static final long serialVersionUID = 1L;
    protected final ISimulation sim;
    protected IBox box;
    protected final IVector angularVelocity;
    protected final RotationTensor3D rotationTensor;
 }