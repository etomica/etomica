package etomica.normalmode;

import etomica.api.IAtomList;
import etomica.api.IMolecule;
import etomica.api.IMoleculeList;
import etomica.api.IPotentialMaster;
import etomica.api.IRandom;
import etomica.api.IVectorMutable;
import etomica.space.ISpace;


/**
 * Class which handles orientation of linear molecules during volume change
 * trials.  During compression, the orientation is scaled back toward the
 * nominal orientation according to latticeScale as calculated by the
 * superClass.
 * 
 * @author Andrew Schultz
 */
public class MCMoveVolumeSolidNPTMolecularOriented extends
        MCMoveVolumeSolidNPTMolecular {

    protected final IVectorMutable orientation, orientationLattice, delta, dr0, dr1, com0, com1;
    
    public MCMoveVolumeSolidNPTMolecularOriented(IPotentialMaster potentialMaster, IRandom random,
            ISpace space, double pressure) {
        super(potentialMaster, random, space, pressure, 5);
        orientation = space.makeVector();
        orientationLattice = space.makeVector();
        delta = space.makeVector();
        dr0 = space.makeVector();
        dr1 = space.makeVector();
        com0 = space.makeVector();
        com1 = space.makeVector();
    }
    
    public void doTransform(double vScaleLocal) {
        super.doTransform(vScaleLocal);

        IMoleculeList moleculeList = box.getMoleculeList();
        IMoleculeList moleculeLatticeList = latticeBox.getMoleculeList();
        int nMolecules = moleculeList.getMoleculeCount();

        for (int i=0; i<nMolecules; i++) {
            IMolecule moleculei = moleculeList.getMolecule(i);
            orientation.Ev1Mv2(moleculei.getChildList().getAtom(1).getPosition(), moleculei.getChildList().getAtom(0).getPosition());
            orientation.normalize();
            IMolecule moleculeLattice = moleculeLatticeList.getMolecule(i);
            orientationLattice.Ev1Mv2(moleculeLattice.getChildList().getAtom(1).getPosition(), moleculeLattice.getChildList().getAtom(0).getPosition());
            orientationLattice.normalize();
            
            delta.E(orientation);
            delta.ME(orientationLattice);
            
            double d2 = delta.squared();
            // d2==0 means the molecule is already in its nominal orientation (can't scale)
            if (d2 == 0) continue;

            double x = d2*0.5;
            double costheta = 1. - x; 
            double sintheta = Math.sqrt((2-x)*x);  // sqrt(1 - (1-x)^2) = sqrt(1 - 1 + 2x - x^2) = sqrt(2x - x^2)
            
            // r0 is lattice orientation, r1 is current orientation
            double rScale = Math.exp(vScaleLocal/boxSize.getD());
            double newX = x*latticeScale*latticeScale*rScale*rScale;
            double cosNewTheta = 1 - newX;
            double sinNewTheta = Math.sqrt((2-newX)*newX);
            
            double a = sinNewTheta/sintheta;
            double b = cosNewTheta - sinNewTheta/sintheta*costheta;

            doRotationTransform(moleculei, moleculeLattice, a, b);
        }
    }
    
    protected void doRotationTransform(IMolecule molecule, IMolecule moleculeLattice, double a, double b) {
        com0.E(moleculeCenter.position(moleculeLattice));
        com1.E(moleculeCenter.position(molecule));
        IAtomList childList = molecule.getChildList();
        IAtomList latticeChildList = moleculeLattice.getChildList();
        for (int i=0; i<childList.getAtomCount(); i++) {
            IVectorMutable p1 = childList.getAtom(i).getPosition();
            IVectorMutable p0 = latticeChildList.getAtom(i).getPosition();
            p1.ME(com1);
            p1.TE(a);
            p1.PEa1Tv1(b, p0);
            p1.PEa1Tv1(-b, com0);
            p1.PE(com1);
        }
    }
}
