package etomica.potential;

import etomica.api.IAtomList;
import etomica.api.IBoundary;
import etomica.api.IBox;
import etomica.api.IMolecule;
import etomica.api.IMoleculeList;
import etomica.api.IVector;
import etomica.api.IVectorMutable;
import etomica.atom.IAtomOriented;
import etomica.exception.MethodNotImplementedException;
import etomica.space.ISpace;
import etomica.space.Tensor;

/**
 * Hard sphere molecule with a dipole sitting at the center.
 *
 * @author Shu
 * date:June 2014
 * 
 */
public class P2HSDipole extends PotentialMolecular  {

    public P2HSDipole(ISpace space) {
        this(space, 1, 1);
    }

    public P2HSDipole(ISpace space, double dipole) {
        this(space, 1,dipole);
    }
    
    public P2HSDipole(ISpace space, double sigma,double dipole) {
        super(2, space);
        setSigma(sigma);
        dr = space.makeVector();
        setDipole(dipole);
    }


    public void setBox(IBox box) {
        boundary = box.getBoundary();
    }

    public double getRange() {
        return Double.POSITIVE_INFINITY;
    }

    public double energy(IMoleculeList pair){
		IMolecule molecule1 = pair.getMolecule(0);
		IMolecule molecule2 = pair.getMolecule(1);
		IAtomList atomList1 = molecule1.getChildList();
		IAtomList atomList2 = molecule2.getChildList();
		
        IAtomOriented atom1 = (IAtomOriented)atomList1.getAtom(0);
        IAtomOriented atom2 = (IAtomOriented)atomList2.getAtom(0);
        // dr is a r1-r2
        dr.Ev1Mv2(atom1.getPosition(), atom2.getPosition());
        double r2 = dr.squared();
        if(r2 < sigma2) { // hard core
            return Double.POSITIVE_INFINITY;
        }
        // normalize dr, the vector between the molecules
        dr.normalize();
        
        // v1 (unit vector) is the orientation of molecule 1: dipole1 direction
        IVector v1 = atom1.getOrientation().getDirection();
        // v2 (unit vector) is the orientation of molecule 2: dipole1 direction
        IVector v2 = atom2.getOrientation().getDirection();

        // cos(dipole 1 and dipole 2)=cos(v1 and v2)
        double cos_D1_D2 = v1.dot(v2);
        //cos(dipole 1 and r12)
        double cos_D1_r = v1.dot(dr);
        //cos(r12 and dipole 2)
        double cos_r_D2=dr.dot(v2);
        double r12Magnitude = Math.sqrt(r2);
        double ener = dipole * dipole *  (cos_D1_D2 - 3.0  * cos_D1_r* cos_r_D2);
        ener = ener/r12Magnitude /r2 ;
        return ener;
    }

    public double getSigma() {return sigma;}

    public final void setSigma(double s) {
        sigma = s;
        sigma2 = s * s;
    }

    public void setDipole(double moment){
        dipole=moment;
    }
    private static final long serialVersionUID = 1L;
    private double sigma , sigma2;
    private IBoundary boundary;
    private final IVectorMutable dr;
    private double dipole;

}
