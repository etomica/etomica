package etomica.data.meter;

import etomica.atom.IAtomList;
import etomica.atom.IAtomOriented;
import etomica.box.Box;
import etomica.data.DataSourceScalar;
import etomica.molecule.IMoleculeList;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Vector3D;
import etomica.units.dimensions.CompoundDimension;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Dipole;
/**
 * meter for (sum dipole)^2, used for dielectric constant calculation
 * 
 * @author shu
 */
public class MeterDipoleSumSquared1site extends DataSourceScalar {
	 
    private Box box;
    private Vector dipoleSum;
    private double dipoleMagnitude;
    
	public MeterDipoleSumSquared1site(Space space, Box box, double dipoleMagnitude) {
		super("dipoleSum^2", new CompoundDimension(new Dimension[]{Dipole.DIMENSION},new double[]{2.0}));
		this.box=box;
		this.dipoleMagnitude = dipoleMagnitude;
		dipoleSum = space.makeVector();
	}

	public double getDataAsScalar() {
		dipoleSum = new Vector3D();
		if (box == null) throw new IllegalStateException("no box");
		IMoleculeList moleculeList = box.getMoleculeList();
		int numMolecule = moleculeList.getMoleculeCount();
		for (int i=0;i<numMolecule; i++){
			IAtomList atomList = moleculeList.getMolecule(i).getChildList();
			IAtomOriented atom = (IAtomOriented) atomList.getAtom(0);
	        Vector v = atom.getOrientation().getDirection();
			dipoleSum.PE(v);
        }
        double squared = dipoleMagnitude*dipoleMagnitude*dipoleSum.squared();
        return squared;
	}

    public Box getBox() {
    	return box;
    }
    public void setBox(Box _box) {
    	box = _box;
    }

}
