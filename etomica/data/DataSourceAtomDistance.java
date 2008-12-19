package etomica.data;

import etomica.api.IAtomPositioned;
import etomica.api.IVectorMutable;
import etomica.data.types.DataDouble;
import etomica.data.types.DataDouble.DataInfoDouble;
import etomica.space.ISpace;
import etomica.units.Length;

public class DataSourceAtomDistance extends DataSourceScalar {
	
	public DataSourceAtomDistance(ISpace space) {
		super("interatomic distance", Length.DIMENSION);
		
		this.space = space;
	
		vector = space.makeVector(); // to avoid making the vector each time getData() is called

	}


	public double getDataAsScalar() {

        vector.Ev1Mv2(atom1.getPosition(), atom2.getPosition());
        
		return Math.sqrt(vector.squared());
	}

	
	public void setAtoms(IAtomPositioned atom1, IAtomPositioned atom2) {
		this.atom1 = atom1;
		this.atom2 = atom2;
	}
	

	protected final ISpace space;
	protected final IVectorMutable vector;
	protected IAtomPositioned atom1;
	protected IAtomPositioned atom2;

}
