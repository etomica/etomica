package etomica.potential;

import etomica.api.IBox;
import etomica.api.IMoleculeList;
import etomica.api.IVector;
import etomica.atom.DipoleSource;
import etomica.space.ISpace;

public class P1ExternalField extends PotentialMolecular {

	protected DipoleSource dipoleSource;
	protected IVector externalField;
	
	
	public P1ExternalField( ISpace space) {
		super(1,space);
	}

	public double getRange() {
		return Double.POSITIVE_INFINITY;
	}
	
	public void setExternalField(IVector newExternalField){
		externalField = newExternalField;
	}
	

	public double energy(IMoleculeList molecules) {
		IVector dr = dipoleSource.getDipole(molecules.getMolecule(0));
		double energy = -externalField.dot(dr);
//		System.out.println("E = " + externalField + " dipole = " + dr );
		return energy;
	}

	public void setBox(IBox box) {
		
	}
	
	public void setDipoleSource(DipoleSource newDipoleSource) {
		dipoleSource = newDipoleSource;
	}

}
