package etomica.paracetamol;

import etomica.atom.AtomType;
import etomica.atom.IAtomPositioned;
import etomica.atom.Molecule;

public class AtomParacetamol extends Molecule {

	public AtomParacetamol(AtomType atomType){
		super(atomType);
	}
	
	public IAtomPositioned O1, O2, C1, C2, C3, C4, C5, C6, C7, C8, N1;
	public IAtomPositioned H1, H2, H3, H4, H5, H6, H7, H8, H9;
	
	public final static int indexC1 = 0;
	public final static int indexC2 = 1;
	public final static int indexH1 = 2;
	public final static int indexC3 = 3;
	public final static int indexH2 = 4;
	public final static int indexC4 = 5;
	public final static int indexO1 = 6;	
	public final static int indexC5 = 7;
	public final static int indexH3 = 8;
	public final static int indexC6 = 9;
	public final static int indexH4 =10;
	public final static int indexH5 =11;
	public final static int indexN1  =12;
	public final static int indexH6 =13;
	public final static int indexC7 =14;
	public final static int indexO2 =15;
	public final static int indexC8 =16;	
	public final static int indexH7 =17;
	public final static int indexH8 =18;
	public final static int indexH9 =19;
	
	public final static double [] Echarge = new double [20];

	private static final long serialVersionUID = 1L;
	
}

