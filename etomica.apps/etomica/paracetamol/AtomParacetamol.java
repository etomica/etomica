package etomica.paracetamol;

import etomica.atom.AtomGroup;
import etomica.atom.AtomLeaf;
import etomica.atom.AtomType;

public class AtomParacetamol extends AtomGroup{

	public AtomParacetamol(AtomType atomType){
		super(atomType);
	}
	
	public AtomLeaf O1, O2, C1, C2, C3, C4, C5, C6, C7, C8, N;
	public final static int indexC1 = 0;
	public final static int indexC2 = 1;
	public final static int indexC3 = 2;
	public final static int indexC4 = 3;
	public final static int indexO1 = 4;
	public final static int indexC5 = 5;
	public final static int indexC6 = 6;
	public final static int indexN  = 7;	
	public final static int indexC7 = 8;
	public final static int indexC8 = 9;
	public final static int indexO2 =10;
	
	public final static double [] Echarge = new double [11];

	private static final long serialVersionUID = 1L;
	
}

