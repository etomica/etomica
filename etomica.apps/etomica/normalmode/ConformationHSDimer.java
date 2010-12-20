package etomica.normalmode;

import etomica.api.IAtom;
import etomica.api.IAtomList;
import etomica.api.IVectorMutable;
import etomica.config.IConformation;
import etomica.space.ISpace;

 /**
  *  Conformation for HS Dimer
  *
  * @author Tai Boon Tan
  *
  */
public class ConformationHSDimer implements IConformation, java.io.Serializable{
	
	public ConformationHSDimer(ISpace space){
		this.space = space;
		vector = space.makeVector();
	}

	public void initializePositions(IAtomList atomList) {
			
		IAtom n1 = atomList.getAtom(SpeciesHSDimer.indexAtom1);
		n1.getPosition().E(new double[] {-L/2, 0, 0});
		
		IAtom n2 = atomList.getAtom(SpeciesHSDimer.indexAtom2);
		n2.getPosition().E(new double[] {L/2, 0, 0});
		
	}
		
	protected final ISpace space;
	protected static final double L = 1.0;
	
	protected IVectorMutable vector;
	
	private static final long serialVersionUID = 1L;
}
