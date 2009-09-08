package etomica.data.meter;

import etomica.EtomicaInfo;
import etomica.api.IAtomPositioned;
import etomica.api.IBox;
import etomica.api.IMolecule;
import etomica.api.IMoleculeList;
import etomica.api.IVector;
import etomica.api.IVectorMutable;
import etomica.data.DataSourceScalar;
import etomica.lattice.BravaisLatticeCrystal;
import etomica.space.ISpace;
import etomica.units.Undefined;

/**
 * Meter for calculation of structure factor of a group of atoms based on a
 * particular wave vector.  The meter is initially constructed the smallest 
 * reciprocal lattice vector as the wave vector and loops over all atoms in
 * the box.  GetData returns the square of the magnitutde of the structure factor.
 *
 * @author Michael Sellers
 */

public class MeterStructureFactor extends DataSourceScalar {
	
	protected BravaisLatticeCrystal lattice;
	protected final ISpace space;
    protected IBox box;
    protected double struct;
    protected IVectorMutable waveVec;
    protected IMoleculeList moleculeList;

    
    
	
	/**
	 * Creates meter with default to compute the structure factor
	 * for all atoms in the box.
	 * @param parent
	 */
	public MeterStructureFactor(ISpace space, BravaisLatticeCrystal aLattice, IBox aBox){
		super("Structure factor", Undefined.DIMENSION);
		this.space = space;
        
        this.lattice = aLattice;
        this.box = aBox;
        
        waveVec = space.makeVector();
        waveVec.E(lattice.getPrimitive().makeReciprocal().vectors()[0]);
        waveVec.PE(lattice.getPrimitive().makeReciprocal().vectors()[1]);
        waveVec.PE(lattice.getPrimitive().makeReciprocal().vectors()[2]);
        
        struct = 0;
        moleculeList = box.getMoleculeList();
	}
	
	public void reset() {
		waveVec = space.makeVector();
		waveVec.E(lattice.getPrimitive().makeReciprocal().vectors()[0]);
        waveVec.PE(lattice.getPrimitive().makeReciprocal().vectors()[1]);
        waveVec.PE(lattice.getPrimitive().makeReciprocal().vectors()[2]);
		struct = 0;
		moleculeList = box.getMoleculeList();
	}
	
	public void actionPerformed() {
		
		int numAtoms = moleculeList.getMoleculeCount();
		//System.out.println(numAtoms);
		IVectorMutable workvector = space.makeVector();
		double term1 = 0;
		double term2 = 0;
		double dotprod = 0;
		//System.out.println(waveVec);
		for(int i=0; i<numAtoms; i++){
			workvector.E(((IAtomPositioned)moleculeList.getMolecule(i).getChildList().getAtom(0)).getPosition());
			dotprod = waveVec.dot(workvector);
			term1 += Math.abs(Math.cos(dotprod)); 
			term2 += Math.abs(Math.sin(dotprod));
		}
		struct = (term1*term1 + term2*term2)/(numAtoms*numAtoms);
	}
	
	/**
	 * @param waveVec Sets a custom wave vector array.
	 */
	public void setWaveVec(IVectorMutable waveVec){
		this.waveVec = waveVec;
		
	}
	
	/**
	 * @param atomList Sets the list of atoms for factor calculation.
	 */
	public void setAtoms(IMoleculeList moleculeList){
		this.moleculeList = moleculeList;
	}
	
    public static EtomicaInfo getEtomicaInfo() {
        EtomicaInfo info = new EtomicaInfo("Structure factor");
        return info;
    }

	//this is really the sqare of the magnitude of the structure factor
	public double getDataAsScalar() {
		if(struct==Double.NaN){
			struct=0.0;
		}
		return struct;
	}
	
}
