package etomica.data.meter;

import etomica.api.IAtom;
import etomica.api.IAtomList;
import etomica.api.IBox;
import etomica.api.IMolecule;
import etomica.api.IMoleculeList;
import etomica.api.ISpecies;
import etomica.api.IVectorMutable;
import etomica.data.DataSourceMolecular;
import etomica.data.DataSourceScalar;
import etomica.data.IData;
import etomica.data.IEtomicaDataInfo;
import etomica.space.ISpace;
import etomica.space3d.Vector3D;
import etomica.units.CompoundDimension;
import etomica.units.Debye;
import etomica.units.Dimension;
import etomica.units.Dipole;
import etomica.units.Electron;
import etomica.units.Energy;
import etomica.units.Length;
import etomica.units.Mass;
import etomica.units.Time;
/**
 * meter for (sum dipole)^2
 * used for dielectric constant calculation
 * 
 * @author shu
 *
 */
public class MeterDipoleSumSquaredTIP4PWater extends DataSourceScalar{
	 
    private IBox box;
    private IVectorMutable dipole, dipoleSum;
    
	public MeterDipoleSumSquaredTIP4PWater(ISpace space, IBox box) {
		super("TIP4P water, dipoleSum^2", new CompoundDimension(new Dimension[]{Dipole.DIMENSION},new double[]{2.0}));
		this.box=box;
		dipole = space.makeVector();
		dipoleSum = space.makeVector();
	}
	public double getDataAsScalar() {
		dipoleSum = new Vector3D();
		if (box == null) throw new IllegalStateException("no box");
		IMoleculeList moleculeList = box.getMoleculeList();
		int numMolecule = moleculeList.getMoleculeCount();
		for (int i=0;i<numMolecule; i++){
			IAtomList childList = moleculeList.getMolecule(i).getChildList();	
			IAtom atomH1 = childList.getAtom(0);
			IAtom atomH2 = childList.getAtom(1);
			IAtom atomO = childList.getAtom(2);
			IAtom atomM = childList.getAtom(3);
			double chargeH = Electron.UNIT.toSim(+0.52);
			double chargeM = Electron.UNIT.toSim(-1.04);
			
			dipole.Ea1Tv1(chargeH, atomH1.getPosition());
			dipole.PEa1Tv1(chargeH, atomH2.getPosition());
			dipole.PEa1Tv1(chargeM, atomM.getPosition());// now this is the current dipole
//			System.out.println("in meter class, dipole:"+dipole);			
			dipoleSum.PE(dipole);
//			System.out.println("in meter class, dipoleSum:"+dipoleSum);
//			System.out.println("in meter class, squared of dipoleSum:"+dipoleSum.squared());
		}
        return dipoleSum.squared();
	}
    public IBox getBox() {
    	return box;
    }
    public void setBox(IBox _box) {
    	box = _box;
    }

}