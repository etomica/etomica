package etomica.paracetamol;

import java.lang.reflect.Constructor;

import etomica.atom.AtomLeaf;
import etomica.atom.AtomLeafDynamic;
import etomica.atom.AtomPositionGeometricCenter;
import etomica.atom.AtomTypeLeaf;
import etomica.atom.AtomTypeMolecule;
import etomica.atom.AtomTypeSphere;
import etomica.atom.IAtomPositioned;
import etomica.atom.IMolecule;
import etomica.chem.elements.Carbon;
import etomica.chem.elements.Hydrogen;
import etomica.chem.elements.Nitrogen;
import etomica.chem.elements.Oxygen;
import etomica.simulation.ISimulation;
import etomica.space.Space;
import etomica.species.Species;
import etomica.species.SpeciesSignature;
import etomica.units.ElectronVolt;

public class SpeciesParacetamol extends Species {

	public SpeciesParacetamol(ISimulation sim) {
		super(new AtomTypeMolecule(new AtomPositionGeometricCenter(sim.getSpace())));
		space = sim.getSpace();
		isDynamic = sim.isDynamic();
        AtomParacetamol.Echarge[ 0] = ElectronVolt.UNIT.toSim( 0.382743);
        AtomParacetamol.Echarge[ 1] = ElectronVolt.UNIT.toSim(-0.227865);
        AtomParacetamol.Echarge[ 2] = ElectronVolt.UNIT.toSim( 0.168364);
        AtomParacetamol.Echarge[ 3] = ElectronVolt.UNIT.toSim(-0.173831);
        AtomParacetamol.Echarge[ 4] = ElectronVolt.UNIT.toSim( 0.140597);
        AtomParacetamol.Echarge[ 5] = ElectronVolt.UNIT.toSim( 0.343891);
        AtomParacetamol.Echarge[ 6] = ElectronVolt.UNIT.toSim(-0.595151);
        AtomParacetamol.Echarge[ 7] = ElectronVolt.UNIT.toSim(-0.241393);
        AtomParacetamol.Echarge[ 8] = ElectronVolt.UNIT.toSim( 0.121492);
        AtomParacetamol.Echarge[ 9] = ElectronVolt.UNIT.toSim(-0.227674);
        AtomParacetamol.Echarge[10] = ElectronVolt.UNIT.toSim( 0.130830);
        AtomParacetamol.Echarge[11] = ElectronVolt.UNIT.toSim( 0.422442);
        AtomParacetamol.Echarge[12] = ElectronVolt.UNIT.toSim(-0.748988);
        AtomParacetamol.Echarge[13] = ElectronVolt.UNIT.toSim( 0.352011);
        AtomParacetamol.Echarge[14] = ElectronVolt.UNIT.toSim( 0.808224);
        AtomParacetamol.Echarge[15] = ElectronVolt.UNIT.toSim(-0.554225);
        AtomParacetamol.Echarge[16] = ElectronVolt.UNIT.toSim(-0.417884);
        AtomParacetamol.Echarge[17] = ElectronVolt.UNIT.toSim( 0.093811);
        AtomParacetamol.Echarge[18] = ElectronVolt.UNIT.toSim( 0.110836);
        AtomParacetamol.Echarge[19] = ElectronVolt.UNIT.toSim( 0.111772);
        
        //atomic Instance Class, atomic diameter
        //L. Pauling, The Nature of the Chemical Bond, Cornell University Press, USA, 1945.
        //1.7, 1.55, 1.52, 1.2
        oType = new AtomTypeSphere(Oxygen.INSTANCE, 2*1.7);
        atomType.addChildType(oType);
        cType = new AtomTypeSphere(Carbon.INSTANCE, 2*1.55);
        atomType.addChildType(cType);
        nType = new AtomTypeSphere(Nitrogen.INSTANCE, 2*1.52);
        atomType.addChildType(nType);
        hpType = new AtomTypeSphere(HydrogenP.INSTANCE, 2*1.20);
        atomType.addChildType(hpType);
        hyType = new AtomTypeSphere(Hydrogen.INSTANCE, 2*1.20);
        atomType.addChildType(hyType);

        //CoordinateFactory leafCoordFactory = new CoordinateFactorySphere(sim);
        atomType.setConformation(new ConformationParacetamolOrthorhombic(sim.getSpace()));
    }
    
    
    public IMolecule makeMolecule() {

        AtomParacetamol moleculeParacetamol = new AtomParacetamol(atomType);
        
        moleculeParacetamol.O1 = makeLeafAtom(oType);
        moleculeParacetamol.O2 = makeLeafAtom(oType);
        moleculeParacetamol.C1 = makeLeafAtom(cType);
        moleculeParacetamol.C2 = makeLeafAtom(cType);
        moleculeParacetamol.C3 = makeLeafAtom(cType);
        moleculeParacetamol.C4 = makeLeafAtom(cType);
        moleculeParacetamol.C5 = makeLeafAtom(cType);
        moleculeParacetamol.C6 = makeLeafAtom(cType);
        moleculeParacetamol.C7 = makeLeafAtom(cType);
        moleculeParacetamol.C8 = makeLeafAtom(cType);
        moleculeParacetamol.N1 = makeLeafAtom(nType);
        moleculeParacetamol.H1 = makeLeafAtom(hyType);
        moleculeParacetamol.H2 = makeLeafAtom(hyType);
        moleculeParacetamol.H3 = makeLeafAtom(hyType);
        moleculeParacetamol.H4 = makeLeafAtom(hyType);
        moleculeParacetamol.H5 = makeLeafAtom(hpType);
        moleculeParacetamol.H6 = makeLeafAtom(hpType);
        moleculeParacetamol.H7 = makeLeafAtom(hyType);
        moleculeParacetamol.H8 = makeLeafAtom(hyType);
        moleculeParacetamol.H9 = makeLeafAtom(hyType);
        
        //set parent
        moleculeParacetamol.addChildAtom(moleculeParacetamol.C1);
        moleculeParacetamol.addChildAtom(moleculeParacetamol.C2);
        moleculeParacetamol.addChildAtom(moleculeParacetamol.H1);
        moleculeParacetamol.addChildAtom(moleculeParacetamol.C3);
        moleculeParacetamol.addChildAtom(moleculeParacetamol.H2);
        moleculeParacetamol.addChildAtom(moleculeParacetamol.C4);
        moleculeParacetamol.addChildAtom(moleculeParacetamol.O1);
        moleculeParacetamol.addChildAtom(moleculeParacetamol.C5);
        moleculeParacetamol.addChildAtom(moleculeParacetamol.H3);
        moleculeParacetamol.addChildAtom(moleculeParacetamol.C6);
        moleculeParacetamol.addChildAtom(moleculeParacetamol.H4);
        moleculeParacetamol.addChildAtom(moleculeParacetamol.H5);
        moleculeParacetamol.addChildAtom(moleculeParacetamol.N1);
        moleculeParacetamol.addChildAtom(moleculeParacetamol.H6);
        moleculeParacetamol.addChildAtom(moleculeParacetamol.C7);
        moleculeParacetamol.addChildAtom(moleculeParacetamol.O2);
        moleculeParacetamol.addChildAtom(moleculeParacetamol.C8);
        moleculeParacetamol.addChildAtom(moleculeParacetamol.H7);
        moleculeParacetamol.addChildAtom(moleculeParacetamol.H8);
        moleculeParacetamol.addChildAtom(moleculeParacetamol.H9);
        
        atomType.getConformation().initializePositions(moleculeParacetamol.getChildList());
        return moleculeParacetamol; 
    }
    
    protected IAtomPositioned makeLeafAtom(AtomTypeLeaf leafType) {
        return isDynamic ? new AtomLeafDynamic(space, leafType)
                         : new AtomLeaf(space, leafType);
    }

    
    public int getNumLeafAtoms() {
        return 20;
    }
    
	public AtomTypeSphere getCType() {
        return cType;
    }


    public AtomTypeSphere getOType() {
        return oType;
    }


    public AtomTypeSphere getNType() {
        return nType;
    }


    public AtomTypeSphere getHpType() {
        return hpType;
    }


    public AtomTypeSphere getHyType() {
        return hyType;
    }


    public SpeciesSignature getSpeciesSignature(){
		Constructor constructor = null;
		try {
			constructor = this.getClass().getConstructor(new Class[] {ISimulation.class});
		}
		catch(NoSuchMethodException e){
			System.err.println("Have NO CONSTRUCTOR");
		}
		return new SpeciesSignature(constructor,new Object[]{});
	}
	
	private static final long serialVersionUID = 1L;
	protected final Space space;
	protected final boolean isDynamic;
    protected final AtomTypeSphere cType, oType, nType, hpType, hyType;
}