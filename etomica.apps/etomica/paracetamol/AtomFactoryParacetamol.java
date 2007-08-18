package etomica.paracetamol;

/*
 *  This class contains all the atoms for paracetamol molecule
 *  Each atoms in the paracetamol molecules is assigned with a specific e-charge
 *  
 *  @author Tai Tan
 */

import etomica.atom.AtomFactory;
import etomica.atom.AtomFactoryMono;
import etomica.atom.AtomFactoryMonoDynamic;
import etomica.atom.AtomLeaf;
import etomica.atom.AtomPositionGeometricCenter;
import etomica.atom.AtomTypeGroup;
import etomica.atom.AtomTypeSphere;
import etomica.atom.IAtom;
import etomica.chem.elements.Carbon;
import etomica.chem.elements.Hydrogen;
import etomica.chem.elements.Nitrogen;
import etomica.chem.elements.Oxygen;
import etomica.simulation.ISimulation;
import etomica.units.ElectronVolt;

public class AtomFactoryParacetamol extends AtomFactory{

	public AtomFactoryParacetamol(ISimulation sim) {
		super (new AtomTypeGroup(new AtomPositionGeometricCenter(sim.getSpace())));
		
		AtomParacetamol.Echarge[ 0] = ElectronVolt.UNIT.toSim( 0.382743);
		AtomParacetamol.Echarge[ 1] = ElectronVolt.UNIT.toSim(-0.227865);
		AtomParacetamol.Echarge[ 2] = ElectronVolt.UNIT.toSim( 0.168364);
		AtomParacetamol.Echarge[ 3] = ElectronVolt.UNIT.toSim(-0.173831);
		AtomParacetamol.Echarge[ 4] = ElectronVolt.UNIT.toSim( 0.140597);
		AtomParacetamol.Echarge[ 5] = ElectronVolt.UNIT.toSim( 0.343891);
		AtomParacetamol.Echarge[ 6] = ElectronVolt.UNIT.toSim(-0.595151);
		AtomParacetamol.Echarge[ 7] = ElectronVolt.UNIT.toSim( 0.422442);
		AtomParacetamol.Echarge[ 8] = ElectronVolt.UNIT.toSim(-0.241393);
		AtomParacetamol.Echarge[ 9] = ElectronVolt.UNIT.toSim( 0.121492);
		AtomParacetamol.Echarge[10] = ElectronVolt.UNIT.toSim(-0.227674);
		AtomParacetamol.Echarge[11] = ElectronVolt.UNIT.toSim( 0.130830);
		AtomParacetamol.Echarge[12] = ElectronVolt.UNIT.toSim(-0.748988);
		AtomParacetamol.Echarge[13] = ElectronVolt.UNIT.toSim( 0.352011);
		AtomParacetamol.Echarge[14] = ElectronVolt.UNIT.toSim( 0.808224);
		AtomParacetamol.Echarge[15] = ElectronVolt.UNIT.toSim(-0.554225);
		AtomParacetamol.Echarge[16] = ElectronVolt.UNIT.toSim(-0.417884);
		AtomParacetamol.Echarge[17] = ElectronVolt.UNIT.toSim( 0.093811);
		AtomParacetamol.Echarge[18] = ElectronVolt.UNIT.toSim( 0.110836);
		AtomParacetamol.Echarge[19] = ElectronVolt.UNIT.toSim( 0.111772);
		
		oType = new AtomTypeSphere(Oxygen.INSTANCE, 2*1.7); //atomic Instance Class, atomic diameter
		oType.setParentType((AtomTypeGroup)atomType);		//L. Pauling, The Nature of the Chemical Bond, Cornell University Press, USA, 1945.
		cType = new AtomTypeSphere(Carbon.INSTANCE, 2*1.55); //1.7, 1.55, 1.52, 1.2
		cType.setParentType((AtomTypeGroup)atomType);
		nType = new AtomTypeSphere(Nitrogen.INSTANCE, 2*1.52);
		nType.setParentType((AtomTypeGroup)atomType);
		hpType = new AtomTypeSphere(Hydrogen.INSTANCE, 2*1.20);
		hpType.setParentType((AtomTypeGroup)atomType);
		hyType = new AtomTypeSphere(Hydrogen.INSTANCE, 2*1.20);
		hyType.setParentType((AtomTypeGroup)atomType);

		//CoordinateFactory leafCoordFactory = new CoordinateFactorySphere(sim);
		oFactory = new AtomFactoryMonoDynamic(sim.getSpace(), oType);
		cFactory = new AtomFactoryMonoDynamic(sim.getSpace(), cType);
		nFactory = new AtomFactoryMonoDynamic(sim.getSpace(), nType);
		hpFactory = new AtomFactoryMonoDynamic(sim.getSpace(), hpType);
		hyFactory = new AtomFactoryMonoDynamic(sim.getSpace(), hyType);
		
		((AtomTypeGroup)atomType).setConformation(new ConformationParacetamolOrthorhombic(sim.getSpace()));
	}
	
	
	public IAtom makeAtom() {

		AtomParacetamol moleculeParacetamol = new AtomParacetamol(atomType);
		
		moleculeParacetamol.O1 = (AtomLeaf)oFactory.makeAtom();
		moleculeParacetamol.O2 = (AtomLeaf)oFactory.makeAtom();
		moleculeParacetamol.C1 = (AtomLeaf)cFactory.makeAtom();
		moleculeParacetamol.C2 = (AtomLeaf)cFactory.makeAtom();
		moleculeParacetamol.C3 = (AtomLeaf)cFactory.makeAtom();
		moleculeParacetamol.C4 = (AtomLeaf)cFactory.makeAtom();
		moleculeParacetamol.C5 = (AtomLeaf)cFactory.makeAtom();
		moleculeParacetamol.C6 = (AtomLeaf)cFactory.makeAtom();
		moleculeParacetamol.C7 = (AtomLeaf)cFactory.makeAtom();
		moleculeParacetamol.C8 = (AtomLeaf)cFactory.makeAtom();
		moleculeParacetamol.N  = (AtomLeaf)nFactory.makeAtom();
		moleculeParacetamol.H1 = (AtomLeaf)hyFactory.makeAtom();
		moleculeParacetamol.H2 = (AtomLeaf)hyFactory.makeAtom();
		moleculeParacetamol.H3 = (AtomLeaf)hyFactory.makeAtom();
		moleculeParacetamol.H4 = (AtomLeaf)hyFactory.makeAtom();
		moleculeParacetamol.H5 = (AtomLeaf)hpFactory.makeAtom();
		moleculeParacetamol.H6 = (AtomLeaf)hpFactory.makeAtom();
		moleculeParacetamol.H7 = (AtomLeaf)hyFactory.makeAtom();
		moleculeParacetamol.H8 = (AtomLeaf)hyFactory.makeAtom();
		moleculeParacetamol.H9 = (AtomLeaf)hyFactory.makeAtom();
		
		//set parent
		moleculeParacetamol.addChildAtom(moleculeParacetamol.C1);
		moleculeParacetamol.addChildAtom(moleculeParacetamol.C2);
		moleculeParacetamol.addChildAtom(moleculeParacetamol.H1);
		moleculeParacetamol.addChildAtom(moleculeParacetamol.C3);
		moleculeParacetamol.addChildAtom(moleculeParacetamol.H2);
		moleculeParacetamol.addChildAtom(moleculeParacetamol.C4);
		moleculeParacetamol.addChildAtom(moleculeParacetamol.O1);
		moleculeParacetamol.addChildAtom(moleculeParacetamol.H5);
		moleculeParacetamol.addChildAtom(moleculeParacetamol.C5);
		moleculeParacetamol.addChildAtom(moleculeParacetamol.H3);
		moleculeParacetamol.addChildAtom(moleculeParacetamol.C6);
		moleculeParacetamol.addChildAtom(moleculeParacetamol.H4);
		moleculeParacetamol.addChildAtom(moleculeParacetamol.N );
		moleculeParacetamol.addChildAtom(moleculeParacetamol.H6);
		moleculeParacetamol.addChildAtom(moleculeParacetamol.C7);
		moleculeParacetamol.addChildAtom(moleculeParacetamol.O2);
		moleculeParacetamol.addChildAtom(moleculeParacetamol.C8);
		moleculeParacetamol.addChildAtom(moleculeParacetamol.H7);
		moleculeParacetamol.addChildAtom(moleculeParacetamol.H8);
		moleculeParacetamol.addChildAtom(moleculeParacetamol.H9);
		
        ((AtomTypeGroup)atomType).getConformation().initializePositions(moleculeParacetamol.getChildList());
		return moleculeParacetamol;	
	}
	
	
	public int getNumTreeAtoms(){
		return 21;
	}
	
	public int getNumChildAtoms(){
		return 20;
	}
	
    public int getNumLeafAtoms() {
        return 20;
    }
    
	public final AtomFactoryMono oFactory, cFactory, nFactory, hpFactory, hyFactory;
	
	public final AtomTypeSphere cType, oType, nType, hpType, hyType;
	private static final long serialVersionUID = 1L;
}
