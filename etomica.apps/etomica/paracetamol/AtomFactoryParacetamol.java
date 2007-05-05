package etomica.paracetamol;

import etomica.atom.AtomFactory;
import etomica.atom.AtomFactoryMono;
import etomica.atom.AtomFactoryMonoDynamic;
import etomica.atom.AtomPositionGeometricCenter;
import etomica.atom.AtomTypeGroup;
import etomica.atom.AtomTypeSphere;
import etomica.atom.IAtom;
import etomica.atom.IAtomPositioned;
import etomica.chem.elements.Carbon;
import etomica.chem.elements.Nitrogen;
import etomica.chem.elements.Oxygen;
import etomica.simulation.Simulation;
import etomica.units.ElectronVolt;

public class AtomFactoryParacetamol extends AtomFactory{

	public AtomFactoryParacetamol(Simulation sim) {
		super (new AtomTypeGroup(new AtomPositionGeometricCenter(sim.getSpace())));
		
		AtomParacetamol.Echarge[0] = ElectronVolt.UNIT.toSim( 0.10);
		AtomParacetamol.Echarge[1] = ElectronVolt.UNIT.toSim(-0.10);
		AtomParacetamol.Echarge[2] = ElectronVolt.UNIT.toSim(-0.10);
		AtomParacetamol.Echarge[3] = ElectronVolt.UNIT.toSim( 0.03);
		AtomParacetamol.Echarge[4] = ElectronVolt.UNIT.toSim(-0.38);
		AtomParacetamol.Echarge[5] = ElectronVolt.UNIT.toSim(-0.10);
		AtomParacetamol.Echarge[6] = ElectronVolt.UNIT.toSim(-0.10);
		AtomParacetamol.Echarge[7] = ElectronVolt.UNIT.toSim(-0.39);
		AtomParacetamol.Echarge[8] = ElectronVolt.UNIT.toSim( 0.38);
		AtomParacetamol.Echarge[9] = ElectronVolt.UNIT.toSim( 0.30);
		AtomParacetamol.Echarge[10]= ElectronVolt.UNIT.toSim(-0.38);
		
		oType = new AtomTypeSphere(Oxygen.INSTANCE, 2*1.7); //atomic Instance Class, atomic diameter
		oType.setParentType((AtomTypeGroup)atomType);		//L. Pauling, The Nature of the Chemical Bond, Cornell University Press, USA, 1945.
		cType = new AtomTypeSphere(Carbon.INSTANCE, 2*1.55); //1.518, 1.430, 1.338
		cType.setParentType((AtomTypeGroup)atomType);
		nType = new AtomTypeSphere(Nitrogen.INSTANCE, 2*1.52);
		nType.setParentType((AtomTypeGroup)atomType);

		//CoordinateFactory leafCoordFactory = new CoordinateFactorySphere(sim);
		oFactory = new AtomFactoryMonoDynamic(sim.getSpace(), oType);
		cFactory = new AtomFactoryMonoDynamic(sim.getSpace(), cType);
		nFactory = new AtomFactoryMonoDynamic(sim.getSpace(), nType);
		
		conformation = new ConformationParacetamolOrthorhombic(sim.getSpace());
	}
	
	
	public IAtom makeAtom() {

		AtomParacetamol moleculeParacetamol = new AtomParacetamol(atomType);
		
		moleculeParacetamol.O1 = (IAtomPositioned)oFactory.makeAtom();
		moleculeParacetamol.O2 = (IAtomPositioned)oFactory.makeAtom();
		moleculeParacetamol.C1 = (IAtomPositioned)cFactory.makeAtom();
		moleculeParacetamol.C2 = (IAtomPositioned)cFactory.makeAtom();
		moleculeParacetamol.C3 = (IAtomPositioned)cFactory.makeAtom();
		moleculeParacetamol.C4 = (IAtomPositioned)cFactory.makeAtom();
		moleculeParacetamol.C5 = (IAtomPositioned)cFactory.makeAtom();
		moleculeParacetamol.C6 = (IAtomPositioned)cFactory.makeAtom();
		moleculeParacetamol.C7 = (IAtomPositioned)cFactory.makeAtom();
		moleculeParacetamol.C8 = (IAtomPositioned)cFactory.makeAtom();
		moleculeParacetamol.N  = (IAtomPositioned)nFactory.makeAtom();
		//set parent
		moleculeParacetamol.addChildAtom(moleculeParacetamol.C1);
		moleculeParacetamol.addChildAtom(moleculeParacetamol.C2);
		moleculeParacetamol.addChildAtom(moleculeParacetamol.C3);
		moleculeParacetamol.addChildAtom(moleculeParacetamol.C4);
		moleculeParacetamol.addChildAtom(moleculeParacetamol.O1);
		moleculeParacetamol.addChildAtom(moleculeParacetamol.C5);
		moleculeParacetamol.addChildAtom(moleculeParacetamol.C6);
		moleculeParacetamol.addChildAtom(moleculeParacetamol.N );
		moleculeParacetamol.addChildAtom(moleculeParacetamol.C7);
		moleculeParacetamol.addChildAtom(moleculeParacetamol.C8);
		moleculeParacetamol.addChildAtom(moleculeParacetamol.O2);
		
		conformation.initializePositions(moleculeParacetamol.getChildList());
		return moleculeParacetamol;	
	}
	
	
	public int getNumTreeAtoms(){
		return 12;
	}
	
	public int getNumChildAtoms(){
		return 11;
	}
	
    public int getNumLeafAtoms() {
        return 11;
    }
    
	public final AtomFactoryMono oFactory, cFactory, nFactory;
	
	public final AtomTypeSphere cType, oType, nType;
	private static final long serialVersionUID = 1L;
}
