/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data.meter;

import etomica.atom.AtomLeafAgentManager;
import etomica.box.Box;
import etomica.data.*;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.integrator.IntegratorRigidIterative.MoleculeAgent;
import etomica.molecule.DipoleSource;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.molecule.MoleculeAgentManager;
import etomica.molecule.MoleculeAgentManager.MoleculeAgentSource;
import etomica.potential.*;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.units.dimensions.Null;

/**
 * meter for AEE use mapping average
 * 
 * @author Weisong
 */
public class MeterDipoleSumSquaredMappedAverage implements IDataSource, MoleculeAgentSource {

	protected final DataDoubleArray data;
	protected final DataInfoDoubleArray dataInfo;
	protected final DataTag tag;
//	private IBoundary boundary;
	protected PotentialCalculationEnergySum energySum;
	protected PotentialCalculationTorqueSum torqueSum;
	protected PotentialCalculationPhiSum secondDerivativeSum;
	protected final Space space;
	private Box box;
	private Vector vectorSum;
//	private IVectorMutable r;
//	private IVectorMutable [] a;
	private double dipoleMagnitude;
//	private double truncation;
	private double temperature;
    protected final PotentialMaster potentialMaster;
    private final IteratorDirective allAtoms;
    protected Vector dr;
    protected Vector work;
    protected MoleculeAgentManager moleculeAgentManager;
    protected DipoleSource dipoleSource;
    protected AtomLeafAgentManager<Vector> atomAgentManager;
    protected PotentialCalculationForceSum pcForce;

	public MeterDipoleSumSquaredMappedAverage(final Space space, Box box, Simulation sim, double dipoleMagnitude, double temperature, PotentialMaster potentialMaster) {
		data = new DataDoubleArray(2);
		dataInfo = new DataInfoDoubleArray("stuff", Null.DIMENSION, new int[]{2});
		tag = new DataTag();
		dataInfo.addTag(tag);
		this.box = box;
		this.dipoleMagnitude = dipoleMagnitude;
		this.temperature = temperature;
		this.space = space;
		this.potentialMaster = potentialMaster;
		vectorSum = space.makeVector();
//		r = space.makeVector();
		vectorSum.setX(2, 1);
		torqueSum = new PotentialCalculationTorqueSum();
		energySum = new PotentialCalculationEnergySum();
		secondDerivativeSum = new  PotentialCalculationPhiSum(space);
		moleculeAgentManager  = new MoleculeAgentManager(sim,box,this);
		torqueSum.setMoleculeAgentManager(moleculeAgentManager);
		allAtoms = new IteratorDirective();
		dr = space.makeVector();
		work = space.makeVector();
		
		pcForce = new PotentialCalculationForceSum();
		atomAgentManager = new AtomLeafAgentManager<Vector>(a -> space.makeVector() , box);
        pcForce.setAgentManager(atomAgentManager);
		
	}

	public IData getData() {		
//		IBoundary boundary = box.getBoundary();// TODO
		double[] x = data.getData();
		double bt = 1/(temperature);//beta
		
		double mu = dipoleMagnitude;//miu
		double mu2 = mu*mu;
		double bt2 = bt*bt;
		double bt3 = bt*bt*bt;
		if (box == null) throw new IllegalStateException("no box");
		IMoleculeList moleculeList = box.getMoleculeList();
		
		
		int nM = moleculeList.size();
		
		//TODO
//		IAtomList atomList0 = moleculeList.getMolecule(0).getChildList();
//		IAtomOriented atom0 = (IAtomOriented) atomList0.getAtom(0);
//		IVectorMutable  v0 =  (IVectorMutable) atom0.getOrientation().getDirection();  
//		IAtomList atomList1 = moleculeList.getMolecule(1).getChildList();
//		IAtomOriented atom1 = (IAtomOriented) atomList1.getAtom(0);
//		IVectorMutable  v1 =  (IVectorMutable) atom1.getOrientation().getDirection();  
//		v1.E(0);
//		v1.setX(0, 1);
//		System.out.println("v0 = " + v0);
//		System.out.println("v1 = " + v1);
		
		 torqueSum.reset();
		 potentialMaster.calculate(box, allAtoms, torqueSum);
		 secondDerivativeSum.zeroSum();
		 potentialMaster.calculate(box, allAtoms, secondDerivativeSum);
		
			
		 IteratorDirective id = new IteratorDirective();
		 potentialMaster.calculate(box, id, pcForce);
		 
		 
		 double A = 0;
		 vectorSum.E(0);
		 for (int i = 0;i < nM; i++){
			 dr.E(dipoleSource.getDipole(moleculeList.get(i)));
			 dr.normalize();
			 
			 A += -2.0/3.0*bt2*mu2*(dr.squared()-1);
			 
			 MoleculeAgent torqueAgent = (MoleculeAgent) moleculeAgentManager.getAgent(moleculeList.get(i));
			 dr.XE(torqueAgent.torque);
			 vectorSum.PE(dr);
		 }//i loop
		 
		 //TODO
//		x[0] = ( -2*nM*bt2*mu2+0.25*bt3*mu2*secondDerivativeSum.getSum())/3.0;
//		x[1] = 0.25*bt2*bt2*mu2*vectorSum.squared() ;
		
		
		
		x[0] = -nM*bt2*mu2 - 0.25*bt2*bt2*mu2*vectorSum.squared()+ 0.25*bt3*mu2*secondDerivativeSum.getSum();//TODO
		
		
//		x[1] = vectorSum.getX(0) +vectorSum.getX(1) + vectorSum.getX(2);
//		x[1] = -nM*bt2*mu2 - 0.25*bt2*bt2*mu2*vectorSum.squared()+ 0.25*bt3*mu2*secondDerivativeSum.getSum() + A;
		return data;
	}
	
	public void setDipoleSource(DipoleSource newDipoleSource) {
		dipoleSource = newDipoleSource;
		secondDerivativeSum.setDipoleSource(newDipoleSource); 
	}
	
	public DataTag getTag() {
		return tag;
	}

	public IDataInfo getDataInfo() {
		return dataInfo;
	}

	public Box getBox() {
		return box;
	}


	public Object makeAgent(IMolecule a) {
		// TODO Auto-generated method stub
		 return new MoleculeAgent(space);
	}

	public void releaseAgent(Object agent, IMolecule a) {
		// TODO Auto-generated method stub
		
	}

}
