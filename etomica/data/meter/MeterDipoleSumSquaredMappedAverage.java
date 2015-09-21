package etomica.data.meter;

import etomica.api.IAtomList;
import etomica.api.IBoundary;
import etomica.api.IBox;
import etomica.api.IMolecule;
import etomica.api.IMoleculeList;
import etomica.api.IPotentialMaster;
import etomica.api.ISimulation;
import etomica.api.IVector;
import etomica.api.IVectorMutable;
import etomica.atom.DipoleSource;
import etomica.atom.IAtomOriented;
import etomica.atom.MoleculeAgentManager;
import etomica.atom.MoleculeAgentManager.MoleculeAgentSource;
import etomica.atom.iterator.IteratorDirective;
import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IEtomicaDataInfo;
import etomica.data.IEtomicaDataSource;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.integrator.IntegratorRigidIterative.MoleculeAgent;
import etomica.potential.PotentialCalculationEnergySum;
import etomica.potential.PotentialCalculationPhiSum;
import etomica.potential.PotentialCalculationTorqueSum;
import etomica.space.ISpace;
import etomica.units.Null;

/**
 * meter for AEE use mapping average
 * 
 * @author Weisong
 */
public class MeterDipoleSumSquaredMappedAverage implements IEtomicaDataSource,MoleculeAgentSource  {

	protected final DataDoubleArray data;
	protected final DataInfoDoubleArray dataInfo;
	protected final DataTag tag;
	private IBoundary boundary;
	protected PotentialCalculationTorqueSum torqueSum;
	protected PotentialCalculationPhiSum secondDerivativeSum;
	protected final ISpace space;
	private IBox box;
	private IVectorMutable fieldE;
	private IVectorMutable r;
	private IVectorMutable [] a;
	private double dipoleMagnitude;
	private double truncation;
	private double temperature;
    protected final IPotentialMaster potentialMaster;
    private final IteratorDirective allAtoms;
    protected IVectorMutable dr;
    protected IVectorMutable work;
    protected MoleculeAgentManager moleculeAgentManager;
    protected DipoleSource dipoleSource;
    


	public MeterDipoleSumSquaredMappedAverage(ISpace space, IBox box,ISimulation sim, double dipoleMagnitude,double truncation,double temperature,IPotentialMaster potentialMaster) {
		data = new DataDoubleArray(2);
		dataInfo = new DataInfoDoubleArray("stuff", Null.DIMENSION, new int[]{2});
		tag = new DataTag();
		dataInfo.addTag(tag);
		this.box = box;
		this.dipoleMagnitude = dipoleMagnitude;
		this.temperature = temperature;
		this.space = space;
		this.potentialMaster = potentialMaster;
		this.truncation = truncation;
		fieldE = space.makeVector();
		r = space.makeVector();
		fieldE.setX(2, 1);
		torqueSum = new PotentialCalculationTorqueSum();
		secondDerivativeSum = new  PotentialCalculationPhiSum(space);
		moleculeAgentManager  = new MoleculeAgentManager(sim,box,this);
		torqueSum.setMoleculeAgentManager(moleculeAgentManager);
		allAtoms = new IteratorDirective();
		dr = space.makeVector();
		work = space.makeVector();
		
	}

	public IData getData() {		
		IBoundary boundary = box.getBoundary();
		double[] x = data.getData();
		double bt = 1/(temperature);//beta
		double mu = dipoleMagnitude;//miu
		if (box == null) throw new IllegalStateException("no box");
		IMoleculeList moleculeList = box.getMoleculeList();
		double rx, ry,rz;
		double mxi,myi,mzi;
		double mxj,myj,mzj;
		double Fi = 0;
		double phi = 0;
		double phii = 0;
		double ci,si,cj,sj;
		double rN;
		double sum0 = 0;
		double sum00 =0;
		double sum1 = 0;
		double testSum = 0;
		int nM = moleculeList.getMoleculeCount();	
		double Ftesti = 0;
		double phtestij = 0;
		double phtestii = 0;
		
		if(false)			//debug only for fixing the dipoles  TODO
		{	IAtomList atomListi = moleculeList.getMolecule(0).getChildList();
    		IAtomOriented atomi = (IAtomOriented) atomListi.getAtom(0);		
    		IVectorMutable  vi = (IVectorMutable) atomi.getOrientation().getDirection();  
    		IVectorMutable pi = atomi.getPosition();
    		pi.setX(0,0);
    		pi.setX(1,0);
    		pi.setX(2,0);
    		
    		vi.setX(0,0);
    		vi.setX(1,1);
    		vi.setX(2,0);
    		vi.normalize();
    		
    		IAtomList atomListj = moleculeList.getMolecule(1).getChildList();
			IAtomOriented atomj = (IAtomOriented) atomListj.getAtom(0);
			IVectorMutable  vj = (IVectorMutable) atomj.getOrientation().getDirection();
			IVectorMutable pj = atomj.getPosition();
			pj.setX(0,1.05);
    		pj.setX(1,0);
    		pj.setX(2,0);
    		
    		vj.setX(0,0);
    		vj.setX(1,1);
    		vj.setX(2,0);
		}
		
		
		 torqueSum.reset();
		 potentialMaster.calculate(box, allAtoms, torqueSum);
		 secondDerivativeSum.zeroSum();
		 potentialMaster.calculate(box, allAtoms, secondDerivativeSum);
		
		 
	        for(int k = 0; k < 3; k++){
//		 for(int k = 0; k < 1; k++){							
	        	fieldE.E(0);
	        	fieldE.setX(k, 1);
	        	
	        	for (int i = 0;i < nM; i++){
	        		IAtomList atomListi = moleculeList.getMolecule(i).getChildList();
	        		
	        		dr.E(dipoleSource.getDipole(moleculeList.getMolecule(i)));
//	        		IAtomOriented atomi = (IAtomOriented) atomListi.getAtom(0);		
//	        		IVector  vi =  atomi.getOrientation().getDirection();  
//	        		IVectorMutable pi = atomi.getPosition();
	        		dr.normalize();
	        		ci = dr.dot(fieldE);//the xi or cos thetai
	        		si = Math.sqrt(1-ci*ci);// sin thetai

	        		MoleculeAgent torqueAgent = (MoleculeAgent) moleculeAgentManager.getAgent(moleculeList.getMolecule(i));
	        		// method 1
	        		dr.XE(fieldE);//dr is not unit vector here so we divid si again in Ftesti
//	        		System.out.println("torque = " + torqueAgent.torque);
	        		Ftesti = torqueAgent.torque.dot(dr)/si/si;
//	        		System.out.println("k = " + k + "Ftesti = " +  Ftesti);
	        		
	        		
	        		if(false){	
	        			dr.E(dipoleSource.getDipole(moleculeList.getMolecule(i)));//TODO!!!!!!!!!!>>????????????????
	        			dr.normalize();
	        		mxi = mu * dr.getX(k%3);
	        		myi = mu * dr.getX((k+1)%3);
	        		mzi = mu * dr.getX((k+2)%3);
	        		
	        		Fi = 0;
	        		sum00 = 0;
	        		phii = 0;
	        		phtestij = 0;
	        		
	        		for(int j = 0; j< nM;j++){
	        			if(j == i) 	continue;      			
	        			IAtomList atomListj = moleculeList.getMolecule(j).getChildList();
	        			IAtomOriented atomj = (IAtomOriented) atomListj.getAtom(0);
	        			IVector  vj = atomj.getOrientation().getDirection();
	        			IVectorMutable pj = atomj.getPosition();
//	        			r.Ev1Mv2(pj, pi);
	        			boundary.nearestImage(r);
	        			rN = Math.sqrt(r.squared());    
	        			
	        			//Debug only
//	        			System.out.println("rN = " + rN + " rCut = " +  truncation);
//	        			System.out.println("r = " + r);
//	        			System.out.println("ei = " + dr);
//	        			System.out.println("ej = " + vj);
	        			
	        			if(rN > truncation ){
	        				Fi += 0;
	        				phi = 0;
	        				phii += 0;
	        				continue;
	        			}
	        			
	        			rx = r.getX(k%3);
	        			ry = r.getX((k+1)%3);
	        			rz = r.getX((k+2)%3);
	        			mxj = mu * vj.getX(k%3);
	        			myj = mu * vj.getX((k+1)%3);
	        			mzj = mu * vj.getX((k+2)%3);
	        			cj = vj.dot(fieldE);
	        			sj = Math.sqrt(1-cj*cj);
	        			
	        			Fi += ( (   3*(rx*mxi+ry*myi)*(rx*mxj+ry*myj) - rN*rN*(mxi*mxj+myi*myj) + 3*rz*mu*(rx*mxi+ry*myi)*cj  )*ci/si 
	        					+ (-3*rz*(rx*mxj+ry*myj) + (rN*rN-3*rz*rz)*mu*cj )*mu*si )/rN/rN/rN/rN/rN/(-si);

	        			phi = (  mu*si*( 3*rz*( rx*mxj+ry*myj )*cj/sj  + (rN*rN-3*rz*rz)*mu*sj )	+
	        					( (-3*(rx*mxi+ry*myi)*(rx*mxj+ry*myj) + rN*rN*(mxi*mxj + myi*myj) )*cj/sj 
	        							+ 3*rz*mu*(rx*mxi+ry*myi)*sj    )*ci/si)/rN/rN/rN/rN/rN/(-si)/(-sj);
	        			
	        			phii +=  ( 3.0*(rx*mxi+ry*myi)*(rx*mxj+ry*myj) - rN*rN*(mxi*mxj+myi*myj) +
	        					3.0*rz*mu*(rx*mxi+ry*myi)*cj + mu*ci*(3.0*rz*(rx*mxj+ry*myj) - (rN*rN - 3.0*rz*rz)*mu*cj)
	        					)/rN/rN/rN/rN/rN/(-si)/(-si);
	        			
	        			sum00 += 0.25*bt*bt*bt*mu*mu*phi*(ci*ci-1)*(cj*cj-1);				//phij and phji
	        			
	        		}//j loop
//	        		System.out.println("phi = \t" + phi + "\t phii = \t" + phii);
//	        		System.out.println("phtestij = \t" + phtestij + "\t phtestii = \t" + phtestii);
//	        		System.out.println("i: " + i + "\tFtesti = " + Ftesti + "\tFi = " + Fi);
//	        		System.exit(2);
	        		sum00 += 0.25*bt*bt*bt*mu*mu*phii*(ci*ci-1)*(ci*ci-1);  			 //phii and phijj 
	        		
//	        		sum0 += sum00;
	        		
	        		
//	        		sum0 += sum00 - Fi*ci*(ci*ci-1)*2.0/3.0 - bt*bt*mu*mu/3.0 
//	        				- 0.25*bt*bt*bt*bt*mu*mu*Fi*Fi*(ci*ci-1)*(ci*ci-1) ;	
//	        		sum0 += - Fi*ci*(ci*ci-1)*2.0/3.0 - bt*bt*mu*mu/3.0 
//	        				- 0.25*bt*bt*bt*bt*mu*mu*Fi*Fi*(ci*ci-1)*(ci*ci-1) ;		//!!!!!! debug only //old method
	        		}else{
	        		sum0 += - Ftesti*ci*(ci*ci-1)*2.0/3.0 - bt*bt*mu*mu/3.0 
	        				- 0.25*bt*bt*bt*bt*mu*mu*Ftesti*Ftesti*(ci*ci-1)*(ci*ci-1) ;
	        		
	        		sum1 += 0.5*bt*bt*mu*mu*(Ftesti*(ci*ci-1));
	        		}
	        	}//i loop
	        }	//k loop
	        
	        
//	        System.out.println("new = \t" + 0.25*bt*bt*bt*mu*mu*secondDerivativeSum.getSum());
//	        System.out.println("old = \t" + sum0 );
//	        System.out.println("old/new = \t" + sum0/(0.25*bt*bt*bt*mu*mu*secondDerivativeSum.getSum()));
//	        System.exit(2);
	        
		x[0] = sum0 + 0.25*bt*bt*bt*mu*mu*secondDerivativeSum.getSum();
		x[1] = sum1;				   
		return data;
	}
	
	public void setDipoleSource(DipoleSource newDipoleSource) {
		dipoleSource = newDipoleSource;
		secondDerivativeSum.setDipoleSource(newDipoleSource); 
	}
	
	public DataTag getTag() {
		return tag;
	}

	public IEtomicaDataInfo getDataInfo() {
		return dataInfo;
	}

	public IBox getBox() {
		return box;
	}


	public Object makeAgent(IMolecule a) {
		// TODO Auto-generated method stub
		 return new MoleculeAgent(space);
	}

	public void releaseAgent(Object agent, IMolecule a) {
		// TODO Auto-generated method stub
		
	}

	public Class getMoleculeAgentClass() {
		// TODO Auto-generated method stub
		return MoleculeAgent.class;
	}
	

}