package etomica.normalmode;

import java.util.Arrays;

import javax.crypto.spec.IvParameterSpec;

import Jama.EigenvalueDecomposition;
import Jama.Matrix;
import etomica.atom.AtomLeafAgentManager;
import etomica.atom.AtomLeafAgentManager.AgentSource;
import etomica.atom.IAtomList;
import etomica.atom.MoleculeAgentManager;
import etomica.atom.iterator.IteratorDirective;
import etomica.box.Box;
import etomica.data.DataSourceScalar;
import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IDataSource;
import etomica.data.IEtomicaDataInfo;
import etomica.data.IEtomicaDataSource;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.integrator.IntegratorVelocityVerlet.MyAgent;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.molecule.MoleculeAgentManager;
import etomica.potential.IteratorDirective;
import etomica.potential.PotentialCalculationForceSum;
import etomica.potential.PotentialMaster;
import etomica.space.ISpace;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.space.Vector;
import etomica.space3d.Orientation3D;
import etomica.space3d.OrientationFull3D;
import etomica.units.Null;
import etomica.units.dimensions.Null;
import etomica.util.Constants;
/**
 @author  Weisong Lin
 */
//TODO I replace IEtomicaDataSource to IDATASource here
public class MeterDADBWaterTIP4P implements IDataSource, AgentSource<MyAgent> {

    protected final DataDoubleArray data;
    protected final DataInfoDoubleArray dataInfo;
    protected final DataTag tag;
    protected final Space space;
    protected final MoleculeAgentManager latticeCoordinates;
    protected final DataSourceScalar meterPE;
    protected final PotentialCalculationForceSum pcForceSum;
    protected final PotentialMaster potentialMaster;
    protected final AtomLeafAgentManager<MyAgent> forceManager;
    protected final IteratorDirective id;
    protected final Vector dr;
    protected double latticeEnergy;
    protected final double temperature;
    public static boolean justDADB = true;
    public static boolean justU = false;
	protected final Vector ph1h2;
	protected final Vector q;
	protected final Vector totalforce;
	private final Vector centermass;
	
	protected final Vector torque;
    //TODO  I use different MoleculeAgentManager
    public MeterDADBWaterTIP4P(Space space, DataSourceScalar meterPE, PotentialMaster potentialMaster, double temperature, MoleculeAgentManager latticeCoordinates ) {
        int nData = justDADB ? 1 : 9;
        //default is 1?
        data = new DataDoubleArray(nData);
        dataInfo = new DataInfoDoubleArray("stuff", Null.DIMENSION, new int[]{nData});
        tag = new DataTag();
        dataInfo.addTag(tag);
        this.space = space;
        this.latticeCoordinates = latticeCoordinates;
        this.meterPE = meterPE;
        this.potentialMaster = potentialMaster;
        id = new IteratorDirective();
        pcForceSum = new PotentialCalculationForceSum();
        forceManager = new AtomLeafAgentManager<MyAgent>(this, latticeCoordinates.getBox(), MyAgent.class);
        pcForceSum.setAgentManager(forceManager);
        dr = space.makeVector();
        MeterPotentialEnergy meterPE2 = new MeterPotentialEnergy(potentialMaster);
        meterPE2.setBox(latticeCoordinates.getBox());
        latticeEnergy = meterPE2.getDataAsScalar();
        this.temperature = temperature;
		ph1h2 = space.makeVector();
		q = space.makeVector();
		totalforce = space.makeVector();
		centermass = space.makeVector();
		torque = space.makeVector();
    }
    
    public void setLatticeEnergy(double newLatticeEnergy) {
        latticeEnergy = newLatticeEnergy;
    }
    
    public IData getData() {
        Box box = latticeCoordinates.getBox();
        pcForceSum.reset();
        double[] x = data.getData();
        double x0 = meterPE.getDataAsScalar() - latticeEnergy;
        potentialMaster.calculate(box, id, pcForceSum);
        IMoleculeList molecules = box.getMoleculeList();
        double ForceSum = 0;
        double orientationSum = 0;
        
        
        
        //test make k3 rotate back to its nominal orientation
        for(int j = 99; j > 0 &&false; j--){
        	for (int i = 0; i<molecules.getMoleculeCount(); i++){
        		IMolecule molecule = molecules.getMolecule(i);
             	IAtomList leafList = molecule.getChildList();
             	Vector h1 = leafList.getAtom(0).getPosition();
            	Vector h2 = leafList.getAtom(1).getPosition();
            	Vector o = leafList.getAtom(2).getPosition();
    			Vector m = leafList.getAtom(3).getPosition();
    		    OrientationFull3D or = ((MoleculeSiteSource.LatticeCoordinate)latticeCoordinates.getAgent(molecule)).orientation;
                Vector a0 = (Vector) or.getDirection();//om
                Vector a1 = (Vector) or.getSecondaryDirection();//h1h2
//                dr.Ev1Mv2(m, o);
//                dr.normalize();
//                System.out.println("om = " + dr);
//                System.out.println("a0 =  " + a0);
//                System.out.println("a1 = " + a1);
//                if(i == 3){System.exit(2);}
            	dr.Ev1Mv2(h2, h1);
    			dr.normalize();
    			ph1h2.Ea1Tv1((dr.dot(a0)), a0);//ph1h2 is the vector perpendicular to oh1h2 plane
    			ph1h2.ME(dr);
    			ph1h2.normalize();
    			double acos = -a1.dot(ph1h2);
    			if(acos > 1) acos = 1;
    			if(acos < -1) acos = -1;
    			ph1h2.XE(a1);
    			double kappa = Math.signum(ph1h2.dot(a0))* Math.acos(acos);
                h1.ME(o);
    			double lenth = Math.sqrt(h1.squared());
    			h1.normalize();
    			Orientation3D orientation = new Orientation3D(space);
    			orientation.setDirection(h1);
    			orientation.rotateBy(-kappa/(j+1), a0);
    			h1.Ea1Tv1(lenth, orientation.getDirection());
    			h1.PE(o);
    			h2.ME(o);
    		    h2.normalize();
    		    orientation.setDirection(h2);
    			orientation.rotateBy(-kappa/(j+1), a0);
    			h2.Ea1Tv1(lenth, orientation.getDirection());
    			h2.PE(o);
//    			if(i == 0){
//    				System.out.println(kappa);
//    			}
        	}
//        	double u = meterPE.getDataAsScalar() - latticeEnergy;
//        	System.out.println( j*0.01+ " " + u );
        }
//        System.exit(2);
        
        
      //test make configuration translate back to its nominal orientation
        for(int j = 99; j > 0 && false ; j--){
        	for (int i = 0; i<molecules.getMoleculeCount(); i++){
        		
        		IMolecule molecule = molecules.getMolecule(i);
             	IAtomList leafList = molecule.getChildList();
             	Vector h1 = leafList.getAtom(0).getPosition();
            	Vector h2 = leafList.getAtom(1).getPosition();
            	Vector o = leafList.getAtom(2).getPosition();
    			Vector m = leafList.getAtom(3).getPosition();
    			double hmass = leafList.getAtom(0).getType().getMass();
				double omass = leafList.getAtom(2).getType().getMass();
    			centermass.Ea1Tv1(hmass, h1);
				centermass.PEa1Tv1(hmass, h2) ;
				centermass.PEa1Tv1(omass, o) ;
				centermass.TE(1/(2*hmass + omass));
				Vector nominalcentermass = ((MoleculeSiteSource.LatticeCoordinate)latticeCoordinates.getAgent(molecule)).position;
                dr.Ev1Mv2(nominalcentermass,centermass);
                dr.TE(1.0/(j+1));
                o.PE( dr);
                h1.PE(dr);
                h2.PE(dr);
                m.PE(dr);              
                
//    			if(i == 0){
//    				System.out.println(latticeEnergy);
//    			}
        	}
        	double u = meterPE.getDataAsScalar()-latticeEnergy;
        	System.out.println( j*0.01+ " " + u );
        }
        
//        System.out.println(latticeEnergy);
//        System.exit(2);
        
        
        //debug for rotation angle and axis
        for(int j = 99; j > 0 && false; j--){
        	for (int i = 0; i<molecules.getMoleculeCount(); i++){
        		IMolecule molecule = molecules.getMolecule(i);
             	IAtomList leafList = molecule.getChildList();
             	Vector h1 = leafList.getAtom(0).getPosition();
            	Vector h2 = leafList.getAtom(1).getPosition();
            	Vector o = leafList.getAtom(2).getPosition();
    			Vector m = leafList.getAtom(3).getPosition();
    		    OrientationFull3D or = ((MoleculeSiteSource.LatticeCoordinate)latticeCoordinates.getAgent(molecule)).orientation;
                Vector a0 = (Vector) or.getDirection();//om
                Vector a1 = (Vector) or.getSecondaryDirection();//h1h2
//                dr.Ev1Mv2(m, o);
//                dr.normalize();
//                System.out.println("om = " + dr);
//                System.out.println("a0 =  " + a0);
//                System.out.println("a1 = " + a1);
//                if(i == 3){System.exit(2);}
    			
    			Vector om = space.makeVector();
              	Vector hh = space.makeVector();
              	Vector axis = space.makeVector();
                Vector a2 = space.makeVector();
                
                a2.E(a0);
                a2.XE(a1);
                a2.normalize();
                double [][] array = new double [3][3]; 
                a0.assignTo(array[0]);
                a1.assignTo(array[1]);
                a2.assignTo(array[2]);
                Matrix a = new Matrix( array).transpose();
                om.Ev1Mv2(m, o);
                om.normalize();
                hh.Ev1Mv2(h2,h1);
              //try to test shake tolerance effect the purpendicular or not
//                Vector p = space.makeVector();
//                p.E(om);
//                p.TE(hh.dot(om));
//                hh.ME(p);
                
                
                Vector p = space.makeVector();
                p.E(om);
              	p.TE(hh.dot(om));
              	hh.ME(p);
                
                
                hh.normalize();
                a2.E(om);
                a2.XE(hh);
                a2.normalize();
                
                double [][] array1 = new double [3][3]; 
                om.assignTo(array1[0]);
                hh.assignTo(array1[1]);
                a2.assignTo(array1[2]);
                Matrix newa = new Matrix(array1).transpose();
                a = a.inverse();
                Matrix matrix = newa.times(a); 
                double beta = 0;
                
                EigenvalueDecomposition  eigenvalueDecomposition= matrix.eig();
                Matrix eigenvectormatrix = eigenvalueDecomposition.getV();
                double [] eigenValueArray = eigenvalueDecomposition.getRealEigenvalues();
            	double[][] eigenVectors = eigenvectormatrix.transpose().getArrayCopy();
               
            	int best = 0;
            	double value = Math.abs(eigenValueArray[0] - 1) ;
            	if (value > Math.abs(eigenValueArray[1] - 1)){
            		best = 1;
            		value = Math.abs(eigenValueArray[1] - 1);
            	}
            	if(value > Math.abs(eigenValueArray[2] - 1)){
            		best = 2;
            		value = Math.abs(eigenValueArray[2] - 1);
            	}
            	
            	
            	axis.E(eigenVectors[best]);
                
                
                om.Ev1Mv2(m, o);
                hh.Ev1Mv2(h2, h1);
                if(om.dot(a0) > hh.dot(a0)){
                	dr.E(om);
                	dr.XE(a0);
                }
                else{
                	dr.E(hh);
                	dr.XE(a1);
                }
                beta =  Math.signum(axis.dot(dr))*Math.acos((matrix.trace()-1)/2.0);//rotation angle
                
                
//                dr.Ev1Mv2(m, o);
//    			ph1h2.Ev1Mv2(h2, h1);
//    			System.out.println("old  =" + dr.dot(ph1h2));
                double hmass = leafList.getAtom(0).getType().getMass();
				double omass = leafList.getAtom(2).getType().getMass();
                centermass.Ea1Tv1(hmass, h1);
				centermass.PEa1Tv1(hmass, h2) ;
				centermass.PEa1Tv1(omass, o) ;
				centermass.TE(1/(2*hmass + omass));
				
//              h1.ME(o);
//    			double lenth = Math.sqrt(h1.squared());
//    			h1.normalize();
//    			Orientation3D orientation = new Orientation3D(space);
//    			orientation.setDirection(h1);
//    			orientation.rotateBy(-beta/(j+1), a0);
//    			h1.Ea1Tv1(lenth, orientation.getDirection());
//    			h1.PE(o);
//    			h2.ME(o);
//    		    h2.normalize();
//    		    orientation.setDirection(h2);
//    			orientation.rotateBy(-beta/(j+1), a0);
//    			h2.Ea1Tv1(lenth, orientation.getDirection());
//    			h2.PE(o);
                h1.ME(centermass);
                h2.ME(centermass);
                o.ME(centermass);
                m.ME(centermass);
    			double hlength = Math.sqrt(h1.squared());
    			double olength = Math.sqrt(o.squared());
    			double mlength = Math.sqrt(m.squared());
    			h1.normalize();
    			h2.normalize();
    			o.normalize();
    			m.normalize();
    			Orientation3D orientation = new Orientation3D(space);
    			orientation.setDirection(h1);
    			orientation.rotateBy(beta/(j+1), axis);
    			h1.Ea1Tv1(hlength, orientation.getDirection());
    			h1.PE(centermass);
    		    orientation.setDirection(h2);
    			orientation.rotateBy(beta/(j+1), axis);
    			h2.Ea1Tv1(hlength, orientation.getDirection());
    			h2.PE(centermass);
    			orientation.setDirection(o);
    			orientation.rotateBy(beta/(j+1), axis);
    			o.Ea1Tv1(olength, orientation.getDirection());
    			o.PE(centermass);
    			orientation.setDirection(m);
    			orientation.rotateBy(beta/(j+1), axis);
    			m.Ea1Tv1(mlength, orientation.getDirection());
    			m.PE(centermass);
//    			
//    			dr.Ev1Mv2(m, o);
//    			ph1h2.Ev1Mv2(h2, h1);
//    			System.out.println("new =" +dr.dot(ph1h2));
    			
    			
//    			if(i == 0){
//    				System.out.println(beta);
//    			}
        	}
        	
        	double u = meterPE.getDataAsScalar() - latticeEnergy + orientationSum;
        	System.out.println( j*0.01+ " " + u );
        }
        
//        System.exit(2);
        
        for (int i = 0; i<molecules.getMoleculeCount(); i++){
            IMolecule molecule = molecules.getMolecule(i);
        	IAtomList leafList = molecule.getChildList();
			Vector h1force = ((MyAgent)forceManager.getAgent(leafList.getAtom(0))).force();
			Vector h2force = ((MyAgent)forceManager.getAgent(leafList.getAtom(1))).force();
			Vector oforce = ((MyAgent)forceManager.getAgent(leafList.getAtom(2))).force();
			Vector mforce = ((MyAgent)forceManager.getAgent(leafList.getAtom(3))).force();
			totalforce.E(h1force);
			totalforce.PE(h2force);
			totalforce.PE(mforce);
			totalforce.PE(oforce);
			Vector h1 = leafList.getAtom(0).getPosition();
        	Vector h2 = leafList.getAtom(1).getPosition();
        	Vector o = leafList.getAtom(2).getPosition();
			Vector m = leafList.getAtom(3).getPosition();
			double hmass = leafList.getAtom(0).getType().getMass();
			double omass = leafList.getAtom(2).getType().getMass();
			
			centermass.Ea1Tv1(hmass, h1);
			centermass.PEa1Tv1(hmass, h2) ;
			centermass.PEa1Tv1(omass, o) ;
			centermass.TE(1/(2*hmass + omass));
			
			dr.Ev1Mv2(m, centermass);
			torque.E(mforce);
			torque.XE(dr);
			q.E(torque);
			dr.Ev1Mv2(h1, centermass);
			torque.E(h1force);
			torque.XE(dr);
			q.PE(torque);
			dr.Ev1Mv2(h2, centermass);
			torque.E(h2force);
			torque.XE(dr);
			q.PE(torque);
			dr.Ev1Mv2(o,centermass);
			torque.E(oforce);
			torque.XE(dr);
			q.PE(torque);
			//for the total torque q
			
            Vector lPos = ((MoleculeSiteSource.LatticeCoordinate)latticeCoordinates.getAgent(molecule)).position;
			dr.Ev1Mv2(centermass, lPos);
            ForceSum += totalforce.dot(dr);
            //get the forcesum!!
            
            
            OrientationFull3D or = ((MoleculeSiteSource.LatticeCoordinate)latticeCoordinates.getAgent(molecule)).orientation;
//            OrientationFull3D newor = new OrientationFull3D(space);
//            newor.E(or);
            Vector a0 = (Vector) or.getDirection();//om
            Vector a1 = (Vector) or.getSecondaryDirection();//h1h2
           
            dr.Ev1Mv2(h2, h1);
            dr.normalize();
          
            Vector axis = space.makeVector();
        	Vector oh1 = space.makeVector();
        	Vector oh2 = space.makeVector();
        	Vector h1lpos = space.makeVector();
          	Vector h2lpos = space.makeVector();
          	Vector om = space.makeVector();
          	Vector hh = space.makeVector();
          	
          	
//        	  Vector p1 = space.makeVector();
//            Vector p2 = space.makeVector();
//            p1.Ev1Mv2(dr, a1);
//            dr.Ev1Mv2(m, o);meterPE.getDataAsScalar() - latticeEnergy
//            dr.normalize();
//            p2.Ev1Mv2(dr, a0);
//            axis.E(p1);
//            axis.XE(p2);
//            double alpha = -1.0*Math.signum(p1.dot(p2))*Math.atan(Math.sqrt(p1.squared()/p2.squared()));
//            if(axis.squared() < 1e-14){
//            	axis.Ea1Tv1(Math.cos(alpha), a1);
//            	axis.PEa1Tv1(Math.sin(alpha), a0);
//            }
//            else{
//            	axis.normalize();
//            }
//            om.Ev1Mv2(m, o);
//            hh.Ev1Mv2(h2,h1);
//            om.normalize();
//            hh.normalize();
//            System.out.println("om = " + om);
//            System.out.println("hh =  " + hh);
//            double rotationAngle = 0;
//            double omAngle = 0;
//            double hhAngle = 0;
//            double distance = 0;
//            omAngle = Math.acos(om.dot(axis));
//            hhAngle = Math.acos(hh.dot(axis));
//            
//            if(omAngle > hhAngle){
//            	distance = Math.sqrt(om.squared())*Math.sin(omAngle);
//            	om.ME(a0);
//            	rotationAngle = 2.0*Math.asin(Math.sqrt(om.squared()/2/distance));
//            	
//            }
//            else{
//            	distance = Math.sqrt(hh.squared())*Math.sin(hhAngle);
//            	hh.ME(a1);
//            	rotationAngle = 2.0*Math.asin(Math.sqrt(hh.squared()/2/distance));
//            }
////            System.out.println("axis = " + axis);
//            
//            double DUDT = q.dot(axis);
//            orientationSum -= 0.5*rotationAngle*DUDT;
//            newor.rotateBy(-1.0*rotationAngle, axis);
//            System.out.println("a0 = " + a0);
//            System.out.println("a1 = " + a1);
//            System.out.println(newor.getDirection());
//            System.out.println(newor.getSecondaryDirection());
            
            
           
            Vector a2 = space.makeVector();
            a2.E(a0);
            a2.XE(a1);
            a2.normalize();
            double [][] array = new double [3][3]; 
            a0.assignTo(array[0]);
            a1.assignTo(array[1]);
            a2.assignTo(array[2]);
            Matrix a = new Matrix( array).transpose();
            
            om.Ev1Mv2(m, o);
            om.normalize();
            hh.Ev1Mv2(h2,h1);
            
            
            //try to test shake tolerance effect the purpendicular or not
            Vector p = space.makeVector();
            p.E(om);
          	p.TE(hh.dot(om));
          	hh.ME(p);
            hh.normalize();
            
            a2.E(om);
            a2.XE(hh);
            a2.normalize();
            
            double [][] array1 = new double [3][3]; 
            om.assignTo(array1[0]);
            hh.assignTo(array1[1]);
            a2.assignTo(array1[2]);
            Matrix newa = new Matrix( array1).transpose();
            
            
            a = a.inverse();
            Matrix matrix = newa.times(a); 
            double beta = 0;
            EigenvalueDecomposition  eigenvalueDecomposition= matrix.eig();
            Matrix eigenvectormatrix = eigenvalueDecomposition.getV();
            double [] eigenValueArray = eigenvalueDecomposition.getRealEigenvalues();
        	double[][] eigenVectors = eigenvectormatrix.transpose().getArrayCopy();
        	
        	int best = 0;
        	double value = Math.abs(eigenValueArray[0] - 1) ;
        	if (value > Math.abs(eigenValueArray[1] - 1)){
        		best = 1;
        		value = Math.abs(eigenValueArray[1] - 1);
        	}
        	if(value > Math.abs(eigenValueArray[2] - 1)){
        		best = 2;
        		value = Math.abs(eigenValueArray[2] - 1);
        	}
        	axis.E(eigenVectors[best]);
        	
        	
        	if(value > 1e-12){
        		System.out.println(Arrays.toString(eigenValueArray));
        		System.out.println("non of the eigenvalue is close to one, use the closest one instead");
        	}
            om.Ev1Mv2(m, o);
            hh.Ev1Mv2(h2, h1);
            if(om.dot(a0) > hh.dot(a0)){
            	dr.E(om);
            	dr.XE(a0);
            }
            else{
            	dr.E(hh);
            	dr.XE(a1);
            }
            beta =  -1.0*Math.signum(axis.dot(dr))*Math.acos((matrix.trace()-1)/2.0);
            
//            newor.rotateBy(-1.0*beta, axis);
//            System.out.println("a0 = " + a0);
//            System.out.println("a1 = " + a1);
            om.Ev1Mv2(m, o);
            hh.Ev1Mv2(h2,h1);
            om.normalize();
            hh.normalize();
//            System.out.println("om = " + om);
//            System.out.println("hh =  " + hh);
//            System.out.println(newor.getDirection());
//            System.out.println(newor.getSecondaryDirection());
            double DUDT = q.dot(axis);
//            if(i == 2){
////            	System.out.println(axis + " " + beta + " " + DUDT);
////            	System.out.println(Math.sqrt(q.squared()));
//            	System.out.println(" " + DUDT);
//            }
            orientationSum -= 1.5*(beta-Math.sin(beta))/(1-Math.cos(beta))*DUDT;
            // find kapa33
//			dr.Ev1Mv2(h2, h1);
//			dr.normalize();
//			ph1h2.Ea1Tv1((dr.dot(a0)), a0);
//			ph1h2.ME(dr);
//			ph1h2.normalize();
//			double acos = -a1.dot(ph1h2);
//			if(acos > 1) acos = 1;
//			if(acos < -1) acos = -1;
//			ph1h2.XE(a1);
//			double kappa = Math.signum(ph1h2.dot(a0))* Math.acos(acos);
//			
//			
//			
//			//kapa3 1 or theta
//			dr.Ev1Mv2(m, o);
//			dr.normalize();
//			double theta = Math.acos(dr.dot(a0));
//			if(theta > 1) theta = 1;
//			if(theta < -1) theta = -1;
//			Vector a = space.makeVector();
//			a.E(a0);
//			a.XE(dr);
//			a.normalize();
////			double DUDT = q.dot(a);
//			double DUDK = q.dot(dr);
//			double DUDK1 = DUDT/Math.cos(theta/2);
////			System.out.println(kappa +" " +DUDK);
//			double kappa1 = 2*Math.sin(theta/2);
//			orientationSum -= (1-Math.cos(theta))*DUDT/Math.sin(theta) + 0.5*kappa*DUDK;
//			orientationSum -= 0.5*kappa1*DUDK1 + 0.5*kappa*DUDK;
//			orientationSum -= 0.5*kappa*DUDK;

//			h1.ME(o);
//			double lenth = Math.sqrt(h1.squared());
//			h1.normalize();
//			Orientation3D orientation = new Orientation3D(space);
//			orientation.setDirection(h1);
//			orientation.rotateBy(0.002, a0);
//			h1.Ea1Tv1(lenth, orientation.getDirection());
//			h1.PE(o);
//			h2.ME(o);
//		    h2.normalize();
//		    orientation.setDirection(h2);
//			orientation.rotateBy(0.002, a0);
//			h2.Ea1Tv1(lenth, orientation.getDirection());
//			h2.PE(o);
//			double xplus = meterPE.getDataAsScalar() - latticeEnergy;
//			h1.ME(o);
//			h1.normalize();
//			orientation.setDirection(h1);
//			orientation.rotateBy(-0.002*2, a0);
//			h1.Ea1Tv1(lenth, orientation.getDirection());
//			h1.PE(o);
//			h2.ME(o);
//		    h2.normalize();
//		    orientation.setDirection(h2);
//			orientation.rotateBy(-0.002*2, a0);
//			h2.Ea1Tv1(lenth, orientation.getDirection());
//			h2.PE(o);
//			double xminus = meterPE.getDataAsScalar() - latticeEnergy;
//			double d = (xplus - xminus)/(2*0.002); 
//			h1.ME(o);
//			h1.normalize();
//			orientation.setDirection(h1);
//			orientation.rotateBy(0.002, a0);
//			h1.Ea1Tv1(lenth, orientation.getDirection());
//			h1.PE(o);
//			h2.ME(o);
//		    h2.normalize();
//		    orientation.setDirection(h2);
//			orientation.rotateBy(0.002, a0);
//			h2.Ea1Tv1(lenth, orientation.getDirection());
//			h2.PE(o);
//			double diff = (d - DUDK)/DUDK;
//			if(Math.abs(diff) > 0.01){
//				System.out.println("d =  " + d);
//				System.out.println("DUDK = " + DUDK);
//			}
			
        }
        
        if (justDADB) {
            if (justU) {
                x[0] = (x0+latticeEnergy) + (molecules.getMoleculeCount()*3*temperature) + 0.5*ForceSum + orientationSum;
            }
            else {
//                System.out.println(x0+" "+(0.5*sum)+" "+(x0+0.5*sum)/atoms.getAtomCount());
                x[0] = x0 + 0.5*ForceSum + orientationSum; //translation and rotation
//                x[0] = x0 + 0.5*ForceSum ;//only translation
//                x[0] = x0 + orientationSum;//only rotation
            }
            
//            if (Math.random() < 0.01) {
//                System.out.println(0+" "+(x0+latticeEnergy));
//                for (int j=0; j<99; j++) {
//                    for (int i=0; i<atoms.getAtomCount(); i++) {
//                        IAtom atom = atoms.getAtom(i);
//                        IVector lPos = coordinateDefinition.getLatticePosition(atom);
//                        Vector pos = atom.getPosition();
//                        dr.Ev1Mv2(pos, lPos);
//                        pos.PEa1Tv1(-1.0/(100-j), dr);
//                    }
//                    double u = meterPE.getDataAsScalar();
//                    System.out.println(j+1+" "+(u-latticeEnergy));
//                }
//                System.exit(1);
//            }
            if(data.isNaN()){
            	throw new RuntimeException();
            }
            return data;
        }
        x[0] = x0;
        x[1] = ForceSum;
        x[2] = x[0]+0.5*x[1];
        // so this is the energy
        x[3] = (2*temperature-x[0])*x[0];   // x[0]*x[0] = 2*T*x[0] - x[3] 
        x[4] = (2*temperature-x[0])*x[2];
//        x[5] = (x[0]*x[0]+2*temperature*temperature)*x[2];
//        x[5] = (-6*temperature*(temperature - x[0]) - x[0]*x[0])*x[2];
        x[5] = x[0]*x[0];
        x[6] = x[0]*x[0]*x[2];
        //what does x[3-7] used for?
        
//        x[7] = x[0]*x[1];
        //dAc/dT = -(uAvg+0.5*fdrAvg))/ T^2
        //       = -<x2>/T^2
        //d2Ac/dT2 = ( (2*temperature*(uAvg + 0.5*fdrAvg) - u2Avg + uAvg*uAvg + 0.5*(uAvg*fdrAvg - ufdrAvg) ) / T^4
        //         = (<x4> + <x0>*<x2>)/T^4;
        //d3Ac/dT3 = ( -6*T*T*(uAvg+0.5fdrAvg) + 6*T*(u2Avg - uAvg + 0.5*(ufdrAvg - uAvg*fdrAvg)) - uuu# - 0.5*uuf#)
        //  xyz# = <xyz> - <x><yz> - <y><xz> - <z><xy> + 2<x><y><z>
        //d3Ac/dT3 = (<x5> - 6T*<x0><x2> + 3*<u2>*<x2> + 2*<x0>*(-<x0>*<x2> + 0.5*<uF>) ) / T^6
        //d3Ac/dT3 = (<x5> - 6T*<x0><x2> + 3*<x6>*<x2> + 2*<x0>*(-<x0>*<x2> + 0.5*<x7>) ) / T^6
        return data;
    }

    private void Ev1Mv2(Vector a0, double d) {
		// TODO Auto-generated method stub
		
	}
    
    public void debug(){
    	IBox box = latticeCoordinates.getBox();
    	IMoleculeList molecules = box.getMoleculeList();
    	for(int j = 99; j > 0 ; j--){
        	for (int i = 0; i<molecules.getMoleculeCount(); i++){
        		IMolecule molecule = molecules.getMolecule(i);
             	IAtomList leafList = molecule.getChildList();
             	Vector h1 = leafList.getAtom(0).getPosition();
            	Vector h2 = leafList.getAtom(1).getPosition();
            	Vector o = leafList.getAtom(2).getPosition();
    			Vector m = leafList.getAtom(3).getPosition();
    		    OrientationFull3D or = ((MoleculeSiteSource.LatticeCoordinate)latticeCoordinates.getAgent(molecule)).orientation;
                Vector a0 = (Vector) or.getDirection();//om
                Vector a1 = (Vector) or.getSecondaryDirection();//h1h2
//                dr.Ev1Mv2(m, o);
//                dr.normalize();
//                System.out.println("om = " + dr);
//                System.out.println("a0 =  " + a0);
//                System.out.println("a1 = " + a1);
//                if(i == 3){System.exit(2);}
    			
    			Vector om = space.makeVector();
              	Vector hh = space.makeVector();
              	Vector axis = space.makeVector();
                Vector a2 = space.makeVector();
                
                a2.E(a0);
                a2.XE(a1);
                a2.normalize();
                double [][] array = new double [3][3]; 
                a0.assignTo(array[0]);
                a1.assignTo(array[1]);
                a2.assignTo(array[2]);
                Matrix a = new Matrix( array).transpose();
                om.Ev1Mv2(m, o);
                om.normalize();
                hh.Ev1Mv2(h2,h1);
              //try to test shake tolerance effect the purpendicular or not
//                Vector p = space.makeVector();
//                p.E(om);
//                p.TE(hh.dot(om));
//                hh.ME(p);
                
                
                Vector p = space.makeVector();
                p.E(om);
              	p.TE(hh.dot(om));
              	hh.ME(p);
                
                
                hh.normalize();
                a2.E(om);
                a2.XE(hh);
                a2.normalize();
                
                double [][] array1 = new double [3][3]; 
                om.assignTo(array1[0]);
                hh.assignTo(array1[1]);
                a2.assignTo(array1[2]);
                Matrix newa = new Matrix(array1).transpose();
                a = a.inverse();
                Matrix matrix = newa.times(a); 
                double beta = 0;
                
                EigenvalueDecomposition  eigenvalueDecomposition= matrix.eig();
                Matrix eigenvectormatrix = eigenvalueDecomposition.getV();
                double [] eigenValueArray = eigenvalueDecomposition.getRealEigenvalues();
            	double[][] eigenVectors = eigenvectormatrix.transpose().getArrayCopy();
               
            	int best = 0;
            	double value = Math.abs(eigenValueArray[0] - 1) ;
            	if (value > Math.abs(eigenValueArray[1] - 1)){
            		best = 1;
            		value = Math.abs(eigenValueArray[1] - 1);
            	}
            	if(value > Math.abs(eigenValueArray[2] - 1)){
            		best = 2;
            		value = Math.abs(eigenValueArray[2] - 1);
            	}
            	
            	
            	axis.E(eigenVectors[best]);
                
                
                om.Ev1Mv2(m, o);
                hh.Ev1Mv2(h2, h1);
                if(om.dot(a0) > hh.dot(a0)){
                	dr.E(om);
                	dr.XE(a0);
                }
                else{
                	dr.E(hh);
                	dr.XE(a1);
                }
                
                beta =  Math.signum(axis.dot(dr))*Math.acos((matrix.trace()-1)/2.0);//rotation angle
                
                
                
                
                
                
//                dr.Ev1Mv2(m, o);
//    			ph1h2.Ev1Mv2(h2, h1);
//    			System.out.println("old  =" + dr.dot(ph1h2));
                double hmass = leafList.getAtom(0).getType().getMass();
				double omass = leafList.getAtom(2).getType().getMass();
                centermass.Ea1Tv1(hmass, h1);
				centermass.PEa1Tv1(hmass, h2) ;
				centermass.PEa1Tv1(omass, o) ;
				centermass.TE(1/(2*hmass + omass));
				
//              h1.ME(o);
//    			double lenth = Math.sqrt(h1.squared());
//    			h1.normalize();
//    			Orientation3D orientation = new Orientation3D(space);
//    			orientation.setDirection(h1);
//    			orientation.rotateBy(-beta/(j+1), a0);
//    			h1.Ea1Tv1(lenth, orientation.getDirection());
//    			h1.PE(o);
//    			h2.ME(o);
//    		    h2.normalize();
//    		    orientation.setDirection(h2);
//    			orientation.rotateBy(-beta/(j+1), a0);
//    			h2.Ea1Tv1(lenth, orientation.getDirection());
//    			h2.PE(o);
                h1.ME(centermass);
                h2.ME(centermass);
                o.ME(centermass);
                m.ME(centermass);
    			double hlength = Math.sqrt(h1.squared());
    			double olength = Math.sqrt(o.squared());
    			double mlength = Math.sqrt(m.squared());
    			h1.normalize();
    			h2.normalize();
    			o.normalize();
    			m.normalize();
    			Orientation3D orientation = new Orientation3D(space);
    			orientation.setDirection(h1);
    			orientation.rotateBy(beta/(j+1), axis);
    			h1.Ea1Tv1(hlength, orientation.getDirection());
    			h1.PE(centermass);
    		    orientation.setDirection(h2);
    			orientation.rotateBy(beta/(j+1), axis);
    			h2.Ea1Tv1(hlength, orientation.getDirection());
    			h2.PE(centermass);
    			orientation.setDirection(o);
    			orientation.rotateBy(beta/(j+1), axis);
    			o.Ea1Tv1(olength, orientation.getDirection());
    			o.PE(centermass);
    			orientation.setDirection(m);
    			orientation.rotateBy(beta/(j+1), axis);
    			m.Ea1Tv1(mlength, orientation.getDirection());
    			m.PE(centermass);
//    			
//    			dr.Ev1Mv2(m, o);
//    			ph1h2.Ev1Mv2(h2, h1);
//    			System.out.println("new =" +dr.dot(ph1h2));
    			
    			
//    			if(i == 0){
//    				System.out.println(beta);
//    			}
        	}
        	System.out.print( j*0.01 + " " ); 
        	double y = getData().getValue(0);
//        	System.out.println( j*0.01 + " " + y );  

        }
        
    }

	public DataTag getTag() {
        return tag;
    }

    public IEtomicaDataInfo getDataInfo() {
        return dataInfo;
    }

    public final MyAgent makeAgent(IAtom a) {
        return new MyAgent(space);
    }
    
    public void releaseAgent(MyAgent agent, IAtom atom) {}
}
