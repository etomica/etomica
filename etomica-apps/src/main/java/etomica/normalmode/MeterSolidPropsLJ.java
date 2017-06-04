package etomica.normalmode;
  
import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.space.Vector;
import etomica.atom.AtomLeafAgentManager;
import etomica.atom.AtomLeafAgentManager.AgentSource;
import etomica.atom.iterator.IteratorDirective;
import etomica.data.DataSourceScalar;
import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IEtomicaDataInfo;
import etomica.data.IEtomicaDataSource;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.integrator.IntegratorVelocityVerlet.MyAgent;
import etomica.potential.PotentialCalculation;
import etomica.potential.PotentialCalculationForceSum;
import etomica.potential.PotentialMaster;
import etomica.space.Space;
import etomica.units.Null;

public class MeterSolidPropsLJ implements IEtomicaDataSource, AgentSource<MyAgent> {

    protected final DataDoubleArray data;
    protected final DataInfoDoubleArray dataInfo;
    protected final DataTag tag;
    protected final Space space;
    protected final CoordinateDefinition coordinateDefinition;
    protected final DataSourceScalar meterPE;
    protected final PotentialMaster potentialMaster;
    protected final AtomLeafAgentManager<MyAgent> forceManager;
    protected final IteratorDirective id;
    protected final Vector dr;
    protected double ULat, PLat;
    protected double volume, density;
    protected int nMol, nRowsA, nColumnsA, nCells, nBasis;
    protected double dP, ddP, f1, f2, dc11;
    protected double fe, fee;
    protected double dU_est;
    protected double dSigma, gamma , ex;
//    protected double[][] A;
//    protected final IVectorMutable[] posCells;
//    protected final IVectorMutable[] kVectors;
//    protected final int[] kDeg;
//    protected final  BasisCell[] cells;
    protected double gruneisen_knu;
    
    protected final double temperature;
    private final PotentialCalculation pcSolidProps;
    protected  PotentialCalculationForceSum pcForceSum;
    
    public MeterSolidPropsLJ(Space space, DataSourceScalar meterPE, PotentialMaster potentialMaster, CoordinateDefinition coordinateDefinition, double temperature, double dP, double dU, double ddP, double ULat, double PLat, double dSigma, double dc11, double gamma, double ex) {
        int nData = 2;
//        int nData = 18;
        data = new DataDoubleArray(nData);
        dataInfo = new DataInfoDoubleArray("stuff", Null.DIMENSION, new int[]{nData});
        tag = new DataTag();
        dataInfo.addTag(tag);
        this.space = space;
        this.coordinateDefinition = coordinateDefinition;
        this.meterPE = meterPE;
        this.potentialMaster = potentialMaster;
        id = new IteratorDirective();
        pcForceSum = new PotentialCalculationForceSum();
        forceManager = new AtomLeafAgentManager<MyAgent>(this, coordinateDefinition.getBox(), MyAgent.class);
        pcForceSum.setAgentManager(forceManager);
        dr = space.makeVector();
        this.temperature = temperature;
        this.dP = dP;
        this.dU_est = dU;
        this.ddP = ddP;
        this.ULat = ULat;
        this.PLat = PLat;
        this.dSigma = dSigma;
        this.dc11 = dc11;
        this.gamma = gamma;
        this.ex = ex;
        
        
        volume = coordinateDefinition.getBox().getBoundary().volume();
        nMol = coordinateDefinition.getBox().getLeafList().getAtomCount();
        density = nMol/volume;
        dP = 9.530521593903*temperature;
    	f1 = (-1.0/volume + dP/temperature)/3.0/(nMol-1.0);
    	f2 = (1.0/volume/volume + ddP/temperature)/3.0/(nMol-1.0)  +  f1*f1;
    	
    	fe  = -(1+dSigma/temperature*volume)/(3.0*(nMol-1));
    	fee = fe*fe + (1-dc11*volume/temperature)/(3.0*(nMol-1));
    	
		pcSolidProps = new PotentialCalculationLJSP(space,coordinateDefinition.getBox(),coordinateDefinition,temperature,dP, f1, fe, fee);


//        pcUP.setAgentManager(forceManager);

//        this.A = A;
//        nRowsA = A.length;
//        nColumnsA = A[0].length;
//        kDeg = new int[nRowsA];
//        kVectors = new IVectorMutable[nRowsA];
//        for(int k=0;k<nRowsA;k++){
//        	kDeg[k] = (int)A[k][0];
//        	kVectors[k] = space.makeVector();
//        	kVectors[k].setX(0, A[k][1]); kVectors[k].setX(1, A[k][2]); kVectors[k].setX(2, A[k][3]);
//        }
//
//        cells  = coordinateDefinition.getBasisCells();
//        nCells = cells.length;
//        int cornerAtomIdx;
//        posCells = new IVectorMutable[nCells];
//        for(int l=0;l<nCells;l++){
//        	posCells[l] = space.makeVector();
//        	cornerAtomIdx = cells[l].molecules.getMolecule(0).getIndex();
//        	posCells[l].E(coordinateDefinition.getBox().getLeafList().getAtom(cornerAtomIdx).getPosition());
//        }
//        
//        nBasis = cells[0].molecules.getMoleculeCount();
    }
    

    //U_direct is computed in the main classes using MeterPotentialEnergyFromIntegrator.
    public IData getData() {
        Box box = coordinateDefinition.getBox();
        double[] sum;
        ((PotentialCalculationLJSP)pcSolidProps).reset();
        potentialMaster.calculate(box, id, pcSolidProps);
        sum = ((PotentialCalculationLJSP)pcSolidProps).getSum();

        potentialMaster.calculate(box, id, pcForceSum);

        double[] x = data.getData();
        //ALL with "ij" pairs.
        double fdr      = sum[0];
//        double fr       = sum[1];
//        double fR       = sum[2];
//        double T1PhiT1  = sum[3];
//        double T1Phidr  = sum[4];
//        double rPhir    = sum[5];
//        double drPhidr  = sum[6];
//        double fxrx     = sum[7];
//        double fxRx     = sum[8];
//        double fxrz     = sum[9];
//        double fxRz     = sum[10];
//        double T2PhiT2  = sum[11];
//        double T3PhiT3  = sum[12];
//        double fxdrx    = sum[13];


        
        double dU    = meterPE.getDataAsScalar() - ULat;
        double Ur   = dU + 0.5*fdr;
//        double Ur2   = dU*(2.0 + 1/(3*(nMol-1.0)*temperature)*fdr);
        
//        double Ur2   = dU + dU_est/(3*(nMol-1.0)*temperature)*fdr;
//
//        double Pvir = fr/3.0/volume;
//        double Pr   = 1.0/3.0/volume*fR + f1*fdr - PLat;
//        double sig11r = - (fxrx + fe*fdr)/volume;
//    	fe  = -(1+dSigma/temperature*volume)/(3.0*(nMol-1));

//        System.out.println(dU/nMol  + "   " + ((3.0/2.0*(nMol-1.0)*temperature +Ur)/nMol));
//        x[10] = Ur2;
        
        
        x[0] =  dU;//U
        x[1] = Ur; //Um   
        
//        x[2] = Pvir;//P	    
//        x[3] = Pr; //Pm
//        
//        x[4] = dU*dU/temperature/temperature;//Cv
//        x[5] = -1.0/2.0/temperature*(0.5*fdr+drPhidr) + Ur*Ur/temperature/temperature;//Cvm
//        
//        x[6] = 2.0/3.0*Pvir + 1.0/9.0/volume*rPhir - volume/temperature*Pvir*Pvir; //B
//        x[7] = -volume*(-2.0/9.0/volume/volume*fR + f2*fdr  - T1PhiT1) - volume/temperature*Pr*Pr; //Bm
//        
//        x[8] = 1.0/temperature/temperature*dU*Pvir; //dPdT
//        x[9] = 1.0/temperature*(f1/2.0*fdr - 1.0/2.0*T1Phidr) + 1.0/temperature/temperature*Ur*Pr; //dPdTm
//
//        
////sigma_11
//        x[10] = -1/volume*fxrx;//conv
//        x[11] = -1/volume*fxrx + (dSigma/temperature*volume+1)/(3.0*(nMol-1)*volume)*fdr;// Map1
//        x[12] = -1/volume*fxRx ;//+ (dSigma/temperature*volume+1)/(3.0*(nMol-1)*volume)*fdr;//Map2
////sigma_12
//        x[13] = -1/volume*fxrz;//conv
//        x[14] = -1/volume*fxrz + dSigma/temperature/(3.0*(nMol-1))*fdr;//Map1
//
//        x[15] = -1/volume*fxRz + dSigma/temperature/(3.0*(nMol-1))*fdr;//Map2
//
//        double dx,dy,dz;
//        double sig11C = -1/volume*fxrx;
//        x[16] =   1/volume*T2PhiT2 - volume/temperature*sig11C*sig11C; //c11C
//        x[17] = -2*fe/volume*fxdrx - fee/volume*fdr + 1/volume*T3PhiT3 - volume/temperature*sig11r*sig11r;//c11M
//        
//        
//        int iBasis;
//        IVectorMutable rBasis, rLatBasis;
//        IVectorMutable  drBasis = space.makeVector();
//        IAtom atomBasis;
//        IVectorMutable[] dr_kR     = new IVectorMutable[nBasis];
//        IVectorMutable[] dr_kI     = new IVectorMutable[nBasis];
//        IVectorMutable[] eVec_knuR = new IVectorMutable[nBasis];
//        IVectorMutable[] eVec_knuI = new IVectorMutable[nBasis];
//        IVectorMutable[] f_kR      = new IVectorMutable[nBasis];
//        IVectorMutable[] f_kI      = new IVectorMutable[nBasis];
//        IVectorMutable[] sumk_R    = new IVectorMutable[nBasis];
//        IVectorMutable[] sumk_I    = new IVectorMutable[nBasis];
//        
//        for(int b=0;b<nBasis;b++){
//          dr_kR[b] = space.makeVector();     dr_kI[b] = space.makeVector();
//      	  eVec_knuR[b] = space.makeVector(); eVec_knuI[b] = space.makeVector();
//	      f_kR[b] = space.makeVector();      f_kI[b] = space.makeVector();
//          sumk_R[b] = space.makeVector();    sumk_I[b] = space.makeVector();
//        }
//        double sumFR = 0;
//        
//        
//        
//        for(int l=0;l<nCells;l++){
//    	  
//          for(int b=0;b<nBasis;b++){
//            sumk_R[b].E(0);
//            sumk_I[b].E(0);
//          }
//        
//          
//          
//  	      for(int k=0;k<nRowsA;k++){
//  	    	  
//	        for(int b=0;b<nBasis;b++){
//	          dr_kR[b].E(0);
//	          dr_kI[b].E(0);
//	        }
//	        for(int lp=0;lp<nCells;lp++){
//	          for(int b=0;b<nBasis;b++){
//	            iBasis  = cells[lp].molecules.getMolecule(b).getIndex();
//	            atomBasis= coordinateDefinition.getBox().getLeafList().getAtom(iBasis);
//	            rLatBasis = coordinateDefinition.getLatticePosition(atomBasis);
//	            rBasis =   atomBasis.getPosition();
//	            drBasis.Ev1Mv2(rBasis, rLatBasis);
//	            
////	            if(lp==l && b==0 && k==0) System.out.println(drBasis);
//                box.getBoundary().nearestImage(drBasis);
//	            dr_kR[b].PEa1Tv1( Math.cos(kVectors[k].dot(posCells[lp])), drBasis);
//	            dr_kI[b].PEa1Tv1(-Math.sin(kVectors[k].dot(posCells[lp])), drBasis);
//	          }//b
//	        }//l_prime
//	        
//	      
//	      
//	        for (int b=0;b<nBasis;b++){
//	          f_kR[b].E(0);
//	          f_kI[b].E(0);
//	        }
//	      
//	//mode
//
//  	        for(int nu=0;nu<3*nBasis;nu++){
//	          gruneisen_knu = A[k][4+nu*25]/(3*volume);
////	          System.out.println(gruneisen_knu);
////	          if(kDeg[k]==2) System.out.println(gruneisen_knu);
//	          
//	          for (int b=0;b<nBasis;b++){
//	            for(int i=0;i<3;i++){
//	              eVec_knuR[b].setX(i, A[k][5 + nu*(2*3*nBasis+1) + b*6 + i*2]);   
//	              eVec_knuI[b].setX(i, A[k][5 + nu*(2*3*nBasis+1) + b*6 + i*2 + 1]);   
//	            }//i
//	          }//n
//	          double krDotR = 0, krDotI = 0;
//	          for (int b=0;b<nBasis;b++){
//	            krDotR += eVec_knuR[b].dot(dr_kR[b]) + eVec_knuI[b].dot(dr_kI[b]);
//	            krDotI += eVec_knuR[b].dot(dr_kI[b]) - eVec_knuI[b].dot(dr_kR[b]);
////	            tmpSum += eVec_knuR[b].dot(eVec_knuR[b]) + eVec_knuI[b].dot(eVec_knuI[b]);//Is is NORMALIZED, yay!
//	          }
//	          for (int b=0;b<nBasis;b++){
//	            f_kR[b].PEa1Tv1( gruneisen_knu*krDotR , eVec_knuR[b] );
//	            f_kR[b].PEa1Tv1(-gruneisen_knu*krDotI , eVec_knuI[b] );
//                f_kI[b].PEa1Tv1( gruneisen_knu*krDotR , eVec_knuI[b] );
//	            f_kI[b].PEa1Tv1( gruneisen_knu*krDotI , eVec_knuR[b] );
//	          }
//	        }//nu
//	      
//
//            for(int b=0;b<nBasis;b++){
//              sumk_R[b].PEa1Tv1( Math.cos(kVectors[k].dot(posCells[l])), f_kR[b]);
//	          sumk_R[b].PEa1Tv1(-Math.sin(kVectors[k].dot(posCells[l])), f_kI[b]);
//	          sumk_I[b].PEa1Tv1(Math.cos(kVectors[k].dot(posCells[l])), f_kI[b]);
//	          sumk_I[b].PEa1Tv1(Math.sin(kVectors[k].dot(posCells[l])), f_kR[b]);
//	          if(kDeg[k] == 2){
//	            sumk_R[b].PEa1Tv1( Math.cos(kVectors[k].dot(posCells[l])), f_kR[b]);
//		        sumk_R[b].PEa1Tv1(-Math.sin(kVectors[k].dot(posCells[l])), f_kI[b]);
//		        sumk_I[b].PEa1Tv1(-Math.cos(kVectors[k].dot(posCells[l])), f_kI[b]);
//		        sumk_I[b].PEa1Tv1(-Math.sin(kVectors[k].dot(posCells[l])), f_kR[b]);
//	          }
//	        }//b
////  	      System.out.println("WWWW");//(0.010067032676144283, 0.0030921968986246107, -0.041397144041727785)
//  	      }//k
//	    
//          for(int b=0;b<nBasis;b++){ //(0.007455269472657289, -0.008813100713411792, -0.007843585115492768)
//            sumk_R[b].TE(1.0/nCells);
//            sumk_I[b].TE(1.0/nCells);
//          }
//          for(int b=0;b<nBasis;b++){
//            iBasis  = cells[l].molecules.getMolecule(b).getIndex();
//            atomBasis= coordinateDefinition.getBox().getLeafList().getAtom(iBasis);
//            sumFR += sumk_R[b].dot(((IntegratorVelocityVerlet.MyAgent)forceManager.getAgent(atomBasis)).force);
//          }
//
//        }//l
//	      System.out.println(f1*fdr +"  "+ -0.5*sumFR);
//	      pcForceSum.reset();
        return data;
    }

    public DataTag getTag() {
        return tag;
    }

    public IEtomicaDataInfo getDataInfo() {
        return dataInfo;
    }

    public final MyAgent makeAgent(IAtom a, Box box) {
        return new MyAgent(space);
    }
    
    public void releaseAgent(MyAgent agent, IAtom atom, Box box) {}
}
