package etomica.virial;

import java.io.FileWriter;
import java.io.IOException;

import javax.management.RuntimeErrorException;

import Jama.EigenvalueDecomposition;
import etomica.api.IAtomList;
import etomica.api.IAtomType;
import etomica.api.IBox;
import etomica.api.IPotentialMaster;
import etomica.api.IRandom;
import etomica.api.IVectorMutable;
import etomica.atom.AtomHydrogen;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorLeafAtoms;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.integrator.mcmove.MCMoveBoxStep;
import etomica.math.linearalgebra.Matrix;
import etomica.space.ISpace;
import etomica.space.Space;
import etomica.units.BohrRadius;
import etomica.units.Hartree;
import etomica.units.Kelvin;
import etomica.util.Arrays;
import etomica.util.Constants;
import etomica.util.DoubleRange;
import etomica.util.HistogramSimple;

public class MCMoveChangeBondLength extends MCMoveBoxStep {
    protected final AtomIteratorLeafAtoms leafIterator;
    protected double[][] prevBondLength;
    protected final IRandom random;
    protected final MeterPotentialEnergy mpe;
    protected double wOld,wNew;
    protected double kHarmonic,pRatio;
    protected final IVectorMutable dr,dr1,dr2;
    protected double kIntra;    
    protected final double re = 1.401;
    protected double t,r0;
    protected double[] B;
    protected double[][] A;
    public double sD = 0, eAvg = 0, uActAvg1 = 0, count1 = 0, uActAvg2 = 0, count2 = 0, uActAvg3 = 0, count3 = 0;

    protected final double c6 = -6.499027;
    protected final double c8 = -124.3991;
    protected final double c10 = -3285.828;
    protected final double beta = 0.2;    
    protected final double alpha = 3.980917850296971;
    protected final double[] xlp = {-1.709663318869429E-1,-6.675286769482015E-1,-1.101876055072129E+0,-1.106460095658282E+0,-7.414724284525231E-1,-3.487882731923523E-1,-1.276736255756598E-1,-5.875965709151867E-2,-4.030128017933840E-2,-2.038653237221007E-2,-7.639198558462706E-4,2.912954920483885E-3,-2.628273116815280E-4,-4.622088855684211E-4,1.278468948126147E-4,-1.157434070240206E-5,-2.609691840882097E-12};
    protected int P;
    protected int[] leftAdjusted,rightAdjusted;
    protected double u0;    
    public HistogramSimple h1;
    public HistogramSimple h2;
    public boolean deleteMode = true, doHist = true, flagAdjusted;
    public int moveCount;
    protected int nBins = 100;
    protected double xMin = 0.35,xMax = 1.25;
    protected double deltaX = (xMax - xMin)/nBins;
    protected double [] hValues;
    protected double[] pCummulative;
    protected double pTotal = 0.0;    

    public MCMoveChangeBondLength(IPotentialMaster potentialMaster, IRandom random, ISpace space, double temperature) {
        super(potentialMaster);    
        leafIterator = new AtomIteratorLeafAtoms();
        this.random = random;
        mpe = new MeterPotentialEnergy(potentialMaster);
        dr = space.makeVector();
        dr1 = space.makeVector();
        dr2 = space.makeVector();
        setTemperature(temperature);
        h1 = new HistogramSimple(nBins, new DoubleRange(xMin, xMax));
        h2 = new HistogramSimple(1000, new DoubleRange(xMin, xMax));
		hValues = h1.getHistogram();
		pCummulative = new double[nBins];
        moveCount = 0;
        sD = 0;
        leftAdjusted = new int[2];
        rightAdjusted = new int[2];
    }
    public void setTemperature(double x){
        t = x;
    }
    public double v(double pr) {
        if (false) {
            return 100000*(pr - 1)*(pr - 1);
        }
        // H2 singlet potential for FCI/CBS
        double r = BohrRadius.UNIT.fromSim(pr);
        double damp=1.0-Math.exp(-beta*r*r);
        double xd = damp/r ;
        double xd2 = xd*xd;
        double xd6 = xd2*xd2*xd2;
        double xd8 = xd6*xd2;
        double xd10 = xd8*xd2;
        double prefac = Math.exp(- alpha*(r-re));
        double eLong = c6*xd6 + c8*xd8 + c10*xd10;
        double eShort = 0.00;        
        for (int i=0; i<17; i++) {
            eShort += xlp[i]*prefac;
            prefac *= (r-re);
        }
        double f = Hartree.UNIT.toSim(eLong + eShort);  
        return f;
    }
    
    public double dvdr(double pr) {
        if (false) {
            return -200000*(pr - 1)/(t*P);
        }
        double r = BohrRadius.UNIT.fromSim(pr);
        double y = Math.exp(-beta*r*r);        
        double xd = (1.0 - y)/r ;
        double xd2 = xd*xd;
        double xd4 = xd2*xd2;
        double xd6 = xd2*xd2*xd2;
        double xd8 = xd6*xd2;
        double xd10 = xd8*xd2;
        double prefac = Math.exp(- alpha*(r-re));        
        double deLongdr = c6*(12*beta*y*xd4*xd - 6*xd6/r) + c8*(16*beta*y*xd6*xd - 8*xd8/r) + c10*(20*beta*y*xd8*xd - 10*xd10/r);
        double deShortdr = 0.00;
        double term0 = 1.0/(r-re);        
        for (int i=0; i<17; i++) {
            deShortdr += xlp[i]*prefac*term0*(i - alpha*(r-re));
            term0 *= (r-re);
        }
        double dv = 2/(P*pr);
//        System.out.println("eLong = "+deLongdr+" eShort = "+deShortdr);
//        double dv = 2/pr;
        dv += -(Hartree.UNIT.toSim(deLongdr + deShortdr))/(t*P*BohrRadius.UNIT.toSim(1));         
        return dv;
    }
    
    public double d2vdr2(double pr) {
        if (false)return -200000/(t*P);
        double r = BohrRadius.UNIT.fromSim(pr);
        double r1 = 1.0/r;
        double r2 = r1*r1;
        double y = Math.exp(-beta*r*r);
        double xd = (1.0 - y)/r ;
        double xd2 = xd*xd;
        double xd4 = xd2*xd2;
        double xd6 = xd2*xd2*xd2;
        double xd8 = xd6*xd2;
        double xd10 = xd8*xd2;
        double prefac = Math.exp(- alpha*(r-re));
        double d2eL = c6*(42*xd6*r2 - 132*xd4*xd*beta*y*r1 + 120*xd4*beta*beta*y*y -24*r*xd4*xd*beta*beta*y) + c8*(72*xd8*r2 - 240*xd6*xd*beta*y*r1 + 224*xd6*beta*beta*y*y - 32*r*xd6*xd*beta*beta*y) + c10*(110*xd10*r2 - 380*xd8*xd*beta*y*r1 + 360*xd8*beta*beta*y*y - 40*r*xd8*xd*beta*beta*y);
        double d2eS = 0.00;
        double term0 = 1.0/((r-re)*(r-re));        
        for (int i=0; i<17; i++) {
            d2eS += xlp[i]*prefac*term0*((i - alpha*(r-re))*(i - alpha*(r-re)) - i);
            term0 *= (r-re);
        }
        double d2v = -2/(P*pr*pr);
//        double d2v = -2/(pr*pr);
        d2v += -Hartree.UNIT.toSim(d2eL + d2eS)/(t*P*BohrRadius.UNIT.toSim(1)*BohrRadius.UNIT.toSim(1));         
        return d2v;
    }
    public void setBox(IBox p) {
        super.setBox(p);
        int nMolecules = box.getMoleculeList().getMoleculeCount();
        P = box.getMoleculeList().getMolecule(0).getChildList().getAtomCount();
        prevBondLength = new double[nMolecules][P];        
        leafIterator.setBox(p);
//        for (int i=0; i<10; i++) {
//            double r = 1.0 + i/10.0;
//            System.out.println("r = "+BohrRadius.UNIT.fromSim(r)+" dv = "+dvdr(r)+" u = "+v(r)+" d2v = "+d2vdr2(r));            
//        }
//        System.exit(0);
//        
    }
        
    public AtomIterator affectedAtoms() {        
        return leafIterator;
    }
    
    public double energyChange() {    
        return 0;
    }    
    public boolean doTrial() {        
        wOld = ((BoxCluster)box).getSampleCluster().value((BoxCluster)box);
        mpe.setBox(box);        
        
        double uA1Old = 0.0;
        
        int nMolecules = box.getMoleculeList().getMoleculeCount();
        double uGenOld = 0.0;
        double uGenNew = 0.0;
        double uA2Old = 0.0;
        double uA2New = 0.0;
        double uActOld2 = 0;
        double uActNew2 = 0;
        double uActOld3 = 0;
        double uActNew3 = 0;
        int iRand = random.nextInt(2);
//        if (box.getIndex() == 1 & moveCount == 2) System.exit(1);
		
		    		
		
//        if (MCMoveChangeBondLengthBruteForce.U0 == 0) MCMoveChangeBondLengthBruteForce.U0 = nMolecules*P*v(r0);
//        {IAtomList atoms = box.getMoleculeList().getMolecule(0).getChildList();
//        ((AtomHydrogen)atoms.getAtom(0)).setBondLength(0.981838106371022);
//        ((AtomHydrogen)atoms.getAtom(1)).setBondLength(0.7651082776025715);
//        atoms = box.getMoleculeList().getMolecule(1).getChildList();
//        ((AtomHydrogen)atoms.getAtom(0)).setBondLength(0.8231323944420721);
//        ((AtomHydrogen)atoms.getAtom(1)).setBondLength(0.7860125242207836);}
        
//        System.out.println(dvdr(0.7479652840490955));
//        System.exit(1);
        boolean flag = moveCount >= 1000;
        flagAdjusted = false;
        
        if (moveCount > 0) {
        	for (int i=0; i<nMolecules; i++) {        		
        		IAtomList atoms = box.getMoleculeList().getMolecule(i).getChildList();
        		for (int j=0; j<P; j++) {
        			sD += (((AtomHydrogen)atoms.getAtom(j)).getBondLength() - 0.75)*(((AtomHydrogen)atoms.getAtom(j)).getBondLength() - 0.75);
        		}        		
        	}        	
        }
        double avgBL = 0;
        if (deleteMode || doHist) {
        	for (int i=0; i<nMolecules; i++) {
        		avgBL = 0;
        		IAtomList atoms = box.getMoleculeList().getMolecule(i).getChildList();
        		for (int j=0; j<P; j++) {
        			avgBL += ((AtomHydrogen)atoms.getAtom(j)).getBondLength();
        		}
        		avgBL /= P;
//        		if (avgBL > 1) System.out.println("avgBL "+avgBL);
//        		if (avgBL > 0.8) throw new RuntimeException();
            	h1.addValue(avgBL);
            	if (flag) h2.addValue(avgBL);
        	}
        }
        
        if (deleteMode && flag) {
        	if (moveCount %100 == 0) {        		
        		
        		hValues = h1.getHistogram();
        		adjustHistogram(h1);
        		pCummulative[0] = hValues[0];
        		for (int k=1; k<nBins; k++) {            	
        			pCummulative[k] = pCummulative[k-1] + hValues[k];
        		}
        		pTotal = pCummulative[nBins-1];
        		for (int k=0; k<nBins; k++) {        			
        			pCummulative[k] /= pTotal;
        			if (Double.isInfinite(pCummulative[k]) || Double.isInfinite(pTotal) || Double.isNaN(pCummulative[k]) || Double.isNaN(pTotal)) throw new RuntimeException("pCummulative ["+k+"] = "+ pCummulative[k]+" pTotal = "+pTotal);
        		}
        	}        	
        	        	
        }
//        if (moveCount == 0) {
//        	for (int i=0; i<nMolecules; i++) {
//
//        		IAtomList atoms = box.getMoleculeList().getMolecule(i).getChildList();
//        		((AtomHydrogen)atoms.getAtom(0)).setBondLength(0.75);
//        		((AtomHydrogen)atoms.getAtom(1)).setBondLength(0.72);
//        		((AtomHydrogen)atoms.getAtom(2)).setBondLength(0.73);
//        		((AtomHydrogen)atoms.getAtom(3)).setBondLength(0.74);
//        		((AtomHydrogen)atoms.getAtom(4)).setBondLength(0.75);
//        		((AtomHydrogen)atoms.getAtom(5)).setBondLength(0.76);
//        		((AtomHydrogen)atoms.getAtom(6)).setBondLength(0.77);
//        		((AtomHydrogen)atoms.getAtom(7)).setBondLength(0.78);
//        	for (int j=0; j<P; j++) {
//        		((AtomHydrogen)atoms.getAtom(j)).setBondLength(BohrRadius.UNIT.toSim(1.448736));
//        	}
//        	}
//        }
        uA1Old = mpe.getDataAsScalar();// - nMolecules*P*v(r0);
        for (int i=0; i<nMolecules; i++) {
//        	if (box.getIndex() == 0 && moveCount == 1007){
//            	System.out.println(moveCount+" "+box.getIndex()+" molecule = "+i);
//            }
            double[] sigma = new double[P];            
            double[] etaOld = new double[P];
            double[] etaNew = new double[P];            
            IAtomList atoms = box.getMoleculeList().getMolecule(i).getChildList();
            
            for (int j=0; j<P; j++) {
                prevBondLength[i][j] = ((AtomHydrogen)atoms.getAtom(j)).getBondLength();                
            }
            

            if (deleteMode && flag) {
            
//              Identifying the uniform mode
            	int choice = -1;
            	boolean done = false;
            	for (int j=0; j<P && !done; j++) {
            		double mag = 0;
                	for (int k=0; k<P; k++) {
                		mag += A[j][k];
                	}
                	double dummy1 = (mag - Math.sqrt(P))*(mag - Math.sqrt(P));
                	double dummy2 = (mag + Math.sqrt(P))*(mag + Math.sqrt(P));
//                	System.out.println(dummy1 + " "+dummy2);
                	if (dummy1  < 1E-16 || dummy2 < 1E-16) {
                		choice = j;
                		done = true;                		
                	}
                }
                if (choice < 0 || !done ) throw new RuntimeException("No uniform mode found!!");
                double oldAvgBL = 0;
                for (int k=0; k<P; k++) {
                	oldAvgBL += prevBondLength[i][k];
                }
                oldAvgBL /= P;
                for (int j=0; j<P; j++) {                                
                    sigma[j] = Math.sqrt(1/B[j]);
                    etaOld[j] = 0.0;
                    etaNew[j] = 0.0;
                    for (int k=0; k<P; k++) {//mode = j, bondlength = k
                    	if (j != choice) etaOld[j] += A[j][k]*(prevBondLength[i][k] - oldAvgBL);
                	}
                    
//                    System.out.println(etaNew[j]);
                    
//                    double newBL = BohrRadius.UNIT.fromSim(((AtomHydrogen)atoms.getAtom(j)).getBondLength());
//                    System.out.println("j ="+j+" oldBL ="+oldBL+" nextBL ="+nextBL+" newBL ="+newBL);
        
//                    System.exit(0);
                    int prev = j-1;
                    if (prev < 0) prev = P-1;

                    AtomHydrogen jAtom = (AtomHydrogen)atoms.getAtom(j);
                    AtomHydrogen jPrev = (AtomHydrogen)atoms.getAtom(prev);
                    
                    
                    if (fixedOrientation) {
                        uA2Old += 2*Math.log(prevBondLength[i][j])/P - kHarmonic*dist(jAtom,jPrev); // to be used when orientation moves are turned off
                    }
                    else {
                        uA2Old += 2*Math.log(prevBondLength[i][j]) - kHarmonic*dist(jAtom,jPrev); // to be used when orientation moves are turned on
                    }
                    uActOld2 += 2*Math.log(prevBondLength[i][j])/P;
                    uActOld3 += - kHarmonic*dist(jAtom,jPrev);
//                    if (i == 1 && moveCount > 1000 && box.getIndex() == 1) System.out.print(prevBondLength[i][j]+" ");
//                    if (i == 1 && moveCount > 1000 && box.getIndex() == 1) System.out.print(dist(jAtom,jPrev)+" ");
                }
                for (int k=0; k<P; k++) {                	
//            		if (moveCount>1000) System.out.println("etaOld["+k+"] "+etaOld[k]);
                	if (k != choice) {
                		uGenOld += -0.5*etaOld[k]*etaOld[k]*B[k];                		
                	}
                }
//                if (i == 1 && moveCount > 1000 && box.getIndex() == 1) System.out.println();
                
                int oldChoice = (int)((oldAvgBL - xMin)/deltaX);
            	if (hValues[oldChoice] == 0 || pTotal == 0) throw new RuntimeException("oldChoice = "+oldChoice+" hValues = "+hValues[oldChoice]+" "+" pTotal = "+pTotal+" moveCount = "+moveCount+ " molecule = "+i+" bl = "+oldAvgBL+" box = "+box.getIndex());
                uGenOld += Math.log(hValues[oldChoice]/pTotal);

//              choose value from avg. bond length distribution
                int newChoice = -1;
                boolean flag1 = false;
                double[] newBL = new double[P];
                double rStar = 0;
                double newAvgBL = 0;
                
                do {
                	newAvgBL = 0;
                	for (int k=0; k<P; k++) {
                    	
                    	if (k != choice) {
                    		etaNew[k] = random.nextGaussian()*sigma[k];
//                    		if (moveCount>1000) System.out.println("etaNew["+k+"] "+etaNew[k]);
                    		uGenNew += -0.5*etaNew[k]*etaNew[k]*B[k];                		
                    	}
                    }
                	double dummy = random.nextDouble();        		
                	done = false;
                	if (dummy <= pCummulative[0]) {
                		newChoice = 0;
                		done = true;
                	}
//        		dummy = 28;
//        		pCummulative[0] = 0;
//        		for (int k=1; k<nBins; k++) {
//        			pCummulative[k] = k;
//        		}
                	if (!done || newChoice < 0) newChoice = getIndex(dummy,pCummulative,1);
//                	newChoice = 8;
                	if (newChoice > 0) done = true;
//        		System.out.print("newChoice = "+newChoice);
//        		System.exit(1);
                	

                	double y1 = xMin + newChoice*deltaX;
                	double y2 = xMin + (newChoice+1)*deltaX;
                	double x1 = 0;
                	if (newChoice > 0) x1 = pCummulative[newChoice-1];                	
                	double x2 = pCummulative[newChoice];
                	if ((newChoice >= leftAdjusted[0] && newChoice <= leftAdjusted[1]) || (newChoice >= rightAdjusted[0] && newChoice <= rightAdjusted[1])) flagAdjusted = true; 
//                	if (box.getIndex() == 0 && moveCount > 4978 && flagAdjusted) {
//                		System.out.print("newChoice = "+newChoice+" left = "+leftAdjusted[0]+" "+leftAdjusted[1]+" right = "+rightAdjusted[0] +" "+rightAdjusted[1]+" moveCount = "+moveCount+" \n");
//                	}
//                	dummy = x1 + (x2-x1)*dummy;
                	rStar = y1 + (dummy - x1)*(y2 - y1)/(x2 - x1);
//                	if (box.getIndex() == 0 && rStar < 0.725) System.out.println("rStar = "+rStar +" "+i +" "+y1+" "+y2);
                	if (rStar>1) foo = true;
                	uGenNew += Math.log(hValues[newChoice]/pTotal);
                	flag1 = true;                	
                	for (int j=0; j<P; j++) {
                		newBL[j] = rStar;
//                    System.out.println(i+" "+ j+" "+newBL[j]);
                		for (int k=0; k<P; k++) {//mode = k, bondlength = j
                			if (k != choice) newBL[j] += A[k][j]*etaNew[k];
                		}
                		newAvgBL += newBL[j];
                		if (newBL[j] < 0) {
                			flag1 = false;
                			break;
                		}                		

//                    h1.addValue(newBL[j]);
//                    System.out.println(i+" "+ j+" "+newBL[j]);
                		((AtomHydrogen)atoms.getAtom(j)).setBondLength(newBL[j]);
                	}
                	newAvgBL /= P;
                } while(!flag1);
//        		System.out.println(oldChoice+" "+Math.log(hValues[oldChoice]/pTotal)+" "+newChoice+" "+Math.log(hValues[newChoice]/pTotal));
//                if (box.getIndex() == 0 && moveCount > 1000) System.out.println(oldAvgBL+" "+newAvgBL);
//                if (newAvgBL < 0.71) throw new RuntimeException(); 
//                for (int k=0; k<P; k++) {
//                	System.out.println("eta ["+k+"] = "+etaOld[k] + " " +etaNew[k]);
//                }
//                if (foo) {
//	        		avgBL = 0;
//	        		for (int j=0; j<P; j++) {
//	        			avgBL += ((AtomHydrogen)atoms.getAtom(j)).getBondLength();
//	        		}
//	        		avgBL /= P;
//	        		System.out.println(rStar+" "+avgBL);
//                }

                for (int j=0; j<P; j++) {
                    int prev = j-1;
                    if (prev <0) prev = P-1;
                    AtomHydrogen jAtom = (AtomHydrogen)atoms.getAtom(j);
                    AtomHydrogen jPrev = (AtomHydrogen)atoms.getAtom(prev);
                    if (fixedOrientation) {
                        uA2New += 2*Math.log(newBL[j])/P - kHarmonic*dist(jAtom,jPrev); // to be used when orientation moves are turned off
                    }
                    else {
                        uA2New += 2*Math.log(newBL[j]) - kHarmonic*dist(jAtom,jPrev); // to be used when orientation moves are turned on
                    }                    
                    uActNew2 += v(jAtom.getBondLength());
                    uActNew3 += - kHarmonic*dist(jAtom,jPrev);
//                    if (i == 1 && moveCount > 1000 && box.getIndex() == 1) System.out.print(newBL[j] +" ");
//                    if (i == 1 && moveCount > 1000 && box.getIndex() == 1) System.out.print(dist(jAtom,jPrev) +" ");
                }
                if (newChoice == 52) {
                	count1++;
                	for (int j=0; j<P; j++) {
                		int prev = j-1;
                        if (prev <0) prev = P-1;
                        AtomHydrogen jAtom = (AtomHydrogen)atoms.getAtom(j);
                        AtomHydrogen jPrev = (AtomHydrogen)atoms.getAtom(prev);
                        uActAvg1 += -v(jAtom.getBondLength())/(t*P);
                        if (fixedOrientation) {
                            uActAvg1 += 2*Math.log(newBL[j])/P - kHarmonic*dist(jAtom,jPrev); // to be used when orientation moves are turned off
                        }
                        else {
                            uActAvg1 += 2*Math.log(newBL[j]) - kHarmonic*dist(jAtom,jPrev); // to be used when orientation moves are turned on
                        }	
                	}
                	  
                }
                if (newChoice == 58) {
                	count2++;
                	for (int j=0; j<P; j++) {
                		int prev = j-1;
                        if (prev <0) prev = P-1;
                        AtomHydrogen jAtom = (AtomHydrogen)atoms.getAtom(j);
                        AtomHydrogen jPrev = (AtomHydrogen)atoms.getAtom(prev);
                        uActAvg2 += -v(jAtom.getBondLength())/(t*P);
                        if (fixedOrientation) {
                            uActAvg2 += 2*Math.log(newBL[j])/P - kHarmonic*dist(jAtom,jPrev); // to be used when orientation moves are turned off
                        }
                        else {
                            uActAvg2 += 2*Math.log(newBL[j]) - kHarmonic*dist(jAtom,jPrev); // to be used when orientation moves are turned on
                        }	
                	}
                	  
                }
                if (newChoice == 64) {
                	count3++;
                	for (int j=0; j<P; j++) {
                		int prev = j-1;
                        if (prev <0) prev = P-1;
                        AtomHydrogen jAtom = (AtomHydrogen)atoms.getAtom(j);
                        AtomHydrogen jPrev = (AtomHydrogen)atoms.getAtom(prev);
                        uActAvg3 += -v(jAtom.getBondLength())/(t*P);
                        if (fixedOrientation) {
                            uActAvg3 += 2*Math.log(newBL[j])/P - kHarmonic*dist(jAtom,jPrev); // to be used when orientation moves are turned off
                        }
                        else {
                            uActAvg3 += 2*Math.log(newBL[j]) - kHarmonic*dist(jAtom,jPrev); // to be used when orientation moves are turned on
                        }	
                	}
                	  
                }
                
//                if (i == 1 && moveCount > 1000 && box.getIndex() == 1) System.out.println();
                
            }
            else{ //verify formula for eta
            	for (int j=0; j<P; j++) {// mode = j, bondlength = k                                
                    sigma[j] = Math.sqrt(1/B[j]);
                    etaOld[j] = 0.0;
                    for (int k=0; k<P; k++) {
                        etaOld[j] += A[j][k]*(prevBondLength[i][k] - r0);
                    }
                    uGenOld += -0.5*etaOld[j]*etaOld[j]*B[j];                    
                    etaNew[j] = random.nextGaussian()*sigma[j]; //(j == 0 ? random.nextGaussian()*sigma[j]:etaOld[j]);
//                    if (iRand == 0){
//                    	if (j == 0) {
//                    		etaNew[j] = etaOld[j];
//                    	}
//                    	else {
//                    		etaNew[j] = random.nextGaussian()*sigma[j];
//                    	}
//                    }
//                    else {
//                    	if (j == 0) {
//                    		etaNew[j] = random.nextGaussian()*sigma[j];
//                    	}
//                    	else {
//                    		etaNew[j] = etaOld[j]; 
//                    	}
//                    }
//                    System.out.println(etaNew[j]);
                    uGenNew += -0.5*etaNew[j]*etaNew[j]*B[j];
//                    double newBL = BohrRadius.UNIT.fromSim(((AtomHydrogen)atoms.getAtom(j)).getBondLength());
//                    System.out.println("j ="+j+" oldBL ="+oldBL+" nextBL ="+nextBL+" newBL ="+newBL);
        
//                    System.exit(0);
                    int prev = j-1;
                    if (prev < 0) prev = P-1;

                    AtomHydrogen jAtom = (AtomHydrogen)atoms.getAtom(j);
                    AtomHydrogen jPrev = (AtomHydrogen)atoms.getAtom(prev);

                    if (fixedOrientation) {
                        uA2Old += 2*Math.log(prevBondLength[i][j])/P - kHarmonic*dist(jAtom,jPrev); // to be used when orientation moves are turned off
                    }
                    else {
                        uA2Old += 2*Math.log(prevBondLength[i][j]) - kHarmonic*dist(jAtom,jPrev); // to be used when orientation moves are turned on
                    }
                    uActOld2 += 2*Math.log(prevBondLength[i][j])/P;
                    uActOld3 += - kHarmonic*dist(jAtom,jPrev);

                }
            	
                double[] newBL = new double[P];
                for (int j=0; j<P; j++) {
                    newBL[j] = r0;
//                    System.out.println(i+" "+ j+" "+newBL[j]);

                    for (int k=0; k<P; k++) {//mode = k, bondlength = j
                        newBL[j] += A[k][j]*etaNew[k];
//                        if (i+j==0) System.out.println(A[j][k]+" "+etaNew[k]);
                    }               
//                    h1.addValue(newBL[j]);
//                    System.out.println(i+" "+ j+" "+newBL[j]);
                    ((AtomHydrogen)atoms.getAtom(j)).setBondLength(newBL[j]);
                }
                
                for (int j=0; j<P; j++) {
                    int prev = j-1;
                    if (prev <0) prev = P-1;
                    AtomHydrogen jAtom = (AtomHydrogen)atoms.getAtom(j);
                    AtomHydrogen jPrev = (AtomHydrogen)atoms.getAtom(prev);
                    if (fixedOrientation) {
                        uA2New += 2*Math.log(newBL[j])/P - kHarmonic*dist(jAtom,jPrev); // to be used when orientation moves are turned off
                    }
                    else {
                        uA2New += 2*Math.log(newBL[j]) - kHarmonic*dist(jAtom,jPrev); // to be used when orientation moves are turned on
                    }
                    uActNew2 += 2*Math.log(newBL[j])/P;
                    uActNew3 += - kHarmonic*dist(jAtom,jPrev); 
                }            	
            }
            eAvg += etaOld[3]*etaOld[3];
        }        
        
        double uA1New = mpe.getDataAsScalar();// - nMolecules*P*v(r0);
//        System.out.println("fancy: "+uActOld2+" "+uActNew2+" "+uActOld3+" "+uActNew3+" "+uA1Old/(t*P)+" "+uA1New/(t*P));
        double uActOld = -uA1Old/(t*P) + uA2Old;        
        double uActNew = -uA1New/(t*P) + uA2New;
        pRatio = Math.exp(-uActOld + uActNew - uGenNew + uGenOld);
        if (Double.isInfinite(uActNew) || Double.isInfinite(uActOld) || Double.isInfinite(uGenNew) || Double.isInfinite(uGenOld) || Double.isInfinite(pRatio)) throw new RuntimeException("uActOld ="+uActOld+" uActNew =" +uActNew+ " uGenOld ="+uGenOld+ " uGenNew ="+uGenNew+" pRatio ="+pRatio+" box = "+box.getIndex()+" moveCount = " +moveCount);
        if (Double.isNaN(uActNew) || Double.isNaN(uActOld) || Double.isNaN(uGenNew) || Double.isNaN(uGenOld) || Double.isNaN(pRatio)) throw new RuntimeException("uActOld ="+uActOld+" uActNew =" +uActNew+ " uGenOld ="+uGenOld+ " uGenNew ="+uGenNew+" pRatio ="+pRatio+" box = "+box.getIndex() +" moveCoun = "+moveCount);
//        if (box.getIndex() == 1 && moveCount > 1000) {
//        	System.out.println("uA1Old = "+uA1Old+" uA2Old ="+uA2Old+" uActOld = "+uActOld);
//        	System.out.println("uA1New = "+uA1New+" uA2New ="+uA2New+" uActNew = "+uActNew);
//        }
//        System.out.println(pRatio);
//        double pActNew = Math.exp(uActNew);        
//        double avgBL = 0;
//        for (int j=0; j<P; j++) {
//            IAtomList atoms = box.getMoleculeList().getMolecule(0).getChildList();            
//            avgBL += ((AtomHydrogen)atoms.getAtom(j)).getBondLength()/P;
//        }
//        double pClassical = avgBL*avgBL*Math.exp(-v(avgBL)/t);
//        System.out.println("Classical = "+pClassical+", PI = "+uActNew+" "+pActNew);
//        if (((AtomHydrogen)box.getMoleculeList().getMolecule(0).getChildList().getAtom(0)).getBondLength() < 0.73 && ((AtomHydrogen)box.getMoleculeList().getMolecule(1).getChildList().getAtom(0)).getBondLength() < 0.73) {
//            System.out.println("uActOld ="+uActOld+" uActNew =" +uActNew+ " uGenNew ="+(2*P*u0+uGenNew)+ " uGenOld ="+(uGenOld+u0*P*2)+" pRatio ="+pRatio);  
//        }      

//        if (box.getIndex() == 0 && moveCount > 4978) {
//        	System.out.println(moveCount+" Old : "+"uActOld ="+uActOld+" uActNew =" +uActNew+ " uGenOld ="+uGenOld+ " uGenNew ="+uGenNew+" pRatio ="+pRatio);
//        	System.out.println("uA1Old = "+uA1Old/(t*P)+" uA2Old ="+uA2Old+" uActOld = "+uActOld);
//        	System.out.println("uA1New = "+uA1New/(t*P)+" uA2New ="+uA2New+" uActNew = "+uActNew);
//            System.exit(1);
//        }
//        if (random.nextDouble() < 0.00001) {
//            double[] y = h1.getHistogram();
//            double[] x = h1.xValues();
//            for (int i=0; i<x.length; i++) {
//                System.out.println(x[i]+" "+y[i]);
//            }
//        }
        ((BoxCluster)box).trialNotify();
        
        moveCount++;
        return true;
    }
    protected double dist(AtomHydrogen a0, AtomHydrogen a1) {
        dr.Ev1Mv2(a0.getPosition(),a1.getPosition());
        dr1.Ea1Tv1(a0.getBondLength()/2, a0.getOrientation().getDirection());
        dr1.PEa1Tv1(-a1.getBondLength()/2, a1.getOrientation().getDirection());
        double r2 = dr.Mv1Squared(dr1);
        dr1.TE(-1);
        r2 += dr.Mv1Squared(dr1);
        return r2;
    }
    protected boolean fixedOrientation = false;
    public void setFixedOrientation(boolean a) {
        fixedOrientation = a;
    }
    protected int getIndex(double a, double [] b, int i) {
    	
    	int start = 0;
    	int end = nBins;
    	
    	
//    	System.out.println("moveCount = "+moveCount+" box = "+box.getIndex());
    	
    	
    	do {
    		int mid = (start+end)/2;
    		if (a < b[mid]) {
    			end = mid;
    		}
    		else {
    			start = mid;
    		}
    	} while (end-start>1);
    	if (end < 0 || end == nBins) throw new RuntimeException("you shall not pass "+end);
    	return end;
    	
    }
    
    protected void adjustHistogram(HistogramSimple h) {
    	long [] count = h.getBinCounts();
    	boolean flag1 = false;
    	boolean flag2 = false;
    	boolean flag3 = false;
    	boolean flag4 = false;
    	int right = h.getNBins()/2 + 1;
    	int left = h.getNBins()/2 - 1;
//    	if (box.getIndex() == 0) {
//    		System.out.println("hi");
//    	}
    	int iCount = 0;
    	double deltaX = (h.getXRange().maximum() - h.getXRange().minimum())/h.getNBins();
    	double sigmaC = h.getCount();
    	double value = 10.0/(deltaX*sigmaC);
    	int nAdjusted = 20;
    	while (!flag1 && !flag2) {
    		if (count[right] == 0) {
    			hValues[right] = value;
    			iCount++;
    			if (iCount == nAdjusted) {
    				flag1 = true;
    				break;
    			}    			
    		}
    		if (count[right] > 0 && count[right] < 10) {    			
//    			count[right] = 10;
    			hValues[right] = value;    		
    		}    		    		
    		right++;
    		if (right >= h.getNBins()) {
				flag2 = true;
				break;
			}
    	}
    	iCount = 0;
    	while (!flag3 && !flag4) {
    		if (count[left] == 0) {
    			hValues[left] = value;
    			iCount++;
    			if (iCount == nAdjusted) {
    				flag3 = true;
    				break;
    			}
    		}
    		if (count[left] > 0 && count[left] < 10) {    			
//    			count[left] = 10;
    			hValues[left] = value;    			
    		}    		    		
    		left--;
    		if (left < 0) {
				flag4 = true;
				break;
			}
    	}
//    	while (!flag1 && !flag2) {
//    		if (right == h.getNBins()) {
//				flag2 = true;
//				break;
//			}
//    		if (count[right] > 0 && count[right] < 10) {    			
////    			count[right] = 10;
//    			hValues[right] = value;
//    		}
//    		if (count[right] == 0) {
//    			right++;
//    			if (right == h.getNBins()) {
//    				flag2 = true;
//    				break;
//    			}
//    			if (count[right] == 0) {
//    				right++;
//    				if (right == h.getNBins()) {
//        				flag2 = true;
//        				break;
//        			}
//    				if (count[right] == 0) {    					
////    					count[right] = 10;
//    					hValues[right] = value;
//    					right --;
////    					count[right] = 10;
//    					hValues[right] = value;
//    					right --;
////    					count[right] = 10;
//    					hValues[right] = value;
//    					flag1 = true;
//    				}
//    			}
//    		}
//    		right++;
//    	}
//    	while (!flag3 && !flag4) {
//    		if (left < 0) {
//				flag3 = true;
//				break;
//			}
//    		if (count[left] > 0 && count[left] < 10) {    			
////    			count[left] = 10;
//    			hValues[left] = value;
//    		}
//    		if (count[left] == 0) {
//    			left--;
//    			if (left < 0) {
//    				flag3 = true;
//    				break;
//    			}
//    			if (count[left] == 0) {
//    				left--;
//    				if (left < 0) {
//        				flag3 = true;
//        				break;
//        			}
//    				if (count[left] == 0) {    					
////    					count[left] = 10;
//    					hValues[left] = value;
//    					left ++;
////    					count[left] = 10;
//    					hValues[left] = value;
//    					left ++;
////    					count[left] = 10;
//    					hValues[left] = value;
//    					flag4 = true;
//    				}
//    			}
//    		}
//    		left--;
//    	}
    }
    
    public void setStiffness(IAtomType a) {        
        double lambda = Constants.PLANCK_H/(Math.sqrt(2*Math.PI*a.getMass()*t));
        kHarmonic = Math.PI*P/(lambda*lambda);        
//        System.out.println("Stiffness = "+kHarmonic);
        double[][] m = new double[P][P];
        int maxIter = 100000;
        double tol = 1E-15;
        double x0 = BohrRadius.UNIT.toSim(1.0);
        double xNew = 0;        
        boolean done = false;
        for (int i=0; i<maxIter && !done; i++) {
            if (d2vdr2(x0) == 0 || d2vdr2(x0) != d2vdr2(x0)) throw new RuntimeException("x0 = "+x0+" d2vdr2 = "+d2vdr2(x0));
            xNew = x0 - dvdr(x0)/d2vdr2(x0);
            if (xNew != xNew) throw new RuntimeException("xNew =" + xNew+ " x0 = "+x0+" dvdr = "+dvdr(x0)+" d2vdr2 ="+d2vdr2(x0));
            if (xNew < 0) xNew = x0;
//            System.out.println("iter = "+i+" ,x0 = "+BohrRadius.UNIT.fromSim(x0)+" ,xNew = "+BohrRadius.UNIT.fromSim(xNew)+ " ,error = "+(x0-xNew)*(x0-xNew));
//            System.out.println("x0 = "+x0+" dvdr ="+dvdr(x0)+" v ="+v(x0)+" d2v = "+d2vdr2(x0));            
            if ((x0 - xNew)*(x0 - xNew) < tol*tol) {
                r0 = xNew;
                done = true;
//                System.out.println("r0 = "+(r0/10));//BohrRadius.UNIT.toSim(1.414772533)
            }
            x0 = xNew;
        }
        kIntra = -d2vdr2(r0);
        u0 = v(r0);
        for (int i=0; i<P; i++) {
            int next = i+1;
            if (next == P) next = 0;            
            int prev = i-1;
            if (prev < 0) prev = P-1;
            for (int j=0; j<P; j++) {
                m[i][j] = 0.0;
                if (j == next || j == prev) m[i][j] = -kHarmonic;
                if (j == i) m[i][j] = kIntra + 2*kHarmonic;
            }            
        }
        Jama.Matrix M = new Jama.Matrix(m); 
        EigenvalueDecomposition E = M.eig();
        A = E.getV().transpose().getArray();
        B = E.getRealEigenvalues();                
        
//        System.out.println("Old: "+r0+" "+kIntra);
//        for (int i=0; i<P; i++) {
//        	double mag = 0;
//        	for (int j=0; j<P; j++) {
//        		System.out.print(A[i][j]+" ");
//        		mag += A[i][j]*A[i][j];
//        	}        	
//        	System.out.println(" "+Math.sqrt(mag));
//        }
//        System.out.println(dvdr(r0));
//        System.exit(1);
        
    }
    public double getStiffness() {
        return kHarmonic;
    }

    
    public double getA() {
        double wNew = ((BoxCluster)box).getSampleCluster().value((BoxCluster)box);    
        return wNew*pRatio/wOld;
    }

    
    public double getB() {
        return 0.00;
    }

    
    public void acceptNotify() {
//        if (box.getIndex() == 0 && moveCount > 1000) {-
//        	System.out.println("Old accepted");
//        }
    	if (foo) {
//    		System.out.println("accepted");
    		foo = false;
    	}
//    	if (flagAdjusted) System.out.print("moveCount = "+moveCount+" accepted \n");
        ((BoxCluster)box).acceptNotify();        
        
    }
    protected boolean foo;
    
    public void rejectNotify() {
//    	if (box.getIndex() == 1 && moveCount > 1000) System.out.println("Old rejected");
		foo = false;
        int nMolecules = box.getMoleculeList().getMoleculeCount();        
        for (int i=0; i<nMolecules; i++) {
            IAtomList atoms = box.getMoleculeList().getMolecule(i).getChildList();
            for (int j=0; j<P; j++) {                
                ((AtomHydrogen)atoms.getAtom(j)).setBondLength(prevBondLength[i][j]);
            }
        }
        ((BoxCluster)box).rejectNotify();
    }

}
