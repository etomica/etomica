/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

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
import etomica.potential.P1HydrogenMielke;
import etomica.potential.P1HydrogenMielke.P1HydrogenMielkeAtomic;
import etomica.potential.P1IntraMolecular;
import etomica.potential.Potential2SoftSpherical;
import etomica.space.ISpace;
import etomica.space.Space;
import etomica.units.BohrRadius;
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
    protected final IVectorMutable dr,dr1,dr2,dummy;
    protected double t,r0;
    protected double[] B;
    protected double[][] A;
    public double uActAvg1 = 0, count1 = 0, uActAvg2 = 0, count2 = 0, uActAvg3 = 0, count3 = 0;
    public double [] e2Avg,eAvg;
    protected P1IntraMolecular p11;
    protected int P;
    protected int[] leftAdjusted,rightAdjusted;
    protected double u0 = -1;    
    public HistogramSimple h1,hEta23;
    public HistogramSimple h2;
    public boolean deleteMode = true, doHist = false, flagAdjusted, doDebug = false;
    public int moveCount;
    protected int nBins = 100;
    protected double xMin = 0.35,xMax = 1.25, uMin = 0.0;
    protected double deltaX = (xMax - xMin)/nBins;
    public double [] hValues;
    protected double[] pCummulative;
    protected double pTotal = 0.0;    
    public double [] sigma;
    protected boolean printOnce = true;


    public MCMoveChangeBondLength(IPotentialMaster potentialMaster, IRandom random, ISpace space, double temperature) {
        super(potentialMaster);    
        leafIterator = new AtomIteratorLeafAtoms();
        this.random = random;
        mpe = new MeterPotentialEnergy(potentialMaster);
        dr = space.makeVector();
        dr1 = space.makeVector();
        dr2 = space.makeVector();
        dummy = space.makeVector();
        setTemperature(temperature);
        h1 = new HistogramSimple(nBins, new DoubleRange(xMin, xMax));
        h2 = new HistogramSimple(1000, new DoubleRange(xMin, xMax));
        hEta23 = new HistogramSimple(1000, new DoubleRange(0, 2));

        hValues = h1.getHistogram();
        pCummulative = new double[nBins];
        moveCount = 0;

        leftAdjusted = new int[2];
        rightAdjusted = new int[2];

    }

    public void setTemperature(double x){
        t = x;
    }  

    public void setBox(IBox p) {
        super.setBox(p);
        int nMolecules = box.getMoleculeList().getMoleculeCount();
        P = box.getMoleculeList().getMolecule(0).getChildList().getAtomCount();
        prevBondLength = new double[nMolecules][P];
        e2Avg = new double [P];
        eAvg = new double [P];
        sigma = new double[P];
        leafIterator.setBox(p);
    }

    public AtomIterator affectedAtoms() {        
        return leafIterator;
    }

    public double energyChange() {    
        return 0;
    }    
    public boolean doTrial() {
        if (doDebug && moveCount == 0) System.out.println("Warning: doDebug is on");
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

        boolean flag = moveCount >= 1000;
        
        flagAdjusted = false;

        double avgBL = 0;
        if (deleteMode || doHist) {
            for (int i=0; i<nMolecules; i++) {
                avgBL = 0;
                IAtomList atoms = box.getMoleculeList().getMolecule(i).getChildList();
                for (int j=0; j<P; j++) {
                    avgBL += ((AtomHydrogen)atoms.getAtom(j)).getBondLength();
                }
                avgBL /= P;
                h1.addValue(avgBL);
                if (flag) h2.addValue(avgBL);
            }
        }

        if (deleteMode && flag) {
            if (moveCount %100 == 0) {
                hValues = h1.getHistogram();
                double [] xValues = h1.xValues();
                adjustHistogram(h1);
//                if (printOnce && box.getIndex() == 0) {
//                    for (int i=0; i< xValues.length; i++) {
//                        System.out.println(xValues[i] + " " + hValues[i]);
//                    }
//                    printOnce = false;
//                }
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
        if (u0 == -1) throw new RuntimeException("u0 has not been set!");

        uA1Old = mpe.getDataAsScalar();// - nMolecules*P*u0;
        for (int i=0; i<nMolecules; i++) {

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
                for (int j=0; j<P && !done; j++) {// j = bond length; k = mode
                    double mag = 0;
                    for (int k=0; k<P; k++) {
                        mag += A[j][k];
                    }
                    double dummy1 = (mag - Math.sqrt(P))*(mag - Math.sqrt(P));
                    double dummy2 = (mag + Math.sqrt(P))*(mag + Math.sqrt(P));

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
                    eAvg[j] = 0.0;
                    e2Avg[j] = 0.0;
                    etaOld[j] = 0.0;
                    etaNew[j] = 0.0;
                    for (int k=0; k<P; k++) {//mode = j, bondlength = k
                        if (j != choice) {
                            etaOld[j] += A[j][k]*(prevBondLength[i][k] - oldAvgBL);                    		
                            eAvg[j] += etaOld[j];
                            e2Avg[j] += etaOld[j]*etaOld[j];                    		
                        }
                    }                    


                    int prev = j-1;
                    if (prev < 0) prev = P-1;

                    AtomHydrogen jAtom = (AtomHydrogen)atoms.getAtom(j);
                    AtomHydrogen jPrev = (AtomHydrogen)atoms.getAtom(prev);


                    if (fixedOrientation) {
                        uA2Old += -2*Math.log(prevBondLength[i][j])/P + kHarmonic*dist(jAtom,jPrev); // to be used when orientation moves are turned off
                    }
                    else {
                        uA2Old += -2*Math.log(prevBondLength[i][j]) + kHarmonic*dist(jAtom,jPrev); // to be used when orientation moves are turned on
                    }
                    uActOld2 += -2*Math.log(prevBondLength[i][j])/P;
                    uActOld3 += + kHarmonic*dist(jAtom,jPrev);
                }
                double hTemp = Math.sqrt(etaOld[1]*etaOld[1] + etaOld[2]*etaOld[2]);
                hEta23.addValue(hTemp);
                for (int k=0; k<P; k++) {
                    if (k != choice ) uGenOld += 0.5*etaOld[k]*etaOld[k]*B[k];
                }

                int oldChoice = (int)((oldAvgBL - xMin)/deltaX);
                if (hValues[oldChoice] == 0 || pTotal == 0) throw new RuntimeException("oldChoice = "+oldChoice+" hValues = "+hValues[oldChoice]+" "+" pTotal = "+pTotal+" moveCount = "+moveCount+ " molecule = "+i+" bl = "+oldAvgBL+" box = "+box.getIndex());
                uGenOld += -Math.log(hValues[oldChoice]/pTotal);

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
                            uGenNew += 0.5*etaNew[k]*etaNew[k]*B[k];                		
                        }
                    }
                    double dummy = random.nextDouble();        		
                    done = false;
                    if (dummy <= pCummulative[0]) {
                        newChoice = 0;
                        done = true;
                    }

                    if (!done || newChoice < 0) newChoice = getIndex(dummy,pCummulative,1);

                    if (newChoice > 0) done = true;

                    double y1 = xMin + newChoice*deltaX;
                    double y2 = xMin + (newChoice+1)*deltaX;
                    double x1 = 0;
                    if (newChoice > 0) x1 = pCummulative[newChoice-1];                	
                    double x2 = pCummulative[newChoice];
                    if ((newChoice >= leftAdjusted[0] && newChoice <= leftAdjusted[1]) || (newChoice >= rightAdjusted[0] && newChoice <= rightAdjusted[1])) flagAdjusted = true;
                    rStar = y1 + (dummy - x1)*(y2 - y1)/(x2 - x1);

                    if (hValues[newChoice] == 0 || pTotal == 0) throw new RuntimeException("newChoice = "+newChoice+" hValues = "+hValues[newChoice]+" "+" pTotal = "+pTotal+" moveCount = "+moveCount+ " molecule = "+i+" box = "+box.getIndex());
                    uGenNew += -Math.log(hValues[newChoice]/pTotal);
                    flag1 = true;                	
                    for (int j=0; j<P; j++) {
                        newBL[j] = rStar;
                        for (int k=0; k<P; k++) {//mode = k, bondlength = j
                            if (k != choice) newBL[j] += A[k][j]*etaNew[k];
                        }                		
                        newAvgBL += newBL[j];
                        if (newBL[j] < 0) {                            
                            flag1 = false;
                            break;
                        }
                        ((AtomHydrogen)atoms.getAtom(j)).setBondLength(newBL[j]);
                    }
                    newAvgBL /= P;
                } while(!flag1);

                for (int j=0; j<P; j++) {
                    int prev = j-1;
                    if (prev <0) prev = P-1;
                    AtomHydrogen jAtom = (AtomHydrogen)atoms.getAtom(j);
                    AtomHydrogen jPrev = (AtomHydrogen)atoms.getAtom(prev);
                    if (fixedOrientation) {
                        uA2New += -2*Math.log(newBL[j])/P + kHarmonic*dist(jAtom,jPrev); // to be used when orientation moves are turned off
                    }
                    else {
                        uA2New += -2*Math.log(newBL[j]) + kHarmonic*dist(jAtom,jPrev); // to be used when orientation moves are turned on
                    }                    

                }
            }
            else{ //verify formula for eta

                for (int j=0; j<P; j++) {// mode = j, bondlength = k                                
                    sigma[j] = Math.sqrt(1/B[j]);
                    etaOld[j] = 0.0;
                    for (int k=0; k<P; k++) {
                        etaOld[j] += A[j][k]*(prevBondLength[i][k] - r0);
                    }

                    uGenOld += 0.5*etaOld[j]*etaOld[j]*B[j];

                    etaNew[j] = random.nextGaussian()*sigma[j]; //(j == 0 ? random.nextGaussian()*sigma[j]:etaOld[j]);

                    uGenNew += 0.5*etaNew[j]*etaNew[j]*B[j];


                    int prev = j-1;
                    if (prev < 0) prev = P-1;

                    AtomHydrogen jAtom = (AtomHydrogen)atoms.getAtom(j);
                    AtomHydrogen jPrev = (AtomHydrogen)atoms.getAtom(prev);

                    if (fixedOrientation) {                    	
                        uA2Old += -2*Math.log(prevBondLength[i][j])/P + kHarmonic*dist(jAtom,jPrev); // to be used when orientation moves are turned off                        
                    }
                    else {
                        uA2Old += -2*Math.log(prevBondLength[i][j]) + kHarmonic*dist(jAtom,jPrev); // to be used when orientation moves are turned on
                    }
                    uActOld2 += -2*Math.log(prevBondLength[i][j])/P;
                    uActOld3 += kHarmonic*dist(jAtom,jPrev);
                }
                double hTemp = Math.sqrt(etaOld[1]*etaOld[1] + etaOld[2]*etaOld[2]);
                hEta23.addValue(hTemp);

                double[] newBL = new double[P];
                for (int j=0; j<P; j++) {
                    newBL[j] = r0;                    

                    for (int k=0; k<P; k++) {//mode = k, bondlength = j
                        newBL[j] += A[k][j]*etaNew[k];
                    }

                    ((AtomHydrogen)atoms.getAtom(j)).setBondLength(newBL[j]);
                }

                for (int j=0; j<P; j++) {

                    int prev = j-1;
                    if (prev <0) prev = P-1;
                    AtomHydrogen jAtom = (AtomHydrogen)atoms.getAtom(j);
                    AtomHydrogen jPrev = (AtomHydrogen)atoms.getAtom(prev);
                    if (fixedOrientation) {
                        double dy1 = -2*Math.log(newBL[j])/P;
                        double dy2 = kHarmonic*dist(jAtom,jPrev);
                        uA2New += dy1 + dy2;                        
                    }
                    else {
                        uA2New += -2*Math.log(newBL[j]) + kHarmonic*dist(jAtom,jPrev); // to be used when orientation moves are turned on
                    }
                    uActNew2 += -2*Math.log(newBL[j])/P;
                    uActNew3 += kHarmonic*dist(jAtom,jPrev); 
                }                
            }            
        }        

        double uA1New = mpe.getDataAsScalar();// - nMolecules*P*u0;

        ((BoxCluster)box).trialNotify();
        
        double uActOld = uA2Old;
        double uActNew = uA2New;      
        if (!doExchange){
            uActOld += uA1Old/(t*P);
            uActNew += uA1New/(t*P);
        }

        pRatio = Math.exp(uActOld - uActNew + uGenNew - uGenOld);


        if (Double.isInfinite(uActNew) || Double.isInfinite(uActOld) || Double.isInfinite(uGenNew) || Double.isInfinite(uGenOld) || Double.isInfinite(pRatio)) {
            throw new RuntimeException("uActOld ="+uActOld+" uActNew =" +uActNew+ " uGenOld ="+uGenOld+ " uGenNew ="+uGenNew+" pRatio ="+pRatio+" box = "+box.getIndex()+" moveCount = " +moveCount);
        }
        if (Double.isNaN(uActNew) || Double.isNaN(uActOld) || Double.isNaN(uGenNew) || Double.isNaN(uGenOld) || Double.isNaN(pRatio)) throw new RuntimeException("uActOld ="+uActOld+" uActNew =" +uActNew+ " uGenOld ="+uGenOld+ " uGenNew ="+uGenNew+" pRatio ="+pRatio+" box = "+box.getIndex() +" moveCoun = "+moveCount);
        moveCount++;

        return true;
    }
    protected boolean doExchange = false;
    public void setDoExchange(boolean b) {
        doExchange = b;
    }
    protected boolean fixedOrientation = false;
    public void setFixedOrientation(boolean a) {
        fixedOrientation = a;
    }
    protected double dist(AtomHydrogen a0, AtomHydrogen a1) {
        dummy.E(a0.getOrientation().getDirection());
        if (a0.getIndex() == 0 && doExchange) {
            dummy.TE(-1);
        }
        dr.Ev1Mv2(a0.getPosition(),a1.getPosition());
        dr1.Ea1Tv1(a0.getBondLength()/2, dummy);
        dr1.PEa1Tv1(-a1.getBondLength()/2, a1.getOrientation().getDirection());
        double r2 = dr.Mv1Squared(dr1);
        dr1.TE(-1);
        r2 += dr.Mv1Squared(dr1);
        dr.normalize();
        return r2;
    }

    protected int getIndex(double a, double [] b, int i) {

        int start = 0;
        int end = nBins;    	

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
        double [] dummyXValues = h.xValues();
        int right = 0;        
        double hMax = Double.NEGATIVE_INFINITY;
        
        for (int i=0; i< dummyXValues.length; i++) {
            if (h.getHistogram()[i] > hMax) {
                hMax = h.getHistogram()[i];
                right = i;
            }
        }
        int left = right - 1;
        
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
                hValues[left] = value;    			
            }    		    		
            left--;
            if (left < 0) {
                flag4 = true;
                break;
            }
        }
    }

    public void setStiffness(IAtomType a, P1IntraMolecular p1) {
        p11 = p1;
        double lambda = Constants.PLANCK_H/(Math.sqrt(2*Math.PI*a.getMass()*t));
        kHarmonic = Math.PI*P/(lambda*lambda);
        double cT = 1.0;
        if (doExchange) cT = Math.cos(Math.PI/P);        
        double[][] m = new double[P][P];
        int maxIter = 100000;
        double tol = 1E-15;
        double x0 = BohrRadius.UNIT.toSim(1.0);
        double xNew = 0;        
        boolean done = false;
        for (int i=0; i<maxIter && !done; i++) {
            cT = 1 - (P - 1)/(P*kHarmonic*x0*x0);
            if (Double.isNaN(cT)) throw new RuntimeException("cT = "+cT);
            if (cT < -1.0) {
                cT = -1.0;
            }
            if (cT >  1.0) {
                cT = 1.0;
            }
            double d2vdr2 = 2*kHarmonic + 2/(x0*x0) + p1.d2u(x0)/(x0*x0*t*P);
            if (x0 < 0) x0 = BohrRadius.UNIT.toSim(random.nextDouble());
            if (d2vdr2 == 0 || d2vdr2 != d2vdr2) throw new RuntimeException("x0 = "+x0+" d2vdr2 = "+d2vdr2);
            double dvdr = 2*kHarmonic*x0*(1 - cT) - 2/x0 + p1.du(x0)/(x0*t*P);            
            xNew = x0 - dvdr/d2vdr2;
            if (xNew != xNew) throw new RuntimeException("xNew =" + xNew+ " x0 = "+x0+" dvdr = "+dvdr+" d2vdr2 ="+d2vdr2);
            
            if ((x0 - xNew)*(x0 - xNew) < tol*tol) {
                r0 = xNew;
                done = true;
            }
            x0 = xNew;
        }
        if (!done) throw new RuntimeException("Newton Raphson method has not converged!!!");
        cT = 1 - (P - 1)/(P*kHarmonic*r0*r0);
        u0 = p1.u(r0);
        double diagTerm = 2*kHarmonic + 2/(r0*r0) + p1.d2u(r0)/(r0*r0*t*P);
        //          System.out.println(r0+" "+kHarmonic+" "+diagTerm);
        //          System.exit(1);
        //          if (box.getIndex() == 0) {
        //        	double la = kHarmonic*r0*r0*0.02*0.02;
        //        	double lb = -2*Math.log(1.02)/P;
        //        	double lc = (p1.u(1.02*r0) - u0)/(t*P);
        //        	double lhs = la + lb +lc;
        //        	double da = 0;
        //        	double d2a = 2*kHarmonic;
        //        	double db = -2/(P*r0);
        //        	double d2b = 2/(P*r0*r0);
        //        	double dc = p1.du(r0)/(t*P*r0);
        //        	double d2c = p1.d2u(r0)/(t*P*r0*r0);
        //        	double ra = 0.02*r0*da + 0.5*0.02*0.02*r0*r0*d2a;
        //        	double rb = 0.02*r0*db + 0.5*0.02*0.02*r0*r0*d2b;
        //        	double rc = 0.02*r0*dc + 0.5*0.02*0.02*r0*r0*d2c;        
        //        	double rhs = ra + rb + rc;
        //
        //        	System.out.println("la = "+la+" ra = "+ra+ " difference = "+(la - ra));
        //        	System.out.println("lb = "+lb+" rb = "+rb+ " difference = "+(lb - rb));
        //        	System.out.println("lc = "+lc+" rc = "+rc+ " difference = "+(lc - rc));
        //
        //        	System.out.println("dS = "+da+" + "+db+" + "+dc+" = "+(da+db+dc));
        //        	System.out.println("lhs = "+lhs+" rhs = "+rhs+" difference = "+(lhs - rhs));
        //          }
        for (int i=0; i<P; i++) {
            int next = i+1;
            if (next == P) next = 0;            
            int prev = i-1;
            if (prev < 0) prev = P-1;
            for (int j=0; j<P; j++) {
                m[i][j] = 0.0;
                if (j == next || j == prev) m[i][j] = -kHarmonic*cT;                
                if (j == i) m[i][j] = diagTerm ;
            }            
        }
        Jama.Matrix M = new Jama.Matrix(m); 
        EigenvalueDecomposition E = M.eig();
        A = E.getV().transpose().getArray();
        B = E.getRealEigenvalues();
    }

    public double getStiffness() {
        return kHarmonic;
    }

    public double getA() {
        int nMolecules = box.getMoleculeList().getMoleculeCount();        
        for (int i=0; i<nMolecules; i++) {
            IAtomList atoms = box.getMoleculeList().getMolecule(i).getChildList();
            for (int j=0; j<P; j++) {                
                if (((AtomHydrogen)atoms.getAtom(j)).getBondLength() < 0 ) return 0;
            }
        }
        double wNew = ((BoxCluster)box).getSampleCluster().value((BoxCluster)box);

        return wNew*pRatio/wOld;
    }


    public double getB() {
        return 0.00;
    }

    public void acceptNotify() {
        ((BoxCluster)box).acceptNotify();
    }

    public void rejectNotify() {
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
