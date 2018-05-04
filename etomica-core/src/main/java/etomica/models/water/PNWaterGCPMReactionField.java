package etomica.models.water;

import Jama.Matrix;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.math.SpecialFunctions;
import etomica.molecule.IMoleculeList;
import etomica.molecule.MoleculePair;
import etomica.potential.PotentialMolecular;
import etomica.potential.PotentialPolarizable;
import etomica.space.Boundary;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.space.Vector;
import etomica.units.*;
import etomica.util.Arrays;

/**
 * GCPM Water potential class including Reaction Field Method.
 * The polarization energy is solved for using matrix
 * inversion rather than iteration, so this class may not be suitable for
 * large systems.
 * 
 * @author Hye Min Kim
 */
public class PNWaterGCPMReactionField extends PotentialMolecular implements PotentialPolarizable {

    public PNWaterGCPMReactionField(Space space) {
	    super(Integer.MAX_VALUE, space);
    	//super(2, space);//ignore many-body interaction
	    pair = new MoleculePair();
        sigma = 3.69;
        epsilon = Kelvin.UNIT.toSim(110);
        gamma = 12.75;
        chargeH = Electron.UNIT.toSim(0.6113);
        chargeM = Electron.UNIT.toSim(-1.2226);
        core = 4.41; //4.41 = 2.1^2; value according to Cummings
        sigmaM = 0.610;
        sigmaH = 0.455;
        sqrtHMsigmas = Math.sqrt(2*(sigmaH*sigmaH+sigmaM*sigmaM));
        massH = 1.00794;
        massO = 15.9995;
        totalMass = 18.01538;
        sqrtPiHMsigmas = Math.sqrt(Math.PI*(sigmaH*sigmaH+sigmaM*sigmaM));
        sqrtPiMMsigmas = Math.sqrt(Math.PI*(2*sigmaM*sigmaM));
        alphaPol = 1.444;
        
        comW1 = space.makeVector();
        comW2 = space.makeVector();
        comWi = space.makeVector();
        comWj = space.makeVector();
        shift = space.makeVector();
        r12Vector = space.makeVector();
        rijVector = space.makeVector();

        work = space.makeVector();
        
        dipoleMoment1 = space.makeVector();
        dipoleMoment2 = space.makeVector();
        iDipoleMoment = space.makeVector();
        jDipoleMoment = space.makeVector();
        myRq = space.makeVector();
        
        Tunit = space.makeTensor();
        Tij = space.makeTensor();

        Eq = new Matrix[0];
        A = new Matrix[0];
	}   

    public double energy(IMoleculeList atoms){
        double volume = box.getBoundary().volume();
        double boxLength = Math.pow(volume, 1.0/3.0);      
    	setCutOffDistance(boxLength*0.49);
    	initRqFactor();
        double sum = 0;
        for (int i=0; i<atoms.getMoleculeCount()-1; i++) {
            pair.atom0 = atoms.getMolecule(i);
            
            IAtomList iLeafAtoms = pair.atom0.getChildList();
            Vector O1r = iLeafAtoms.getAtom(SpeciesWater4P.indexO).getPosition();
            Vector H11r = iLeafAtoms.getAtom(SpeciesWater4P.indexH1).getPosition();
            Vector H12r = iLeafAtoms.getAtom(SpeciesWater4P.indexH2).getPosition();
            Vector M1r = iLeafAtoms.getAtom(SpeciesWater4P.indexM).getPosition();

            comWi.Ea1Tv1(massH, H11r);
            comWi.PEa1Tv1(massO, O1r);
            comWi.PEa1Tv1(massH, H12r);
            comWi.TE(1.0/totalMass);//c.o.m of molecule i
            
            comHi1r = space.makeVector();
            comHi2r = space.makeVector();
            comMir = space.makeVector();
            comHi1r.Ev1Mv2(comWi, H11r);
            comHi2r.Ev1Mv2(comWi, H12r);
            comMir.Ev1Mv2(comWi, M1r);
            iDipoleMoment.Ea1Tv1(chargeH, comHi1r);
            iDipoleMoment.PEa1Tv1(chargeH, comHi2r);
            iDipoleMoment.PEa1Tv1(chargeM, comMir);
            
            //if (i ==0)System.out.println("1 dipole moment "+iDipoleMoment);
            
            //System.out.println("fixd dipole "+Debye.UNIT.fromSim(iDipoleMoment.getX(2)));
            sum -= myRqFactor*0.5*iDipoleMoment.squared();//reaction field contribution
            for (int j=i+1; j<atoms.getMoleculeCount(); j++) {
                pair.atom1 = atoms.getMolecule(j);
                double nonPolE = getNonPolarizationEnergy(pair);               
                sum += nonPolE;
                if (nonPolE !=0){
	            	IAtomList jLeafAtoms = pair.atom1.getChildList();
	            	Vector Ojr = jLeafAtoms.getAtom(SpeciesWater4P.indexO).getPosition();
	            	Vector Hj1r = jLeafAtoms.getAtom(SpeciesWater4P.indexH1).getPosition();
	            	Vector Hj2r = jLeafAtoms.getAtom(SpeciesWater4P.indexH2).getPosition();
	            	Vector Mjr = jLeafAtoms.getAtom(SpeciesWater4P.indexM).getPosition();
	
	                comWj.Ea1Tv1(massH, Hj1r);
	                comWj.PEa1Tv1(massO, Ojr);
	                comWj.PEa1Tv1(massH, Hj2r);
	                comWj.TE(1.0/totalMass);//c.o.m of molecule j
	                
	                comHj1r = space.makeVector();
	                comHj2r = space.makeVector();
	                comMjr = space.makeVector();
	                comHj1r.Ev1Mv2(comWj, Hj1r);
	                comHj2r.Ev1Mv2(comWj, Hj2r);
	                comMjr.Ev1Mv2(comWj, Mjr);
	                jDipoleMoment.Ea1Tv1(chargeH, comHj1r);
	                jDipoleMoment.PEa1Tv1(chargeH, comHj2r);
	                jDipoleMoment.PEa1Tv1(chargeM, comMjr);
	                sum -= myRqFactor*iDipoleMoment.dot(jDipoleMoment);//reaction field contribution
	                if (Double.isInfinite(sum)) {
	                    return sum;
	                }
                }
                
            }
        }
        sum += getPolarizationEnergy(atoms);
        return sum;
    }
    
    /**
     * This returns the pairwise-additive portion of the GCPM potential for a
     * pair of atoms (dispersion + fixed-charge electrostatics)
     */
    public double getNonPolarizationEnergy(IMoleculeList atoms) {
        double volume = box.getBoundary().volume();
        double boxLength = Math.pow(volume, 1.0/3.0);      
    	setCutOffDistance(boxLength*0.49);
    	
        IAtomList water1Atoms = atoms.getMolecule(0).getChildList();
        IAtomList water2Atoms = atoms.getMolecule(1).getChildList();

        Vector O1r = water1Atoms.getAtom(SpeciesWater4P.indexO).getPosition();
        Vector O2r = water2Atoms.getAtom(SpeciesWater4P.indexO).getPosition();
        
        work.Ev1Mv2(O1r, O2r);
        shift.Ea1Tv1(-1,work);
		box.getBoundary().nearestImage(work);
        shift.PE(work);
        final boolean zeroShift = shift.squared() < 0.1;
        
        double r2 = work.squared();
        
        if(r2<=core) {
            return Double.POSITIVE_INFINITY;
        }

        Vector H11r = water1Atoms.getAtom(SpeciesWater4P.indexH1).getPosition();
        Vector H12r = water1Atoms.getAtom(SpeciesWater4P.indexH2).getPosition();
        Vector H21r = water2Atoms.getAtom(SpeciesWater4P.indexH1).getPosition();
        Vector H22r = water2Atoms.getAtom(SpeciesWater4P.indexH2).getPosition();

        Vector M1r = water1Atoms.getAtom(SpeciesWater4P.indexM).getPosition();
        Vector M2r = water2Atoms.getAtom(SpeciesWater4P.indexM).getPosition();
        
        comW1.Ea1Tv1(massH, H11r);
        comW1.PEa1Tv1(massO, O1r);
        comW1.PEa1Tv1(massH, H12r);
        comW1.TE(1.0/totalMass);//c.o.m of molecule 1
        
//        comH11r = space.makeVector();
//        comH12r = space.makeVector();
//        comM1r = space.makeVector();
//        comH11r.Ev1Mv2(comW1, H11r);
//        comH12r.Ev1Mv2(comW1, H12r);
//        comM1r.Ev1Mv2(comW1, M1r);
//        dipoleMoment1.Ea1Tv1(chargeH, comH11r);
//        dipoleMoment1.PEa1Tv1(chargeH, comH12r);
//        dipoleMoment1.PEa1Tv1(chargeM, comM1r);
        
        comW2.Ea1Tv1(massH, H21r);
        comW2.PEa1Tv1(massO, O2r);
        comW2.PEa1Tv1(massH, H22r);
        comW2.TE(1.0/totalMass);//c.o.m of molecule 2
        
//        comH21r = space.makeVector();
//        comH22r = space.makeVector();
//        comM2r = space.makeVector();
//        comH21r.Ev1Mv2(comW2, H21r);
//        comH22r.Ev1Mv2(comW2, H22r);
//        comM2r.Ev1Mv2(comW2, M2r);
//        dipoleMoment2.Ea1Tv1(chargeH, comH21r);
//        dipoleMoment2.PEa1Tv1(chargeH, comH22r);
//        dipoleMoment2.PEa1Tv1(chargeM, comM2r);
        
        rijVector.E(comW1);
        rijVector.ME(comW2);
        rijVector.PE(shift);
        
        double r12 = rijVector.squared();
        if (r12 > cutOffDistance*cutOffDistance){
        	return 0.0;
        }


        double r = Math.sqrt(r2);
        double rOverSigma = r/sigma;
        double sigma2OverR2 = 1/(rOverSigma*rOverSigma);
        double sixOverGamma = 6/gamma;
   
        double sum = epsilon/(1 - sixOverGamma)*(sixOverGamma*Math.exp(gamma*(1 - rOverSigma)) - sigma2OverR2*sigma2OverR2*sigma2OverR2);//exp-6 potential(Udisp)
        if (dodebug){
        	//System.out.println("sum in nonPol "+sum);
        }
        if (zeroShift){//Uqq, no pbc
	        r2 = H11r.Mv1Squared(H21r);
	        sum += chargeH*chargeH/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaH)));

	
	        r2 = H11r.Mv1Squared(H22r);
	        sum += chargeH*chargeH/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaH)));

	
	        r2 = H12r.Mv1Squared(H21r);
	        sum += chargeH*chargeH/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaH)));

	
	        r2 = H12r.Mv1Squared(H22r);
	        sum += chargeH*chargeH/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaH)));

	        
	        r2 = M1r.Mv1Squared(H21r);
	        sum += chargeH*chargeM/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/sqrtHMsigmas));

	
	        r2 = M1r.Mv1Squared(H22r);
	        sum += chargeH*chargeM/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/sqrtHMsigmas));

	
	        r2 = M2r.Mv1Squared(H11r);
	        sum += chargeH*chargeM/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/sqrtHMsigmas));

	
	        r2 = M2r.Mv1Squared(H12r);
	        sum += chargeH*chargeM/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/sqrtHMsigmas));

	
	        r2 = M1r.Mv1Squared(M2r);
	        sum += chargeM*chargeM/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaM)));

	        
        }
        else {//with pbc
        	shift.PE(H11r);
        	r2 = H21r.Mv1Squared(shift);
        	shift.ME(H11r);
	        sum += chargeH*chargeH/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaH)));
	        
	        shift.PE(H11r);
        	r2 = H22r.Mv1Squared(shift);
        	shift.ME(H11r);
	        sum += chargeH*chargeH/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaH)));
	
	        shift.PE(H12r);
        	r2 = H21r.Mv1Squared(shift);
        	shift.ME(H12r);
	        sum += chargeH*chargeH/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaH)));
	
	        shift.PE(H12r);
        	r2 = H22r.Mv1Squared(shift);
        	shift.ME(H12r);
	        sum += chargeH*chargeH/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaH)));
	        
	        shift.PE(M1r);
        	r2 = H21r.Mv1Squared(shift);
        	shift.ME(M1r);
	        sum += chargeH*chargeM/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/sqrtHMsigmas));
	
	        shift.PE(M1r);
        	r2 = H22r.Mv1Squared(shift);
        	shift.ME(M1r);
	        sum += chargeH*chargeM/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/sqrtHMsigmas));
	
	        shift.PE(H11r);
        	r2 = M2r.Mv1Squared(shift);
        	shift.ME(H11r);
	        sum += chargeH*chargeM/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/sqrtHMsigmas));
	
	        shift.PE(H12r);
        	r2 = M2r.Mv1Squared(shift);
        	shift.ME(H12r);
	        sum += chargeH*chargeM/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/sqrtHMsigmas));
	
	        shift.PE(M1r);
        	r2 = M2r.Mv1Squared(shift);
        	shift.ME(M1r);
	        sum += chargeM*chargeM/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaM)));

        	
        }
        return sum;
    }

    /**
     * This returns the polarizable portion of the GCPM potential for any
     * number of atoms.
     */
    public double getPolarizationEnergy(IMoleculeList atoms) {
    	double volume = box.getBoundary().volume();
        double boxLength = Math.pow(volume, 1.0/3.0);  
    	setCutOffDistance(boxLength*0.49);
        
        final int atomCount = atoms.getMoleculeCount();
        if (Eq.length < atomCount+1) {
            Eq = (Matrix[])Arrays.resizeArray(Eq, atomCount+1);
            A = (Matrix[])Arrays.resizeArray(A, atomCount+1);
        }
        
        if (Eq[atomCount] == null) {
            Eq[atomCount] = new Matrix(3*atomCount, 1);
            A[atomCount] = new Matrix(3*atomCount, 3*atomCount);
            
            for (int i=0; i<3*atomCount; i++) {
                A[atomCount].set(i, i, 1);
            }
        }
        final Matrix myEq = Eq[atomCount];
        final Matrix myA = A[atomCount];
        for (int i=0; i<3*atomCount; i++) {
            myEq.set(i, 0, 0);
        }
        
        /*
         * Finding the Electric fields at the center of mass of each molecule, Eqi
         * kmb, 8/7/06
         */

        for (int i=0; i<atoms.getMoleculeCount(); i++) {
            IAtomList iLeafAtoms = atoms.getMolecule(i).getChildList();
            Vector O1r = iLeafAtoms.getAtom(SpeciesWater4P.indexO).getPosition();
            Vector H11r = iLeafAtoms.getAtom(SpeciesWater4P.indexH1).getPosition();
            Vector H12r = iLeafAtoms.getAtom(SpeciesWater4P.indexH2).getPosition();
            Vector M1r = iLeafAtoms.getAtom(SpeciesWater4P.indexM).getPosition();

            comWi.Ea1Tv1(massH, H11r);
            comWi.PEa1Tv1(massO, O1r);
            comWi.PEa1Tv1(massH, H12r);
            comWi.TE(1.0/totalMass);//c.o.m of molecule i
            
            comHi1r = space.makeVector();
            comHi2r = space.makeVector();
            comMir = space.makeVector();
            comHi1r.Ev1Mv2(comWi, H11r);
            comHi2r.Ev1Mv2(comWi, H12r);
            comMir.Ev1Mv2(comWi, M1r);
            iDipoleMoment.Ea1Tv1(chargeH, comHi1r);
            iDipoleMoment.PEa1Tv1(chargeH, comHi2r);
            iDipoleMoment.PEa1Tv1(chargeM, comMir);
            
            myRq.Ea1Tv1(myRqFactor,iDipoleMoment);
            
            int neighborCount = 0;
            
            for (int j=0; j<atoms.getMoleculeCount(); j++) {
                if  (i == j) continue;
                IAtomList jLeafAtoms = atoms.getMolecule(j).getChildList();
                Vector Mjr = jLeafAtoms.getAtom(SpeciesWater4P.indexM).getPosition();
                Vector Ojr = jLeafAtoms.getAtom(SpeciesWater4P.indexO).getPosition();
                Vector Hj1r = jLeafAtoms.getAtom(SpeciesWater4P.indexH1).getPosition();
                Vector Hj2r = jLeafAtoms.getAtom(SpeciesWater4P.indexH2).getPosition();
                
                comWj.Ea1Tv1(massH, Hj1r);
                comWj.PEa1Tv1(massO, Ojr);
                comWj.PEa1Tv1(massH, Hj2r);
                comWj.TE(1.0/totalMass);
                

                
                rijVector.Ev1Mv2(comWi,comWj);
                shift.Ea1Tv1(-1,rijVector);
        		box.getBoundary().nearestImage(rijVector);
                shift.PE(rijVector);
                final boolean zeroShift = shift.squared() < 0.1;
                
                double r12 = Math.sqrt(rijVector.squared());
                if (i==0 && j==3){
                	//System.out.println("r12 "+r12);
                }
                if (r12 > cutOffDistance) {
                	if (i < j){
	                    int mOffset = i*3;
	                    int nOffset = j*3;
	                    for (int m=0; m<3; m++) {
	                        for (int n=0; n<3; n++) {
	                            myA.set(mOffset+m, nOffset+n, 0.0);
	                            myA.set(nOffset+n, mOffset+m, 0.0);
	                        }
	                    }
                	}
                	continue;
                	
                }
                	
                
                neighborCount++;
                comHj1r = space.makeVector();
                comHj2r = space.makeVector();
                comMjr = space.makeVector();
                comHj1r.Ev1Mv2(comWj, Hj1r);
                comHj2r.Ev1Mv2(comWj, Hj2r);
                comMjr.Ev1Mv2(comWj, Mjr);
                jDipoleMoment.Ea1Tv1(chargeH, comHj1r);
                jDipoleMoment.PEa1Tv1(chargeH, comHj2r);
                jDipoleMoment.PEa1Tv1(chargeM, comMjr);
                
            	myRq.PEa1Tv1(myRqFactor,jDipoleMoment);//summation of Mu_j
                
                double comWtoH1,comWtoH2,comWtoM;
                
                if (zeroShift){
                
                	comWtoH1 = Math.sqrt(comWi.Mv1Squared(Hj1r));
                	comWtoH2 = Math.sqrt(comWi.Mv1Squared(Hj2r));
                	comWtoM = Math.sqrt(comWi.Mv1Squared(Mjr));
                }
                else {
                	shift.PE(comWi);
                	comWtoH1 = Math.sqrt(Hj1r.Mv1Squared(shift));
                	shift.ME(comWi);
                	
                	shift.PE(comWi);
                	comWtoH2 = Math.sqrt(Hj2r.Mv1Squared(shift));
                	shift.ME(comWi);
                	
                	shift.PE(comWi);
                	comWtoM = Math.sqrt(Mjr.Mv1Squared(shift));
                	shift.ME(comWi);
                	
                }
                               

                // For molecules that are far apart, fac=chargeX/comWtoX^3, but we add up
                // facs for H and M, which mostly cancel each other out, so we lose quite 
                // a bit of precision (~2-3 digits).
                double fac = chargeH/(comWtoH1*comWtoH1*comWtoH1)*((1-SpecialFunctions.erfc(comWtoH1/sqrtHMsigmas))
                        -Math.sqrt(2)*comWtoH1/sqrtPiHMsigmas*Math.exp(-comWtoH1*comWtoH1/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
                work.Ev1Mv2(comWi, Hj1r);
                work.PE(shift);
                work.TE(fac);
                myEq.set(i*3+0, 0, myEq.get(i*3+0, 0)+work.getX(0));
                myEq.set(i*3+1, 0, myEq.get(i*3+1, 0)+work.getX(1));
                myEq.set(i*3+2, 0, myEq.get(i*3+2, 0)+work.getX(2));
    
                fac = chargeH/(comWtoH2*comWtoH2*comWtoH2)*((1-SpecialFunctions.erfc(comWtoH2/sqrtHMsigmas))
                        -Math.sqrt(2)*comWtoH2/sqrtPiHMsigmas*Math.exp(-comWtoH2*comWtoH2/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
                work.Ev1Mv2(comWi, Hj2r);
                work.PE(shift);
                work.TE(fac);
                myEq.set(i*3+0, 0, myEq.get(i*3+0, 0)+work.getX(0));
                myEq.set(i*3+1, 0, myEq.get(i*3+1, 0)+work.getX(1));
                myEq.set(i*3+2, 0, myEq.get(i*3+2, 0)+work.getX(2));
    
                fac = chargeM/(comWtoM*comWtoM*comWtoM)*((1-SpecialFunctions.erfc(comWtoM/(2*sigmaM)))
                        -Math.sqrt(2)*comWtoM/sqrtPiMMsigmas*Math.exp(-comWtoM*comWtoM/(4*sigmaM*sigmaM)));
                work.Ev1Mv2(comWi, Mjr);
                work.PE(shift);
                work.TE(fac);
                myEq.set(i*3+0, 0, myEq.get(i*3+0, 0)+work.getX(0));
                myEq.set(i*3+1, 0, myEq.get(i*3+1, 0)+work.getX(1));
                myEq.set(i*3+2, 0, myEq.get(i*3+2, 0)+work.getX(2));
                     
                if (i<j) {
                    double OOr2;
                    if(zeroShift){
                    	OOr2 = O1r.Mv1Squared(Ojr);
                    }else{
                    	shift.PE(O1r);
                    	OOr2 = Ojr.Mv1Squared(shift);
                    	shift.ME(O1r);
                    }
                    if (OOr2 < core) {
                        UpolAtkins = Double.NaN;
                        return UpolAtkins;
                    }
                    comWj.Ea1Tv1(massH, Hj1r);
                    comWj.PEa1Tv1(massO, Ojr);
                    comWj.PEa1Tv1(massH, Hj2r);
                    comWj.TE(1.0/totalMass);
                    
                    rijVector.Ev1Mv2(comWi,comWj);
                    rijVector.PE(shift);
                    
                    r12 = Math.sqrt(rijVector.squared());



                    double f = (1-SpecialFunctions.erfc(r12/(2*sigmaM)))-(r12/(sigmaM*Math.sqrt(Math.PI)) + (r12*r12*r12)/(6*Math.sqrt(Math.PI)*sigmaM*sigmaM*sigmaM))*Math.exp(-r12*r12/(4*sigmaM*sigmaM));
                    
                    double g = (1-SpecialFunctions.erfc(r12/(2*sigmaM)))-(r12/(sigmaM*Math.sqrt(Math.PI)))*Math.exp(-r12*r12/(4*sigmaM*sigmaM));
                    
                    // Filling the unit matrix I
                    Tij.Ev1v2(rijVector,rijVector);//Each tensor Tij is a 3X3 matrix
                    
                    Tij.TE(3*f/(r12*r12));
                    
                    Tunit.E(g);
                    
                    Tij.ME(Tunit);
                    Tij.TE(1/(r12*r12*r12));
                    if (i==0&&j==3){
                    	//System.out.println("Tij "+Tij);
                    }
                    
                    //Try matrix inversion solution with Jama library
                            
                    Tij.TE(alphaPol);
                    
                    int mOffset = i*3;
                    int nOffset = j*3;
                    for (int m=0; m<3; m++) {
                        for (int n=0; n<3; n++) {
                            myA.set(mOffset+m, nOffset+n, -Tij.component(m, n));
                            myA.set(nOffset+n, mOffset+m, -Tij.component(n, m));
                        }
                    }
                }
            }
            //System.out.println("neighbor "+neighborCount);
            //add reaction field to myEq
            myEq.set(i*3+0, 0, myEq.get(i*3+0, 0)+myRq.getX(0));
            myEq.set(i*3+1, 0, myEq.get(i*3+1, 0)+myRq.getX(1));
            myEq.set(i*3+2, 0, myEq.get(i*3+2, 0)+myRq.getX(2));
            
        }
        
        //x here represents P (almost).
        //For x to be P, the A of the Ax=b actually needs an extra factor of
        //alphaPol.  We'll add that bit in when we calculate UpolAtkins.  

        Matrix x = myA.solve(myEq);//myA*x=myEq
        //myA.print(10, 5);

        if (false) {
            // this is (mathematically) what we want.  But Jama is slow.
            UpolAtkins = -0.5*(x.transpose().times(myEq)).get(0,0)*alphaPol;
        }
        else {
            UpolAtkins = 0;
            Matrix dipolePol = x.times(alphaPol);            
            //System.out.println("induced dipole "+Debye.UNIT.fromSim(dipolePol.get(2,0)));
            //dipolePol.print(10, 5);
            //System.out.println("x "+x.get);
            for (int i=0; i<3*atomCount; i++) {
                UpolAtkins += x.get(i,0)*myEq.get(i,0);
                //System.out.println("i "+i+" myA "+myA.get(i, 0));
                //System.out.println("UpolAtkins "+UpolAtkins+" i "+i+" x "+x.get(i,0)+" myEq "+myEq.get(i,0));
            }

            UpolAtkins *= -0.5*alphaPol;
        }
        // only needed for more complicated Eq8 from Cummings paper 
        if (false) {
            
            // for the sake of clarity (over perf), just multiply x by alphaPol
            // (see comment above about A lacking alphaPol)
            x.timesEquals(alphaPol);
            Matrix Ep = myA.times(x).minus(x);
            Ep.timesEquals(-1/alphaPol);

            double x2NormF = x.normF();
            double UpolEquation8 = 2*UpolAtkins -0.5*(x.transpose().times(Ep).get(0,0))+(0.5/alphaPol)*(x2NormF*x2NormF);

            if (Math.abs(UpolAtkins-UpolEquation8) > 1.e-6) {
                throw new RuntimeException("oops "+UpolAtkins+" "+UpolEquation8);
            }
        }
        
        return UpolAtkins;
    }
    
    public double getLastPolarizationEnergy() {
        return UpolAtkins;
    }

    public final double getRange() {
        return Double.POSITIVE_INFINITY;
    }
    
    public void setCutOffDistance(double a){
    	cutOffDistance = a;
    	initRqFactor();
    }
    
    public void setTemperature(double t){    	
    	temperature = t;
    	initRqFactor();
    }
    
    public void setRho(double c){
    	rho = c;
    	CompoundUnit rhoUnit = new CompoundUnit(new Unit[]{Mole.UNIT,Liter.UNIT},new double[]{1,-1});
        rho = rhoUnit.fromSim(rho);
    	initRqFactor();
    }
    
    public void initRqFactor(){
    	rho =  box.getMoleculeList().getMoleculeCount()/ box.getBoundary().volume();
    	CompoundUnit rhoUnit = new CompoundUnit(new Unit[]{Mole.UNIT,Liter.UNIT},new double[]{1,-1});
        rho = rhoUnit.fromSim(rho);
        double reducedT = temperature/298.15;
        //System.out.println("reduced T "+reducedT);
        double reducedT2 = reducedT*reducedT;
        //double reducedRho = rho/6.02214*18.02*1000;//molecule/A3
        double reducedRho = rho*18.02/1000;//mol/L
        //System.out.println("volume "+box.getBoundary().volume()+"rho "+rho+" reduced rho "+reducedRho);
        double reducedRho2 = reducedRho*reducedRho;
        double reducedRho3 = reducedRho*reducedRho2;
        double reducedRho4 = reducedRho2*reducedRho2;
        double epsilonRF = 1+7.62571/reducedT*reducedRho+(244.003/reducedT-140.569+27.7841*reducedT)*reducedRho2+
        					(-96.2805/reducedT+41.7909*reducedT-10.2099*reducedT2)*reducedRho3+(-45.2059/reducedT2+84.6395/reducedT-35.8644)*reducedRho4;
        myRqFactor = (epsilonRF-1)*2/(2*epsilonRF+1)/Math.pow(cutOffDistance, 3);  
    }
    
    public void setBox(Box box) {
    	this.box = box;
    }
    public void setDodebug(boolean a){
    	dodebug = a;
    }

    private static final long serialVersionUID = 1L;
    protected final MoleculePair pair;
    protected Boundary boundary;
    protected final double sigma;
    protected final double epsilon, gamma;
    protected final double chargeH, chargeM;
    protected final double core; //4.41 = 2.1^2; value according to Cummings
    protected Matrix[] Eq, Rq, A;
    protected Vector comW1, comW2, comWi, comWj;
    protected Vector comH11r, comH12r, comM1r, comH21r, comH22r, comM2r, comHi1r, comHi2r, comMir, comHj1r, comHj2r, comMjr;
    protected final Vector r12Vector, rijVector;
    protected final Vector work, shift, dipoleMoment1,dipoleMoment2, iDipoleMoment, jDipoleMoment, myRq;
    protected final Tensor Tunit, Tij;
    protected final double sigmaM;
    protected final double sigmaH;
    protected final double sqrtHMsigmas;
    protected final double massH;
    protected final double massO;
    protected final double totalMass;
    protected final double sqrtPiHMsigmas;
    protected final double sqrtPiMMsigmas;
    protected final double alphaPol;
    private double cutOffDistance;
    private double temperature, rho;
    protected double UpolAtkins;
    protected double myRqFactor;
    protected Box box;
    public static boolean dodebug=false;

}