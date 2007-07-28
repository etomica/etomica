package etomica.models.water;

import etomica.atom.AtomSet;
import etomica.box.Box;
import etomica.exception.MethodNotImplementedException;
import etomica.math.SpecialFunctions;
import etomica.potential.Potential2;
import etomica.potential.Potential2Soft;
import etomica.potential.PotentialPolarizable;
import etomica.space.IVector;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.space3d.Tensor3D;
import etomica.space3d.Vector3D;
import etomica.units.Debye;
import etomica.units.Electron;
import etomica.units.Kelvin;

/** 
 * 
 * @author kofke
 *
 * Special version of water model that permits no truncation of the potential
 * and no periodic boundary. Used for cluster calculations.
 */
public class PotentialWaterGCPM3forB5 extends Potential2 implements Potential2Soft, PotentialPolarizable {

	public PotentialWaterGCPM3forB5(Space space) {
	    super(space);
	    setSigma(3.69);
	    setEpsilon(Kelvin.UNIT.toSim(110));
	    //setGamma(12.75);  Do I want to introduce a method for this, to be addressed in contructor of superclass (or this one)? kmb, 8/7/06
        chargeH11 = Electron.UNIT.toSim(0.6113);
        chargeH12 = Electron.UNIT.toSim(0.6113);
        chargeM1 = Electron.UNIT.toSim(-1.2226);
        chargeH21 = Electron.UNIT.toSim(0.6113);
        chargeH22 = Electron.UNIT.toSim(0.6113);
        chargeM2 = Electron.UNIT.toSim(-1.2226);
        chargeH31 = Electron.UNIT.toSim(0.6113);
        chargeH32 = Electron.UNIT.toSim(0.6113);
        chargeM3 = Electron.UNIT.toSim(-1.2226);
        chargeH41 = Electron.UNIT.toSim(0.6113);
        chargeH42 = Electron.UNIT.toSim(0.6113);
        chargeM4 = Electron.UNIT.toSim(-1.2226);
        core = 4.41; //4.41 = 2.1^2; value according to Cummings
        sigmaM = 0.610;
        sigmaH = 0.455;
        sqrtHMsigmas = Math.sqrt(2*(sigmaH*sigmaH+sigmaM*sigmaM));
        massH = 1.01;
        massO = 16.0;
        totalMass = 18.02;
        sqrtPiHMsigmas = Math.sqrt(Math.PI*(sigmaH*sigmaH+sigmaM*sigmaM));
        sqrtPiMMsigmas = Math.sqrt(Math.PI*(2*sigmaM*sigmaM));

        Eq1 = space.makeVector();
        Eq2 = space.makeVector();
        Eq3 = space.makeVector();
        Ep1 = space.makeVector();
        Ep2 = space.makeVector();
        Ep3 = space.makeVector();
        
        P1 = space.makeVector();
        P2 = space.makeVector();
        P3 = space.makeVector();
        P1old = space.makeVector();
        P2old = space.makeVector();
        P3old = space.makeVector();

        comW1 = space.makeVector();
        comW2 = space.makeVector();
        comW3 = space.makeVector();

        r12Vector = space.makeVector();
        T12row1 = space.makeVector();
        T12row2 = space.makeVector();
        T12row3 = space.makeVector();
        T12P1 = space.makeVector();
        T12P2 = space.makeVector();
        B1 = space.makeVector();
        T12Eq2 = space.makeVector();
        
        Tunit = space.makeTensor();
        T12 = space.makeTensor();
        A1 = space.makeTensor();
        T12T12 = space.makeTensor();
        inverseA1 = space.makeTensor();
    }   

    /* This energy method returns the pair energies with SCF solution for
	 * the third virial coefficient.
	 * 
	 * No guarantee that one of the water molecules is at the origin.; kmb 7/12/06
	 */
    public double energy(AtomSet atoms){
        if (atoms.getAtomCount() == 3) return energySCF(atoms);
        if (atoms.getAtomCount() == 4) return energySCF4(atoms);
        //if (atoms.getAtomCount() == 5) return energySCF5(atoms);
    		double sum = 0.0;
        double r2 = 0.0;
        
        
        AtomWater4P node1 = (AtomWater4P)atoms.getAtom(0);
        AtomWater4P node2 = (AtomWater4P)atoms.getAtom(1);

        IVector O1r = node1.O.getPosition();
        IVector O2r = node2.O.getPosition();
        IVector H11r = node1.H1.getPosition();
        IVector H12r = node1.H2.getPosition();
        IVector H21r = node2.H1.getPosition();
        IVector H22r = node2.H2.getPosition();

        IVector M1r = node1.M.getPosition();
        IVector M2r = node2.M.getPosition();
                
// C2v geometry for dimer        
/*        O2r.setX(2,2.2);
        O2r.setX(0,0.0);
        O2r.setX(1,0.0);
        H21r.setX(2,2.785882276618295);
        H21r.setX(0,0.0);
        H21r.setX(1,0.756950327263661);
        H22r.setX(2,2.785882276618295);
        H22r.setX(0,0.0);
        H22r.setX(1,-0.756950327263661);
        M2r.setX(2,2.47);
        M2r.setX(0,0.0);
        M2r.setX(1,0.0);
*/
        
/*// Cs geometry for dimer        
        O2r.setX(2,1.5921074871731);
        O2r.setX(0,2.39631440952296);
        O2r.setX(1,0.0);
        H21r.setX(2,1.37452965421812);
        H21r.setX(0,2.94029779958377);
        H21r.setX(1,0.756950327263661);
        H22r.setX(2,1.37452965421812);
        H22r.setX(0,2.94029779958377);
        H22r.setX(1,-0.756950327263661);
        M2r.setX(2,1.49183817157454);
        M2r.setX(0,2.64700558278081);
        M2r.setX(1,0.0);
*/
        
	// moved to the constructor; KMB, 7/20/07	
//        final double core = 4.41; //4.41 = 2.1^2; value according to Cummings
  

        // Initializing vectors below moved to the constructor; KMB, 7/20/07
/*        Vector Eq1 = new Vector3D();
        Vector Eq2 = new Vector3D();
        Vector Ep1 = new Vector3D();
        Vector Ep2 = new Vector3D();


        Vector P1 = new Vector3D();
        Vector P2 = new Vector3D();
        Vector P1old = new Vector3D();
        Vector P2old = new Vector3D();
*/


        // Need loop to check for configuration overlap between charged particles
        
//      compute O-O distance to consider bypassing the SCF loop   
//      compute O-O distance to consider truncation   
        r2 = O1r.Mv1Squared(O2r);
        
        if(r2<=core) {
        		return Double.POSITIVE_INFINITY;
        }
        
        gamma = 12.75;
	    double r = Math.sqrt(r2);
        double rOverSigma = r/sigma;
        double sigma2OverR2 = 1/(rOverSigma*rOverSigma);
        double sixOverGamma = 6/gamma;
   
        sum = epsilon/(1 - sixOverGamma)*(sixOverGamma*Math.exp(gamma*(1 - rOverSigma)) - sigma2OverR2*sigma2OverR2*sigma2OverR2);

// 		Moved to constructor; KMB, 7/20/07        
/*        double sigmaM = 0.610;
        double sigmaH = 0.455;
        double sqrtHMsigmas = Math.sqrt(2*(sigmaH*sigmaH+sigmaM*sigmaM));
*/        
        // MUST INCLUDE ERF FUNCTION STUFF TO COULOMBIC ENERGY PART!
        // KMB 8/3/06
        
        r2 = H11r.Mv1Squared(H21r);
        sum += chargeH11*chargeH21/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaH)));

        r2 = H11r.Mv1Squared(H22r);
        sum += chargeH11*chargeH22/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaH)));

        r2 = H12r.Mv1Squared(H21r);
        sum += chargeH12*chargeH21/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaH)));

        r2 = H12r.Mv1Squared(H22r);
        sum += chargeH12*chargeH22/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaH)));
        
        r2 = M1r.Mv1Squared(H21r);
        sum += chargeH21*chargeM1/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/sqrtHMsigmas));

        r2 = M1r.Mv1Squared(H22r);
        sum += chargeH22*chargeM1/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/sqrtHMsigmas));

        r2 = M2r.Mv1Squared(H11r);
        sum += chargeH11*chargeM2/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/sqrtHMsigmas));

        r2 = M2r.Mv1Squared(H12r);
        sum += chargeH12*chargeM2/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/sqrtHMsigmas));

        r2 = M1r.Mv1Squared(M2r);
        sum += chargeM1*chargeM2/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaM)));
        
        /*
         * Finding the Electric fields at the center of mass of each molecule, Eqi
         * kmb, 8/7/06
         */

        // Moved to the constructor; KMB, 7/20/07
/*        Vector3D comW1 = new Vector3D();
        Vector3D comW2 = new Vector3D();
*/        
        // How to use COM etomica code correctly? kmb, 8/7/06
/*        DataSourceCOM com = new DataSourceCOM(space);
        com.actionPerformed(atoms.getAtom(0));
        com.getCOM();
*/        

        // Moved to the constructor; KMB, 7/20/07
/*        double massH = 1.01;
        double massO = 16.0;
        double totalMass = 18.02;
*/ 
        double comW1Xcomp = 0.0;
        double comW1Ycomp = 0.0;
        double comW1Zcomp = 0.0;
        
        comW1Xcomp = massH*H11r.x(0) + massO*O1r.x(0) + massH*H12r.x(0);
        comW1Ycomp = massH*H11r.x(1) + massO*O1r.x(1) + massH*H12r.x(1);
        comW1Zcomp = massH*H11r.x(2) + massO*O1r.x(2) + massH*H12r.x(2);
        
        comW1.setX(0,comW1Xcomp);
        comW1.setX(1,comW1Ycomp);
        comW1.setX(2,comW1Zcomp);
        
        comW1.Ea1Tv1(1/totalMass,comW1);

        double comW2Xcomp = 0.0;
        double comW2Ycomp = 0.0;
        double comW2Zcomp = 0.0;
        
        comW2Xcomp = massH*H21r.x(0) + massO*O2r.x(0) + massH*H22r.x(0);
        comW2Ycomp = massH*H21r.x(1) + massO*O2r.x(1) + massH*H22r.x(1);
        comW2Zcomp = massH*H21r.x(2) + massO*O2r.x(2) + massH*H22r.x(2);
        
        comW2.setX(0,comW2Xcomp);
        comW2.setX(1,comW2Ycomp);
        comW2.setX(2,comW2Zcomp);
        
        comW2.Ea1Tv1(1/totalMass,comW2);

//      Moved to constructor; KMB, 7/20/07
/*        double sqrtPiHMsigmas = Math.sqrt(Math.PI*(sigmaH*sigmaH+sigmaM*sigmaM));
        double sqrtPiMMsigmas = Math.sqrt(Math.PI*(2*sigmaM*sigmaM));
*/        
        double Eq1Xcomp = 0.0;
        double Eq1Ycomp = 0.0;
        double Eq1Zcomp = 0.0;
        
        double comW1toH21 = Math.sqrt(comW1.Mv1Squared(H21r));
        double comW1toH22 = Math.sqrt(comW1.Mv1Squared(H22r));
        double comW1toM2 = Math.sqrt(comW1.Mv1Squared(M2r));

        Eq1Xcomp += chargeH21*(comW1.x(0)-H21r.x(0))/(comW1toH21*comW1toH21*comW1toH21)*((1-SpecialFunctions.erfc(comW1toH21/sqrtHMsigmas))-Math.sqrt(2)*comW1toH21/sqrtPiHMsigmas*Math.exp(-comW1toH21*comW1toH21/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq1Xcomp += chargeH22*(comW1.x(0)-H22r.x(0))/(comW1toH22*comW1toH22*comW1toH22)*((1-SpecialFunctions.erfc(comW1toH22/sqrtHMsigmas))-Math.sqrt(2)*comW1toH22/sqrtPiHMsigmas*Math.exp(-comW1toH22*comW1toH22/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq1Xcomp += chargeM2*(comW1.x(0)-M2r.x(0))/(comW1toM2*comW1toM2*comW1toM2)*((1-SpecialFunctions.erfc(comW1toM2/(2*sigmaM)))-Math.sqrt(2)*comW1toM2/sqrtPiMMsigmas*Math.exp(-comW1toM2*comW1toM2/(4*sigmaM*sigmaM)));

        Eq1Ycomp += chargeH21*(comW1.x(1)-H21r.x(1))/(comW1toH21*comW1toH21*comW1toH21)*((1-SpecialFunctions.erfc(comW1toH21/sqrtHMsigmas))-Math.sqrt(2)*comW1toH21/sqrtPiHMsigmas*Math.exp(-comW1toH21*comW1toH21/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq1Ycomp += chargeH22*(comW1.x(1)-H22r.x(1))/(comW1toH22*comW1toH22*comW1toH22)*((1-SpecialFunctions.erfc(comW1toH22/sqrtHMsigmas))-Math.sqrt(2)*comW1toH22/sqrtPiHMsigmas*Math.exp(-comW1toH22*comW1toH22/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq1Ycomp += chargeM2*(comW1.x(1)-M2r.x(1))/(comW1toM2*comW1toM2*comW1toM2)*((1-SpecialFunctions.erfc(comW1toM2/(2*sigmaM)))-Math.sqrt(2)*comW1toM2/sqrtPiMMsigmas*Math.exp(-comW1toM2*comW1toM2/(4*sigmaM*sigmaM)));

        Eq1Zcomp += chargeH21*(comW1.x(2)-H21r.x(2))/(comW1toH21*comW1toH21*comW1toH21)*((1-SpecialFunctions.erfc(comW1toH21/sqrtHMsigmas))-Math.sqrt(2)*comW1toH21/sqrtPiHMsigmas*Math.exp(-comW1toH21*comW1toH21/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq1Zcomp += chargeH22*(comW1.x(2)-H22r.x(2))/(comW1toH22*comW1toH22*comW1toH22)*((1-SpecialFunctions.erfc(comW1toH22/sqrtHMsigmas))-Math.sqrt(2)*comW1toH22/sqrtPiHMsigmas*Math.exp(-comW1toH22*comW1toH22/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq1Zcomp += chargeM2*(comW1.x(2)-M2r.x(2))/(comW1toM2*comW1toM2*comW1toM2)*((1-SpecialFunctions.erfc(comW1toM2/(2*sigmaM)))-Math.sqrt(2)*comW1toM2/sqrtPiMMsigmas*Math.exp(-comW1toM2*comW1toM2/(4*sigmaM*sigmaM)));

        Eq1.setX(0,Eq1Xcomp);
        Eq1.setX(1,Eq1Ycomp);
        Eq1.setX(2,Eq1Zcomp);

                
        double Eq2Xcomp = 0.0;
        double Eq2Ycomp = 0.0;
        double Eq2Zcomp = 0.0;
        
        double comW2toH11 = Math.sqrt(comW2.Mv1Squared(H11r));
        double comW2toH12 = Math.sqrt(comW2.Mv1Squared(H12r));
        double comW2toM1 = Math.sqrt(comW2.Mv1Squared(M1r));

        Eq2Xcomp += chargeH11*(comW2.x(0)-H11r.x(0))/(comW2toH11*comW2toH11*comW2toH11)*((1-SpecialFunctions.erfc(comW2toH11/sqrtHMsigmas))-Math.sqrt(2)*comW2toH11/sqrtPiHMsigmas*Math.exp(-comW2toH11*comW2toH11/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq2Xcomp += chargeH12*(comW2.x(0)-H12r.x(0))/(comW2toH12*comW2toH12*comW2toH12)*((1-SpecialFunctions.erfc(comW2toH12/sqrtHMsigmas))-Math.sqrt(2)*comW2toH12/sqrtPiHMsigmas*Math.exp(-comW2toH12*comW2toH12/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq2Xcomp += chargeM1*(comW2.x(0)-M1r.x(0))/(comW2toM1*comW2toM1*comW2toM1)*((1-SpecialFunctions.erfc(comW2toM1/(2*sigmaM)))-Math.sqrt(2)*comW2toM1/sqrtPiMMsigmas*Math.exp(-comW2toM1*comW2toM1/(4*sigmaM*sigmaM)));

        Eq2Ycomp += chargeH11*(comW2.x(1)-H11r.x(1))/(comW2toH11*comW2toH11*comW2toH11)*((1-SpecialFunctions.erfc(comW2toH11/sqrtHMsigmas))-Math.sqrt(2)*comW2toH11/sqrtPiHMsigmas*Math.exp(-comW2toH11*comW2toH11/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq2Ycomp += chargeH12*(comW2.x(1)-H12r.x(1))/(comW2toH12*comW2toH12*comW2toH12)*((1-SpecialFunctions.erfc(comW2toH12/sqrtHMsigmas))-Math.sqrt(2)*comW2toH12/sqrtPiHMsigmas*Math.exp(-comW2toH12*comW2toH12/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq2Ycomp += chargeM1*(comW2.x(1)-M1r.x(1))/(comW2toM1*comW2toM1*comW2toM1)*((1-SpecialFunctions.erfc(comW2toM1/(2*sigmaM)))-Math.sqrt(2)*comW2toM1/sqrtPiMMsigmas*Math.exp(-comW2toM1*comW2toM1/(4*sigmaM*sigmaM)));

        Eq2Zcomp += chargeH11*(comW2.x(2)-H11r.x(2))/(comW2toH11*comW2toH11*comW2toH11)*((1-SpecialFunctions.erfc(comW2toH11/sqrtHMsigmas))-Math.sqrt(2)*comW2toH11/sqrtPiHMsigmas*Math.exp(-comW2toH11*comW2toH11/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq2Zcomp += chargeH12*(comW2.x(2)-H12r.x(2))/(comW2toH12*comW2toH12*comW2toH12)*((1-SpecialFunctions.erfc(comW2toH12/sqrtHMsigmas))-Math.sqrt(2)*comW2toH12/sqrtPiHMsigmas*Math.exp(-comW2toH12*comW2toH12/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq2Zcomp += chargeM1*(comW2.x(2)-M1r.x(2))/(comW2toM1*comW2toM1*comW2toM1)*((1-SpecialFunctions.erfc(comW2toM1/(2*sigmaM)))-Math.sqrt(2)*comW2toM1/sqrtPiMMsigmas*Math.exp(-comW2toM1*comW2toM1/(4*sigmaM*sigmaM)));

        Eq2.setX(0,Eq2Xcomp);
        Eq2.setX(1,Eq2Ycomp);
        Eq2.setX(2,Eq2Zcomp);

        
        /*
         * Finding the tensor used to relate the induced dipole moment Pi with the induced electric field Epi.
         * kmb, 8/7/06
         */
        
        double r12 = Math.sqrt(comW1.Mv1Squared(comW2));

        r12Vector.Ev1Mv2(comW2,comW1);  // is this the correct direction? kmb, 8/7/06
        

        double f = (1-SpecialFunctions.erfc(r12/(2*sigmaM)))-(r12/(sigmaM*Math.sqrt(Math.PI)) + (r12*r12*r12)/(6*Math.sqrt(Math.PI)*sigmaM*sigmaM*sigmaM))*Math.exp(-r12*r12/(4*sigmaM*sigmaM));
        
        double g = (1-SpecialFunctions.erfc(r12/(2*sigmaM)))-(r12/(sigmaM*Math.sqrt(Math.PI)))*Math.exp(-r12*r12/(4*sigmaM*sigmaM));
        
        // Filling the unit matrix I
        
        T12.PEv1v2(r12Vector,r12Vector);
        
        // just added these to test whether T12 = T21; IT DOES!
/*        Tensor3D T21 = new Tensor3D();
        
        T21.E(r21Vector,r21Vector);
*/        
        T12.TE(3*f/(r12*r12));
        
        Tunit.E(g);
        
        T12.ME(Tunit);
        T12.TE(1/(r12*r12*r12));
        
        // T12 = T21, so I can get by for now in the case of B2!
        
        // Now distribute the elements of the tensor into 3 separate "row" vectors
        // so I can do dot products with etomica math methods
        // kmb, 8/7/06
        
        T12row1.setX(0,T12.component(0,0));
        T12row1.setX(1,T12.component(0,1));
        T12row1.setX(2,T12.component(0,2));
        T12row2.setX(0,T12.component(1,0));
        T12row2.setX(1,T12.component(1,1));
        T12row2.setX(2,T12.component(1,2));
        T12row3.setX(0,T12.component(2,0));
        T12row3.setX(1,T12.component(2,1));
        T12row3.setX(2,T12.component(2,2));
        
        
        //Begin analytical solution for Upol for pair energy
        
        double alphaPol = 1.444;
        double alphaPol2 = alphaPol*alphaPol;
        
        A1.E(0);
        A1.setComponent(0,0,1);
        A1.setComponent(1,1,1);
        A1.setComponent(2,2,1);
        T12T12.E(T12);
        T12T12.TE(T12);
        T12T12.TE(alphaPol2);
        A1.ME(T12T12);
        
        T12Eq2.setX(0,T12row1.dot(Eq2));
        T12Eq2.setX(1,T12row2.dot(Eq2));
        T12Eq2.setX(2,T12row3.dot(Eq2));
        T12Eq2.TE(alphaPol2);
        
        B1.E(Eq1);
        B1.TE(alphaPol);
        B1.PE(T12Eq2);

        inverseA1.E(A1);
        inverseA1.inverse();
        
        //use transform method in Tensor3D!
        //inverseA3.transform(P1);
        
        P1.E(B1);
        inverseA1.transform(P1);

        // Now find P2
        T12P1.E(P1);
        T12P1.TE(alphaPol);
        T12.transform(T12P1);
        
        T12P2.E(P2);
        T12P2.TE(alphaPol);
        T12.transform(T12P2);
        
        P2.E(Eq2);
        P2.TE(alphaPol);
        P2.PE(T12P1);
        
        
        //System.out.println("From analytical: P1 = " + P1 + ", P2 = " + P2);
        
        // Now find Ep1, Ep2
        
        Ep1.E(T12P2);
        Ep2.E(T12P1);
        
//        double UpolAtkinsAnalytical = -0.5*(P1.dot(Eq1)+P2.dot(Eq2));

        
        
/*        // Set the induced dipole moments equal to 10% of the permanent dipole value
        // kmb, 8/7/06
        P1.setX(0,0); // 14.3952507082236);
        P1.setX(1,0); //14.3952507082236);
        P1.setX(2,0); //14.3952507082236);
        P2.setX(0,0); //14.3952507082236);
        P2.setX(1,0); // 14.3952507082236);
        P2.setX(2,0); //14.3952507082236);
        
        P1old.E(P1);
        P2old.E(P2);
        
        // kmb add these here and initialize, 8/22/06
        // trying to solve convergence problem between dimer and trimer
        
        double deltaP1 = 1.0;
        double deltaP2 = 1.0;
        
        
        while (noSCFforP1 || noSCFforP2) {
        	        	
        		// First calculate Ep1, based upon guess for P2
        	
        		Ep1.setX(0,T12row1.dot(P2));
        		Ep1.setX(1,T12row2.dot(P2));
        		Ep1.setX(2,T12row3.dot(P2));
        		
        		// Now calculate P1 from the value of Ep1
        		
//        		double alphaPol = 1.444;
        		
        		P1.setX(0,alphaPol*(Eq1.x(0) + Ep1.x(0)));
        		P1.setX(1,alphaPol*(Eq1.x(1) + Ep1.x(1)));
        		P1.setX(2,alphaPol*(Eq1.x(2) + Ep1.x(2)));

        		// Next calculate Ep2
        		
        		Ep2.setX(0,T12row1.dot(P1));
        		Ep2.setX(1,T12row2.dot(P1));
        		Ep2.setX(2,T12row3.dot(P1));
        		
        		// Now calculate new P2
        		
        		P2.setX(0,alphaPol*(Eq2.x(0) + Ep2.x(0)));
        		P2.setX(1,alphaPol*(Eq2.x(1) + Ep2.x(1)));
        		P2.setX(2,alphaPol*(Eq2.x(2) + Ep2.x(2)));

        		// Evaluate the criteria
        		
	        	deltaP1 = Math.sqrt(P1.Mv1Squared(P1old));
	        	deltaP2 = Math.sqrt(P2.Mv1Squared(P2old));
	        	
	        	counterSCFloop = counterSCFloop + 1;
	    
	        	
	        	if (deltaP1<1e-15) { //0.00338) { // 0.00338 orig value
	        		noSCFforP1 = false;
	        	}
	        	else {
	        		P1old.E(P1);
	        	}
	        	if (deltaP2<1e-15) { //0.00338) {
	        		noSCFforP2 = false;
	        	}
	        	else {
	        		P2old.E(P2);
	        	}

	        	if (counterSCFloop >= 1000) {
	        		noSCFforP1 = false;
	        		noSCFforP2 = false;
    		        chargeH11 = Electron.UNIT.toSim(0.6113);
    		        chargeH12 = Electron.UNIT.toSim(0.6113);
    		        chargeM1 = Electron.UNIT.toSim(-1.2226);
    		        chargeH21 = Electron.UNIT.toSim(0.6113);
    		        chargeH22 = Electron.UNIT.toSim(0.6113);
    		        chargeM2 = Electron.UNIT.toSim(-1.2226);
        		
	        	}
	        	
        }
        
	        // REPEAT LOOP HERE.
*/    
        
/*
 * Here is where I need to add the polarization term to the energy sum.
 * kmb 5/4/06
 */        
        
  //      double alpha = 1.444;
        
        UpolAtkins = -0.5*(P1.dot(Eq1)+P2.dot(Eq2));
    
//        double UpolEquation8 = -(P1.dot(Eq1)+P2.dot(Eq2))-0.5*(P1.dot(Ep1)+P2.dot(Ep2))+(1/2/alpha)*(P1.squared()+P2.squared());
        
//        if (rOO <= 2.1) UpolAtkins = 0; // value from Cummings
        
  //      if (rOO <= 2.1) UpolEquation8 = 0; // value from Cummings
        
        sum += UpolAtkins;
        
        //System.out.println("From iteration: P1 = " + P1 + ", P2 = " + P2);
        //double PolEnergyDiff = UpolAtkinsAnalytical - UpolAtkins;
        //System.out.println("Difference in Analytical and Iterative Energies = " + PolEnergyDiff);

//        System.out.println("energy = " + sum);
        return sum;
    }//end of energy

	

	
    public double energySCF(AtomSet atomsSCF) { //Atom atom1, Atom atom2, Atom atom3) {//(AtomSet atoms){
        double sumSCF = 0.0;
        double r2 = 0.0;
        double rO1O2 = 0.0;
        double rO1O3 = 0.0;
        double rO2O3 = 0.0;
      
        AtomWater4P node1 = (AtomWater4P)atomsSCF.getAtom(0);
        AtomWater4P node2 = (AtomWater4P)atomsSCF.getAtom(1);//(AtomTreeNodeWaterPPC)pair.atom1.node;
        AtomWater4P node3 = (AtomWater4P)atomsSCF.getAtom(2);

        
        IVector O1r = node1.O.getPosition();
        IVector O2r = node2.O.getPosition();
        IVector O3r = node3.O.getPosition();
        IVector H11r = node1.H1.getPosition();
        IVector H12r = node1.H2.getPosition();
        IVector H21r = node2.H1.getPosition();
        IVector H22r = node2.H2.getPosition();
        IVector H31r = node3.H1.getPosition();
        IVector H32r = node3.H2.getPosition();

        IVector M1r = node1.M.getPosition();
        IVector M2r = node2.M.getPosition();
        IVector M3r = node3.M.getPosition();
        
        
//        System.out.println("O1r coordinates before: " + O1r.x(0) + ", " + O1r.x(1) + ", " + O1r.x(2));
//        System.out.println("O2r coordinates before: " + O2r.x(0) + ", " + O2r.x(1) + ", " + O2r.x(2));
//        System.out.println("H11r coordinates before: " + H11r.x(0) + ", " + H11r.x(1) + ", " + H11r.x(2));
//        System.out.println("H12r coordinates before: " + H12r.x(0) + ", " + H12r.x(1) + ", " + H12r.x(2));
//        System.out.println("H21r coordinates before: " + H21r.x(0) + ", " + H21r.x(1) + ", " + H21r.x(2));
//        System.out.println("H22r coordinates before: " + H22r.x(0) + ", " + H22r.x(1) + ", " + H22r.x(2));
//        System.out.println("M1r coordinates before: " + M1r.x(0) + ", " + M1r.x(1) + ", " + M1r.x(2));
//        System.out.println("M2r coordinates before: " + M2r.x(0) + ", " + M2r.x(1) + ", " + M2r.x(2));
        
/*        O2r.PE(0,4); //PE(5);
        H21r.PE(0,4);
        H22r.PE(0,4);
        M2rOriginal.PE(0,4);
        
        O2r.setX(2,2); //PE(5);
        H21r.setX(2,2);
        H22r.setX(2,2);
        M2rOriginal.setX(2,2);
        //O2r.setX(2,2);
*/        
//		Dimer configuration
/*        O2r.setX(2,1.7683903);
        O2r.setX(0,2.18378015);
        O2r.setX(1,0.0);
        H21r.setX(2,1.90087324);
        H21r.setX(0,2.73561133);
        H21r.setX(1,-0.7531133);
        H22r.setX(2,1.90087324);
        H22r.setX(0,2.73561133);
        H22r.setX(1,0.7531133);
        M2rOriginal.setX(2,1.79406929);
        M2rOriginal.setX(0,2.29074084);
        M2rOriginal.setX(1,0.0);
*/

/*        O2r.setX(2,2.2);
        O2r.setX(0,0.0);
        O2r.setX(1,0.0);
        H21r.setX(2,2.767511567);
        H21r.setX(0,-0.75311329);
        H21r.setX(1,0.0);
        H22r.setX(2,2.767511567);
        H22r.setX(0,0.75311329);
        H22r.setX(1,0.0);
        M2rOriginal.setX(2,2.31);
        M2rOriginal.setX(0,0.0);
        M2rOriginal.setX(1,0.0);


        O1r.setX(2,24600.0);
        O1r.setX(0,0.0);
        O1r.setX(1,0.0);
        H11r.setX(2,24600.567511567);
        H11r.setX(0,-0.75311329);
        H11r.setX(1,0.0);
        H12r.setX(2,24600.567511567);
        H12r.setX(0,0.75311329);
        H12r.setX(1,0.0);
        M1rOriginal.setX(2,24600.11);
        M1rOriginal.setX(0,0.0);
        M1rOriginal.setX(1,0.0);

        O2r.setX(2,0.0);
        O2r.setX(0,0.0);
        O2r.setX(1,0.0);
        H21r.setX(2,0.567511567);
        H21r.setX(0,-0.75311329);
        H21r.setX(1,0.0);
        H22r.setX(2,0.567511567);
        H22r.setX(0,0.75311329);
        H22r.setX(1,0.0);
        M2rOriginal.setX(2,0.11);
        M2rOriginal.setX(0,0.0);
        M2rOriginal.setX(1,0.0);
*/
/*        O3r.setX(2,100.0);
        O3r.setX(0,0.0);
        O3r.setX(1,0.0);
        H31r.setX(2,-24600.567511567);
        H31r.setX(0,-0.75311329);
        H31r.setX(1,0.0);
        H32r.setX(2,-24600.567511567);
        H32r.setX(0,0.75311329);
        H32r.setX(1,0.0);
        M3rOriginal.setX(2,-24600.11);
        M3rOriginal.setX(0,0.0);
        M3rOriginal.setX(1,0.0);
*/
		
//      C2v geometry for dimer        
/*        O2r.setX(2,2.5);
        O2r.setX(0,0.0);
        O2r.setX(1,0.0);
        H21r.setX(2,3.085882276618295);
        H21r.setX(0,0.0);
        H21r.setX(1,0.756950327263661);
        H22r.setX(2,3.085882276618295);
        H22r.setX(0,0.0);
        H22r.setX(1,-0.756950327263661);
        M2r.setX(2,2.77);
        M2r.setX(0,0.0);
        M2r.setX(1,0.0);
*/
// Cs geometry for dimer        
/*        O2r.setX(2,1.5921074871731);
        O2r.setX(0,2.39631440952296);
        O2r.setX(1,0.0);
        H21r.setX(2,1.37452965421812);
        H21r.setX(0,2.94029779958377);
        H21r.setX(1,0.756950327263661);
        H22r.setX(2,1.37452965421812);
        H22r.setX(0,2.94029779958377);
        H22r.setX(1,-0.756950327263661);
        M2r.setX(2,1.49183817157454);
        M2r.setX(0,2.64700558278081);
        M2r.setX(1,0.0);
        
        O3r.setX(2,1000000.0);
        O3r.setX(0,0.0);
        O3r.setX(1,0.0);
        H31r.setX(2,1000000.585882276618295);
        H31r.setX(0,0.0);
        H31r.setX(1,0.756950327263661);
        H32r.setX(2,1000000.585882276618295);
        H32r.setX(0,0.0);
        H32r.setX(1,-0.756950327263661);
        M3r.setX(2,1000000.27);
        M3r.setX(0,0.0);
        M3r.setX(1,0.0);
*/        

        // moved to constructor; KMB, 7/23/07
//        final double core = 4.41; //4.41 = 2.1^2; value according to Cummings
  
        // Here is the beginning of my self-consistent solution algorithm
        // 12/26/05
   
        // moved to constructor; KMB, 7/23/07
/*        Vector Eq1old = new Vector3D();
        Vector Eq2old = new Vector3D();
        Vector Eq3old = new Vector3D();
        Vector Eq1 = new Vector3D();
        Vector Eq2 = new Vector3D();
        Vector Eq3 = new Vector3D();
        Vector Ep1old = new Vector3D();
        Vector Ep2old = new Vector3D();
        Vector Ep3old = new Vector3D();
        Vector Ep1 = new Vector3D();
        Vector Ep2 = new Vector3D();
        Vector Ep3 = new Vector3D();
        
        Vector Eq2LabFrame = new Vector3D();
        Vector Eq1LabFrame = new Vector3D();
        Vector Eq3LabFrame = new Vector3D();
        Vector Ep3LabFrame = new Vector3D();
        Vector Ep2LabFrame = new Vector3D();
        Vector Ep1LabFrame = new Vector3D();

        Vector Eon1Total = new Vector3D();
        Vector Eon2Total = new Vector3D();
        Vector Eon3Total = new Vector3D();
        
        Eon1Total.Ev1Pv2(Eq1, Ep1);
        Eon2Total.Ev1Pv2(Eq2, Ep2);
        Eon3Total.Ev1Pv2(Eq3, Ep3);

        Vector P1 = new Vector3D();
        Vector P2 = new Vector3D();
        Vector P3 = new Vector3D();
        Vector P1old = new Vector3D();
        Vector P2old = new Vector3D();
        Vector P3old = new Vector3D();
*/

        int counterSCFloop = 0; //1;
        
        // moved to constructor; KMB, 7/23/07
        // boolean counterSCFloopOK = true;
        
        int loopFailures = 0;
        
        
        // These booleans monitor whether the SCF algorithm has iteratively found
        // a solution given the criterion of Cummings, J.Chem.Phys. 2005
        boolean noSCFforP1 = true;
        boolean noSCFforP2 = true;
        boolean noSCFforP3 = true;
		
        // Need loop to check for configuration overlap between charged particles
        
//      compute O-O distance to consider bypassing the SCF loop   
        rO1O2 = Math.sqrt(O1r.Mv1Squared(O2r));
        rO1O3 = Math.sqrt(O1r.Mv1Squared(O3r));
        rO2O3 = Math.sqrt(O3r.Mv1Squared(O2r));
        
//        if (rO1O2 >= 2.1 && rO1O3 >= 2.1 && rO2O3 >= 2.1) {
//        		System.out.println("O1-O2 distance is " + rO1O2 + ", O1-O3 distance is " + rO1O3 + ", O2-O3 distance is " + rO2O3);
//        }

//        if (rO1O2 > 100 || rO1O3 > 100 || rO2O3 > 100) {
//    			return -123456789.0;
//        }
        		
        		
        if(rO1O2 <= 2.1 || rO1O3 <= 2.1 || rO2O3 <= 2.1) { // use to be 2.1 for Cummings cutoff
 
	    //    return 0; // this is for = UpolAtkins
	        return Double.POSITIVE_INFINITY;  // this is for += UpolAtkins

        }
                
                
        gamma = 12.75;
//	    double r = Math.sqrt(r2);
        double rO1O2OverSigma = rO1O2/sigma;
        double sigma2OverRO1O2sq = 1/(rO1O2OverSigma*rO1O2OverSigma);
        double rO1O3OverSigma = rO1O3/sigma;
        double sigma2OverRO1O3sq = 1/(rO1O3OverSigma*rO1O3OverSigma);
        double rO2O3OverSigma = rO2O3/sigma;
        double sigma2OverRO2O3sq = 1/(rO2O3OverSigma*rO2O3OverSigma);

        double sixOverGamma = 6/gamma;
   
        sumSCF += epsilon/(1 - sixOverGamma)*(sixOverGamma*Math.exp(gamma*(1 - rO1O2OverSigma)) - sigma2OverRO1O2sq*sigma2OverRO1O2sq*sigma2OverRO1O2sq);
        sumSCF += epsilon/(1 - sixOverGamma)*(sixOverGamma*Math.exp(gamma*(1 - rO1O3OverSigma)) - sigma2OverRO1O3sq*sigma2OverRO1O3sq*sigma2OverRO1O3sq);
        sumSCF += epsilon/(1 - sixOverGamma)*(sixOverGamma*Math.exp(gamma*(1 - rO2O3OverSigma)) - sigma2OverRO2O3sq*sigma2OverRO2O3sq*sigma2OverRO2O3sq);

        
/*        if (Math.abs(Math.sqrt(Eon1.squared())) >= 0.2083 || Math.abs(Math.sqrt(Eon2.squared())) >= 0.2083 ) {
    			System.out.println("About to return Double.Positive_Infinity");
        		System.out.println("Outside core, with bad electric field values: rOO = " + rOO + ", sum = " + sum + ", Eon1(V/A) = " + Math.sqrt(Eon1.squared())*EFconverttoVperA + ", Eon2(V/A) = " + Math.sqrt(Eon2.squared())*EFconverttoVperA);
        		return Double.POSITIVE_INFINITY;
        }
*/

        // moved to constructor; KMB, 7/23/07
/*        double sigmaM = 0.610;
        double sigmaH = 0.455;
        double sqrtHMsigmas = Math.sqrt(2*(sigmaH*sigmaH+sigmaM*sigmaM));
*/        
        // MUST INCLUDE ERF FUNCTION STUFF TO COULOMBIC ENERGY PART!
        // KMB 8/3/06
        
        double sumElecO1O2 = 0.0;
        
        r2 = H11r.Mv1Squared(H21r);  // COUNTED-2
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        	sumElecO1O2 += chargeH11*chargeH21/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaH)));
        	sumSCF += chargeH11*chargeH21/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaH)));
        	
        r2 = H11r.Mv1Squared(H22r); // COUNTED-2
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumElecO1O2 += chargeH11*chargeH22/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaH)));
        sumSCF += chargeH11*chargeH22/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaH)));
        
        r2 = H12r.Mv1Squared(H21r); // COUNTED-2
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumElecO1O2 += chargeH12*chargeH21/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaH)));
        sumSCF += chargeH12*chargeH21/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaH)));

        r2 = H12r.Mv1Squared(H22r); // COUNTED-2
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumElecO1O2 += chargeH12*chargeH22/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaH)));
        sumSCF += chargeH12*chargeH22/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaH)));

//        System.out.println("sum of all O-H terms is " + sum);
        
        r2 = M1r.Mv1Squared(H21r); // COUNTED-2
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumElecO1O2 += chargeH21*chargeM1/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/sqrtHMsigmas));
        sumSCF += chargeH21*chargeM1/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/sqrtHMsigmas));
        //System.out.println("sum is " + sum);

        r2 = M1r.Mv1Squared(H22r); // COUNTED-2
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumElecO1O2 += chargeH22*chargeM1/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/sqrtHMsigmas));
        sumSCF += chargeH22*chargeM1/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/sqrtHMsigmas));
        //System.out.println("sum is " + sum);

        r2 = M2r.Mv1Squared(H11r); // COUNTED-2
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumElecO1O2 += chargeH11*chargeM2/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/sqrtHMsigmas));
        sumSCF += chargeH11*chargeM2/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/sqrtHMsigmas));

        ///System.out.println("sum is " + sum);
        r2 = M2r.Mv1Squared(H12r); // COUNTED-2
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumElecO1O2 += chargeH12*chargeM2/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/sqrtHMsigmas));
        sumSCF += chargeH12*chargeM2/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/sqrtHMsigmas));
        //System.out.println("sum is " + sum);

        r2 = M1r.Mv1Squared(M2r); // COUNTED-2
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumElecO1O2 += chargeM1*chargeM2/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaM)));
        sumSCF += chargeM1*chargeM2/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaM)));
        //System.out.println("sum is " + sum);


        
        // EXTRA TERMS FOR 3 MOLECULE PERMUTATIONS AND COMBINATIONS FOR ELECTROSTATICS HERE
        


        r2 = H11r.Mv1Squared(H31r);  // COUNTED-2
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        	sumSCF += chargeH11*chargeH31/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaH)));

        r2 = H11r.Mv1Squared(H32r);  // COUNTED-2
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH11*chargeH32/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaH)));



        r2 = H12r.Mv1Squared(H31r);  // COUNTED-2
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH12*chargeH31/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaH)));

        r2 = H12r.Mv1Squared(H32r);  // COUNTED-2
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH12*chargeH32/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaH)));

        r2 = H21r.Mv1Squared(H31r);  // COUNTED-2
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH21*chargeH31/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaH)));

        r2 = H21r.Mv1Squared(H32r);  // COUNTED-2
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH21*chargeH32/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaH)));

        r2 = H22r.Mv1Squared(H31r);  // COUNTED-2
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH22*chargeH31/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaH)));

        r2 = H22r.Mv1Squared(H32r);  // COUNTED-2
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH22*chargeH32/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaH)));

//        System.out.println("sum of all O-H terms is " + sum);
        


        r2 = M1r.Mv1Squared(H31r);  // COUNTED-2
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH31*chargeM1/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/sqrtHMsigmas));
        //System.out.println("sum is " + sum);

        r2 = M1r.Mv1Squared(H32r);  // COUNTED-2
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH32*chargeM1/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/sqrtHMsigmas));
        //System.out.println("sum is " + sum);



        r2 = M2r.Mv1Squared(H31r);  // COUNTED-2
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH31*chargeM2/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/sqrtHMsigmas));

        ///System.out.println("sum is " + sum);
        r2 = M2r.Mv1Squared(H32r);  // COUNTED-2
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH32*chargeM2/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/sqrtHMsigmas));
        //System.out.println("sum is " + sum);

        r2 = M3r.Mv1Squared(H11r);  // COUNTED-2
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH11*chargeM3/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/sqrtHMsigmas));

        ///System.out.println("sum is " + sum);
        r2 = M3r.Mv1Squared(H12r);  // COUNTED-2
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH12*chargeM3/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/sqrtHMsigmas));
        //System.out.println("sum is " + sum);

        r2 = M3r.Mv1Squared(H21r);  // COUNTED-2
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH21*chargeM3/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/sqrtHMsigmas));

        ///System.out.println("sum is " + sum);
        r2 = M3r.Mv1Squared(H22r);  // COUNTED-2
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH22*chargeM3/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/sqrtHMsigmas));
        //System.out.println("sum is " + sum);


        r2 = M1r.Mv1Squared(M3r);  // COUNTED-2
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeM1*chargeM3/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaM)));
        //System.out.println("sum is " + sum);

        r2 = M3r.Mv1Squared(M2r);  // COUNTED-2
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeM3*chargeM2/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaM)));
        //System.out.println("sum is " + sum);

        
        /*
         * Finding the Electric fields at the center of mass of each molecule, Eqi
         * kmb, 8/7/06
         */

        // moved to constructor; KMB, 7/23/07
/*        Vector3D comW1 = new Vector3D();
        Vector3D comW2 = new Vector3D();
        Vector3D comW3 = new Vector3D();
*/        
        // How to use COM etomica code correctly? kmb, 8/7/06
/*        DataSourceCOM com = new DataSourceCOM(space);
        com.actionPerformed(atoms.getAtom(0));
        com.getCOM();
*/        
        
        // moved to constructor; KMB, 7/23/07
/*        double massH = 1.01;
        double massO = 16.0;
        double totalMass = 18.02;*/
 
        double comW1Xcomp = 0.0;
        double comW1Ycomp = 0.0;
        double comW1Zcomp = 0.0;
        
        comW1Xcomp = massH*H11r.x(0) + massO*O1r.x(0) + massH*H12r.x(0);
        comW1Ycomp = massH*H11r.x(1) + massO*O1r.x(1) + massH*H12r.x(1);
        comW1Zcomp = massH*H11r.x(2) + massO*O1r.x(2) + massH*H12r.x(2);
        
        comW1.setX(0,comW1Xcomp);
        comW1.setX(1,comW1Ycomp);
        comW1.setX(2,comW1Zcomp);
        
        comW1.Ea1Tv1(1/totalMass,comW1);

        double comW2Xcomp = 0.0;
        double comW2Ycomp = 0.0;
        double comW2Zcomp = 0.0;
        
        comW2Xcomp = massH*H21r.x(0) + massO*O2r.x(0) + massH*H22r.x(0);
        comW2Ycomp = massH*H21r.x(1) + massO*O2r.x(1) + massH*H22r.x(1);
        comW2Zcomp = massH*H21r.x(2) + massO*O2r.x(2) + massH*H22r.x(2);
        
        comW2.setX(0,comW2Xcomp);
        comW2.setX(1,comW2Ycomp);
        comW2.setX(2,comW2Zcomp);
        
        comW2.Ea1Tv1(1/totalMass,comW2);

        double comW3Xcomp = 0.0;
        double comW3Ycomp = 0.0;
        double comW3Zcomp = 0.0;
        
        comW3Xcomp = massH*H31r.x(0) + massO*O3r.x(0) + massH*H32r.x(0);
        comW3Ycomp = massH*H31r.x(1) + massO*O3r.x(1) + massH*H32r.x(1);
        comW3Zcomp = massH*H31r.x(2) + massO*O3r.x(2) + massH*H32r.x(2);
        
        comW3.setX(0,comW3Xcomp);
        comW3.setX(1,comW3Ycomp);
        comW3.setX(2,comW3Zcomp);
        
        comW3.Ea1Tv1(1/totalMass,comW3);

        
        /*
         * DOUBLE SUMMING NOW COMPLETE; KMB 8/9/06
         * 
         * THESE ELECTRIC FIELDS ARE WRONG! KMB 8/8/06
         * I HAVE NOT DONE THE DOUBLE SUM REQUIRED FOR 3 WATER
         * MOLECULES AS PER CUMMINGS EQUATION 4!
         */
        
//        double sqrtHMsigmas = Math.sqrt(2*(sigmaH*sigmaH+sigmaM*sigmaM));
        
        // moved to constructor; KMB, 7/23/07
/*        double sqrtPiHMsigmas = Math.sqrt(Math.PI*(sigmaH*sigmaH+sigmaM*sigmaM));
        double sqrtPiMMsigmas = Math.sqrt(Math.PI*(2*sigmaM*sigmaM));
*/        
        double Eq1XcompW2 = 0.0;
        double Eq1YcompW2 = 0.0;
        double Eq1ZcompW2 = 0.0;
        double Eq1XcompW3 = 0.0;
        double Eq1YcompW3 = 0.0;
        double Eq1ZcompW3 = 0.0;

        
        double comW1toH21 = Math.sqrt(comW1.Mv1Squared(H21r));
        double comW1toH22 = Math.sqrt(comW1.Mv1Squared(H22r));
        double comW1toM2 = Math.sqrt(comW1.Mv1Squared(M2r));

        double comW1toH31 = Math.sqrt(comW1.Mv1Squared(H31r));
        double comW1toH32 = Math.sqrt(comW1.Mv1Squared(H32r));
        double comW1toM3 = Math.sqrt(comW1.Mv1Squared(M3r));

        // Contributions to sum from water#2
        Eq1XcompW2 += chargeH21*(comW1.x(0)-H21r.x(0))/(comW1toH21*comW1toH21*comW1toH21)*((1-SpecialFunctions.erfc(comW1toH21/sqrtHMsigmas))-Math.sqrt(2)*comW1toH21/sqrtPiHMsigmas*Math.exp(-comW1toH21*comW1toH21/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq1XcompW2 += chargeH22*(comW1.x(0)-H22r.x(0))/(comW1toH22*comW1toH22*comW1toH22)*((1-SpecialFunctions.erfc(comW1toH22/sqrtHMsigmas))-Math.sqrt(2)*comW1toH22/sqrtPiHMsigmas*Math.exp(-comW1toH22*comW1toH22/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq1XcompW2 += chargeM2*(comW1.x(0)-M2r.x(0))/(comW1toM2*comW1toM2*comW1toM2)*((1-SpecialFunctions.erfc(comW1toM2/(2*sigmaM)))-Math.sqrt(2)*comW1toM2/sqrtPiMMsigmas*Math.exp(-comW1toM2*comW1toM2/(4*sigmaM*sigmaM)));

        Eq1YcompW2 += chargeH21*(comW1.x(1)-H21r.x(1))/(comW1toH21*comW1toH21*comW1toH21)*((1-SpecialFunctions.erfc(comW1toH21/sqrtHMsigmas))-Math.sqrt(2)*comW1toH21/sqrtPiHMsigmas*Math.exp(-comW1toH21*comW1toH21/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq1YcompW2 += chargeH22*(comW1.x(1)-H22r.x(1))/(comW1toH22*comW1toH22*comW1toH22)*((1-SpecialFunctions.erfc(comW1toH22/sqrtHMsigmas))-Math.sqrt(2)*comW1toH22/sqrtPiHMsigmas*Math.exp(-comW1toH22*comW1toH22/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq1YcompW2 += chargeM2*(comW1.x(1)-M2r.x(1))/(comW1toM2*comW1toM2*comW1toM2)*((1-SpecialFunctions.erfc(comW1toM2/(2*sigmaM)))-Math.sqrt(2)*comW1toM2/sqrtPiMMsigmas*Math.exp(-comW1toM2*comW1toM2/(4*sigmaM*sigmaM)));

        Eq1ZcompW2 += chargeH21*(comW1.x(2)-H21r.x(2))/(comW1toH21*comW1toH21*comW1toH21)*((1-SpecialFunctions.erfc(comW1toH21/sqrtHMsigmas))-Math.sqrt(2)*comW1toH21/sqrtPiHMsigmas*Math.exp(-comW1toH21*comW1toH21/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq1ZcompW2 += chargeH22*(comW1.x(2)-H22r.x(2))/(comW1toH22*comW1toH22*comW1toH22)*((1-SpecialFunctions.erfc(comW1toH22/sqrtHMsigmas))-Math.sqrt(2)*comW1toH22/sqrtPiHMsigmas*Math.exp(-comW1toH22*comW1toH22/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq1ZcompW2 += chargeM2*(comW1.x(2)-M2r.x(2))/(comW1toM2*comW1toM2*comW1toM2)*((1-SpecialFunctions.erfc(comW1toM2/(2*sigmaM)))-Math.sqrt(2)*comW1toM2/sqrtPiMMsigmas*Math.exp(-comW1toM2*comW1toM2/(4*sigmaM*sigmaM)));

        
        // Contributions to sum from water#3
        Eq1XcompW3 += chargeH31*(comW1.x(0)-H31r.x(0))/(comW1toH31*comW1toH31*comW1toH31)*((1-SpecialFunctions.erfc(comW1toH31/sqrtHMsigmas))-Math.sqrt(2)*comW1toH31/sqrtPiHMsigmas*Math.exp(-comW1toH31*comW1toH31/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq1XcompW3 += chargeH32*(comW1.x(0)-H32r.x(0))/(comW1toH32*comW1toH32*comW1toH32)*((1-SpecialFunctions.erfc(comW1toH32/sqrtHMsigmas))-Math.sqrt(2)*comW1toH32/sqrtPiHMsigmas*Math.exp(-comW1toH32*comW1toH32/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq1XcompW3 += chargeM3*(comW1.x(0)-M3r.x(0))/(comW1toM3*comW1toM3*comW1toM3)*((1-SpecialFunctions.erfc(comW1toM3/(2*sigmaM)))-Math.sqrt(2)*comW1toM3/sqrtPiMMsigmas*Math.exp(-comW1toM3*comW1toM3/(4*sigmaM*sigmaM)));

        Eq1YcompW3 += chargeH31*(comW1.x(1)-H31r.x(1))/(comW1toH31*comW1toH31*comW1toH31)*((1-SpecialFunctions.erfc(comW1toH31/sqrtHMsigmas))-Math.sqrt(2)*comW1toH31/sqrtPiHMsigmas*Math.exp(-comW1toH31*comW1toH31/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq1YcompW3 += chargeH32*(comW1.x(1)-H32r.x(1))/(comW1toH32*comW1toH32*comW1toH32)*((1-SpecialFunctions.erfc(comW1toH32/sqrtHMsigmas))-Math.sqrt(2)*comW1toH32/sqrtPiHMsigmas*Math.exp(-comW1toH32*comW1toH32/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq1YcompW3 += chargeM3*(comW1.x(1)-M3r.x(1))/(comW1toM3*comW1toM3*comW1toM3)*((1-SpecialFunctions.erfc(comW1toM3/(2*sigmaM)))-Math.sqrt(2)*comW1toM3/sqrtPiMMsigmas*Math.exp(-comW1toM3*comW1toM3/(4*sigmaM*sigmaM)));

        Eq1ZcompW3 += chargeH31*(comW1.x(2)-H31r.x(2))/(comW1toH31*comW1toH31*comW1toH31)*((1-SpecialFunctions.erfc(comW1toH31/sqrtHMsigmas))-Math.sqrt(2)*comW1toH31/sqrtPiHMsigmas*Math.exp(-comW1toH31*comW1toH31/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq1ZcompW3 += chargeH32*(comW1.x(2)-H32r.x(2))/(comW1toH32*comW1toH32*comW1toH32)*((1-SpecialFunctions.erfc(comW1toH32/sqrtHMsigmas))-Math.sqrt(2)*comW1toH32/sqrtPiHMsigmas*Math.exp(-comW1toH32*comW1toH32/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq1ZcompW3 += chargeM3*(comW1.x(2)-M3r.x(2))/(comW1toM3*comW1toM3*comW1toM3)*((1-SpecialFunctions.erfc(comW1toM3/(2*sigmaM)))-Math.sqrt(2)*comW1toM3/sqrtPiMMsigmas*Math.exp(-comW1toM3*comW1toM3/(4*sigmaM*sigmaM)));
        
        
        Eq1.setX(0,Eq1XcompW2+Eq1XcompW3);
        Eq1.setX(1,Eq1YcompW2+Eq1YcompW3);
        Eq1.setX(2,Eq1ZcompW2+Eq1ZcompW3);

                
        double Eq2XcompW1 = 0.0;
        double Eq2YcompW1 = 0.0;
        double Eq2ZcompW1 = 0.0;
        double Eq2XcompW3 = 0.0;
        double Eq2YcompW3 = 0.0;
        double Eq2ZcompW3 = 0.0;

        
        double comW2toH11 = Math.sqrt(comW2.Mv1Squared(H11r));
        double comW2toH12 = Math.sqrt(comW2.Mv1Squared(H12r));
        double comW2toM1 = Math.sqrt(comW2.Mv1Squared(M1r));
        double comW2toH31 = Math.sqrt(comW2.Mv1Squared(H31r));
        double comW2toH32 = Math.sqrt(comW2.Mv1Squared(H32r));
        double comW2toM3 = Math.sqrt(comW2.Mv1Squared(M3r));

        // Contributions to sum from water molecule#1
        Eq2XcompW1 += chargeH11*(comW2.x(0)-H11r.x(0))/(comW2toH11*comW2toH11*comW2toH11)*((1-SpecialFunctions.erfc(comW2toH11/sqrtHMsigmas))-Math.sqrt(2)*comW2toH11/sqrtPiHMsigmas*Math.exp(-comW2toH11*comW2toH11/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq2XcompW1 += chargeH12*(comW2.x(0)-H12r.x(0))/(comW2toH12*comW2toH12*comW2toH12)*((1-SpecialFunctions.erfc(comW2toH12/sqrtHMsigmas))-Math.sqrt(2)*comW2toH12/sqrtPiHMsigmas*Math.exp(-comW2toH12*comW2toH12/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq2XcompW1 += chargeM1*(comW2.x(0)-M1r.x(0))/(comW2toM1*comW2toM1*comW2toM1)*((1-SpecialFunctions.erfc(comW2toM1/(2*sigmaM)))-Math.sqrt(2)*comW2toM1/sqrtPiMMsigmas*Math.exp(-comW2toM1*comW2toM1/(4*sigmaM*sigmaM)));

        Eq2YcompW1 += chargeH11*(comW2.x(1)-H11r.x(1))/(comW2toH11*comW2toH11*comW2toH11)*((1-SpecialFunctions.erfc(comW2toH11/sqrtHMsigmas))-Math.sqrt(2)*comW2toH11/sqrtPiHMsigmas*Math.exp(-comW2toH11*comW2toH11/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq2YcompW1 += chargeH12*(comW2.x(1)-H12r.x(1))/(comW2toH12*comW2toH12*comW2toH12)*((1-SpecialFunctions.erfc(comW2toH12/sqrtHMsigmas))-Math.sqrt(2)*comW2toH12/sqrtPiHMsigmas*Math.exp(-comW2toH12*comW2toH12/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq2YcompW1 += chargeM1*(comW2.x(1)-M1r.x(1))/(comW2toM1*comW2toM1*comW2toM1)*((1-SpecialFunctions.erfc(comW2toM1/(2*sigmaM)))-Math.sqrt(2)*comW2toM1/sqrtPiMMsigmas*Math.exp(-comW2toM1*comW2toM1/(4*sigmaM*sigmaM)));

        Eq2ZcompW1 += chargeH11*(comW2.x(2)-H11r.x(2))/(comW2toH11*comW2toH11*comW2toH11)*((1-SpecialFunctions.erfc(comW2toH11/sqrtHMsigmas))-Math.sqrt(2)*comW2toH11/sqrtPiHMsigmas*Math.exp(-comW2toH11*comW2toH11/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq2ZcompW1 += chargeH12*(comW2.x(2)-H12r.x(2))/(comW2toH12*comW2toH12*comW2toH12)*((1-SpecialFunctions.erfc(comW2toH12/sqrtHMsigmas))-Math.sqrt(2)*comW2toH12/sqrtPiHMsigmas*Math.exp(-comW2toH12*comW2toH12/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq2ZcompW1 += chargeM1*(comW2.x(2)-M1r.x(2))/(comW2toM1*comW2toM1*comW2toM1)*((1-SpecialFunctions.erfc(comW2toM1/(2*sigmaM)))-Math.sqrt(2)*comW2toM1/sqrtPiMMsigmas*Math.exp(-comW2toM1*comW2toM1/(4*sigmaM*sigmaM)));


        // Contributions to sum from water molecule#3
        Eq2XcompW3 += chargeH31*(comW2.x(0)-H31r.x(0))/(comW2toH31*comW2toH31*comW2toH31)*((1-SpecialFunctions.erfc(comW2toH31/sqrtHMsigmas))-Math.sqrt(2)*comW2toH31/sqrtPiHMsigmas*Math.exp(-comW2toH31*comW2toH31/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq2XcompW3 += chargeH32*(comW2.x(0)-H32r.x(0))/(comW2toH32*comW2toH32*comW2toH32)*((1-SpecialFunctions.erfc(comW2toH32/sqrtHMsigmas))-Math.sqrt(2)*comW2toH32/sqrtPiHMsigmas*Math.exp(-comW2toH32*comW2toH32/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq2XcompW3 += chargeM3*(comW2.x(0)-M3r.x(0))/(comW2toM3*comW2toM3*comW2toM3)*((1-SpecialFunctions.erfc(comW2toM3/(2*sigmaM)))-Math.sqrt(2)*comW2toM3/sqrtPiMMsigmas*Math.exp(-comW2toM3*comW2toM3/(4*sigmaM*sigmaM)));

        Eq2YcompW3 += chargeH31*(comW2.x(1)-H31r.x(1))/(comW2toH31*comW2toH31*comW2toH31)*((1-SpecialFunctions.erfc(comW2toH31/sqrtHMsigmas))-Math.sqrt(2)*comW2toH31/sqrtPiHMsigmas*Math.exp(-comW2toH31*comW2toH31/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq2YcompW3 += chargeH32*(comW2.x(1)-H32r.x(1))/(comW2toH32*comW2toH32*comW2toH32)*((1-SpecialFunctions.erfc(comW2toH32/sqrtHMsigmas))-Math.sqrt(2)*comW2toH32/sqrtPiHMsigmas*Math.exp(-comW2toH32*comW2toH32/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq2YcompW3 += chargeM3*(comW2.x(1)-M3r.x(1))/(comW2toM3*comW2toM3*comW2toM3)*((1-SpecialFunctions.erfc(comW2toM3/(2*sigmaM)))-Math.sqrt(2)*comW2toM3/sqrtPiMMsigmas*Math.exp(-comW2toM3*comW2toM3/(4*sigmaM*sigmaM)));

        Eq2ZcompW3 += chargeH31*(comW2.x(2)-H31r.x(2))/(comW2toH31*comW2toH31*comW2toH31)*((1-SpecialFunctions.erfc(comW2toH31/sqrtHMsigmas))-Math.sqrt(2)*comW2toH31/sqrtPiHMsigmas*Math.exp(-comW2toH31*comW2toH31/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq2ZcompW3 += chargeH32*(comW2.x(2)-H32r.x(2))/(comW2toH32*comW2toH32*comW2toH32)*((1-SpecialFunctions.erfc(comW2toH32/sqrtHMsigmas))-Math.sqrt(2)*comW2toH32/sqrtPiHMsigmas*Math.exp(-comW2toH32*comW2toH32/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq2ZcompW3 += chargeM3*(comW2.x(2)-M3r.x(2))/(comW2toM3*comW2toM3*comW2toM3)*((1-SpecialFunctions.erfc(comW2toM3/(2*sigmaM)))-Math.sqrt(2)*comW2toM3/sqrtPiMMsigmas*Math.exp(-comW2toM3*comW2toM3/(4*sigmaM*sigmaM)));
        
        
        Eq2.setX(0,Eq2XcompW1+Eq2XcompW3);
        Eq2.setX(1,Eq2YcompW1+Eq2YcompW3);
        Eq2.setX(2,Eq2ZcompW1+Eq2ZcompW3);


        // Find Eq3
        double Eq3XcompW1	 = 0.0;
        double Eq3YcompW1 = 0.0;
        double Eq3ZcompW1 = 0.0;
        double Eq3XcompW2 = 0.0;
        double Eq3YcompW2 = 0.0;
        double Eq3ZcompW2 = 0.0;

        
        double comW3toH11 = Math.sqrt(comW3.Mv1Squared(H11r));
        double comW3toH12 = Math.sqrt(comW3.Mv1Squared(H12r));
        double comW3toM1 = Math.sqrt(comW3.Mv1Squared(M1r));

        double comW3toH21 = Math.sqrt(comW3.Mv1Squared(H21r));
        double comW3toH22 = Math.sqrt(comW3.Mv1Squared(H22r));
        double comW3toM2 = Math.sqrt(comW3.Mv1Squared(M2r));

        // Contributions to sum from water molecule#1       
        Eq3XcompW1 += chargeH11*(comW3.x(0)-H11r.x(0))/(comW3toH11*comW3toH11*comW3toH11)*((1-SpecialFunctions.erfc(comW3toH11/sqrtHMsigmas))-Math.sqrt(2)*comW3toH11/sqrtPiHMsigmas*Math.exp(-comW3toH11*comW3toH11/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq3XcompW1 += chargeH12*(comW3.x(0)-H12r.x(0))/(comW3toH12*comW3toH12*comW3toH12)*((1-SpecialFunctions.erfc(comW3toH12/sqrtHMsigmas))-Math.sqrt(2)*comW3toH12/sqrtPiHMsigmas*Math.exp(-comW3toH12*comW3toH12/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq3XcompW1 += chargeM1*(comW3.x(0)-M1r.x(0))/(comW3toM1*comW3toM1*comW3toM1)*((1-SpecialFunctions.erfc(comW3toM1/(2*sigmaM)))-Math.sqrt(2)*comW3toM1/sqrtPiMMsigmas*Math.exp(-comW3toM1*comW3toM1/(4*sigmaM*sigmaM)));

        Eq3YcompW1 += chargeH11*(comW3.x(1)-H11r.x(1))/(comW3toH11*comW3toH11*comW3toH11)*((1-SpecialFunctions.erfc(comW3toH11/sqrtHMsigmas))-Math.sqrt(2)*comW3toH11/sqrtPiHMsigmas*Math.exp(-comW3toH11*comW3toH11/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq3YcompW1 += chargeH12*(comW3.x(1)-H12r.x(1))/(comW3toH12*comW3toH12*comW3toH12)*((1-SpecialFunctions.erfc(comW3toH12/sqrtHMsigmas))-Math.sqrt(2)*comW3toH12/sqrtPiHMsigmas*Math.exp(-comW3toH12*comW3toH12/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq3YcompW1 += chargeM1*(comW3.x(1)-M1r.x(1))/(comW3toM1*comW3toM1*comW3toM1)*((1-SpecialFunctions.erfc(comW3toM1/(2*sigmaM)))-Math.sqrt(2)*comW3toM1/sqrtPiMMsigmas*Math.exp(-comW3toM1*comW3toM1/(4*sigmaM*sigmaM)));

        Eq3ZcompW1 += chargeH11*(comW3.x(2)-H11r.x(2))/(comW3toH11*comW3toH11*comW3toH11)*((1-SpecialFunctions.erfc(comW3toH11/sqrtHMsigmas))-Math.sqrt(2)*comW3toH11/sqrtPiHMsigmas*Math.exp(-comW3toH11*comW3toH11/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq3ZcompW1 += chargeH12*(comW3.x(2)-H12r.x(2))/(comW3toH12*comW3toH12*comW3toH12)*((1-SpecialFunctions.erfc(comW3toH12/sqrtHMsigmas))-Math.sqrt(2)*comW3toH12/sqrtPiHMsigmas*Math.exp(-comW3toH12*comW3toH12/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq3ZcompW1 += chargeM1*(comW3.x(2)-M1r.x(2))/(comW3toM1*comW3toM1*comW3toM1)*((1-SpecialFunctions.erfc(comW3toM1/(2*sigmaM)))-Math.sqrt(2)*comW3toM1/sqrtPiMMsigmas*Math.exp(-comW3toM1*comW3toM1/(4*sigmaM*sigmaM)));

        
        // Contributions to sum from water molecule#2
        Eq3XcompW2 += chargeH21*(comW3.x(0)-H21r.x(0))/(comW3toH21*comW3toH21*comW3toH21)*((1-SpecialFunctions.erfc(comW3toH21/sqrtHMsigmas))-Math.sqrt(2)*comW3toH21/sqrtPiHMsigmas*Math.exp(-comW3toH21*comW3toH21/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq3XcompW2 += chargeH22*(comW3.x(0)-H22r.x(0))/(comW3toH22*comW3toH22*comW3toH22)*((1-SpecialFunctions.erfc(comW3toH22/sqrtHMsigmas))-Math.sqrt(2)*comW3toH22/sqrtPiHMsigmas*Math.exp(-comW3toH22*comW3toH22/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq3XcompW2 += chargeM2*(comW3.x(0)-M2r.x(0))/(comW3toM2*comW3toM2*comW3toM2)*((1-SpecialFunctions.erfc(comW3toM2/(2*sigmaM)))-Math.sqrt(2)*comW3toM2/sqrtPiMMsigmas*Math.exp(-comW3toM2*comW3toM2/(4*sigmaM*sigmaM)));

        Eq3YcompW2 += chargeH21*(comW3.x(1)-H21r.x(1))/(comW3toH21*comW3toH21*comW3toH21)*((1-SpecialFunctions.erfc(comW3toH21/sqrtHMsigmas))-Math.sqrt(2)*comW3toH21/sqrtPiHMsigmas*Math.exp(-comW3toH21*comW3toH21/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq3YcompW2 += chargeH22*(comW3.x(1)-H22r.x(1))/(comW3toH22*comW3toH22*comW3toH22)*((1-SpecialFunctions.erfc(comW3toH22/sqrtHMsigmas))-Math.sqrt(2)*comW3toH22/sqrtPiHMsigmas*Math.exp(-comW3toH22*comW3toH22/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq3YcompW2 += chargeM2*(comW3.x(1)-M2r.x(1))/(comW3toM2*comW3toM2*comW3toM2)*((1-SpecialFunctions.erfc(comW3toM2/(2*sigmaM)))-Math.sqrt(2)*comW3toM2/sqrtPiMMsigmas*Math.exp(-comW3toM2*comW3toM2/(4*sigmaM*sigmaM)));

        Eq3ZcompW2 += chargeH21*(comW3.x(2)-H21r.x(2))/(comW3toH21*comW3toH21*comW3toH21)*((1-SpecialFunctions.erfc(comW3toH21/sqrtHMsigmas))-Math.sqrt(2)*comW3toH21/sqrtPiHMsigmas*Math.exp(-comW3toH21*comW3toH21/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq3ZcompW2 += chargeH22*(comW3.x(2)-H22r.x(2))/(comW3toH22*comW3toH22*comW3toH22)*((1-SpecialFunctions.erfc(comW3toH22/sqrtHMsigmas))-Math.sqrt(2)*comW3toH22/sqrtPiHMsigmas*Math.exp(-comW3toH22*comW3toH22/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq3ZcompW2 += chargeM2*(comW3.x(2)-M2r.x(2))/(comW3toM2*comW3toM2*comW3toM2)*((1-SpecialFunctions.erfc(comW3toM2/(2*sigmaM)))-Math.sqrt(2)*comW3toM2/sqrtPiMMsigmas*Math.exp(-comW3toM2*comW3toM2/(4*sigmaM*sigmaM)));
        
        
        Eq3.setX(0,Eq3XcompW1+Eq3XcompW2);
        Eq3.setX(1,Eq3YcompW1+Eq3YcompW2);
        Eq3.setX(2,Eq3ZcompW1+Eq3ZcompW2);

        
        /*
         * Finding the tensor used to relate the induced dipole moment Pi with the induced electric field Epi.
         * kmb, 8/9/06
         */
        
        double r12 = Math.sqrt(comW1.Mv1Squared(comW2));
        double r13 = Math.sqrt(comW1.Mv1Squared(comW3));
        double r23 = Math.sqrt(comW2.Mv1Squared(comW3));
        
        Vector3D r12Vector = new Vector3D();
        Vector3D r13Vector = new Vector3D();
        Vector3D r23Vector = new Vector3D();
        
        r12Vector.Ev1Mv2(comW1,comW2);  // is this the correct direction? kmb, 8/7/06 / Direction doesn't matter; kmb, 8/10/06
        r13Vector.Ev1Mv2(comW1,comW3);  // is this the correct direction? kmb, 8/7/06 / Direction doesn't matter; kmb, 8/10/06
        r23Vector.Ev1Mv2(comW2,comW3);  // is this the correct direction? kmb, 8/7/06 / Direction doesn't matter; kmb, 8/10/06
        
        double f12 = (1-SpecialFunctions.erfc(r12/(2*sigmaM)))-(r12/(sigmaM*Math.sqrt(Math.PI)) + (r12*r12*r12)/(6*Math.sqrt(Math.PI)*sigmaM*sigmaM*sigmaM))*Math.exp(-r12*r12/(4*sigmaM*sigmaM));
        double f13 = (1-SpecialFunctions.erfc(r13/(2*sigmaM)))-(r13/(sigmaM*Math.sqrt(Math.PI)) + (r13*r13*r13)/(6*Math.sqrt(Math.PI)*sigmaM*sigmaM*sigmaM))*Math.exp(-r13*r13/(4*sigmaM*sigmaM));
        double f23 = (1-SpecialFunctions.erfc(r23/(2*sigmaM)))-(r23/(sigmaM*Math.sqrt(Math.PI)) + (r23*r23*r23)/(6*Math.sqrt(Math.PI)*sigmaM*sigmaM*sigmaM))*Math.exp(-r23*r23/(4*sigmaM*sigmaM));
        
        double g12 = (1-SpecialFunctions.erfc(r12/(2*sigmaM)))-(r12/(sigmaM*Math.sqrt(Math.PI)))*Math.exp(-r12*r12/(4*sigmaM*sigmaM));
        double g13 = (1-SpecialFunctions.erfc(r13/(2*sigmaM)))-(r13/(sigmaM*Math.sqrt(Math.PI)))*Math.exp(-r13*r13/(4*sigmaM*sigmaM));
        double g23 = (1-SpecialFunctions.erfc(r23/(2*sigmaM)))-(r23/(sigmaM*Math.sqrt(Math.PI)))*Math.exp(-r23*r23/(4*sigmaM*sigmaM));
        
        // Filling the unit matrix I
        
/*        double[][] I = new double[3][3];
                
        int i = 0;
        int j = 0;
        
        while (i < 3) {
        		while (j < 3) {
        			I[i][j] = 1;
            		j = j + 1;
        		}
        		i = i + 1;
        }
*/        
        Tensor3D I = new Tensor3D();
        
        I.E(1);
        
//        double[][] T12 = new double[3][3];

        Tensor3D T12 = new Tensor3D();
        
        T12.PEv1v2(r12Vector,r12Vector);
        T12.TE(3*f12/(r12*r12));
        
        I.TE(g12);
        
        T12.ME(I);
        T12.TE(1/(r12*r12*r12));
        
        // T12 = T21, so I can get by for now in the case of B2!

        
        I.E(1);
        
//        double[][] T12 = new double[3][3];

        Tensor3D T13 = new Tensor3D();
        
        T13.PEv1v2(r13Vector,r13Vector);
        T13.TE(3*f13/(r13*r13));
        
        I.TE(g13);
        
        T13.ME(I);
        T13.TE(1/(r13*r13*r13));
        
        
        I.E(1);
        
//        double[][] T12 = new double[3][3];

        Tensor3D T23 = new Tensor3D();
        
        T23.PEv1v2(r23Vector,r23Vector);
        T23.TE(3*f23/(r23*r23));
        
        I.TE(g23);
        
        T23.ME(I);
        T23.TE(1/(r23*r23*r23));


        // Now distribute the elements of the tensor into 3 separate "row" vectors
        // so I can do dot products with etomica math methods
        // kmb, 8/7/06
        
        Vector3D T12row1 = new Vector3D();
        Vector3D T12row2 = new Vector3D();
        Vector3D T12row3 = new Vector3D();
        
        T12row1.setX(0,T12.component(0,0));
        T12row1.setX(1,T12.component(0,1));
        T12row1.setX(2,T12.component(0,2));
        T12row2.setX(0,T12.component(1,0));
        T12row2.setX(1,T12.component(1,1));
        T12row2.setX(2,T12.component(1,2));
        T12row3.setX(0,T12.component(2,0));
        T12row3.setX(1,T12.component(2,1));
        T12row3.setX(2,T12.component(2,2));
        

        Vector3D T13row1 = new Vector3D();
        Vector3D T13row2 = new Vector3D();
        Vector3D T13row3 = new Vector3D();
        
        T13row1.setX(0,T13.component(0,0));
        T13row1.setX(1,T13.component(0,1));
        T13row1.setX(2,T13.component(0,2));
        T13row2.setX(0,T13.component(1,0));
        T13row2.setX(1,T13.component(1,1));
        T13row2.setX(2,T13.component(1,2));
        T13row3.setX(0,T13.component(2,0));
        T13row3.setX(1,T13.component(2,1));
        T13row3.setX(2,T13.component(2,2));

        
        Vector3D T23row1 = new Vector3D();
        Vector3D T23row2 = new Vector3D();
        Vector3D T23row3 = new Vector3D();
        
        T23row1.setX(0,T23.component(0,0));
        T23row1.setX(1,T23.component(0,1));
        T23row1.setX(2,T23.component(0,2));
        T23row2.setX(0,T23.component(1,0));
        T23row2.setX(1,T23.component(1,1));
        T23row2.setX(2,T23.component(1,2));
        T23row3.setX(0,T23.component(2,0));
        T23row3.setX(1,T23.component(2,1));
        T23row3.setX(2,T23.component(2,2));

        
        //Begin analytical solution for Upol
        
        double alphaPol = 1.444;
        double alphaPol2 = alphaPol*alphaPol;
        
        I.E(0);
        I.setComponent(0,0,1);
        I.setComponent(1,1,1);
        I.setComponent(2,2,1);
        Tensor3D A1 = new Tensor3D();
        Tensor3D A2 = new Tensor3D();
        Tensor3D C1 = new Tensor3D();
        Tensor3D C2 = new Tensor3D();
//        Tensor3D II = new Tensor3D();
        Tensor3D T23T23 = new Tensor3D();
        Tensor3D T13T13 = new Tensor3D();
        Tensor3D T13T23 = new Tensor3D();
        
        A1.E(I);
//        A1.TE(I);
        T23T23.E(T23);
        T23T23.TE(T23);
        T23T23.TE(alphaPol2);
        T13T13.E(T13);
        T13T13.TE(T13);
        T13T13.TE(alphaPol2);
        T13T23.E(T13);
        T13T23.TE(T23);
        T13T23.TE(alphaPol);
        A2.E(I);
  //      A2.TE(I);
        
        A1.ME(T13T13);
        A2.ME(T23T23);
        
        C1.E(T12);
        C1.TE(alphaPol);
        C1.PE(T13T23);
        C2.E(C1);
        
        Vector3D B1 = new Vector3D();
        Vector3D B2 = new Vector3D();
                
        Vector3D T23Eq3 = new Vector3D();
        Vector3D T13Eq3 = new Vector3D();
        
        T23Eq3.setX(0,T23row1.dot(Eq3));
        T23Eq3.setX(1,T23row2.dot(Eq3));
        T23Eq3.setX(2,T23row3.dot(Eq3));
        T23Eq3.TE(alphaPol2);
                        
        T13Eq3.setX(0,T13row1.dot(Eq3));
        T13Eq3.setX(1,T13row2.dot(Eq3));
        T13Eq3.setX(2,T13row3.dot(Eq3));
        T13Eq3.TE(alphaPol2);
        
        B1.E(Eq1);
        B1.TE(alphaPol);
        B1.PE(T13Eq3);

        B2.E(Eq2);
        B2.TE(alphaPol);
        B2.PE(T23Eq3);
        
        Tensor3D inverseA1 = new Tensor3D();
        Tensor3D inverseA2 = new Tensor3D();
        
        inverseA1.E(A1);
        inverseA1.inverse();
        
        inverseA2.E(A2);
        inverseA2.inverse();
        
        //use transform method in Tensor3D!
        
        Tensor3D C1A2inverseC2 = new Tensor3D();
        Tensor3D C1A2inverse = new Tensor3D();
        
        C1A2inverseC2.E(C2);
        C1A2inverseC2.TE(inverseA2);
        C1A2inverseC2.TE(C1);
        //C2A1inverseC1.TE(alphaPol);
        
        Vector3D C1A2inverseB2 = new Vector3D();
        C1A2inverseB2.E(B2);
//        C2A1inverseB1.TE(alphaPol);
        
        C1A2inverse.E(C1);
        C1A2inverse.TE(inverseA2);
        C1A2inverse.transform(C1A2inverseB2);
        
        Tensor3D inverseA4 = new Tensor3D();
        Vector3D D4 = new Vector3D();
        
        inverseA4.E(A1);
        inverseA4.ME(C1A2inverseC2);
        inverseA4.inverse();
        
        D4.E(B1);
        D4.PE(C1A2inverseB2);
        
        P1.E(D4);
        inverseA4.transform(P1);
        
//        double alphaPol = 1.444;
//        P1.TE(alphaPol);
        
        // Now find P2
        
        Tensor3D A2inverseC2 = new Tensor3D();
        A2inverseC2.E(C2);
        A2inverseC2.TE(inverseA2);
        
        Vector3D B4 = new Vector3D();
        
        B4.E(P1);
        A2inverseC2.transform(B4);
        
        Vector3D A4 = new Vector3D();
        
        A4.E(B2);
        inverseA2.transform(A4);
        
        P2.E(A4);
        P2.PE(B4);
        
//        P2.TE(alphaPol);
        
        // Now get P3
        
        Vector3D T12P1 = new Vector3D();
        Vector3D T12P2 = new Vector3D();
        Vector3D T13P1 = new Vector3D();
        Vector3D T13P3 = new Vector3D();
        Vector3D T23P2 = new Vector3D();
        Vector3D T23P3 = new Vector3D();
        
        T13P1.E(P1);
        T13.transform(T13P1);
        
        T23P2.E(P2);
        T23.transform(T23P2);
        
        P3.E(Eq3);
        P3.PE(T13P1);
        P3.PE(T23P2);
        
        P3.TE(alphaPol);
        
//        System.out.println("From analytical: P1 = " + P1 + ", P2 = " + P2 + ", P3 = " + P3);
        
        // Now find Ep1, Ep2, Ep3
        
        Ep1.E(T12P2);
        Ep1.PE(T13P3);
        Ep2.E(T12P1);
        Ep2.PE(T23P3);
        Ep3.E(T13P1);
        Ep3.PE(T23P2);
        
//        double UpolAtkinsAnalytical = -0.5*(P1.dot(Eq1)+P2.dot(Eq2)+P3.dot(Eq3));
        
        
/*        // Set the induced dipole moments equal to 10% of the permanent dipole value
        // kmb, 8/7/06
        P1.setX(0,14.3952507082236);
        P1.setX(1,14.3952507082236);
        P1.setX(2,14.3952507082236);
        P2.setX(0,14.3952507082236);
        P2.setX(1,14.3952507082236);
        P2.setX(2,14.3952507082236);
        P3.setX(0,14.3952507082236);
        P3.setX(1,14.3952507082236);
        P3.setX(2,14.3952507082236);
        
        P1old.E(P1);
        P2old.E(P2);
        P3old.E(P3);

        
        // kmb add these here and initialize, 8/22/06
        // trying to solve convergence problem between dimer and trimer
        
        double deltaP1 = 1.0;
        double deltaP2 = 1.0;
        double deltaP3 = 1.0;

        
        while (noSCFforP1 || noSCFforP2 || noSCFforP3) {

    			// First calculate Ep1, based upon guess for P2 and P3
        	
        		Ep1.setX(0,T12row1.dot(P2)+T13row1.dot(P3));
        		Ep1.setX(1,T12row2.dot(P2)+T13row2.dot(P3));
        		Ep1.setX(2,T12row3.dot(P2)+T13row3.dot(P3));
        		
        		// Now calculate new P1 from the value of Ep1
        		
        		//double alphaPol = 1.444;
        		
        		P1.setX(0,alphaPol*(Eq1.x(0) + Ep1.x(0)));
        		P1.setX(1,alphaPol*(Eq1.x(1) + Ep1.x(1)));
        		P1.setX(2,alphaPol*(Eq1.x(2) + Ep1.x(2)));

        		// Next calculate Ep2
        		
        		Ep2.setX(0,T12row1.dot(P1)+T23row1.dot(P3));
        		Ep2.setX(1,T12row2.dot(P1)+T23row2.dot(P3));
        		Ep2.setX(2,T12row3.dot(P1)+T23row3.dot(P3));
        		
        		// Now calculate new P2
        		
        		P2.setX(0,alphaPol*(Eq2.x(0) + Ep2.x(0)));
        		P2.setX(1,alphaPol*(Eq2.x(1) + Ep2.x(1)));
        		P2.setX(2,alphaPol*(Eq2.x(2) + Ep2.x(2)));
        		
        		// Next calculate Ep3
        		
        		Ep3.setX(0,T13row1.dot(P1)+T23row1.dot(P2));
        		Ep3.setX(1,T13row2.dot(P1)+T23row2.dot(P2));
        		Ep3.setX(2,T13row3.dot(P1)+T23row3.dot(P2));
        		
        		// Now calculate new P3
        		
        		P3.setX(0,alphaPol*(Eq3.x(0) + Ep3.x(0)));
        		P3.setX(1,alphaPol*(Eq3.x(1) + Ep3.x(1)));
        		P3.setX(2,alphaPol*(Eq3.x(2) + Ep3.x(2)));

        		
        		// Evaluate the criteria
        		
	        	deltaP1 = Math.sqrt(P1.Mv1Squared(P1old));
	        	deltaP2 = Math.sqrt(P2.Mv1Squared(P2old));
	        	deltaP3 = Math.sqrt(P3.Mv1Squared(P3old));
	        	
	        	counterSCFloop = counterSCFloop + 1;
	    
	        	
	        	if (deltaP1<1e-15) { //0.00338) {
	        		noSCFforP1 = false;
	        	}
	        	else {
	        		P1old.E(P1);
	        	}
	        	if (deltaP2<1e-15) { //0.00338) {
	        		noSCFforP2 = false;
	        	}
	        	else {
	        		P2old.E(P2);
	        	}
	        	if (deltaP3<1e-15) { //0.00338) {
	        		noSCFforP3 = false;
	        	}
	        	else {
	        		P3old.E(P3);
	        	}

	        	
	        	if (counterSCFloop >= 1000) {
//	        		counterSCFloopOK = false;
	        		noSCFforP1 = false;
	        		noSCFforP2 = false;
	        		noSCFforP3 = false;
	        		loopFailures = loopFailures + 1;
//	        		System.out.println("counterSCFloop = " + counterSCFloop + ", exiting SCF loop due to likely large repulsion");
    		        chargeH11 = Electron.UNIT.toSim(0.6113);
    		        chargeH12 = Electron.UNIT.toSim(0.6113);
    		        chargeM1 = Electron.UNIT.toSim(-1.2226);
    		        chargeH21 = Electron.UNIT.toSim(0.6113);
    		        chargeH22 = Electron.UNIT.toSim(0.6113);
    		        chargeM2 = Electron.UNIT.toSim(-1.2226);
    		        chargeH31 = Electron.UNIT.toSim(0.6113);
    		        chargeH32 = Electron.UNIT.toSim(0.6113);
    		        chargeM3 = Electron.UNIT.toSim(-1.2226);
    		        
	        	}
	     
	        // REPEAT LOOP HERE.  RE-EVALUATE CHARGES ON W1 AND THEN REPEAT.
        }
*/        
        
        /*
         * Here is where I need to add the polarization term to the energy sum.
         * kmb 5/4/06
         */        
                
//                double alpha = 1.444;
                
                double uW1, uW2, uW3; //UpolAtkins;
                        
                UpolAtkins = -0.5*(P1.dot(Eq1)+P2.dot(Eq2)+P3.dot(Eq3));
                       
                double UpolW1W2 = -0.5*(P1.dot(Eq1)+P2.dot(Eq2));
                
//                double UpolEquation8 = -(P1.dot(Eq1)+P2.dot(Eq2)+P3.dot(Eq3))-0.5*(P1.dot(Ep1)+P2.dot(Ep2)+P3.dot(Ep3))+(1/2/alpha)*(P1.squared()+P2.squared()+P3.squared());
                
                if (rO1O2 <= 2.1 || rO1O3 <= 2.1 || rO2O3 <= 2.1) UpolAtkins = 0; // value from Cummings
                
                sumSCF += UpolAtkins; // used to be += UpolAtkins
                
                
                uW1 = 1.855 + Debye.UNIT.fromSim(Math.sqrt(P1.squared()));
                uW2 = 1.855 + Debye.UNIT.fromSim(Math.sqrt(P2.squared()));
                uW3 = 1.855 + Debye.UNIT.fromSim(Math.sqrt(P3.squared()));

                
            //System.out.println("Actually going to return an energy not equal to infinity" + sum);
            //if (rOO >= 500) System.out.println("O-O distance = " + rOO + ", sum = " + sum);
                
                if (sumSCF <= -50000) {
/*            		System.out.println("O1-O2 distance = " + rO1O2 + ", sumSCF = " + sumSCF);
            		System.out.println("O1-O3 distance = " + rO1O3 + ", O2-O3 distance = " + rO2O3);
            	    System.out.println("O1x = " + O1r.x(0) + ", O1y = " + O1r.x(1) + ", O1z = " + O1r.x(2));
            	    System.out.println("O2x = " + O2r.x(0) + ", O2y = " + O2r.x(1) + ", O2z = " + O2r.x(2));
            	    System.out.println("O3x = " + O3r.x(0) + ", O3y = " + O3r.x(1) + ", O3z = " + O3r.x(2));
            	    System.out.println("M1x = " + M1r.x(0) + ", M1y = " + M1r.x(1) + ", M1z = " + M1r.x(2));
            	    System.out.println("M2x = " + M2r.x(0) + ", M2y = " + M2r.x(1) + ", M2z = " + M2r.x(2));
            	    System.out.println("M3x = " + M3r.x(0) + ", M3y = " + M3r.x(1) + ", M3z = " + M3r.x(2));
            	    System.out.println("H11x = " + H11r.x(0) + ", H11y = " + H11r.x(1) + ", H11z = " + H11r.x(2));
            	    System.out.println("H12x = " + H12r.x(0) + ", H12y = " + H12r.x(1) + ", H12z = " + H12r.x(2));
            	    System.out.println("H21x = " + H21r.x(0) + ", H21y = " + H21r.x(1) + ", H21z = " + H21r.x(2));
            	    System.out.println("H22x = " + H22r.x(0) + ", H22y = " + H22r.x(1) + ", H22z = " + H22r.x(2));
            	    System.out.println("H31x = " + H31r.x(0) + ", H31y = " + H31r.x(1) + ", H31z = " + H31r.x(2));
            	    System.out.println("H32x = " + H32r.x(0) + ", H32y = " + H32r.x(1) + ", H32z = " + H32r.x(2));
*/
            }

        // call energy method, removing the third atom from the picture
                
               /* AtomPair twoAtoms = new AtomPair(atomsSCF.get(0), atomsSCF.get(1));
                AtomSet setTwoAtoms = (AtomSet)twoAtoms;
                
                energy(setTwoAtoms);
*/              
                //System.out.println("From iteration: P1 = " + P1 + ", P2 = " + P2 + ", P3 = " + P3);
                //double PolEnergyDiff = UpolAtkinsAnalytical - UpolAtkins;
                //System.out.println("Difference in Analytical and Iterative Energies = " + PolEnergyDiff);
                
        return sumSCF;
        
//        return scfSums;
    }//end of energy

    
    public double energySCF4(AtomSet atomsSCF) { //Atom atom1, Atom atom2, Atom atom3) {//(AtomSet atoms){
        double sumSCF = 0.0;
        double r2 = 0.0;
        double rO1O2 = 0.0;
        double rO1O3 = 0.0;
        double rO1O4 = 0.0;
        double rO2O3 = 0.0;
        double rO2O4 = 0.0;
        double rO3O4 = 0.0;
        
        AtomWater4P node1 = (AtomWater4P)atomsSCF.getAtom(0);
        AtomWater4P node2 = (AtomWater4P)atomsSCF.getAtom(1);//(AtomTreeNodeWaterPPC)pair.atom1.node;
        AtomWater4P node3 = (AtomWater4P)atomsSCF.getAtom(2);
        AtomWater4P node4 = (AtomWater4P)atomsSCF.getAtom(3);

        
        IVector O1r = node1.O.getPosition();
        IVector O2r = node2.O.getPosition();
        IVector O3r = node3.O.getPosition();
        IVector O4r = node4.O.getPosition();
        IVector H11r = node1.H1.getPosition();
        IVector H12r = node1.H2.getPosition();
        IVector H21r = node2.H1.getPosition();
        IVector H22r = node2.H2.getPosition();
        IVector H31r = node3.H1.getPosition();
        IVector H32r = node3.H2.getPosition();
        IVector H41r = node4.H1.getPosition();
        IVector H42r = node4.H2.getPosition();

        IVector M1r = node1.M.getPosition();
        IVector M2r = node2.M.getPosition();
        IVector M3r = node3.M.getPosition();
        IVector M4r = node4.M.getPosition();
        
        
//        System.out.println("O1r coordinates before: " + O1r.x(0) + ", " + O1r.x(1) + ", " + O1r.x(2));
//        System.out.println("O2r coordinates before: " + O2r.x(0) + ", " + O2r.x(1) + ", " + O2r.x(2));
//        System.out.println("H11r coordinates before: " + H11r.x(0) + ", " + H11r.x(1) + ", " + H11r.x(2));
//        System.out.println("H12r coordinates before: " + H12r.x(0) + ", " + H12r.x(1) + ", " + H12r.x(2));
//        System.out.println("H21r coordinates before: " + H21r.x(0) + ", " + H21r.x(1) + ", " + H21r.x(2));
//        System.out.println("H22r coordinates before: " + H22r.x(0) + ", " + H22r.x(1) + ", " + H22r.x(2));
//        System.out.println("M1r coordinates before: " + M1r.x(0) + ", " + M1r.x(1) + ", " + M1r.x(2));
//        System.out.println("M2r coordinates before: " + M2r.x(0) + ", " + M2r.x(1) + ", " + M2r.x(2));
        
/*        O2r.PE(0,4); //PE(5);
        H21r.PE(0,4);
        H22r.PE(0,4);
        M2rOriginal.PE(0,4);
        
        O2r.setX(2,2); //PE(5);
        H21r.setX(2,2);
        H22r.setX(2,2);
        M2rOriginal.setX(2,2);
        //O2r.setX(2,2);
*/        
//		Dimer configuration
/*        O2r.setX(2,1.7683903);
        O2r.setX(0,2.18378015);
        O2r.setX(1,0.0);
        H21r.setX(2,1.90087324);
        H21r.setX(0,2.73561133);
        H21r.setX(1,-0.7531133);
        H22r.setX(2,1.90087324);
        H22r.setX(0,2.73561133);
        H22r.setX(1,0.7531133);
        M2rOriginal.setX(2,1.79406929);
        M2rOriginal.setX(0,2.29074084);
        M2rOriginal.setX(1,0.0);
*/

/*        O2r.setX(2,2.2);
        O2r.setX(0,0.0);
        O2r.setX(1,0.0);
        H21r.setX(2,2.767511567);
        H21r.setX(0,-0.75311329);
        H21r.setX(1,0.0);
        H22r.setX(2,2.767511567);
        H22r.setX(0,0.75311329);
        H22r.setX(1,0.0);
        M2rOriginal.setX(2,2.31);
        M2rOriginal.setX(0,0.0);
        M2rOriginal.setX(1,0.0);


        O1r.setX(2,24600.0);
        O1r.setX(0,0.0);
        O1r.setX(1,0.0);
        H11r.setX(2,24600.567511567);
        H11r.setX(0,-0.75311329);
        H11r.setX(1,0.0);
        H12r.setX(2,24600.567511567);
        H12r.setX(0,0.75311329);
        H12r.setX(1,0.0);
        M1rOriginal.setX(2,24600.11);
        M1rOriginal.setX(0,0.0);
        M1rOriginal.setX(1,0.0);

        O2r.setX(2,0.0);
        O2r.setX(0,0.0);
        O2r.setX(1,0.0);
        H21r.setX(2,0.567511567);
        H21r.setX(0,-0.75311329);
        H21r.setX(1,0.0);
        H22r.setX(2,0.567511567);
        H22r.setX(0,0.75311329);
        H22r.setX(1,0.0);
        M2rOriginal.setX(2,0.11);
        M2rOriginal.setX(0,0.0);
        M2rOriginal.setX(1,0.0);

        O3r.setX(2,-24600.0);
        O3r.setX(0,0.0);
        O3r.setX(1,0.0);
        H31r.setX(2,-24600.567511567);
        H31r.setX(0,-0.75311329);
        H31r.setX(1,0.0);
        H32r.setX(2,-24600.567511567);
        H32r.setX(0,0.75311329);
        H32r.setX(1,0.0);
        M3rOriginal.setX(2,-24600.11);
        M3rOriginal.setX(0,0.0);
        M3rOriginal.setX(1,0.0);
*/
		
        
		// moved to constructor; KMB, 7/23/07
//        final double core = 4.41; //4.41 = 2.1^2; value according to Cummings
  
        // Here is the beginning of my self-consistent solution algorithm
        // 12/26/05
   
        // moved to constructor; KMB, 7/23/07
/*        Vector Eq1old = new Vector3D();
        Vector Eq2old = new Vector3D();
        Vector Eq3old = new Vector3D();
        Vector Eq1 = new Vector3D();
        Vector Eq2 = new Vector3D();
        Vector Eq3 = new Vector3D();
        Vector Eq4 = new Vector3D();
        Vector Ep1old = new Vector3D();
        Vector Ep2old = new Vector3D();
        Vector Ep3old = new Vector3D();
        Vector Ep1 = new Vector3D();
        Vector Ep2 = new Vector3D();
        Vector Ep3 = new Vector3D();
        Vector Ep4 = new Vector3D();
*/        
/*        Vector Eq2LabFrame = new Vector3D();
        Vector Eq1LabFrame = new Vector3D();
        Vector Eq3LabFrame = new Vector3D();
        Vector Ep3LabFrame = new Vector3D();
        Vector Ep2LabFrame = new Vector3D();
        Vector Ep1LabFrame = new Vector3D();
*/
/*        Vector Eon1Total = new Vector3D();
        Vector Eon2Total = new Vector3D();
        Vector Eon3Total = new Vector3D();
        Vector Eon4Total = new Vector3D();
        
        Eon1Total.Ev1Pv2(Eq1, Ep1);
        Eon2Total.Ev1Pv2(Eq2, Ep2);
        Eon3Total.Ev1Pv2(Eq3, Ep3);
        Eon4Total.Ev1Pv2(Eq4, Ep4);

        Vector P1 = new Vector3D();
        Vector P2 = new Vector3D();
        Vector P3 = new Vector3D();
        Vector P4 = new Vector3D();
        Vector P1old = new Vector3D();
        Vector P2old = new Vector3D();
        Vector P3old = new Vector3D();
        Vector P4old = new Vector3D();
*/        
        int counterSCFloop = 0; //1;
//        boolean counterSCFloopOK = true;
        int loopFailures = 0;
        
        
        // These booleans monitor whether the SCF algorithm has iteratively found
        // a solution given the criterion of Cummings, J.Chem.Phys. 2005
        boolean noSCFforP1 = true;
        boolean noSCFforP2 = true;
        boolean noSCFforP3 = true;
        boolean noSCFforP4 = true;
		
        // Need loop to check for configuration overlap between charged particles
        
//      compute O-O distance to consider bypassing the SCF loop   
        rO1O2 = Math.sqrt(O1r.Mv1Squared(O2r));
        rO1O3 = Math.sqrt(O1r.Mv1Squared(O3r));
        rO1O4 = Math.sqrt(O1r.Mv1Squared(O4r));
        rO2O3 = Math.sqrt(O3r.Mv1Squared(O2r));
        rO2O4 = Math.sqrt(O4r.Mv1Squared(O2r));
        rO3O4 = Math.sqrt(O4r.Mv1Squared(O3r));
        
        if (rO1O2 > 100 || rO1O3 > 100 || rO1O4 > 100 || rO2O3 > 100 || rO2O4 > 100 || rO3O4 > 100) {
			return -123456789.0;
    }

        
        //System.out.println("O-O distance is " + r2);
        
		//System.out.println(rO1O2 + ", " + rO1O3 + ", " + rO1O4 + ", " + rO2O3 + ", " + rO2O4 + ", " + rO3O4 );
        
        if(rO1O2 <= 2.1 || rO1O3 <= 2.1 || rO1O4 <= 2.1 || rO2O3 <= 2.1 || rO2O4 <= 2.1 || rO3O4 <= 2.1) { // use to be 2.1 for Cummings cutoff
/*        		noSCFforP1 = false;
        		noSCFforP2 = false;
        		noSCFforP3 = false;
        		noSCFforP4 = false;
        		
        		// need to give some value to these or energy method will return garbage
        		
	        chargeH11 = Electron.UNIT.toSim(0.6113);
	        chargeH12 = Electron.UNIT.toSim(0.6113);
	        chargeM1 = Electron.UNIT.toSim(-1.2226);
	        chargeH21 = Electron.UNIT.toSim(0.6113);
	        chargeH22 = Electron.UNIT.toSim(0.6113);
	        chargeM2 = Electron.UNIT.toSim(-1.2226);
	        chargeH31 = Electron.UNIT.toSim(0.6113);
	        chargeH32 = Electron.UNIT.toSim(0.6113);
	        chargeM3 = Electron.UNIT.toSim(-1.2226);
	        chargeH41 = Electron.UNIT.toSim(0.6113);
	        chargeH42 = Electron.UNIT.toSim(0.6113);
	        chargeM4 = Electron.UNIT.toSim(-1.2226);*/
	        
	        return Double.POSITIVE_INFINITY;

        }
                
                
        gamma = 12.75;
//	    double r = Math.sqrt(r2);
        double rO1O2OverSigma = rO1O2/sigma;
        double sigma2OverRO1O2sq = 1/(rO1O2OverSigma*rO1O2OverSigma);
        double rO1O3OverSigma = rO1O3/sigma;
        double sigma2OverRO1O3sq = 1/(rO1O3OverSigma*rO1O3OverSigma);
        double rO1O4OverSigma = rO1O4/sigma;
        double sigma2OverRO1O4sq = 1/(rO1O4OverSigma*rO1O4OverSigma);
        double rO2O3OverSigma = rO2O3/sigma;
        double sigma2OverRO2O3sq = 1/(rO2O3OverSigma*rO2O3OverSigma);
        double rO2O4OverSigma = rO2O4/sigma;
        double sigma2OverRO2O4sq = 1/(rO2O4OverSigma*rO2O4OverSigma);
        double rO3O4OverSigma = rO3O4/sigma;
        double sigma2OverRO3O4sq = 1/(rO3O4OverSigma*rO3O4OverSigma);
        
        double sixOverGamma = 6/gamma;
   
        sumSCF += epsilon/(1 - sixOverGamma)*(sixOverGamma*Math.exp(gamma*(1 - rO1O2OverSigma)) - sigma2OverRO1O2sq*sigma2OverRO1O2sq*sigma2OverRO1O2sq);
        sumSCF += epsilon/(1 - sixOverGamma)*(sixOverGamma*Math.exp(gamma*(1 - rO1O3OverSigma)) - sigma2OverRO1O3sq*sigma2OverRO1O3sq*sigma2OverRO1O3sq);
        sumSCF += epsilon/(1 - sixOverGamma)*(sixOverGamma*Math.exp(gamma*(1 - rO1O4OverSigma)) - sigma2OverRO1O4sq*sigma2OverRO1O4sq*sigma2OverRO1O4sq);
        sumSCF += epsilon/(1 - sixOverGamma)*(sixOverGamma*Math.exp(gamma*(1 - rO2O3OverSigma)) - sigma2OverRO2O3sq*sigma2OverRO2O3sq*sigma2OverRO2O3sq);
        sumSCF += epsilon/(1 - sixOverGamma)*(sixOverGamma*Math.exp(gamma*(1 - rO2O4OverSigma)) - sigma2OverRO2O4sq*sigma2OverRO2O4sq*sigma2OverRO2O4sq);
        sumSCF += epsilon/(1 - sixOverGamma)*(sixOverGamma*Math.exp(gamma*(1 - rO3O4OverSigma)) - sigma2OverRO3O4sq*sigma2OverRO3O4sq*sigma2OverRO3O4sq);
        
//        double sumO2Exp6 = sumSCF;
        
/*        if (Math.abs(Math.sqrt(Eon1.squared())) >= 0.2083 || Math.abs(Math.sqrt(Eon2.squared())) >= 0.2083 ) {
    			System.out.println("About to return Double.Positive_Infinity");
        		System.out.println("Outside core, with bad electric field values: rOO = " + rOO + ", sum = " + sum + ", Eon1(V/A) = " + Math.sqrt(Eon1.squared())*EFconverttoVperA + ", Eon2(V/A) = " + Math.sqrt(Eon2.squared())*EFconverttoVperA);
        		return Double.POSITIVE_INFINITY;
        }
*/

        // moved to constructor; KMB, 7/24/07
/*        double sigmaM = 0.610;
        double sigmaH = 0.455;
        double sqrtHMsigmas = Math.sqrt(2*(sigmaH*sigmaH+sigmaM*sigmaM));
*/        
        // MUST INCLUDE ERF FUNCTION STUFF TO COULOMBIC ENERGY PART!
        // KMB 8/3/06
        
        r2 = H11r.Mv1Squared(H21r);  // COUNTED-2
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        	sumSCF += chargeH11*chargeH21/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaH)));

        r2 = H11r.Mv1Squared(H22r); // COUNTED-2
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH11*chargeH22/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaH)));

        r2 = H12r.Mv1Squared(H21r); // COUNTED-2
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH12*chargeH21/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaH)));

        r2 = H12r.Mv1Squared(H22r); // COUNTED-2
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH12*chargeH22/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaH)));

//        System.out.println("sum of all O-H terms is " + sum);
        
        r2 = M1r.Mv1Squared(H21r); // COUNTED-2
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH21*chargeM1/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/sqrtHMsigmas));
        //System.out.println("sum is " + sum);

        r2 = M1r.Mv1Squared(H22r); // COUNTED-2
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH22*chargeM1/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/sqrtHMsigmas));
        //System.out.println("sum is " + sum);

        r2 = M2r.Mv1Squared(H11r); // COUNTED-2
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH11*chargeM2/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/sqrtHMsigmas));

        ///System.out.println("sum is " + sum);
        r2 = M2r.Mv1Squared(H12r); // COUNTED-2
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH12*chargeM2/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/sqrtHMsigmas));
        //System.out.println("sum is " + sum);

        r2 = M1r.Mv1Squared(M2r); // COUNTED-2
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeM1*chargeM2/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaM)));
        //System.out.println("sum is " + sum);

        
        
        // EXTRA TERMS FOR 3 MOLECULE PERMUTATIONS AND COMBINATIONS FOR ELECTROSTATICS HERE
        


        r2 = H11r.Mv1Squared(H31r);  // COUNTED-2
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        	sumSCF += chargeH11*chargeH31/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaH)));

        r2 = H11r.Mv1Squared(H32r);  // COUNTED-2
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH11*chargeH32/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaH)));



        r2 = H12r.Mv1Squared(H31r);  // COUNTED-2
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH12*chargeH31/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaH)));

        r2 = H12r.Mv1Squared(H32r);  // COUNTED-2
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH12*chargeH32/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaH)));

        r2 = H21r.Mv1Squared(H31r);  // COUNTED-2
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH21*chargeH31/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaH)));

        r2 = H21r.Mv1Squared(H32r);  // COUNTED-2
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH21*chargeH32/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaH)));

        r2 = H22r.Mv1Squared(H31r);  // COUNTED-2
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH22*chargeH31/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaH)));

        r2 = H22r.Mv1Squared(H32r);  // COUNTED-2
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH22*chargeH32/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaH)));

//        System.out.println("sum of all O-H terms is " + sum);
        


        r2 = M1r.Mv1Squared(H31r);  // COUNTED-2
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH31*chargeM1/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/sqrtHMsigmas));
        //System.out.println("sum is " + sum);

        r2 = M1r.Mv1Squared(H32r);  // COUNTED-2
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH32*chargeM1/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/sqrtHMsigmas));
        //System.out.println("sum is " + sum);



        r2 = M2r.Mv1Squared(H31r);  // COUNTED-2
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH31*chargeM2/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/sqrtHMsigmas));

        ///System.out.println("sum is " + sum);
        r2 = M2r.Mv1Squared(H32r);  // COUNTED-2
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH32*chargeM2/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/sqrtHMsigmas));
        //System.out.println("sum is " + sum);

        r2 = M3r.Mv1Squared(H11r);  // COUNTED-2
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH11*chargeM3/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/sqrtHMsigmas));

        ///System.out.println("sum is " + sum);
        r2 = M3r.Mv1Squared(H12r);  // COUNTED-2
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH12*chargeM3/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/sqrtHMsigmas));
        //System.out.println("sum is " + sum);

        r2 = M3r.Mv1Squared(H21r);  // COUNTED-2
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH21*chargeM3/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/sqrtHMsigmas));

        ///System.out.println("sum is " + sum);
        r2 = M3r.Mv1Squared(H22r);  // COUNTED-2
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH22*chargeM3/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/sqrtHMsigmas));
        //System.out.println("sum is " + sum);


        r2 = M1r.Mv1Squared(M3r);  // COUNTED-2
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeM1*chargeM3/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaM)));
        //System.out.println("sum is " + sum);

        r2 = M3r.Mv1Squared(M2r);  // COUNTED-2
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeM3*chargeM2/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaM)));
        //System.out.println("sum is " + sum);

        
        
        /*
         * NOW PUTTING IN THE EXTRA ELECTROSTATIC TERMS FOR THE FOURTH WATER MOLECULE
         * KMB, 8/14/06
         */
        
        
        r2 = H11r.Mv1Squared(H41r);  // OK FOR B4
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH11*chargeH41/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaH)));
        
        r2 = H11r.Mv1Squared(H42r);  // OK FOR B4
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH11*chargeH42/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaH)));
        
        r2 = H12r.Mv1Squared(H41r);  // OK FOR B4
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH12*chargeH41/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaH)));
        
        r2 = H12r.Mv1Squared(H42r);  // OK FOR B4
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH12*chargeH42/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaH)));
        
        r2 = H21r.Mv1Squared(H41r);  // OK FOR B4
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH21*chargeH41/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaH)));
        
        r2 = H21r.Mv1Squared(H42r);  // OK FOR B4
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH21*chargeH42/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaH)));
        
        r2 = H22r.Mv1Squared(H41r);  // OK FOR B4
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH22*chargeH41/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaH)));
        
        r2 = H22r.Mv1Squared(H42r);  // OK FOR B4
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH22*chargeH42/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaH)));
        
        r2 = H31r.Mv1Squared(H41r);  // OK FOR B4
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH31*chargeH41/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaH)));
        
        r2 = H31r.Mv1Squared(H42r);  // OK FOR B4
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH31*chargeH42/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaH)));
        
        r2 = H32r.Mv1Squared(H41r);  // OK FOR B4
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH32*chargeH41/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaH)));
        
        r2 = H32r.Mv1Squared(H42r);  // OK FOR B4
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH32*chargeH42/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaH)));
        
        r2 = M4r.Mv1Squared(H11r);  // OK FOR B4
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH11*chargeM4/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/sqrtHMsigmas));

        r2 = M4r.Mv1Squared(H12r);  // OK FOR B4
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH12*chargeM4/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/sqrtHMsigmas));

        r2 = M1r.Mv1Squared(H41r);  // OK FOR B4
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH41*chargeM1/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/sqrtHMsigmas));

        r2 = M1r.Mv1Squared(H42r);  // OK FOR B4
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH42*chargeM1/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/sqrtHMsigmas));

        r2 = M4r.Mv1Squared(H21r);  // OK FOR B4
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH21*chargeM4/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/sqrtHMsigmas));

        r2 = M4r.Mv1Squared(H22r);  // OK FOR B4
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH22*chargeM4/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/sqrtHMsigmas));

        r2 = M2r.Mv1Squared(H41r);  // OK FOR B4
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH41*chargeM2/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/sqrtHMsigmas));

        r2 = M2r.Mv1Squared(H42r);  // OK FOR B4
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH42*chargeM2/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/sqrtHMsigmas));

        r2 = M4r.Mv1Squared(H31r);  // OK FOR B4
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH31*chargeM4/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/sqrtHMsigmas));

        r2 = M4r.Mv1Squared(H32r);  // OK FOR B4
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH32*chargeM4/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/sqrtHMsigmas));

        r2 = M3r.Mv1Squared(H41r);  // OK FOR B4
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH41*chargeM3/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/sqrtHMsigmas));

        r2 = M3r.Mv1Squared(H42r);  // OK FOR B4
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH42*chargeM3/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/sqrtHMsigmas));

        r2 = M1r.Mv1Squared(M4r);  // OK FOR B4
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeM1*chargeM4/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaM)));
        //System.out.println("sum is " + sum);

        r2 = M4r.Mv1Squared(M2r);  // OK FOR B4
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeM4*chargeM2/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaM)));
        //System.out.println("sum is " + sum);

        r2 = M3r.Mv1Squared(M4r);  // OK FOR B4
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeM3*chargeM4/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaM)));
        //System.out.println("sum is " + sum);

        
        /*
         * NOW, HOW DO I WANT TO HANDLE THE ELECTRIC FIELDS?
         * KEEP GOING WITH ADDING TERMS, OR START USING ARRAYS?
         * IS IT EVEN POSSIBLE TO USE ARRAYS?
         * 
         * JUST GOING TO KEEP ADDING TERMS
         * WILL HAVE TO ADD ARRAYS LATER
         * ATOM WILL HAVE TO KNOW ITS POSITION AND CHARGE, FOR EXAMPLE
         * 
         * KMB, 8/15/06
         */
		
        
        /*
         * Finding the Electric fields at the center of mass of each molecule, Eqi
         * kmb, 8/7/06
         */


        // moved to constructor; KMB, 7/24/07
/*        Vector3D comW1 = new Vector3D();
        Vector3D comW2 = new Vector3D();
        Vector3D comW3 = new Vector3D();
        Vector3D comW4 = new Vector3D();
*/        
        // How to use COM etomica code correctly? kmb, 8/7/06
/*        DataSourceCOM com = new DataSourceCOM(space);
        com.actionPerformed(atoms.getAtom(0));
        com.getCOM();
*/        
        
        // moved to constructor; KMB, 7/24/07
/*        double massH = 1.01;
        double massO = 16.0;
        double totalMass = 18.02;
*/ 
        double comW1Xcomp = 0.0;
        double comW1Ycomp = 0.0;
        double comW1Zcomp = 0.0;
        
        comW1Xcomp = massH*H11r.x(0) + massO*O1r.x(0) + massH*H12r.x(0);
        comW1Ycomp = massH*H11r.x(1) + massO*O1r.x(1) + massH*H12r.x(1);
        comW1Zcomp = massH*H11r.x(2) + massO*O1r.x(2) + massH*H12r.x(2);
        
        comW1.setX(0,comW1Xcomp);
        comW1.setX(1,comW1Ycomp);
        comW1.setX(2,comW1Zcomp);
        
        comW1.Ea1Tv1(1/totalMass,comW1);

        double comW2Xcomp = 0.0;
        double comW2Ycomp = 0.0;
        double comW2Zcomp = 0.0;
        
        comW2Xcomp = massH*H21r.x(0) + massO*O2r.x(0) + massH*H22r.x(0);
        comW2Ycomp = massH*H21r.x(1) + massO*O2r.x(1) + massH*H22r.x(1);
        comW2Zcomp = massH*H21r.x(2) + massO*O2r.x(2) + massH*H22r.x(2);
        
        comW2.setX(0,comW2Xcomp);
        comW2.setX(1,comW2Ycomp);
        comW2.setX(2,comW2Zcomp);
        
        comW2.Ea1Tv1(1/totalMass,comW2);

        double comW3Xcomp = 0.0;
        double comW3Ycomp = 0.0;
        double comW3Zcomp = 0.0;
        
        comW3Xcomp = massH*H31r.x(0) + massO*O3r.x(0) + massH*H32r.x(0);
        comW3Ycomp = massH*H31r.x(1) + massO*O3r.x(1) + massH*H32r.x(1);
        comW3Zcomp = massH*H31r.x(2) + massO*O3r.x(2) + massH*H32r.x(2);
        
        comW3.setX(0,comW3Xcomp);
        comW3.setX(1,comW3Ycomp);
        comW3.setX(2,comW3Zcomp);
        
        comW3.Ea1Tv1(1/totalMass,comW3);


        double comW4Xcomp = 0.0;
        double comW4Ycomp = 0.0;
        double comW4Zcomp = 0.0;
        
        comW4Xcomp = massH*H41r.x(0) + massO*O4r.x(0) + massH*H42r.x(0);
        comW4Ycomp = massH*H41r.x(1) + massO*O4r.x(1) + massH*H42r.x(1);
        comW4Zcomp = massH*H41r.x(2) + massO*O4r.x(2) + massH*H42r.x(2);
        
        comW4.setX(0,comW4Xcomp);
        comW4.setX(1,comW4Ycomp);
        comW4.setX(2,comW4Zcomp);
        
        comW4.Ea1Tv1(1/totalMass,comW4);

        
        /*
         * DOUBLE SUMMING NOW COMPLETE; KMB 8/9/06
         * 
         * THESE ELECTRIC FIELDS ARE WRONG! KMB 8/8/06
         * I HAVE NOT DONE THE DOUBLE SUM REQUIRED FOR 3 WATER
         * MOLECULES AS PER CUMMINGS EQUATION 4!
         */
        
//        double sqrtHMsigmas = Math.sqrt(2*(sigmaH*sigmaH+sigmaM*sigmaM));
        // moved to constructor; KMB, 7/24/07
/*        double sqrtPiHMsigmas = Math.sqrt(Math.PI*(sigmaH*sigmaH+sigmaM*sigmaM));
        double sqrtPiMMsigmas = Math.sqrt(Math.PI*(2*sigmaM*sigmaM));
*/        
        double Eq1XcompW2 = 0.0;
        double Eq1YcompW2 = 0.0;
        double Eq1ZcompW2 = 0.0;
        double Eq1XcompW3 = 0.0;
        double Eq1YcompW3 = 0.0;
        double Eq1ZcompW3 = 0.0;
        double Eq1XcompW4 = 0.0;
        double Eq1YcompW4 = 0.0;
        double Eq1ZcompW4 = 0.0;

        
        double comW1toH21 = Math.sqrt(comW1.Mv1Squared(H21r));
        double comW1toH22 = Math.sqrt(comW1.Mv1Squared(H22r));
        double comW1toM2 = Math.sqrt(comW1.Mv1Squared(M2r));

        double comW1toH31 = Math.sqrt(comW1.Mv1Squared(H31r));
        double comW1toH32 = Math.sqrt(comW1.Mv1Squared(H32r));
        double comW1toM3 = Math.sqrt(comW1.Mv1Squared(M3r));

        double comW1toH41 = Math.sqrt(comW1.Mv1Squared(H41r));
        double comW1toH42 = Math.sqrt(comW1.Mv1Squared(H42r));
        double comW1toM4 = Math.sqrt(comW1.Mv1Squared(M4r));

        
        // Contributions to sum from water#2
        Eq1XcompW2 += chargeH21*(comW1.x(0)-H21r.x(0))/(comW1toH21*comW1toH21*comW1toH21)*((1-SpecialFunctions.erfc(comW1toH21/sqrtHMsigmas))-Math.sqrt(2)*comW1toH21/sqrtPiHMsigmas*Math.exp(-comW1toH21*comW1toH21/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq1XcompW2 += chargeH22*(comW1.x(0)-H22r.x(0))/(comW1toH22*comW1toH22*comW1toH22)*((1-SpecialFunctions.erfc(comW1toH22/sqrtHMsigmas))-Math.sqrt(2)*comW1toH22/sqrtPiHMsigmas*Math.exp(-comW1toH22*comW1toH22/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq1XcompW2 += chargeM2*(comW1.x(0)-M2r.x(0))/(comW1toM2*comW1toM2*comW1toM2)*((1-SpecialFunctions.erfc(comW1toM2/(2*sigmaM)))-Math.sqrt(2)*comW1toM2/sqrtPiMMsigmas*Math.exp(-comW1toM2*comW1toM2/(4*sigmaM*sigmaM)));

        Eq1YcompW2 += chargeH21*(comW1.x(1)-H21r.x(1))/(comW1toH21*comW1toH21*comW1toH21)*((1-SpecialFunctions.erfc(comW1toH21/sqrtHMsigmas))-Math.sqrt(2)*comW1toH21/sqrtPiHMsigmas*Math.exp(-comW1toH21*comW1toH21/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq1YcompW2 += chargeH22*(comW1.x(1)-H22r.x(1))/(comW1toH22*comW1toH22*comW1toH22)*((1-SpecialFunctions.erfc(comW1toH22/sqrtHMsigmas))-Math.sqrt(2)*comW1toH22/sqrtPiHMsigmas*Math.exp(-comW1toH22*comW1toH22/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq1YcompW2 += chargeM2*(comW1.x(1)-M2r.x(1))/(comW1toM2*comW1toM2*comW1toM2)*((1-SpecialFunctions.erfc(comW1toM2/(2*sigmaM)))-Math.sqrt(2)*comW1toM2/sqrtPiMMsigmas*Math.exp(-comW1toM2*comW1toM2/(4*sigmaM*sigmaM)));

        Eq1ZcompW2 += chargeH21*(comW1.x(2)-H21r.x(2))/(comW1toH21*comW1toH21*comW1toH21)*((1-SpecialFunctions.erfc(comW1toH21/sqrtHMsigmas))-Math.sqrt(2)*comW1toH21/sqrtPiHMsigmas*Math.exp(-comW1toH21*comW1toH21/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq1ZcompW2 += chargeH22*(comW1.x(2)-H22r.x(2))/(comW1toH22*comW1toH22*comW1toH22)*((1-SpecialFunctions.erfc(comW1toH22/sqrtHMsigmas))-Math.sqrt(2)*comW1toH22/sqrtPiHMsigmas*Math.exp(-comW1toH22*comW1toH22/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq1ZcompW2 += chargeM2*(comW1.x(2)-M2r.x(2))/(comW1toM2*comW1toM2*comW1toM2)*((1-SpecialFunctions.erfc(comW1toM2/(2*sigmaM)))-Math.sqrt(2)*comW1toM2/sqrtPiMMsigmas*Math.exp(-comW1toM2*comW1toM2/(4*sigmaM*sigmaM)));

        
        // Contributions to sum from water#3
        Eq1XcompW3 += chargeH31*(comW1.x(0)-H31r.x(0))/(comW1toH31*comW1toH31*comW1toH31)*((1-SpecialFunctions.erfc(comW1toH31/sqrtHMsigmas))-Math.sqrt(2)*comW1toH31/sqrtPiHMsigmas*Math.exp(-comW1toH31*comW1toH31/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq1XcompW3 += chargeH32*(comW1.x(0)-H32r.x(0))/(comW1toH32*comW1toH32*comW1toH32)*((1-SpecialFunctions.erfc(comW1toH32/sqrtHMsigmas))-Math.sqrt(2)*comW1toH32/sqrtPiHMsigmas*Math.exp(-comW1toH32*comW1toH32/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq1XcompW3 += chargeM3*(comW1.x(0)-M3r.x(0))/(comW1toM3*comW1toM3*comW1toM3)*((1-SpecialFunctions.erfc(comW1toM3/(2*sigmaM)))-Math.sqrt(2)*comW1toM3/sqrtPiMMsigmas*Math.exp(-comW1toM3*comW1toM3/(4*sigmaM*sigmaM)));

        Eq1YcompW3 += chargeH31*(comW1.x(1)-H31r.x(1))/(comW1toH31*comW1toH31*comW1toH31)*((1-SpecialFunctions.erfc(comW1toH31/sqrtHMsigmas))-Math.sqrt(2)*comW1toH31/sqrtPiHMsigmas*Math.exp(-comW1toH31*comW1toH31/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq1YcompW3 += chargeH32*(comW1.x(1)-H32r.x(1))/(comW1toH32*comW1toH32*comW1toH32)*((1-SpecialFunctions.erfc(comW1toH32/sqrtHMsigmas))-Math.sqrt(2)*comW1toH32/sqrtPiHMsigmas*Math.exp(-comW1toH32*comW1toH32/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq1YcompW3 += chargeM3*(comW1.x(1)-M3r.x(1))/(comW1toM3*comW1toM3*comW1toM3)*((1-SpecialFunctions.erfc(comW1toM3/(2*sigmaM)))-Math.sqrt(2)*comW1toM3/sqrtPiMMsigmas*Math.exp(-comW1toM3*comW1toM3/(4*sigmaM*sigmaM)));

        Eq1ZcompW3 += chargeH31*(comW1.x(2)-H31r.x(2))/(comW1toH31*comW1toH31*comW1toH31)*((1-SpecialFunctions.erfc(comW1toH31/sqrtHMsigmas))-Math.sqrt(2)*comW1toH31/sqrtPiHMsigmas*Math.exp(-comW1toH31*comW1toH31/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq1ZcompW3 += chargeH32*(comW1.x(2)-H32r.x(2))/(comW1toH32*comW1toH32*comW1toH32)*((1-SpecialFunctions.erfc(comW1toH32/sqrtHMsigmas))-Math.sqrt(2)*comW1toH32/sqrtPiHMsigmas*Math.exp(-comW1toH32*comW1toH32/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq1ZcompW3 += chargeM3*(comW1.x(2)-M3r.x(2))/(comW1toM3*comW1toM3*comW1toM3)*((1-SpecialFunctions.erfc(comW1toM3/(2*sigmaM)))-Math.sqrt(2)*comW1toM3/sqrtPiMMsigmas*Math.exp(-comW1toM3*comW1toM3/(4*sigmaM*sigmaM)));
        

        // Contributions to sum from water#4
        Eq1XcompW4 += chargeH41*(comW1.x(0)-H41r.x(0))/(comW1toH41*comW1toH41*comW1toH41)*((1-SpecialFunctions.erfc(comW1toH41/sqrtHMsigmas))-Math.sqrt(2)*comW1toH41/sqrtPiHMsigmas*Math.exp(-comW1toH41*comW1toH41/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq1XcompW4 += chargeH42*(comW1.x(0)-H42r.x(0))/(comW1toH42*comW1toH42*comW1toH42)*((1-SpecialFunctions.erfc(comW1toH42/sqrtHMsigmas))-Math.sqrt(2)*comW1toH42/sqrtPiHMsigmas*Math.exp(-comW1toH42*comW1toH42/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq1XcompW4 += chargeM4*(comW1.x(0)-M4r.x(0))/(comW1toM4*comW1toM4*comW1toM4)*((1-SpecialFunctions.erfc(comW1toM4/(2*sigmaM)))-Math.sqrt(2)*comW1toM4/sqrtPiMMsigmas*Math.exp(-comW1toM4*comW1toM4/(4*sigmaM*sigmaM)));

        Eq1YcompW4 += chargeH41*(comW1.x(1)-H41r.x(1))/(comW1toH41*comW1toH41*comW1toH41)*((1-SpecialFunctions.erfc(comW1toH41/sqrtHMsigmas))-Math.sqrt(2)*comW1toH41/sqrtPiHMsigmas*Math.exp(-comW1toH41*comW1toH41/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq1YcompW4 += chargeH42*(comW1.x(1)-H42r.x(1))/(comW1toH42*comW1toH42*comW1toH42)*((1-SpecialFunctions.erfc(comW1toH42/sqrtHMsigmas))-Math.sqrt(2)*comW1toH42/sqrtPiHMsigmas*Math.exp(-comW1toH42*comW1toH42/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq1YcompW4 += chargeM4*(comW1.x(1)-M4r.x(1))/(comW1toM4*comW1toM4*comW1toM4)*((1-SpecialFunctions.erfc(comW1toM4/(2*sigmaM)))-Math.sqrt(2)*comW1toM4/sqrtPiMMsigmas*Math.exp(-comW1toM4*comW1toM4/(4*sigmaM*sigmaM)));

        Eq1ZcompW4 += chargeH41*(comW1.x(2)-H41r.x(2))/(comW1toH41*comW1toH41*comW1toH41)*((1-SpecialFunctions.erfc(comW1toH41/sqrtHMsigmas))-Math.sqrt(2)*comW1toH41/sqrtPiHMsigmas*Math.exp(-comW1toH41*comW1toH41/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq1ZcompW4 += chargeH42*(comW1.x(2)-H42r.x(2))/(comW1toH42*comW1toH42*comW1toH42)*((1-SpecialFunctions.erfc(comW1toH42/sqrtHMsigmas))-Math.sqrt(2)*comW1toH42/sqrtPiHMsigmas*Math.exp(-comW1toH42*comW1toH42/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq1ZcompW4 += chargeM4*(comW1.x(2)-M4r.x(2))/(comW1toM4*comW1toM4*comW1toM4)*((1-SpecialFunctions.erfc(comW1toM4/(2*sigmaM)))-Math.sqrt(2)*comW1toM4/sqrtPiMMsigmas*Math.exp(-comW1toM4*comW1toM4/(4*sigmaM*sigmaM)));

        
        Eq1.setX(0,Eq1XcompW2+Eq1XcompW3+Eq1XcompW4);
        Eq1.setX(1,Eq1YcompW2+Eq1YcompW3+Eq1YcompW4);
        Eq1.setX(2,Eq1ZcompW2+Eq1ZcompW3+Eq1ZcompW4);

                
        double Eq2XcompW1 = 0.0;
        double Eq2YcompW1 = 0.0;
        double Eq2ZcompW1 = 0.0;
        double Eq2XcompW3 = 0.0;
        double Eq2YcompW3 = 0.0;
        double Eq2ZcompW3 = 0.0;
        double Eq2XcompW4 = 0.0;
        double Eq2YcompW4 = 0.0;
        double Eq2ZcompW4 = 0.0;

        
        double comW2toH11 = Math.sqrt(comW2.Mv1Squared(H11r));
        double comW2toH12 = Math.sqrt(comW2.Mv1Squared(H12r));
        double comW2toM1 = Math.sqrt(comW2.Mv1Squared(M1r));
        double comW2toH31 = Math.sqrt(comW2.Mv1Squared(H31r));
        double comW2toH32 = Math.sqrt(comW2.Mv1Squared(H32r));
        double comW2toM3 = Math.sqrt(comW2.Mv1Squared(M3r));
        double comW2toH41 = Math.sqrt(comW2.Mv1Squared(H41r));
        double comW2toH42 = Math.sqrt(comW2.Mv1Squared(H42r));
        double comW2toM4 = Math.sqrt(comW2.Mv1Squared(M4r));

        
        // Contributions to sum from water molecule#1
        Eq2XcompW1 += chargeH11*(comW2.x(0)-H11r.x(0))/(comW2toH11*comW2toH11*comW2toH11)*((1-SpecialFunctions.erfc(comW2toH11/sqrtHMsigmas))-Math.sqrt(2)*comW2toH11/sqrtPiHMsigmas*Math.exp(-comW2toH11*comW2toH11/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq2XcompW1 += chargeH12*(comW2.x(0)-H12r.x(0))/(comW2toH12*comW2toH12*comW2toH12)*((1-SpecialFunctions.erfc(comW2toH12/sqrtHMsigmas))-Math.sqrt(2)*comW2toH12/sqrtPiHMsigmas*Math.exp(-comW2toH12*comW2toH12/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq2XcompW1 += chargeM1*(comW2.x(0)-M1r.x(0))/(comW2toM1*comW2toM1*comW2toM1)*((1-SpecialFunctions.erfc(comW2toM1/(2*sigmaM)))-Math.sqrt(2)*comW2toM1/sqrtPiMMsigmas*Math.exp(-comW2toM1*comW2toM1/(4*sigmaM*sigmaM)));

        Eq2YcompW1 += chargeH11*(comW2.x(1)-H11r.x(1))/(comW2toH11*comW2toH11*comW2toH11)*((1-SpecialFunctions.erfc(comW2toH11/sqrtHMsigmas))-Math.sqrt(2)*comW2toH11/sqrtPiHMsigmas*Math.exp(-comW2toH11*comW2toH11/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq2YcompW1 += chargeH12*(comW2.x(1)-H12r.x(1))/(comW2toH12*comW2toH12*comW2toH12)*((1-SpecialFunctions.erfc(comW2toH12/sqrtHMsigmas))-Math.sqrt(2)*comW2toH12/sqrtPiHMsigmas*Math.exp(-comW2toH12*comW2toH12/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq2YcompW1 += chargeM1*(comW2.x(1)-M1r.x(1))/(comW2toM1*comW2toM1*comW2toM1)*((1-SpecialFunctions.erfc(comW2toM1/(2*sigmaM)))-Math.sqrt(2)*comW2toM1/sqrtPiMMsigmas*Math.exp(-comW2toM1*comW2toM1/(4*sigmaM*sigmaM)));

        Eq2ZcompW1 += chargeH11*(comW2.x(2)-H11r.x(2))/(comW2toH11*comW2toH11*comW2toH11)*((1-SpecialFunctions.erfc(comW2toH11/sqrtHMsigmas))-Math.sqrt(2)*comW2toH11/sqrtPiHMsigmas*Math.exp(-comW2toH11*comW2toH11/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq2ZcompW1 += chargeH12*(comW2.x(2)-H12r.x(2))/(comW2toH12*comW2toH12*comW2toH12)*((1-SpecialFunctions.erfc(comW2toH12/sqrtHMsigmas))-Math.sqrt(2)*comW2toH12/sqrtPiHMsigmas*Math.exp(-comW2toH12*comW2toH12/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq2ZcompW1 += chargeM1*(comW2.x(2)-M1r.x(2))/(comW2toM1*comW2toM1*comW2toM1)*((1-SpecialFunctions.erfc(comW2toM1/(2*sigmaM)))-Math.sqrt(2)*comW2toM1/sqrtPiMMsigmas*Math.exp(-comW2toM1*comW2toM1/(4*sigmaM*sigmaM)));


        // Contributions to sum from water molecule#3
        Eq2XcompW3 += chargeH31*(comW2.x(0)-H31r.x(0))/(comW2toH31*comW2toH31*comW2toH31)*((1-SpecialFunctions.erfc(comW2toH31/sqrtHMsigmas))-Math.sqrt(2)*comW2toH31/sqrtPiHMsigmas*Math.exp(-comW2toH31*comW2toH31/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq2XcompW3 += chargeH32*(comW2.x(0)-H32r.x(0))/(comW2toH32*comW2toH32*comW2toH32)*((1-SpecialFunctions.erfc(comW2toH32/sqrtHMsigmas))-Math.sqrt(2)*comW2toH32/sqrtPiHMsigmas*Math.exp(-comW2toH32*comW2toH32/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq2XcompW3 += chargeM3*(comW2.x(0)-M3r.x(0))/(comW2toM3*comW2toM3*comW2toM3)*((1-SpecialFunctions.erfc(comW2toM3/(2*sigmaM)))-Math.sqrt(2)*comW2toM3/sqrtPiMMsigmas*Math.exp(-comW2toM3*comW2toM3/(4*sigmaM*sigmaM)));

        Eq2YcompW3 += chargeH31*(comW2.x(1)-H31r.x(1))/(comW2toH31*comW2toH31*comW2toH31)*((1-SpecialFunctions.erfc(comW2toH31/sqrtHMsigmas))-Math.sqrt(2)*comW2toH31/sqrtPiHMsigmas*Math.exp(-comW2toH31*comW2toH31/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq2YcompW3 += chargeH32*(comW2.x(1)-H32r.x(1))/(comW2toH32*comW2toH32*comW2toH32)*((1-SpecialFunctions.erfc(comW2toH32/sqrtHMsigmas))-Math.sqrt(2)*comW2toH32/sqrtPiHMsigmas*Math.exp(-comW2toH32*comW2toH32/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq2YcompW3 += chargeM3*(comW2.x(1)-M3r.x(1))/(comW2toM3*comW2toM3*comW2toM3)*((1-SpecialFunctions.erfc(comW2toM3/(2*sigmaM)))-Math.sqrt(2)*comW2toM3/sqrtPiMMsigmas*Math.exp(-comW2toM3*comW2toM3/(4*sigmaM*sigmaM)));

        Eq2ZcompW3 += chargeH31*(comW2.x(2)-H31r.x(2))/(comW2toH31*comW2toH31*comW2toH31)*((1-SpecialFunctions.erfc(comW2toH31/sqrtHMsigmas))-Math.sqrt(2)*comW2toH31/sqrtPiHMsigmas*Math.exp(-comW2toH31*comW2toH31/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq2ZcompW3 += chargeH32*(comW2.x(2)-H32r.x(2))/(comW2toH32*comW2toH32*comW2toH32)*((1-SpecialFunctions.erfc(comW2toH32/sqrtHMsigmas))-Math.sqrt(2)*comW2toH32/sqrtPiHMsigmas*Math.exp(-comW2toH32*comW2toH32/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq2ZcompW3 += chargeM3*(comW2.x(2)-M3r.x(2))/(comW2toM3*comW2toM3*comW2toM3)*((1-SpecialFunctions.erfc(comW2toM3/(2*sigmaM)))-Math.sqrt(2)*comW2toM3/sqrtPiMMsigmas*Math.exp(-comW2toM3*comW2toM3/(4*sigmaM*sigmaM)));
        

        // Contributions to sum from water molecule#4
        Eq2XcompW4 += chargeH41*(comW2.x(0)-H41r.x(0))/(comW2toH41*comW2toH41*comW2toH41)*((1-SpecialFunctions.erfc(comW2toH41/sqrtHMsigmas))-Math.sqrt(2)*comW2toH41/sqrtPiHMsigmas*Math.exp(-comW2toH41*comW2toH41/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq2XcompW4 += chargeH42*(comW2.x(0)-H42r.x(0))/(comW2toH42*comW2toH42*comW2toH42)*((1-SpecialFunctions.erfc(comW2toH42/sqrtHMsigmas))-Math.sqrt(2)*comW2toH42/sqrtPiHMsigmas*Math.exp(-comW2toH42*comW2toH42/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq2XcompW4 += chargeM4*(comW2.x(0)-M4r.x(0))/(comW2toM4*comW2toM4*comW2toM4)*((1-SpecialFunctions.erfc(comW2toM4/(2*sigmaM)))-Math.sqrt(2)*comW2toM4/sqrtPiMMsigmas*Math.exp(-comW2toM4*comW2toM4/(4*sigmaM*sigmaM)));

        Eq2YcompW4 += chargeH41*(comW2.x(1)-H41r.x(1))/(comW2toH41*comW2toH41*comW2toH41)*((1-SpecialFunctions.erfc(comW2toH41/sqrtHMsigmas))-Math.sqrt(2)*comW2toH41/sqrtPiHMsigmas*Math.exp(-comW2toH41*comW2toH41/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq2YcompW4 += chargeH42*(comW2.x(1)-H42r.x(1))/(comW2toH42*comW2toH42*comW2toH42)*((1-SpecialFunctions.erfc(comW2toH42/sqrtHMsigmas))-Math.sqrt(2)*comW2toH42/sqrtPiHMsigmas*Math.exp(-comW2toH42*comW2toH42/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq2YcompW4 += chargeM4*(comW2.x(1)-M4r.x(1))/(comW2toM4*comW2toM4*comW2toM4)*((1-SpecialFunctions.erfc(comW2toM4/(2*sigmaM)))-Math.sqrt(2)*comW2toM4/sqrtPiMMsigmas*Math.exp(-comW2toM4*comW2toM4/(4*sigmaM*sigmaM)));

        Eq2ZcompW4 += chargeH41*(comW2.x(2)-H41r.x(2))/(comW2toH41*comW2toH41*comW2toH41)*((1-SpecialFunctions.erfc(comW2toH41/sqrtHMsigmas))-Math.sqrt(2)*comW2toH41/sqrtPiHMsigmas*Math.exp(-comW2toH41*comW2toH41/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq2ZcompW4 += chargeH42*(comW2.x(2)-H42r.x(2))/(comW2toH42*comW2toH42*comW2toH42)*((1-SpecialFunctions.erfc(comW2toH42/sqrtHMsigmas))-Math.sqrt(2)*comW2toH42/sqrtPiHMsigmas*Math.exp(-comW2toH42*comW2toH42/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq2ZcompW4 += chargeM4*(comW2.x(2)-M4r.x(2))/(comW2toM4*comW2toM4*comW2toM4)*((1-SpecialFunctions.erfc(comW2toM4/(2*sigmaM)))-Math.sqrt(2)*comW2toM4/sqrtPiMMsigmas*Math.exp(-comW2toM4*comW2toM4/(4*sigmaM*sigmaM)));

        
        Eq2.setX(0,Eq2XcompW1+Eq2XcompW3+Eq2XcompW4);
        Eq2.setX(1,Eq2YcompW1+Eq2YcompW3+Eq2YcompW4);
        Eq2.setX(2,Eq2ZcompW1+Eq2ZcompW3+Eq2ZcompW4);


        // Find Eq3
        double Eq3XcompW1	 = 0.0;
        double Eq3YcompW1 = 0.0;
        double Eq3ZcompW1 = 0.0;
        double Eq3XcompW2 = 0.0;
        double Eq3YcompW2 = 0.0;
        double Eq3ZcompW2 = 0.0;
        double Eq3XcompW4 = 0.0;
        double Eq3YcompW4 = 0.0;
        double Eq3ZcompW4 = 0.0;

        
        
        
        double comW3toH11 = Math.sqrt(comW3.Mv1Squared(H11r));
        double comW3toH12 = Math.sqrt(comW3.Mv1Squared(H12r));
        double comW3toM1 = Math.sqrt(comW3.Mv1Squared(M1r));

        double comW3toH21 = Math.sqrt(comW3.Mv1Squared(H21r));
        double comW3toH22 = Math.sqrt(comW3.Mv1Squared(H22r));
        double comW3toM2 = Math.sqrt(comW3.Mv1Squared(M2r));

        double comW3toH41 = Math.sqrt(comW3.Mv1Squared(H41r));
        double comW3toH42 = Math.sqrt(comW3.Mv1Squared(H42r));
        double comW3toM4 = Math.sqrt(comW3.Mv1Squared(M4r));

        
        // Contributions to sum from water molecule#1       
        Eq3XcompW1 += chargeH11*(comW3.x(0)-H11r.x(0))/(comW3toH11*comW3toH11*comW3toH11)*((1-SpecialFunctions.erfc(comW3toH11/sqrtHMsigmas))-Math.sqrt(2)*comW3toH11/sqrtPiHMsigmas*Math.exp(-comW3toH11*comW3toH11/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq3XcompW1 += chargeH12*(comW3.x(0)-H12r.x(0))/(comW3toH12*comW3toH12*comW3toH12)*((1-SpecialFunctions.erfc(comW3toH12/sqrtHMsigmas))-Math.sqrt(2)*comW3toH12/sqrtPiHMsigmas*Math.exp(-comW3toH12*comW3toH12/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq3XcompW1 += chargeM1*(comW3.x(0)-M1r.x(0))/(comW3toM1*comW3toM1*comW3toM1)*((1-SpecialFunctions.erfc(comW3toM1/(2*sigmaM)))-Math.sqrt(2)*comW3toM1/sqrtPiMMsigmas*Math.exp(-comW3toM1*comW3toM1/(4*sigmaM*sigmaM)));

        Eq3YcompW1 += chargeH11*(comW3.x(1)-H11r.x(1))/(comW3toH11*comW3toH11*comW3toH11)*((1-SpecialFunctions.erfc(comW3toH11/sqrtHMsigmas))-Math.sqrt(2)*comW3toH11/sqrtPiHMsigmas*Math.exp(-comW3toH11*comW3toH11/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq3YcompW1 += chargeH12*(comW3.x(1)-H12r.x(1))/(comW3toH12*comW3toH12*comW3toH12)*((1-SpecialFunctions.erfc(comW3toH12/sqrtHMsigmas))-Math.sqrt(2)*comW3toH12/sqrtPiHMsigmas*Math.exp(-comW3toH12*comW3toH12/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq3YcompW1 += chargeM1*(comW3.x(1)-M1r.x(1))/(comW3toM1*comW3toM1*comW3toM1)*((1-SpecialFunctions.erfc(comW3toM1/(2*sigmaM)))-Math.sqrt(2)*comW3toM1/sqrtPiMMsigmas*Math.exp(-comW3toM1*comW3toM1/(4*sigmaM*sigmaM)));

        Eq3ZcompW1 += chargeH11*(comW3.x(2)-H11r.x(2))/(comW3toH11*comW3toH11*comW3toH11)*((1-SpecialFunctions.erfc(comW3toH11/sqrtHMsigmas))-Math.sqrt(2)*comW3toH11/sqrtPiHMsigmas*Math.exp(-comW3toH11*comW3toH11/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq3ZcompW1 += chargeH12*(comW3.x(2)-H12r.x(2))/(comW3toH12*comW3toH12*comW3toH12)*((1-SpecialFunctions.erfc(comW3toH12/sqrtHMsigmas))-Math.sqrt(2)*comW3toH12/sqrtPiHMsigmas*Math.exp(-comW3toH12*comW3toH12/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq3ZcompW1 += chargeM1*(comW3.x(2)-M1r.x(2))/(comW3toM1*comW3toM1*comW3toM1)*((1-SpecialFunctions.erfc(comW3toM1/(2*sigmaM)))-Math.sqrt(2)*comW3toM1/sqrtPiMMsigmas*Math.exp(-comW3toM1*comW3toM1/(4*sigmaM*sigmaM)));

        
        // Contributions to sum from water molecule#2
        Eq3XcompW2 += chargeH21*(comW3.x(0)-H21r.x(0))/(comW3toH21*comW3toH21*comW3toH21)*((1-SpecialFunctions.erfc(comW3toH21/sqrtHMsigmas))-Math.sqrt(2)*comW3toH21/sqrtPiHMsigmas*Math.exp(-comW3toH21*comW3toH21/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq3XcompW2 += chargeH22*(comW3.x(0)-H22r.x(0))/(comW3toH22*comW3toH22*comW3toH22)*((1-SpecialFunctions.erfc(comW3toH22/sqrtHMsigmas))-Math.sqrt(2)*comW3toH22/sqrtPiHMsigmas*Math.exp(-comW3toH22*comW3toH22/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq3XcompW2 += chargeM2*(comW3.x(0)-M2r.x(0))/(comW3toM2*comW3toM2*comW3toM2)*((1-SpecialFunctions.erfc(comW3toM2/(2*sigmaM)))-Math.sqrt(2)*comW3toM2/sqrtPiMMsigmas*Math.exp(-comW3toM2*comW3toM2/(4*sigmaM*sigmaM)));

        Eq3YcompW2 += chargeH21*(comW3.x(1)-H21r.x(1))/(comW3toH21*comW3toH21*comW3toH21)*((1-SpecialFunctions.erfc(comW3toH21/sqrtHMsigmas))-Math.sqrt(2)*comW3toH21/sqrtPiHMsigmas*Math.exp(-comW3toH21*comW3toH21/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq3YcompW2 += chargeH22*(comW3.x(1)-H22r.x(1))/(comW3toH22*comW3toH22*comW3toH22)*((1-SpecialFunctions.erfc(comW3toH22/sqrtHMsigmas))-Math.sqrt(2)*comW3toH22/sqrtPiHMsigmas*Math.exp(-comW3toH22*comW3toH22/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq3YcompW2 += chargeM2*(comW3.x(1)-M2r.x(1))/(comW3toM2*comW3toM2*comW3toM2)*((1-SpecialFunctions.erfc(comW3toM2/(2*sigmaM)))-Math.sqrt(2)*comW3toM2/sqrtPiMMsigmas*Math.exp(-comW3toM2*comW3toM2/(4*sigmaM*sigmaM)));

        Eq3ZcompW2 += chargeH21*(comW3.x(2)-H21r.x(2))/(comW3toH21*comW3toH21*comW3toH21)*((1-SpecialFunctions.erfc(comW3toH21/sqrtHMsigmas))-Math.sqrt(2)*comW3toH21/sqrtPiHMsigmas*Math.exp(-comW3toH21*comW3toH21/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq3ZcompW2 += chargeH22*(comW3.x(2)-H22r.x(2))/(comW3toH22*comW3toH22*comW3toH22)*((1-SpecialFunctions.erfc(comW3toH22/sqrtHMsigmas))-Math.sqrt(2)*comW3toH22/sqrtPiHMsigmas*Math.exp(-comW3toH22*comW3toH22/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq3ZcompW2 += chargeM2*(comW3.x(2)-M2r.x(2))/(comW3toM2*comW3toM2*comW3toM2)*((1-SpecialFunctions.erfc(comW3toM2/(2*sigmaM)))-Math.sqrt(2)*comW3toM2/sqrtPiMMsigmas*Math.exp(-comW3toM2*comW3toM2/(4*sigmaM*sigmaM)));
        

        // Contributions to sum from water molecule#4
        Eq3XcompW4 += chargeH41*(comW3.x(0)-H41r.x(0))/(comW3toH41*comW3toH41*comW3toH41)*((1-SpecialFunctions.erfc(comW3toH41/sqrtHMsigmas))-Math.sqrt(2)*comW3toH41/sqrtPiHMsigmas*Math.exp(-comW3toH41*comW3toH41/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq3XcompW4 += chargeH42*(comW3.x(0)-H42r.x(0))/(comW3toH42*comW3toH42*comW3toH42)*((1-SpecialFunctions.erfc(comW3toH42/sqrtHMsigmas))-Math.sqrt(2)*comW3toH42/sqrtPiHMsigmas*Math.exp(-comW3toH42*comW3toH42/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq3XcompW4 += chargeM4*(comW3.x(0)-M4r.x(0))/(comW3toM4*comW3toM4*comW3toM4)*((1-SpecialFunctions.erfc(comW3toM4/(2*sigmaM)))-Math.sqrt(2)*comW3toM4/sqrtPiMMsigmas*Math.exp(-comW3toM4*comW3toM4/(4*sigmaM*sigmaM)));

        Eq3YcompW4 += chargeH41*(comW3.x(1)-H41r.x(1))/(comW3toH41*comW3toH41*comW3toH41)*((1-SpecialFunctions.erfc(comW3toH41/sqrtHMsigmas))-Math.sqrt(2)*comW3toH41/sqrtPiHMsigmas*Math.exp(-comW3toH41*comW3toH41/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq3YcompW4 += chargeH42*(comW3.x(1)-H42r.x(1))/(comW3toH42*comW3toH42*comW3toH42)*((1-SpecialFunctions.erfc(comW3toH42/sqrtHMsigmas))-Math.sqrt(2)*comW3toH42/sqrtPiHMsigmas*Math.exp(-comW3toH42*comW3toH42/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq3YcompW4 += chargeM4*(comW3.x(1)-M4r.x(1))/(comW3toM4*comW3toM4*comW3toM4)*((1-SpecialFunctions.erfc(comW3toM4/(2*sigmaM)))-Math.sqrt(2)*comW3toM4/sqrtPiMMsigmas*Math.exp(-comW3toM4*comW3toM4/(4*sigmaM*sigmaM)));

        Eq3ZcompW4 += chargeH41*(comW3.x(2)-H41r.x(2))/(comW3toH41*comW3toH41*comW3toH41)*((1-SpecialFunctions.erfc(comW3toH41/sqrtHMsigmas))-Math.sqrt(2)*comW3toH41/sqrtPiHMsigmas*Math.exp(-comW3toH41*comW3toH41/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq3ZcompW4 += chargeH42*(comW3.x(2)-H42r.x(2))/(comW3toH42*comW3toH42*comW3toH42)*((1-SpecialFunctions.erfc(comW3toH42/sqrtHMsigmas))-Math.sqrt(2)*comW3toH42/sqrtPiHMsigmas*Math.exp(-comW3toH42*comW3toH42/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq3ZcompW4 += chargeM4*(comW3.x(2)-M4r.x(2))/(comW3toM4*comW3toM4*comW3toM4)*((1-SpecialFunctions.erfc(comW3toM4/(2*sigmaM)))-Math.sqrt(2)*comW3toM4/sqrtPiMMsigmas*Math.exp(-comW3toM4*comW3toM4/(4*sigmaM*sigmaM)));

        
        Eq3.setX(0,Eq3XcompW1+Eq3XcompW2+Eq3XcompW4);
        Eq3.setX(1,Eq3YcompW1+Eq3YcompW2+Eq3YcompW4);
        Eq3.setX(2,Eq3ZcompW1+Eq3ZcompW2+Eq3ZcompW4);

        
        
        
        
        // Find Eq4
        double Eq4XcompW1	 = 0.0;
        double Eq4YcompW1 = 0.0;
        double Eq4ZcompW1 = 0.0;
        double Eq4XcompW2 = 0.0;
        double Eq4YcompW2 = 0.0;
        double Eq4ZcompW2 = 0.0;
        double Eq4XcompW3 = 0.0;
        double Eq4YcompW3 = 0.0;
        double Eq4ZcompW3 = 0.0;

        
        
        
        double comW4toH11 = Math.sqrt(comW4.Mv1Squared(H11r));
        double comW4toH12 = Math.sqrt(comW4.Mv1Squared(H12r));
        double comW4toM1 = Math.sqrt(comW4.Mv1Squared(M1r));

        double comW4toH21 = Math.sqrt(comW4.Mv1Squared(H21r));
        double comW4toH22 = Math.sqrt(comW4.Mv1Squared(H22r));
        double comW4toM2 = Math.sqrt(comW4.Mv1Squared(M2r));

        double comW4toH31 = Math.sqrt(comW4.Mv1Squared(H31r));
        double comW4toH32 = Math.sqrt(comW4.Mv1Squared(H32r));
        double comW4toM3 = Math.sqrt(comW4.Mv1Squared(M3r));

        
        // Contributions to sum from water molecule#1       
        Eq4XcompW1 += chargeH11*(comW4.x(0)-H11r.x(0))/(comW4toH11*comW4toH11*comW4toH11)*((1-SpecialFunctions.erfc(comW4toH11/sqrtHMsigmas))-Math.sqrt(2)*comW4toH11/sqrtPiHMsigmas*Math.exp(-comW4toH11*comW4toH11/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq4XcompW1 += chargeH12*(comW4.x(0)-H12r.x(0))/(comW4toH12*comW4toH12*comW4toH12)*((1-SpecialFunctions.erfc(comW4toH12/sqrtHMsigmas))-Math.sqrt(2)*comW4toH12/sqrtPiHMsigmas*Math.exp(-comW4toH12*comW4toH12/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq4XcompW1 += chargeM1*(comW4.x(0)-M1r.x(0))/(comW4toM1*comW4toM1*comW4toM1)*((1-SpecialFunctions.erfc(comW4toM1/(2*sigmaM)))-Math.sqrt(2)*comW4toM1/sqrtPiMMsigmas*Math.exp(-comW4toM1*comW4toM1/(4*sigmaM*sigmaM)));

        Eq4YcompW1 += chargeH11*(comW4.x(1)-H11r.x(1))/(comW4toH11*comW4toH11*comW4toH11)*((1-SpecialFunctions.erfc(comW4toH11/sqrtHMsigmas))-Math.sqrt(2)*comW4toH11/sqrtPiHMsigmas*Math.exp(-comW4toH11*comW4toH11/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq4YcompW1 += chargeH12*(comW4.x(1)-H12r.x(1))/(comW4toH12*comW4toH12*comW4toH12)*((1-SpecialFunctions.erfc(comW4toH12/sqrtHMsigmas))-Math.sqrt(2)*comW4toH12/sqrtPiHMsigmas*Math.exp(-comW4toH12*comW4toH12/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq4YcompW1 += chargeM1*(comW4.x(1)-M1r.x(1))/(comW4toM1*comW4toM1*comW4toM1)*((1-SpecialFunctions.erfc(comW4toM1/(2*sigmaM)))-Math.sqrt(2)*comW4toM1/sqrtPiMMsigmas*Math.exp(-comW4toM1*comW4toM1/(4*sigmaM*sigmaM)));

        Eq4ZcompW1 += chargeH11*(comW4.x(2)-H11r.x(2))/(comW4toH11*comW4toH11*comW4toH11)*((1-SpecialFunctions.erfc(comW4toH11/sqrtHMsigmas))-Math.sqrt(2)*comW4toH11/sqrtPiHMsigmas*Math.exp(-comW4toH11*comW4toH11/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq4ZcompW1 += chargeH12*(comW4.x(2)-H12r.x(2))/(comW4toH12*comW4toH12*comW4toH12)*((1-SpecialFunctions.erfc(comW4toH12/sqrtHMsigmas))-Math.sqrt(2)*comW4toH12/sqrtPiHMsigmas*Math.exp(-comW4toH12*comW4toH12/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq4ZcompW1 += chargeM1*(comW4.x(2)-M1r.x(2))/(comW4toM1*comW4toM1*comW4toM1)*((1-SpecialFunctions.erfc(comW4toM1/(2*sigmaM)))-Math.sqrt(2)*comW4toM1/sqrtPiMMsigmas*Math.exp(-comW4toM1*comW4toM1/(4*sigmaM*sigmaM)));

        
        // Contributions to sum from water molecule#2
        Eq4XcompW2 += chargeH21*(comW4.x(0)-H21r.x(0))/(comW4toH21*comW4toH21*comW4toH21)*((1-SpecialFunctions.erfc(comW4toH21/sqrtHMsigmas))-Math.sqrt(2)*comW4toH21/sqrtPiHMsigmas*Math.exp(-comW4toH21*comW4toH21/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq4XcompW2 += chargeH22*(comW4.x(0)-H22r.x(0))/(comW4toH22*comW4toH22*comW4toH22)*((1-SpecialFunctions.erfc(comW4toH22/sqrtHMsigmas))-Math.sqrt(2)*comW4toH22/sqrtPiHMsigmas*Math.exp(-comW4toH22*comW4toH22/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq4XcompW2 += chargeM2*(comW4.x(0)-M2r.x(0))/(comW4toM2*comW4toM2*comW4toM2)*((1-SpecialFunctions.erfc(comW4toM2/(2*sigmaM)))-Math.sqrt(2)*comW4toM2/sqrtPiMMsigmas*Math.exp(-comW4toM2*comW4toM2/(4*sigmaM*sigmaM)));

        Eq4YcompW2 += chargeH21*(comW4.x(1)-H21r.x(1))/(comW4toH21*comW4toH21*comW4toH21)*((1-SpecialFunctions.erfc(comW4toH21/sqrtHMsigmas))-Math.sqrt(2)*comW4toH21/sqrtPiHMsigmas*Math.exp(-comW4toH21*comW4toH21/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq4YcompW2 += chargeH22*(comW4.x(1)-H22r.x(1))/(comW4toH22*comW4toH22*comW4toH22)*((1-SpecialFunctions.erfc(comW4toH22/sqrtHMsigmas))-Math.sqrt(2)*comW4toH22/sqrtPiHMsigmas*Math.exp(-comW4toH22*comW4toH22/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq4YcompW2 += chargeM2*(comW4.x(1)-M2r.x(1))/(comW4toM2*comW4toM2*comW4toM2)*((1-SpecialFunctions.erfc(comW4toM2/(2*sigmaM)))-Math.sqrt(2)*comW4toM2/sqrtPiMMsigmas*Math.exp(-comW4toM2*comW4toM2/(4*sigmaM*sigmaM)));

        Eq4ZcompW2 += chargeH21*(comW4.x(2)-H21r.x(2))/(comW4toH21*comW4toH21*comW4toH21)*((1-SpecialFunctions.erfc(comW4toH21/sqrtHMsigmas))-Math.sqrt(2)*comW4toH21/sqrtPiHMsigmas*Math.exp(-comW4toH21*comW4toH21/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq4ZcompW2 += chargeH22*(comW4.x(2)-H22r.x(2))/(comW4toH22*comW4toH22*comW4toH22)*((1-SpecialFunctions.erfc(comW4toH22/sqrtHMsigmas))-Math.sqrt(2)*comW4toH22/sqrtPiHMsigmas*Math.exp(-comW4toH22*comW4toH22/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq4ZcompW2 += chargeM2*(comW4.x(2)-M2r.x(2))/(comW4toM2*comW4toM2*comW4toM2)*((1-SpecialFunctions.erfc(comW4toM2/(2*sigmaM)))-Math.sqrt(2)*comW4toM2/sqrtPiMMsigmas*Math.exp(-comW4toM2*comW4toM2/(4*sigmaM*sigmaM)));
        

        // Contributions to sum from water molecule#3
        Eq4XcompW3 += chargeH31*(comW4.x(0)-H31r.x(0))/(comW4toH31*comW4toH31*comW4toH31)*((1-SpecialFunctions.erfc(comW4toH31/sqrtHMsigmas))-Math.sqrt(2)*comW4toH31/sqrtPiHMsigmas*Math.exp(-comW4toH31*comW4toH31/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq4XcompW3 += chargeH32*(comW4.x(0)-H32r.x(0))/(comW4toH32*comW4toH32*comW4toH32)*((1-SpecialFunctions.erfc(comW4toH32/sqrtHMsigmas))-Math.sqrt(2)*comW4toH32/sqrtPiHMsigmas*Math.exp(-comW4toH32*comW4toH32/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq4XcompW3 += chargeM3*(comW4.x(0)-M3r.x(0))/(comW4toM3*comW4toM3*comW4toM3)*((1-SpecialFunctions.erfc(comW4toM3/(2*sigmaM)))-Math.sqrt(2)*comW4toM3/sqrtPiMMsigmas*Math.exp(-comW4toM3*comW4toM3/(4*sigmaM*sigmaM)));

        Eq4YcompW3 += chargeH31*(comW4.x(1)-H31r.x(1))/(comW4toH31*comW4toH31*comW4toH31)*((1-SpecialFunctions.erfc(comW4toH31/sqrtHMsigmas))-Math.sqrt(2)*comW4toH31/sqrtPiHMsigmas*Math.exp(-comW4toH31*comW4toH31/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq4YcompW3 += chargeH32*(comW4.x(1)-H32r.x(1))/(comW4toH32*comW4toH32*comW4toH32)*((1-SpecialFunctions.erfc(comW4toH32/sqrtHMsigmas))-Math.sqrt(2)*comW4toH32/sqrtPiHMsigmas*Math.exp(-comW4toH32*comW4toH32/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq4YcompW3 += chargeM3*(comW4.x(1)-M3r.x(1))/(comW4toM3*comW4toM3*comW4toM3)*((1-SpecialFunctions.erfc(comW4toM3/(2*sigmaM)))-Math.sqrt(2)*comW4toM3/sqrtPiMMsigmas*Math.exp(-comW4toM3*comW4toM3/(4*sigmaM*sigmaM)));

        Eq4ZcompW3 += chargeH31*(comW4.x(2)-H31r.x(2))/(comW4toH31*comW4toH31*comW4toH31)*((1-SpecialFunctions.erfc(comW4toH31/sqrtHMsigmas))-Math.sqrt(2)*comW4toH31/sqrtPiHMsigmas*Math.exp(-comW4toH31*comW4toH31/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq4ZcompW3 += chargeH32*(comW4.x(2)-H32r.x(2))/(comW4toH32*comW4toH32*comW4toH32)*((1-SpecialFunctions.erfc(comW4toH32/sqrtHMsigmas))-Math.sqrt(2)*comW4toH32/sqrtPiHMsigmas*Math.exp(-comW4toH32*comW4toH32/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq4ZcompW3 += chargeM3*(comW4.x(2)-M3r.x(2))/(comW4toM3*comW4toM3*comW4toM3)*((1-SpecialFunctions.erfc(comW4toM3/(2*sigmaM)))-Math.sqrt(2)*comW4toM3/sqrtPiMMsigmas*Math.exp(-comW4toM3*comW4toM3/(4*sigmaM*sigmaM)));

        
        Eq4.setX(0,Eq4XcompW1+Eq4XcompW2+Eq4XcompW3);
        Eq4.setX(1,Eq4YcompW1+Eq4YcompW2+Eq4YcompW3);
        Eq4.setX(2,Eq4ZcompW1+Eq4ZcompW2+Eq4ZcompW3);

        
        
        /*
         * Finding the tensor used to relate the induced dipole moment Pi with the induced electric field Epi.
         * kmb, 8/9/06
         */
        
        double r12 = Math.sqrt(comW1.Mv1Squared(comW2));
        double r13 = Math.sqrt(comW1.Mv1Squared(comW3));
        double r14 = Math.sqrt(comW1.Mv1Squared(comW4));
        double r23 = Math.sqrt(comW2.Mv1Squared(comW3));
        double r24 = Math.sqrt(comW2.Mv1Squared(comW4));
        double r34 = Math.sqrt(comW3.Mv1Squared(comW4));
        
        Vector3D r12Vector = new Vector3D();
        Vector3D r13Vector = new Vector3D();
        Vector3D r14Vector = new Vector3D();
        Vector3D r23Vector = new Vector3D();
        Vector3D r24Vector = new Vector3D();
        Vector3D r34Vector = new Vector3D();
        
        r12Vector.Ev1Mv2(comW1,comW2);  // is this the correct direction? kmb, 8/7/06 / Direction doesn't matter; kmb, 8/10/06
        r13Vector.Ev1Mv2(comW1,comW3);  // is this the correct direction? kmb, 8/7/06 / Direction doesn't matter; kmb, 8/10/06
        r14Vector.Ev1Mv2(comW1,comW4);  // is this the correct direction? kmb, 8/7/06 / Direction doesn't matter; kmb, 8/10/06
        r23Vector.Ev1Mv2(comW2,comW3);  // is this the correct direction? kmb, 8/7/06 / Direction doesn't matter; kmb, 8/10/06
        r24Vector.Ev1Mv2(comW2,comW4);  // is this the correct direction? kmb, 8/7/06 / Direction doesn't matter; kmb, 8/10/06
        r34Vector.Ev1Mv2(comW3,comW4);  // is this the correct direction? kmb, 8/7/06 / Direction doesn't matter; kmb, 8/10/06
        
        double f12 = (1-SpecialFunctions.erfc(r12/(2*sigmaM)))-(r12/(sigmaM*Math.sqrt(Math.PI)) + (r12*r12*r12)/(6*Math.sqrt(Math.PI)*sigmaM*sigmaM*sigmaM))*Math.exp(-r12*r12/(4*sigmaM*sigmaM));
        double f13 = (1-SpecialFunctions.erfc(r13/(2*sigmaM)))-(r13/(sigmaM*Math.sqrt(Math.PI)) + (r13*r13*r13)/(6*Math.sqrt(Math.PI)*sigmaM*sigmaM*sigmaM))*Math.exp(-r13*r13/(4*sigmaM*sigmaM));
        double f14 = (1-SpecialFunctions.erfc(r14/(2*sigmaM)))-(r14/(sigmaM*Math.sqrt(Math.PI)) + (r14*r14*r14)/(6*Math.sqrt(Math.PI)*sigmaM*sigmaM*sigmaM))*Math.exp(-r14*r14/(4*sigmaM*sigmaM));
        double f23 = (1-SpecialFunctions.erfc(r23/(2*sigmaM)))-(r23/(sigmaM*Math.sqrt(Math.PI)) + (r23*r23*r23)/(6*Math.sqrt(Math.PI)*sigmaM*sigmaM*sigmaM))*Math.exp(-r23*r23/(4*sigmaM*sigmaM));
        double f24 = (1-SpecialFunctions.erfc(r24/(2*sigmaM)))-(r24/(sigmaM*Math.sqrt(Math.PI)) + (r24*r24*r24)/(6*Math.sqrt(Math.PI)*sigmaM*sigmaM*sigmaM))*Math.exp(-r24*r24/(4*sigmaM*sigmaM));
        double f34 = (1-SpecialFunctions.erfc(r34/(2*sigmaM)))-(r34/(sigmaM*Math.sqrt(Math.PI)) + (r34*r34*r34)/(6*Math.sqrt(Math.PI)*sigmaM*sigmaM*sigmaM))*Math.exp(-r34*r34/(4*sigmaM*sigmaM));
        
        double g12 = (1-SpecialFunctions.erfc(r12/(2*sigmaM)))-(r12/(sigmaM*Math.sqrt(Math.PI)))*Math.exp(-r12*r12/(4*sigmaM*sigmaM));
        double g13 = (1-SpecialFunctions.erfc(r13/(2*sigmaM)))-(r13/(sigmaM*Math.sqrt(Math.PI)))*Math.exp(-r13*r13/(4*sigmaM*sigmaM));
        double g14 = (1-SpecialFunctions.erfc(r14/(2*sigmaM)))-(r14/(sigmaM*Math.sqrt(Math.PI)))*Math.exp(-r14*r14/(4*sigmaM*sigmaM));
        double g23 = (1-SpecialFunctions.erfc(r23/(2*sigmaM)))-(r23/(sigmaM*Math.sqrt(Math.PI)))*Math.exp(-r23*r23/(4*sigmaM*sigmaM));
        double g24 = (1-SpecialFunctions.erfc(r24/(2*sigmaM)))-(r24/(sigmaM*Math.sqrt(Math.PI)))*Math.exp(-r24*r24/(4*sigmaM*sigmaM));
        double g34 = (1-SpecialFunctions.erfc(r34/(2*sigmaM)))-(r34/(sigmaM*Math.sqrt(Math.PI)))*Math.exp(-r34*r34/(4*sigmaM*sigmaM));
        
        // Filling the unit matrix I
        
/*        double[][] I = new double[3][3];
                
        int i = 0;
        int j = 0;
        
        while (i < 3) {
        		while (j < 3) {
        			I[i][j] = 1;
            		j = j + 1;
        		}
        		i = i + 1;
        }
*/        
        Tensor3D I = new Tensor3D();
        
        I.E(1);
        
//        double[][] T12 = new double[3][3];

        Tensor3D T12 = new Tensor3D();
        
        T12.PEv1v2(r12Vector,r12Vector);
        T12.TE(3*f12/(r12*r12));
        
        I.TE(g12);
        
        T12.ME(I);
        T12.TE(1/(r12*r12*r12));
        
        // T12 = T21, so I can get by for now in the case of B2!

        
        I.E(1);
        
//        double[][] T12 = new double[3][3];

        Tensor3D T13 = new Tensor3D();
        
        T13.PEv1v2(r13Vector,r13Vector);
        T13.TE(3*f13/(r13*r13));
        
        I.TE(g13);
        
        T13.ME(I);
        T13.TE(1/(r13*r13*r13));
        

        I.E(1);
        
//        double[][] T12 = new double[3][3];

        Tensor3D T14 = new Tensor3D();
        
        T14.PEv1v2(r14Vector,r14Vector);
        T14.TE(3*f14/(r14*r14));
        
        I.TE(g14);
        
        T14.ME(I);
        T14.TE(1/(r14*r14*r14));

        
        
        
        I.E(1);
        
//        double[][] T12 = new double[3][3];

        Tensor3D T23 = new Tensor3D();
        
        T23.PEv1v2(r23Vector,r23Vector);
        T23.TE(3*f23/(r23*r23));
        
        I.TE(g23);
        
        T23.ME(I);
        T23.TE(1/(r23*r23*r23));


        I.E(1);
        
//        double[][] T12 = new double[3][3];

        Tensor3D T24 = new Tensor3D();
        
        T24.PEv1v2(r24Vector,r24Vector);
        T24.TE(3*f24/(r24*r24));
        
        I.TE(g24);
        
        T24.ME(I);
        T24.TE(1/(r24*r24*r24));

        
        I.E(1);
        
//        double[][] T12 = new double[3][3];

        Tensor3D T34 = new Tensor3D();
        
        T34.PEv1v2(r34Vector,r34Vector);
        T34.TE(3*f34/(r34*r34));
        
        I.TE(g34);
        
        T34.ME(I);
        T34.TE(1/(r34*r34*r34));

        
        
        // Now distribute the elements of the tensor into 3 separate "row" vectors
        // so I can do dot products with etomica math methods
        // kmb, 8/7/06
        
        Vector3D T12row1 = new Vector3D();
        Vector3D T12row2 = new Vector3D();
        Vector3D T12row3 = new Vector3D();
        
        T12row1.setX(0,T12.component(0,0));
        T12row1.setX(1,T12.component(0,1));
        T12row1.setX(2,T12.component(0,2));
        T12row2.setX(0,T12.component(1,0));
        T12row2.setX(1,T12.component(1,1));
        T12row2.setX(2,T12.component(1,2));
        T12row3.setX(0,T12.component(2,0));
        T12row3.setX(1,T12.component(2,1));
        T12row3.setX(2,T12.component(2,2));
        

        Vector3D T13row1 = new Vector3D();
        Vector3D T13row2 = new Vector3D();
        Vector3D T13row3 = new Vector3D();
        
        T13row1.setX(0,T13.component(0,0));
        T13row1.setX(1,T13.component(0,1));
        T13row1.setX(2,T13.component(0,2));
        T13row2.setX(0,T13.component(1,0));
        T13row2.setX(1,T13.component(1,1));
        T13row2.setX(2,T13.component(1,2));
        T13row3.setX(0,T13.component(2,0));
        T13row3.setX(1,T13.component(2,1));
        T13row3.setX(2,T13.component(2,2));


        Vector3D T14row1 = new Vector3D();
        Vector3D T14row2 = new Vector3D();
        Vector3D T14row3 = new Vector3D();
        
        T14row1.setX(0,T14.component(0,0));
        T14row1.setX(1,T14.component(0,1));
        T14row1.setX(2,T14.component(0,2));
        T14row2.setX(0,T14.component(1,0));
        T14row2.setX(1,T14.component(1,1));
        T14row2.setX(2,T14.component(1,2));
        T14row3.setX(0,T14.component(2,0));
        T14row3.setX(1,T14.component(2,1));
        T14row3.setX(2,T14.component(2,2));

        
        
        Vector3D T23row1 = new Vector3D();
        Vector3D T23row2 = new Vector3D();
        Vector3D T23row3 = new Vector3D();
        
        T23row1.setX(0,T23.component(0,0));
        T23row1.setX(1,T23.component(0,1));
        T23row1.setX(2,T23.component(0,2));
        T23row2.setX(0,T23.component(1,0));
        T23row2.setX(1,T23.component(1,1));
        T23row2.setX(2,T23.component(1,2));
        T23row3.setX(0,T23.component(2,0));
        T23row3.setX(1,T23.component(2,1));
        T23row3.setX(2,T23.component(2,2));


        Vector3D T24row1 = new Vector3D();
        Vector3D T24row2 = new Vector3D();
        Vector3D T24row3 = new Vector3D();
        
        T24row1.setX(0,T24.component(0,0));
        T24row1.setX(1,T24.component(0,1));
        T24row1.setX(2,T24.component(0,2));
        T24row2.setX(0,T24.component(1,0));
        T24row2.setX(1,T24.component(1,1));
        T24row2.setX(2,T24.component(1,2));
        T24row3.setX(0,T24.component(2,0));
        T24row3.setX(1,T24.component(2,1));
        T24row3.setX(2,T24.component(2,2));

        
        Vector3D T34row1 = new Vector3D();
        Vector3D T34row2 = new Vector3D();
        Vector3D T34row3 = new Vector3D();
        
        T34row1.setX(0,T34.component(0,0));
        T34row1.setX(1,T34.component(0,1));
        T34row1.setX(2,T34.component(0,2));
        T34row2.setX(0,T34.component(1,0));
        T34row2.setX(1,T34.component(1,1));
        T34row2.setX(2,T34.component(1,2));
        T34row3.setX(0,T34.component(2,0));
        T34row3.setX(1,T34.component(2,1));
        T34row3.setX(2,T34.component(2,2));

        
        // Set the induced dipole moments equal to 10% of the permanent dipole value
        // kmb, 8/7/06
        P1.setX(0,14.3952507082236);
        P1.setX(1,14.3952507082236);
        P1.setX(2,14.3952507082236);
        P2.setX(0,14.3952507082236);
        P2.setX(1,14.3952507082236);
        P2.setX(2,14.3952507082236);
        P3.setX(0,14.3952507082236);
        P3.setX(1,14.3952507082236);
        P3.setX(2,14.3952507082236);
        P4.setX(0,14.3952507082236);
        P4.setX(1,14.3952507082236);
        P4.setX(2,14.3952507082236);

        
        
        P1old.E(P1);
        P2old.E(P2);
        P3old.E(P3);
        P4old.E(P4);

        double deltaP1 = 1.0;
        double deltaP2 = 1.0;
        double deltaP3 = 1.0;
        double deltaP4 = 1.0;
        
        while (noSCFforP1 || noSCFforP2 || noSCFforP3 || noSCFforP4) {

    			// First calculate Ep1, based upon guess for P2 and P3 and P4
        	
        		Ep1.setX(0,T12row1.dot(P2)+T13row1.dot(P3)+T14row1.dot(P4));
        		Ep1.setX(1,T12row2.dot(P2)+T13row2.dot(P3)+T14row2.dot(P4));
        		Ep1.setX(2,T12row3.dot(P2)+T13row3.dot(P3)+T14row3.dot(P4));
        		
        		// Now calculate new P1 from the value of Ep1
        		
        		double alphaPol = 1.444;
        		
        		P1.setX(0,alphaPol*(Eq1.x(0) + Ep1.x(0)));
        		P1.setX(1,alphaPol*(Eq1.x(1) + Ep1.x(1)));
        		P1.setX(2,alphaPol*(Eq1.x(2) + Ep1.x(2)));

        		// Next calculate Ep2
        		
        		Ep2.setX(0,T12row1.dot(P1)+T23row1.dot(P3)+T24row1.dot(P4));
        		Ep2.setX(1,T12row2.dot(P1)+T23row2.dot(P3)+T24row2.dot(P4));
        		Ep2.setX(2,T12row3.dot(P1)+T23row3.dot(P3)+T24row3.dot(P4));
        		
        		// Now calculate new P2
        		
        		P2.setX(0,alphaPol*(Eq2.x(0) + Ep2.x(0)));
        		P2.setX(1,alphaPol*(Eq2.x(1) + Ep2.x(1)));
        		P2.setX(2,alphaPol*(Eq2.x(2) + Ep2.x(2)));
        		
        		// Next calculate Ep3
        		
        		Ep3.setX(0,T13row1.dot(P1)+T23row1.dot(P2)+T34row1.dot(P4));
        		Ep3.setX(1,T13row2.dot(P1)+T23row2.dot(P2)+T34row2.dot(P4));
        		Ep3.setX(2,T13row3.dot(P1)+T23row3.dot(P2)+T34row3.dot(P4));
        		
        		// Now calculate new P3
        		
        		P3.setX(0,alphaPol*(Eq3.x(0) + Ep3.x(0)));
        		P3.setX(1,alphaPol*(Eq3.x(1) + Ep3.x(1)));
        		P3.setX(2,alphaPol*(Eq3.x(2) + Ep3.x(2)));

        		
        		// Next calculate Ep4
        		
        		Ep4.setX(0,T14row1.dot(P1)+T24row1.dot(P2)+T34row1.dot(P3));
        		Ep4.setX(1,T14row2.dot(P1)+T24row2.dot(P2)+T34row2.dot(P3));
        		Ep4.setX(2,T14row3.dot(P1)+T24row3.dot(P2)+T34row3.dot(P3));
        		
        		// Now calculate new P4
        		
        		P4.setX(0,alphaPol*(Eq4.x(0) + Ep4.x(0)));
        		P4.setX(1,alphaPol*(Eq4.x(1) + Ep4.x(1)));
        		P4.setX(2,alphaPol*(Eq4.x(2) + Ep4.x(2)));

        		
        		// Evaluate the criteria
        		
	        	deltaP1 = Math.sqrt(P1.Mv1Squared(P1old));
	        	deltaP2 = Math.sqrt(P2.Mv1Squared(P2old));
	        	deltaP3 = Math.sqrt(P3.Mv1Squared(P3old));
	        	deltaP4 = Math.sqrt(P4.Mv1Squared(P4old));
	        	
	        	counterSCFloop = counterSCFloop + 1;
	    
	        	
	        	if (deltaP1<1e-15) {
	        		noSCFforP1 = false;
	        	}
	        	else {
	        		P1old.E(P1);
	        	}
	        	if (deltaP2<1e-15) {
	        		noSCFforP2 = false;
	        	}
	        	else {
	        		P2old.E(P2);
	        	}
	        	if (deltaP3<1e-15) {
	        		noSCFforP3 = false;
	        	}
	        	else {
	        		P3old.E(P3);
	        	}
	        	if (deltaP4<1e-15) {
	        		noSCFforP4 = false;
	        	}
	        	else {
	        		P4old.E(P4);
	        	}

	        	
	        	if (counterSCFloop >= 1000) {
//	        		counterSCFloopOK = false;
	        		noSCFforP1 = false;
	        		noSCFforP2 = false;
	        		noSCFforP3 = false;
	        		noSCFforP4 = false;
	        		loopFailures = loopFailures + 1;
//	        		System.out.println("counterSCFloop = " + counterSCFloop + ", exiting SCF loop due to likely large repulsion");
    		        chargeH11 = Electron.UNIT.toSim(0.6113);
    		        chargeH12 = Electron.UNIT.toSim(0.6113);
    		        chargeM1 = Electron.UNIT.toSim(-1.2226);
    		        chargeH21 = Electron.UNIT.toSim(0.6113);
    		        chargeH22 = Electron.UNIT.toSim(0.6113);
    		        chargeM2 = Electron.UNIT.toSim(-1.2226);
    		        chargeH31 = Electron.UNIT.toSim(0.6113);
    		        chargeH32 = Electron.UNIT.toSim(0.6113);
    		        chargeM3 = Electron.UNIT.toSim(-1.2226);
    		        chargeH41 = Electron.UNIT.toSim(0.6113);
    		        chargeH42 = Electron.UNIT.toSim(0.6113);
    		        chargeM4 = Electron.UNIT.toSim(-1.2226);
	        	}
	     
	        // REPEAT LOOP HERE.  RE-EVALUATE CHARGES ON W1 AND THEN REPEAT.
        }
        
        
        /*
         * Here is where I need to add the polarization term to the energy sum.
         * kmb 5/4/06
         */        
                
                double alpha = 1.444;
                
                double uW1, uW2, uW3, uW4, UpolAtkins;
                        
                UpolAtkins = -0.5*(P1.dot(Eq1)+P2.dot(Eq2)+P3.dot(Eq3)+P4.dot(Eq4));
                       
                
                if (rO1O2 <= 2.1 || rO1O3 <= 2.1 || rO1O4 <= 2.1 || rO2O3 <= 2.1 || rO2O4 <= 2.1 || rO3O4 <= 2.1) UpolAtkins = 0; // value from Cummings
                
                sumSCF += UpolAtkins;
                
                
                uW1 = 1.855 + Debye.UNIT.fromSim(Math.sqrt(P1.squared()));
                uW2 = 1.855 + Debye.UNIT.fromSim(Math.sqrt(P2.squared()));
                uW3 = 1.855 + Debye.UNIT.fromSim(Math.sqrt(P3.squared()));
                uW3 = 1.855 + Debye.UNIT.fromSim(Math.sqrt(P4.squared()));
                
            //System.out.println("Actually going to return an energy not equal to infinity" + sum);
            //if (rOO >= 500) System.out.println("O-O distance = " + rOO + ", sum = " + sum);

        return sumSCF;
        
//        return scfSums;
    }//end of energy
    
    
/*    public double energySCF5(AtomList atomsSCF) { //Atom atom1, Atom atom2, Atom atom3) {//(AtomSet atoms){
        double sumSCF = 0.0;
        double r2 = 0.0;
        double rO1O2 = 0.0;
        double rO1O3 = 0.0;
        double rO1O4 = 0.0;
        double rO1O5 = 0.0;
        double rO2O3 = 0.0;
        double rO2O4 = 0.0;
        double rO2O5 = 0.0;
        double rO3O4 = 0.0;
        double rO3O5 = 0.0;
        double rO4O5 = 0.0;
        
        AtomTreeNodeWaterGCPM node1 = (AtomTreeNodeWaterGCPM)atomsSCF.get(0).node;
        AtomTreeNodeWaterGCPM node2 = (AtomTreeNodeWaterGCPM)atomsSCF.get(1).node;//(AtomTreeNodeWaterPPC)pair.atom1.node;
        AtomTreeNodeWaterGCPM node3 = (AtomTreeNodeWaterGCPM)atomsSCF.get(2).node;
        AtomTreeNodeWaterGCPM node4 = (AtomTreeNodeWaterGCPM)atomsSCF.get(3).node;
        AtomTreeNodeWaterGCPM node5 = (AtomTreeNodeWaterGCPM)atomsSCF.get(4).node;

        
        Vector O1r = node1.O.coord.position();
        Vector O2r = node2.O.coord.position();
        Vector O3r = node3.O.coord.position();
        Vector O4r = node4.O.coord.position();
        Vector O5r = node5.O.coord.position();
        Vector H11r = node1.H1.coord.position();
        Vector H12r = node1.H2.coord.position();
        Vector H21r = node2.H1.coord.position();
        Vector H22r = node2.H2.coord.position();
        Vector H31r = node3.H1.coord.position();
        Vector H32r = node3.H2.coord.position();
        Vector H41r = node4.H1.coord.position();
        Vector H42r = node4.H2.coord.position();
        Vector H51r = node5.H1.coord.position();
        Vector H52r = node5.H2.coord.position();


        Vector M1r = node1.M.coord.position();
        Vector M2r = node2.M.coord.position();
        Vector M3r = node3.M.coord.position();
        Vector M4r = node4.M.coord.position();
        Vector M5r = node5.M.coord.position();
        
        
//        System.out.println("O1r coordinates before: " + O1r.x(0) + ", " + O1r.x(1) + ", " + O1r.x(2));
//        System.out.println("O2r coordinates before: " + O2r.x(0) + ", " + O2r.x(1) + ", " + O2r.x(2));
//        System.out.println("H11r coordinates before: " + H11r.x(0) + ", " + H11r.x(1) + ", " + H11r.x(2));
//        System.out.println("H12r coordinates before: " + H12r.x(0) + ", " + H12r.x(1) + ", " + H12r.x(2));
//        System.out.println("H21r coordinates before: " + H21r.x(0) + ", " + H21r.x(1) + ", " + H21r.x(2));
//        System.out.println("H22r coordinates before: " + H22r.x(0) + ", " + H22r.x(1) + ", " + H22r.x(2));
//        System.out.println("M1r coordinates before: " + M1r.x(0) + ", " + M1r.x(1) + ", " + M1r.x(2));
//        System.out.println("M2r coordinates before: " + M2r.x(0) + ", " + M2r.x(1) + ", " + M2r.x(2));
        
        O2r.PE(0,4); //PE(5);
        H21r.PE(0,4);
        H22r.PE(0,4);
        M2rOriginal.PE(0,4);
        
        O2r.setX(2,2); //PE(5);
        H21r.setX(2,2);
        H22r.setX(2,2);
        M2rOriginal.setX(2,2);
        //O2r.setX(2,2);
        
//		Dimer configuration
        O2r.setX(2,1.7683903);
        O2r.setX(0,2.18378015);
        O2r.setX(1,0.0);
        H21r.setX(2,1.90087324);
        H21r.setX(0,2.73561133);
        H21r.setX(1,-0.7531133);
        H22r.setX(2,1.90087324);
        H22r.setX(0,2.73561133);
        H22r.setX(1,0.7531133);
        M2rOriginal.setX(2,1.79406929);
        M2rOriginal.setX(0,2.29074084);
        M2rOriginal.setX(1,0.0);


        O2r.setX(2,2.2);
        O2r.setX(0,0.0);
        O2r.setX(1,0.0);
        H21r.setX(2,2.767511567);
        H21r.setX(0,-0.75311329);
        H21r.setX(1,0.0);
        H22r.setX(2,2.767511567);
        H22r.setX(0,0.75311329);
        H22r.setX(1,0.0);
        M2rOriginal.setX(2,2.31);
        M2rOriginal.setX(0,0.0);
        M2rOriginal.setX(1,0.0);


        O1r.setX(2,24600.0);
        O1r.setX(0,0.0);
        O1r.setX(1,0.0);
        H11r.setX(2,24600.567511567);
        H11r.setX(0,-0.75311329);
        H11r.setX(1,0.0);
        H12r.setX(2,24600.567511567);
        H12r.setX(0,0.75311329);
        H12r.setX(1,0.0);
        M1rOriginal.setX(2,24600.11);
        M1rOriginal.setX(0,0.0);
        M1rOriginal.setX(1,0.0);

        O2r.setX(2,0.0);
        O2r.setX(0,0.0);
        O2r.setX(1,0.0);
        H21r.setX(2,0.567511567);
        H21r.setX(0,-0.75311329);
        H21r.setX(1,0.0);
        H22r.setX(2,0.567511567);
        H22r.setX(0,0.75311329);
        H22r.setX(1,0.0);
        M2rOriginal.setX(2,0.11);
        M2rOriginal.setX(0,0.0);
        M2rOriginal.setX(1,0.0);

        O3r.setX(2,-24600.0);
        O3r.setX(0,0.0);
        O3r.setX(1,0.0);
        H31r.setX(2,-24600.567511567);
        H31r.setX(0,-0.75311329);
        H31r.setX(1,0.0);
        H32r.setX(2,-24600.567511567);
        H32r.setX(0,0.75311329);
        H32r.setX(1,0.0);
        M3rOriginal.setX(2,-24600.11);
        M3rOriginal.setX(0,0.0);
        M3rOriginal.setX(1,0.0);

		
        
		
        final double core = 4.41; //4.41 = 2.1^2; value according to Cummings
  
        // Here is the beginning of my self-consistent solution algorithm
        // 12/26/05
        
        Vector Eq1old = new Vector3D();
        Vector Eq2old = new Vector3D();
        Vector Eq3old = new Vector3D();
        Vector Eq1 = new Vector3D();
        Vector Eq2 = new Vector3D();
        Vector Eq3 = new Vector3D();
        Vector Eq4 = new Vector3D();
        Vector Eq5 = new Vector3D();
        Vector Ep1old = new Vector3D();
        Vector Ep2old = new Vector3D();
        Vector Ep3old = new Vector3D();
        Vector Ep1 = new Vector3D();
        Vector Ep2 = new Vector3D();
        Vector Ep3 = new Vector3D();
        Vector Ep4 = new Vector3D();
        Vector Ep5 = new Vector3D();
        
        Vector Eq2LabFrame = new Vector3D();
        Vector Eq1LabFrame = new Vector3D();
        Vector Eq3LabFrame = new Vector3D();
        Vector Ep3LabFrame = new Vector3D();
        Vector Ep2LabFrame = new Vector3D();
        Vector Ep1LabFrame = new Vector3D();

        Vector Eon1Total = new Vector3D();
        Vector Eon2Total = new Vector3D();
        Vector Eon3Total = new Vector3D();
        Vector Eon4Total = new Vector3D();
        Vector Eon5Total = new Vector3D();
        
        Eon1Total.Ev1Pv2(Eq1, Ep1);
        Eon2Total.Ev1Pv2(Eq2, Ep2);
        Eon3Total.Ev1Pv2(Eq3, Ep3);
        Eon4Total.Ev1Pv2(Eq4, Ep4);
        Eon5Total.Ev1Pv2(Eq5, Ep5);

        Vector P1 = new Vector3D();
        Vector P2 = new Vector3D();
        Vector P3 = new Vector3D();
        Vector P4 = new Vector3D();
        Vector P5 = new Vector3D();
        Vector P1old = new Vector3D();
        Vector P2old = new Vector3D();
        Vector P3old = new Vector3D();
        Vector P4old = new Vector3D();
        Vector P5old = new Vector3D();
        
        int counterSCFloop = 0; //1;
        boolean counterSCFloopOK = true;
        int loopFailures = 0;
        
        
        // These booleans monitor whether the SCF algorithm has iteratively found
        // a solution given the criterion of Cummings, J.Chem.Phys. 2005
        boolean noSCFforP1 = true;
        boolean noSCFforP2 = true;
        boolean noSCFforP3 = true;
        boolean noSCFforP4 = true;
        boolean noSCFforP5 = true;
		
        // Need loop to check for configuration overlap between charged particles
        
//      compute O-O distance to consider bypassing the SCF loop   
        rO1O2 = Math.sqrt(O1r.Mv1Squared(O2r));
        rO1O3 = Math.sqrt(O1r.Mv1Squared(O3r));
        rO1O4 = Math.sqrt(O1r.Mv1Squared(O4r));
        rO1O5 = Math.sqrt(O1r.Mv1Squared(O5r));
        rO2O3 = Math.sqrt(O3r.Mv1Squared(O2r));
        rO2O4 = Math.sqrt(O4r.Mv1Squared(O2r));
        rO2O5 = Math.sqrt(O5r.Mv1Squared(O2r));
        rO3O4 = Math.sqrt(O4r.Mv1Squared(O3r));
        rO3O5 = Math.sqrt(O5r.Mv1Squared(O3r));
        rO4O5 = Math.sqrt(O4r.Mv1Squared(O5r));
        
        if (rO1O2 > 10 || rO1O3 > 10 || rO1O4 > 10 || rO1O5 > 10 || rO2O3 > 10 || rO2O4 > 10 || rO2O5 > 10 || rO3O4 > 10 || rO3O5 > 10 || rO4O5 > 10) {
			return -123456789.0;
    }

        
        //System.out.println("O-O distance is " + r2);
        
//		System.out.println(rO1O2 + ", " + rO1O3 + ", " + rO1O4 + ", " + rO2O3 + ", " + rO2O4 + ", " + rO3O4 );
        
        if(rO1O2 <= 2.1 || rO1O3 <= 2.1 || rO1O4 <= 2.1 || rO1O5 <= 2.1 || rO2O3 <= 2.1 || rO2O4 <= 2.1 || rO2O5 <= 2.1 || rO3O4 <= 2.1 || rO3O5 <= 2.1 || rO4O5 <= 2.1) { // use to be 2.1 for Cummings cutoff
        		noSCFforP1 = false;
        		noSCFforP2 = false;
        		noSCFforP3 = false;
        		noSCFforP4 = false;
        		noSCFforP5 = false;
        		
        		// need to give some value to these or energy method will return garbage
        		
	        chargeH11 = Electron.UNIT.toSim(0.6113);
	        chargeH12 = Electron.UNIT.toSim(0.6113);
	        chargeM1 = Electron.UNIT.toSim(-1.2226);
	        chargeH21 = Electron.UNIT.toSim(0.6113);
	        chargeH22 = Electron.UNIT.toSim(0.6113);
	        chargeM2 = Electron.UNIT.toSim(-1.2226);
	        chargeH31 = Electron.UNIT.toSim(0.6113);
	        chargeH32 = Electron.UNIT.toSim(0.6113);
	        chargeM3 = Electron.UNIT.toSim(-1.2226);
	        chargeH41 = Electron.UNIT.toSim(0.6113);
	        chargeH42 = Electron.UNIT.toSim(0.6113);
	        chargeM4 = Electron.UNIT.toSim(-1.2226);
	        chargeH51 = Electron.UNIT.toSim(0.6113);
	        chargeH52 = Electron.UNIT.toSim(0.6113);
	        chargeM5 = Electron.UNIT.toSim(-1.2226);
	        
	        return Double.POSITIVE_INFINITY;

        }
                
                
        gamma = 12.75;
//	    double r = Math.sqrt(r2);
        double rO1O2OverSigma = rO1O2/sigma;
        double sigma2OverRO1O2sq = 1/(rO1O2OverSigma*rO1O2OverSigma);
        double rO1O3OverSigma = rO1O3/sigma;
        double sigma2OverRO1O3sq = 1/(rO1O3OverSigma*rO1O3OverSigma);
        double rO1O4OverSigma = rO1O4/sigma;
        double sigma2OverRO1O4sq = 1/(rO1O4OverSigma*rO1O4OverSigma);
        double rO1O5OverSigma = rO1O5/sigma;
        double sigma2OverRO1O5sq = 1/(rO1O5OverSigma*rO1O5OverSigma);
        double rO2O3OverSigma = rO2O3/sigma;
        double sigma2OverRO2O3sq = 1/(rO2O3OverSigma*rO2O3OverSigma);
        double rO2O4OverSigma = rO2O4/sigma;
        double sigma2OverRO2O4sq = 1/(rO2O4OverSigma*rO2O4OverSigma);
        double rO2O5OverSigma = rO2O5/sigma;
        double sigma2OverRO2O5sq = 1/(rO2O5OverSigma*rO2O5OverSigma);
        double rO3O4OverSigma = rO3O4/sigma;
        double sigma2OverRO3O4sq = 1/(rO3O4OverSigma*rO3O4OverSigma);
        double rO3O5OverSigma = rO3O5/sigma;
        double sigma2OverRO3O5sq = 1/(rO3O5OverSigma*rO3O5OverSigma);
        double rO4O5OverSigma = rO4O5/sigma;
        double sigma2OverRO4O5sq = 1/(rO4O5OverSigma*rO4O5OverSigma);
        
        double sixOverGamma = 6/gamma;
   
        sumSCF += epsilon/(1 - sixOverGamma)*(sixOverGamma*Math.exp(gamma*(1 - rO1O2OverSigma)) - sigma2OverRO1O2sq*sigma2OverRO1O2sq*sigma2OverRO1O2sq);
        sumSCF += epsilon/(1 - sixOverGamma)*(sixOverGamma*Math.exp(gamma*(1 - rO1O3OverSigma)) - sigma2OverRO1O3sq*sigma2OverRO1O3sq*sigma2OverRO1O3sq);
        sumSCF += epsilon/(1 - sixOverGamma)*(sixOverGamma*Math.exp(gamma*(1 - rO1O4OverSigma)) - sigma2OverRO1O4sq*sigma2OverRO1O4sq*sigma2OverRO1O4sq);
        sumSCF += epsilon/(1 - sixOverGamma)*(sixOverGamma*Math.exp(gamma*(1 - rO1O5OverSigma)) - sigma2OverRO1O5sq*sigma2OverRO1O5sq*sigma2OverRO1O5sq);
        sumSCF += epsilon/(1 - sixOverGamma)*(sixOverGamma*Math.exp(gamma*(1 - rO2O3OverSigma)) - sigma2OverRO2O3sq*sigma2OverRO2O3sq*sigma2OverRO2O3sq);
        sumSCF += epsilon/(1 - sixOverGamma)*(sixOverGamma*Math.exp(gamma*(1 - rO2O4OverSigma)) - sigma2OverRO2O4sq*sigma2OverRO2O4sq*sigma2OverRO2O4sq);
        sumSCF += epsilon/(1 - sixOverGamma)*(sixOverGamma*Math.exp(gamma*(1 - rO2O5OverSigma)) - sigma2OverRO2O5sq*sigma2OverRO2O5sq*sigma2OverRO2O5sq);
        sumSCF += epsilon/(1 - sixOverGamma)*(sixOverGamma*Math.exp(gamma*(1 - rO3O4OverSigma)) - sigma2OverRO3O4sq*sigma2OverRO3O4sq*sigma2OverRO3O4sq);
        sumSCF += epsilon/(1 - sixOverGamma)*(sixOverGamma*Math.exp(gamma*(1 - rO3O5OverSigma)) - sigma2OverRO3O5sq*sigma2OverRO3O5sq*sigma2OverRO3O5sq);
        sumSCF += epsilon/(1 - sixOverGamma)*(sixOverGamma*Math.exp(gamma*(1 - rO4O5OverSigma)) - sigma2OverRO4O5sq*sigma2OverRO4O5sq*sigma2OverRO4O5sq);
        
        double sumO2Exp6 = sumSCF;
        
        if (Math.abs(Math.sqrt(Eon1.squared())) >= 0.2083 || Math.abs(Math.sqrt(Eon2.squared())) >= 0.2083 ) {
    			System.out.println("About to return Double.Positive_Infinity");
        		System.out.println("Outside core, with bad electric field values: rOO = " + rOO + ", sum = " + sum + ", Eon1(V/A) = " + Math.sqrt(Eon1.squared())*EFconverttoVperA + ", Eon2(V/A) = " + Math.sqrt(Eon2.squared())*EFconverttoVperA);
        		return Double.POSITIVE_INFINITY;
        }


        double sigmaM = 0.610;
        double sigmaH = 0.455;
        double sqrtHMsigmas = Math.sqrt(2*(sigmaH*sigmaH+sigmaM*sigmaM));
        
        // MUST INCLUDE ERF FUNCTION STUFF TO COULOMBIC ENERGY PART!
        // KMB 8/3/06
        
        r2 = H11r.Mv1Squared(H21r);  // COUNTED-2
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        	sumSCF += chargeH11*chargeH21/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaH)));

        r2 = H11r.Mv1Squared(H22r); // COUNTED-2
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH11*chargeH22/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaH)));

        r2 = H12r.Mv1Squared(H21r); // COUNTED-2
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH12*chargeH21/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaH)));

        r2 = H12r.Mv1Squared(H22r); // COUNTED-2
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH12*chargeH22/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaH)));

//        System.out.println("sum of all O-H terms is " + sum);
        
        r2 = M1r.Mv1Squared(H21r); // COUNTED-2
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH21*chargeM1/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/sqrtHMsigmas));
        //System.out.println("sum is " + sum);

        r2 = M1r.Mv1Squared(H22r); // COUNTED-2
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH22*chargeM1/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/sqrtHMsigmas));
        //System.out.println("sum is " + sum);

        r2 = M2r.Mv1Squared(H11r); // COUNTED-2
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH11*chargeM2/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/sqrtHMsigmas));

        ///System.out.println("sum is " + sum);
        r2 = M2r.Mv1Squared(H12r); // COUNTED-2
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH12*chargeM2/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/sqrtHMsigmas));
        //System.out.println("sum is " + sum);

        r2 = M1r.Mv1Squared(M2r); // COUNTED-2
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeM1*chargeM2/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaM)));
        //System.out.println("sum is " + sum);

        
        
        // EXTRA TERMS FOR 3 MOLECULE PERMUTATIONS AND COMBINATIONS FOR ELECTROSTATICS HERE
        


        r2 = H11r.Mv1Squared(H31r);  // COUNTED-2
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        	sumSCF += chargeH11*chargeH31/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaH)));

        r2 = H11r.Mv1Squared(H32r);  // COUNTED-2
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH11*chargeH32/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaH)));



        r2 = H12r.Mv1Squared(H31r);  // COUNTED-2
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH12*chargeH31/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaH)));

        r2 = H12r.Mv1Squared(H32r);  // COUNTED-2
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH12*chargeH32/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaH)));

        r2 = H21r.Mv1Squared(H31r);  // COUNTED-2
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH21*chargeH31/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaH)));

        r2 = H21r.Mv1Squared(H32r);  // COUNTED-2
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH21*chargeH32/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaH)));

        r2 = H22r.Mv1Squared(H31r);  // COUNTED-2
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH22*chargeH31/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaH)));

        r2 = H22r.Mv1Squared(H32r);  // COUNTED-2
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH22*chargeH32/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaH)));

//        System.out.println("sum of all O-H terms is " + sum);
        


        r2 = M1r.Mv1Squared(H31r);  // COUNTED-2
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH31*chargeM1/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/sqrtHMsigmas));
        //System.out.println("sum is " + sum);

        r2 = M1r.Mv1Squared(H32r);  // COUNTED-2
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH32*chargeM1/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/sqrtHMsigmas));
        //System.out.println("sum is " + sum);



        r2 = M2r.Mv1Squared(H31r);  // COUNTED-2
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH31*chargeM2/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/sqrtHMsigmas));

        ///System.out.println("sum is " + sum);
        r2 = M2r.Mv1Squared(H32r);  // COUNTED-2
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH32*chargeM2/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/sqrtHMsigmas));
        //System.out.println("sum is " + sum);

        r2 = M3r.Mv1Squared(H11r);  // COUNTED-2
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH11*chargeM3/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/sqrtHMsigmas));

        ///System.out.println("sum is " + sum);
        r2 = M3r.Mv1Squared(H12r);  // COUNTED-2
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH12*chargeM3/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/sqrtHMsigmas));
        //System.out.println("sum is " + sum);

        r2 = M3r.Mv1Squared(H21r);  // COUNTED-2
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH21*chargeM3/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/sqrtHMsigmas));

        ///System.out.println("sum is " + sum);
        r2 = M3r.Mv1Squared(H22r);  // COUNTED-2
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH22*chargeM3/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/sqrtHMsigmas));
        //System.out.println("sum is " + sum);


        r2 = M1r.Mv1Squared(M3r);  // COUNTED-2
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeM1*chargeM3/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaM)));
        //System.out.println("sum is " + sum);

        r2 = M3r.Mv1Squared(M2r);  // COUNTED-2
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeM3*chargeM2/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaM)));
        //System.out.println("sum is " + sum);

        
        
        
         * NOW PUTTING IN THE EXTRA ELECTROSTATIC TERMS FOR THE FOURTH WATER MOLECULE
         * KMB, 8/14/06
         
        
        
        r2 = H11r.Mv1Squared(H41r);  // OK FOR B4
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH11*chargeH41/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaH)));
        
        r2 = H11r.Mv1Squared(H42r);  // OK FOR B4
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH11*chargeH42/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaH)));
        
        r2 = H12r.Mv1Squared(H41r);  // OK FOR B4
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH12*chargeH41/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaH)));
        
        r2 = H12r.Mv1Squared(H42r);  // OK FOR B4
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH12*chargeH42/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaH)));
        
        r2 = H21r.Mv1Squared(H41r);  // OK FOR B4
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH21*chargeH41/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaH)));
        
        r2 = H21r.Mv1Squared(H42r);  // OK FOR B4
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH21*chargeH42/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaH)));
        
        r2 = H22r.Mv1Squared(H41r);  // OK FOR B4
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH22*chargeH41/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaH)));
        
        r2 = H22r.Mv1Squared(H42r);  // OK FOR B4
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH22*chargeH42/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaH)));
        
        r2 = H31r.Mv1Squared(H41r);  // OK FOR B4
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH31*chargeH41/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaH)));
        
        r2 = H31r.Mv1Squared(H42r);  // OK FOR B4
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH31*chargeH42/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaH)));
        
        r2 = H32r.Mv1Squared(H41r);  // OK FOR B4
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH32*chargeH41/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaH)));
        
        r2 = H32r.Mv1Squared(H42r);  // OK FOR B4
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH32*chargeH42/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaH)));
        
        r2 = M4r.Mv1Squared(H11r);  // OK FOR B4
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH11*chargeM4/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/sqrtHMsigmas));

        r2 = M4r.Mv1Squared(H12r);  // OK FOR B4
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH12*chargeM4/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/sqrtHMsigmas));

        r2 = M1r.Mv1Squared(H41r);  // OK FOR B4
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH41*chargeM1/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/sqrtHMsigmas));

        r2 = M1r.Mv1Squared(H42r);  // OK FOR B4
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH42*chargeM1/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/sqrtHMsigmas));

        r2 = M4r.Mv1Squared(H21r);  // OK FOR B4
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH21*chargeM4/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/sqrtHMsigmas));

        r2 = M4r.Mv1Squared(H22r);  // OK FOR B4
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH22*chargeM4/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/sqrtHMsigmas));

        r2 = M2r.Mv1Squared(H41r);  // OK FOR B4
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH41*chargeM2/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/sqrtHMsigmas));

        r2 = M2r.Mv1Squared(H42r);  // OK FOR B4
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH42*chargeM2/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/sqrtHMsigmas));

        r2 = M4r.Mv1Squared(H31r);  // OK FOR B4
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH31*chargeM4/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/sqrtHMsigmas));

        r2 = M4r.Mv1Squared(H32r);  // OK FOR B4
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH32*chargeM4/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/sqrtHMsigmas));

        r2 = M3r.Mv1Squared(H41r);  // OK FOR B4
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH41*chargeM3/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/sqrtHMsigmas));

        r2 = M3r.Mv1Squared(H42r);  // OK FOR B4
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH42*chargeM3/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/sqrtHMsigmas));

        r2 = M1r.Mv1Squared(M4r);  // OK FOR B4
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeM1*chargeM4/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaM)));
        //System.out.println("sum is " + sum);

        r2 = M4r.Mv1Squared(M2r);  // OK FOR B4
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeM4*chargeM2/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaM)));
        //System.out.println("sum is " + sum);

        r2 = M3r.Mv1Squared(M4r);  // OK FOR B4
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeM3*chargeM4/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaM)));
        //System.out.println("sum is " + sum);

        
        
        
         * NOW PUTTING IN THE EXTRA ELECTROSTATIC TERMS FOR THE FIFTH WATER MOLECULE
         * KMB, 4/30/07
         
        
        
        r2 = H11r.Mv1Squared(H51r);  // OK FOR B5
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH11*chargeH51/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaH)));
        
        r2 = H11r.Mv1Squared(H52r);  // OK FOR B5
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH11*chargeH52/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaH)));
        
        r2 = H12r.Mv1Squared(H51r);  // OK FOR B5
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH12*chargeH51/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaH)));
        
        r2 = H12r.Mv1Squared(H52r);  // OK FOR B5
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH12*chargeH52/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaH)));
        
        r2 = H21r.Mv1Squared(H51r);  // OK FOR B5
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH21*chargeH51/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaH)));
        
        r2 = H21r.Mv1Squared(H52r);  // OK FOR B5
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH21*chargeH52/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaH)));
        
        r2 = H22r.Mv1Squared(H51r);  // OK FOR B5
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH22*chargeH51/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaH)));
        
        r2 = H22r.Mv1Squared(H52r);  // OK FOR B5
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH22*chargeH52/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaH)));
        
        r2 = H31r.Mv1Squared(H51r);  // OK FOR B5
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH31*chargeH51/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaH)));
        
        r2 = H31r.Mv1Squared(H52r);  // OK FOR B5
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH31*chargeH52/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaH)));
        
        r2 = H32r.Mv1Squared(H51r);  // OK FOR B5
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH32*chargeH51/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaH)));
        
        r2 = H32r.Mv1Squared(H52r);  // OK FOR B5
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH32*chargeH52/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaH)));
        
        r2 = M5r.Mv1Squared(H11r);  // OK FOR B5
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH11*chargeM5/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/sqrtHMsigmas));

        r2 = M5r.Mv1Squared(H12r);  // OK FOR B5
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH12*chargeM5/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/sqrtHMsigmas));

        r2 = M1r.Mv1Squared(H51r);  // OK FOR B5
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH51*chargeM1/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/sqrtHMsigmas));

        r2 = M1r.Mv1Squared(H52r);  // OK FOR B5
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH52*chargeM1/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/sqrtHMsigmas));

        r2 = M5r.Mv1Squared(H21r);  // OK FOR B5
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH21*chargeM5/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/sqrtHMsigmas));

        r2 = M5r.Mv1Squared(H22r);  // OK FOR B5
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH22*chargeM5/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/sqrtHMsigmas));

        r2 = M2r.Mv1Squared(H51r);  // OK FOR B5
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH51*chargeM2/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/sqrtHMsigmas));

        r2 = M2r.Mv1Squared(H52r);  // OK FOR B5
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH52*chargeM2/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/sqrtHMsigmas));

        r2 = M5r.Mv1Squared(H31r);  // OK FOR B5
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH31*chargeM5/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/sqrtHMsigmas));

        r2 = M5r.Mv1Squared(H32r);  // OK FOR B5
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH32*chargeM5/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/sqrtHMsigmas));

        r2 = M3r.Mv1Squared(H51r);  // OK FOR B5
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH51*chargeM3/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/sqrtHMsigmas));

        r2 = M3r.Mv1Squared(H52r);  // OK FOR B5
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH52*chargeM3/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/sqrtHMsigmas));

        r2 = M1r.Mv1Squared(M5r);  // OK FOR B5
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeM1*chargeM5/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaM)));
        //System.out.println("sum is " + sum);

        r2 = M5r.Mv1Squared(M2r);  // OK FOR B5
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeM5*chargeM2/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaM)));
        //System.out.println("sum is " + sum);

        r2 = M3r.Mv1Squared(M5r);  // OK FOR B5
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeM3*chargeM5/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaM)));
        //System.out.println("sum is " + sum);

        r2 = H41r.Mv1Squared(H51r);  // OK FOR B5
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH41*chargeH51/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaH)));
        
        r2 = H41r.Mv1Squared(H52r);  // OK FOR B5
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH41*chargeH52/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaH)));
        
        r2 = M5r.Mv1Squared(H41r);  // OK FOR B5
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH41*chargeM5/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/sqrtHMsigmas));

        r2 = H42r.Mv1Squared(H51r);  // OK FOR B5
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH42*chargeH51/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaH)));
        
        r2 = H42r.Mv1Squared(H52r);  // OK FOR B5
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH42*chargeH52/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaH)));

        r2 = M5r.Mv1Squared(H42r);  // OK FOR B5
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH42*chargeM5/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/sqrtHMsigmas));

        r2 = M4r.Mv1Squared(H51r);  // OK FOR B5
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH51*chargeM4/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/sqrtHMsigmas));

        r2 = M4r.Mv1Squared(H52r);  // OK FOR B5
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeH52*chargeM4/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/sqrtHMsigmas));

        r2 = M4r.Mv1Squared(M5r);  // OK FOR B5
        //if(r2<=core) return Double.POSITIVE_INFINITY;
        sumSCF += chargeM4*chargeM5/Math.sqrt(r2)*(1-SpecialFunctions.erfc(Math.sqrt(r2)/(2*sigmaM)));
        //System.out.println("sum is " + sum);

        
         * NOW, HOW DO I WANT TO HANDLE THE ELECTRIC FIELDS?
         * KEEP GOING WITH ADDING TERMS, OR START USING ARRAYS?
         * IS IT EVEN POSSIBLE TO USE ARRAYS?
         * 
         * JUST GOING TO KEEP ADDING TERMS
         * WILL HAVE TO ADD ARRAYS LATER
         * ATOM WILL HAVE TO KNOW ITS POSITION AND CHARGE, FOR EXAMPLE
         * 
         * KMB, 8/15/06
         
		
        
        
         * Finding the Electric fields at the center of mass of each molecule, Eqi
         * kmb, 8/7/06
         

        
        Vector3D comW1 = new Vector3D();
        Vector3D comW2 = new Vector3D();
        Vector3D comW3 = new Vector3D();
        Vector3D comW4 = new Vector3D();
        Vector3D comW5 = new Vector3D();
        
        // How to use COM etomica code correctly? kmb, 8/7/06
        DataSourceCOM com = new DataSourceCOM(space);
        com.actionPerformed(atoms.getAtom(0));
        com.getCOM();
        
        double massH = 1.01;
        double massO = 16.0;
        double totalMass = 18.02;
 
        double comW1Xcomp = 0.0;
        double comW1Ycomp = 0.0;
        double comW1Zcomp = 0.0;
        
        comW1Xcomp = massH*H11r.x(0) + massO*O1r.x(0) + massH*H12r.x(0);
        comW1Ycomp = massH*H11r.x(1) + massO*O1r.x(1) + massH*H12r.x(1);
        comW1Zcomp = massH*H11r.x(2) + massO*O1r.x(2) + massH*H12r.x(2);
        
        comW1.setX(0,comW1Xcomp);
        comW1.setX(1,comW1Ycomp);
        comW1.setX(2,comW1Zcomp);
        
        comW1.Ea1Tv1(1/totalMass,comW1);

        double comW2Xcomp = 0.0;
        double comW2Ycomp = 0.0;
        double comW2Zcomp = 0.0;
        
        comW2Xcomp = massH*H21r.x(0) + massO*O2r.x(0) + massH*H22r.x(0);
        comW2Ycomp = massH*H21r.x(1) + massO*O2r.x(1) + massH*H22r.x(1);
        comW2Zcomp = massH*H21r.x(2) + massO*O2r.x(2) + massH*H22r.x(2);
        
        comW2.setX(0,comW2Xcomp);
        comW2.setX(1,comW2Ycomp);
        comW2.setX(2,comW2Zcomp);
        
        comW2.Ea1Tv1(1/totalMass,comW2);

        double comW3Xcomp = 0.0;
        double comW3Ycomp = 0.0;
        double comW3Zcomp = 0.0;
        
        comW3Xcomp = massH*H31r.x(0) + massO*O3r.x(0) + massH*H32r.x(0);
        comW3Ycomp = massH*H31r.x(1) + massO*O3r.x(1) + massH*H32r.x(1);
        comW3Zcomp = massH*H31r.x(2) + massO*O3r.x(2) + massH*H32r.x(2);
        
        comW3.setX(0,comW3Xcomp);
        comW3.setX(1,comW3Ycomp);
        comW3.setX(2,comW3Zcomp);
        
        comW3.Ea1Tv1(1/totalMass,comW3);


        double comW4Xcomp = 0.0;
        double comW4Ycomp = 0.0;
        double comW4Zcomp = 0.0;
        
        comW4Xcomp = massH*H41r.x(0) + massO*O4r.x(0) + massH*H42r.x(0);
        comW4Ycomp = massH*H41r.x(1) + massO*O4r.x(1) + massH*H42r.x(1);
        comW4Zcomp = massH*H41r.x(2) + massO*O4r.x(2) + massH*H42r.x(2);
        
        comW4.setX(0,comW4Xcomp);
        comW4.setX(1,comW4Ycomp);
        comW4.setX(2,comW4Zcomp);
        
        comW4.Ea1Tv1(1/totalMass,comW4);


        double comW5Xcomp = 0.0;
        double comW5Ycomp = 0.0;
        double comW5Zcomp = 0.0;
        
        comW5Xcomp = massH*H51r.x(0) + massO*O5r.x(0) + massH*H52r.x(0);
        comW5Ycomp = massH*H51r.x(1) + massO*O5r.x(1) + massH*H52r.x(1);
        comW5Zcomp = massH*H51r.x(2) + massO*O5r.x(2) + massH*H52r.x(2);
        
        comW5.setX(0,comW5Xcomp);
        comW5.setX(1,comW5Ycomp);
        comW5.setX(2,comW5Zcomp);
        
        comW5.Ea1Tv1(1/totalMass,comW5);

        
         * DOUBLE SUMMING NOW COMPLETE; KMB 8/9/06
         * 
         * THESE ELECTRIC FIELDS ARE WRONG! KMB 8/8/06
         * I HAVE NOT DONE THE DOUBLE SUM REQUIRED FOR 3 WATER
         * MOLECULES AS PER CUMMINGS EQUATION 4!
         
        
//        double sqrtHMsigmas = Math.sqrt(2*(sigmaH*sigmaH+sigmaM*sigmaM));
        double sqrtPiHMsigmas = Math.sqrt(Math.PI*(sigmaH*sigmaH+sigmaM*sigmaM));
        double sqrtPiMMsigmas = Math.sqrt(Math.PI*(2*sigmaM*sigmaM));
        
        double Eq1XcompW2 = 0.0;
        double Eq1YcompW2 = 0.0;
        double Eq1ZcompW2 = 0.0;
        double Eq1XcompW3 = 0.0;
        double Eq1YcompW3 = 0.0;
        double Eq1ZcompW3 = 0.0;
        double Eq1XcompW4 = 0.0;
        double Eq1YcompW4 = 0.0;
        double Eq1ZcompW4 = 0.0;
        double Eq1XcompW5 = 0.0;
        double Eq1YcompW5 = 0.0;
        double Eq1ZcompW5 = 0.0;


        double comW1toH21 = Math.sqrt(comW1.Mv1Squared(H21r));
        double comW1toH22 = Math.sqrt(comW1.Mv1Squared(H22r));
        double comW1toM2 = Math.sqrt(comW1.Mv1Squared(M2r));

        double comW1toH31 = Math.sqrt(comW1.Mv1Squared(H31r));
        double comW1toH32 = Math.sqrt(comW1.Mv1Squared(H32r));
        double comW1toM3 = Math.sqrt(comW1.Mv1Squared(M3r));

        double comW1toH41 = Math.sqrt(comW1.Mv1Squared(H41r));
        double comW1toH42 = Math.sqrt(comW1.Mv1Squared(H42r));
        double comW1toM4 = Math.sqrt(comW1.Mv1Squared(M4r));

        double comW1toH51 = Math.sqrt(comW1.Mv1Squared(H51r));
        double comW1toH52 = Math.sqrt(comW1.Mv1Squared(H52r));
        double comW1toM5 = Math.sqrt(comW1.Mv1Squared(M5r));

        
        // Contributions to sum from water#2
        Eq1XcompW2 += chargeH21*(comW1.x(0)-H21r.x(0))/(comW1toH21*comW1toH21*comW1toH21)*((1-SpecialFunctions.erfc(comW1toH21/sqrtHMsigmas))-Math.sqrt(2)*comW1toH21/sqrtPiHMsigmas*Math.exp(-comW1toH21*comW1toH21/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq1XcompW2 += chargeH22*(comW1.x(0)-H22r.x(0))/(comW1toH22*comW1toH22*comW1toH22)*((1-SpecialFunctions.erfc(comW1toH22/sqrtHMsigmas))-Math.sqrt(2)*comW1toH22/sqrtPiHMsigmas*Math.exp(-comW1toH22*comW1toH22/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq1XcompW2 += chargeM2*(comW1.x(0)-M2r.x(0))/(comW1toM2*comW1toM2*comW1toM2)*((1-SpecialFunctions.erfc(comW1toM2/(2*sigmaM)))-Math.sqrt(2)*comW1toM2/sqrtPiMMsigmas*Math.exp(-comW1toM2*comW1toM2/(4*sigmaM*sigmaM)));

        Eq1YcompW2 += chargeH21*(comW1.x(1)-H21r.x(1))/(comW1toH21*comW1toH21*comW1toH21)*((1-SpecialFunctions.erfc(comW1toH21/sqrtHMsigmas))-Math.sqrt(2)*comW1toH21/sqrtPiHMsigmas*Math.exp(-comW1toH21*comW1toH21/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq1YcompW2 += chargeH22*(comW1.x(1)-H22r.x(1))/(comW1toH22*comW1toH22*comW1toH22)*((1-SpecialFunctions.erfc(comW1toH22/sqrtHMsigmas))-Math.sqrt(2)*comW1toH22/sqrtPiHMsigmas*Math.exp(-comW1toH22*comW1toH22/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq1YcompW2 += chargeM2*(comW1.x(1)-M2r.x(1))/(comW1toM2*comW1toM2*comW1toM2)*((1-SpecialFunctions.erfc(comW1toM2/(2*sigmaM)))-Math.sqrt(2)*comW1toM2/sqrtPiMMsigmas*Math.exp(-comW1toM2*comW1toM2/(4*sigmaM*sigmaM)));

        Eq1ZcompW2 += chargeH21*(comW1.x(2)-H21r.x(2))/(comW1toH21*comW1toH21*comW1toH21)*((1-SpecialFunctions.erfc(comW1toH21/sqrtHMsigmas))-Math.sqrt(2)*comW1toH21/sqrtPiHMsigmas*Math.exp(-comW1toH21*comW1toH21/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq1ZcompW2 += chargeH22*(comW1.x(2)-H22r.x(2))/(comW1toH22*comW1toH22*comW1toH22)*((1-SpecialFunctions.erfc(comW1toH22/sqrtHMsigmas))-Math.sqrt(2)*comW1toH22/sqrtPiHMsigmas*Math.exp(-comW1toH22*comW1toH22/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq1ZcompW2 += chargeM2*(comW1.x(2)-M2r.x(2))/(comW1toM2*comW1toM2*comW1toM2)*((1-SpecialFunctions.erfc(comW1toM2/(2*sigmaM)))-Math.sqrt(2)*comW1toM2/sqrtPiMMsigmas*Math.exp(-comW1toM2*comW1toM2/(4*sigmaM*sigmaM)));

        
        // Contributions to sum from water#3
        Eq1XcompW3 += chargeH31*(comW1.x(0)-H31r.x(0))/(comW1toH31*comW1toH31*comW1toH31)*((1-SpecialFunctions.erfc(comW1toH31/sqrtHMsigmas))-Math.sqrt(2)*comW1toH31/sqrtPiHMsigmas*Math.exp(-comW1toH31*comW1toH31/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq1XcompW3 += chargeH32*(comW1.x(0)-H32r.x(0))/(comW1toH32*comW1toH32*comW1toH32)*((1-SpecialFunctions.erfc(comW1toH32/sqrtHMsigmas))-Math.sqrt(2)*comW1toH32/sqrtPiHMsigmas*Math.exp(-comW1toH32*comW1toH32/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq1XcompW3 += chargeM3*(comW1.x(0)-M3r.x(0))/(comW1toM3*comW1toM3*comW1toM3)*((1-SpecialFunctions.erfc(comW1toM3/(2*sigmaM)))-Math.sqrt(2)*comW1toM3/sqrtPiMMsigmas*Math.exp(-comW1toM3*comW1toM3/(4*sigmaM*sigmaM)));

        Eq1YcompW3 += chargeH31*(comW1.x(1)-H31r.x(1))/(comW1toH31*comW1toH31*comW1toH31)*((1-SpecialFunctions.erfc(comW1toH31/sqrtHMsigmas))-Math.sqrt(2)*comW1toH31/sqrtPiHMsigmas*Math.exp(-comW1toH31*comW1toH31/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq1YcompW3 += chargeH32*(comW1.x(1)-H32r.x(1))/(comW1toH32*comW1toH32*comW1toH32)*((1-SpecialFunctions.erfc(comW1toH32/sqrtHMsigmas))-Math.sqrt(2)*comW1toH32/sqrtPiHMsigmas*Math.exp(-comW1toH32*comW1toH32/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq1YcompW3 += chargeM3*(comW1.x(1)-M3r.x(1))/(comW1toM3*comW1toM3*comW1toM3)*((1-SpecialFunctions.erfc(comW1toM3/(2*sigmaM)))-Math.sqrt(2)*comW1toM3/sqrtPiMMsigmas*Math.exp(-comW1toM3*comW1toM3/(4*sigmaM*sigmaM)));

        Eq1ZcompW3 += chargeH31*(comW1.x(2)-H31r.x(2))/(comW1toH31*comW1toH31*comW1toH31)*((1-SpecialFunctions.erfc(comW1toH31/sqrtHMsigmas))-Math.sqrt(2)*comW1toH31/sqrtPiHMsigmas*Math.exp(-comW1toH31*comW1toH31/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq1ZcompW3 += chargeH32*(comW1.x(2)-H32r.x(2))/(comW1toH32*comW1toH32*comW1toH32)*((1-SpecialFunctions.erfc(comW1toH32/sqrtHMsigmas))-Math.sqrt(2)*comW1toH32/sqrtPiHMsigmas*Math.exp(-comW1toH32*comW1toH32/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq1ZcompW3 += chargeM3*(comW1.x(2)-M3r.x(2))/(comW1toM3*comW1toM3*comW1toM3)*((1-SpecialFunctions.erfc(comW1toM3/(2*sigmaM)))-Math.sqrt(2)*comW1toM3/sqrtPiMMsigmas*Math.exp(-comW1toM3*comW1toM3/(4*sigmaM*sigmaM)));
        

        // Contributions to sum from water#4
        Eq1XcompW4 += chargeH41*(comW1.x(0)-H41r.x(0))/(comW1toH41*comW1toH41*comW1toH41)*((1-SpecialFunctions.erfc(comW1toH41/sqrtHMsigmas))-Math.sqrt(2)*comW1toH41/sqrtPiHMsigmas*Math.exp(-comW1toH41*comW1toH41/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq1XcompW4 += chargeH42*(comW1.x(0)-H42r.x(0))/(comW1toH42*comW1toH42*comW1toH42)*((1-SpecialFunctions.erfc(comW1toH42/sqrtHMsigmas))-Math.sqrt(2)*comW1toH42/sqrtPiHMsigmas*Math.exp(-comW1toH42*comW1toH42/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq1XcompW4 += chargeM4*(comW1.x(0)-M4r.x(0))/(comW1toM4*comW1toM4*comW1toM4)*((1-SpecialFunctions.erfc(comW1toM4/(2*sigmaM)))-Math.sqrt(2)*comW1toM4/sqrtPiMMsigmas*Math.exp(-comW1toM4*comW1toM4/(4*sigmaM*sigmaM)));

        Eq1YcompW4 += chargeH41*(comW1.x(1)-H41r.x(1))/(comW1toH41*comW1toH41*comW1toH41)*((1-SpecialFunctions.erfc(comW1toH41/sqrtHMsigmas))-Math.sqrt(2)*comW1toH41/sqrtPiHMsigmas*Math.exp(-comW1toH41*comW1toH41/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq1YcompW4 += chargeH42*(comW1.x(1)-H42r.x(1))/(comW1toH42*comW1toH42*comW1toH42)*((1-SpecialFunctions.erfc(comW1toH42/sqrtHMsigmas))-Math.sqrt(2)*comW1toH42/sqrtPiHMsigmas*Math.exp(-comW1toH42*comW1toH42/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq1YcompW4 += chargeM4*(comW1.x(1)-M4r.x(1))/(comW1toM4*comW1toM4*comW1toM4)*((1-SpecialFunctions.erfc(comW1toM4/(2*sigmaM)))-Math.sqrt(2)*comW1toM4/sqrtPiMMsigmas*Math.exp(-comW1toM4*comW1toM4/(4*sigmaM*sigmaM)));

        Eq1ZcompW4 += chargeH41*(comW1.x(2)-H41r.x(2))/(comW1toH41*comW1toH41*comW1toH41)*((1-SpecialFunctions.erfc(comW1toH41/sqrtHMsigmas))-Math.sqrt(2)*comW1toH41/sqrtPiHMsigmas*Math.exp(-comW1toH41*comW1toH41/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq1ZcompW4 += chargeH42*(comW1.x(2)-H42r.x(2))/(comW1toH42*comW1toH42*comW1toH42)*((1-SpecialFunctions.erfc(comW1toH42/sqrtHMsigmas))-Math.sqrt(2)*comW1toH42/sqrtPiHMsigmas*Math.exp(-comW1toH42*comW1toH42/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq1ZcompW4 += chargeM4*(comW1.x(2)-M4r.x(2))/(comW1toM4*comW1toM4*comW1toM4)*((1-SpecialFunctions.erfc(comW1toM4/(2*sigmaM)))-Math.sqrt(2)*comW1toM4/sqrtPiMMsigmas*Math.exp(-comW1toM4*comW1toM4/(4*sigmaM*sigmaM)));


        // Contributions to sum from water#5
        Eq1XcompW5 += chargeH51*(comW1.x(0)-H51r.x(0))/(comW1toH51*comW1toH51*comW1toH51)*((1-SpecialFunctions.erfc(comW1toH51/sqrtHMsigmas))-Math.sqrt(2)*comW1toH51/sqrtPiHMsigmas*Math.exp(-comW1toH51*comW1toH51/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq1XcompW5 += chargeH52*(comW1.x(0)-H52r.x(0))/(comW1toH52*comW1toH52*comW1toH52)*((1-SpecialFunctions.erfc(comW1toH52/sqrtHMsigmas))-Math.sqrt(2)*comW1toH52/sqrtPiHMsigmas*Math.exp(-comW1toH52*comW1toH52/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq1XcompW5 += chargeM5*(comW1.x(0)-M5r.x(0))/(comW1toM5*comW1toM5*comW1toM5)*((1-SpecialFunctions.erfc(comW1toM5/(2*sigmaM)))-Math.sqrt(2)*comW1toM5/sqrtPiMMsigmas*Math.exp(-comW1toM5*comW1toM5/(4*sigmaM*sigmaM)));

        Eq1YcompW5 += chargeH51*(comW1.x(1)-H51r.x(1))/(comW1toH51*comW1toH51*comW1toH51)*((1-SpecialFunctions.erfc(comW1toH51/sqrtHMsigmas))-Math.sqrt(2)*comW1toH51/sqrtPiHMsigmas*Math.exp(-comW1toH51*comW1toH51/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq1YcompW5 += chargeH52*(comW1.x(1)-H52r.x(1))/(comW1toH52*comW1toH52*comW1toH52)*((1-SpecialFunctions.erfc(comW1toH52/sqrtHMsigmas))-Math.sqrt(2)*comW1toH52/sqrtPiHMsigmas*Math.exp(-comW1toH52*comW1toH52/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq1YcompW5 += chargeM5*(comW1.x(1)-M5r.x(1))/(comW1toM5*comW1toM5*comW1toM5)*((1-SpecialFunctions.erfc(comW1toM5/(2*sigmaM)))-Math.sqrt(2)*comW1toM5/sqrtPiMMsigmas*Math.exp(-comW1toM5*comW1toM5/(4*sigmaM*sigmaM)));

        Eq1ZcompW5 += chargeH51*(comW1.x(2)-H51r.x(2))/(comW1toH51*comW1toH51*comW1toH51)*((1-SpecialFunctions.erfc(comW1toH51/sqrtHMsigmas))-Math.sqrt(2)*comW1toH51/sqrtPiHMsigmas*Math.exp(-comW1toH51*comW1toH51/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq1ZcompW5 += chargeH52*(comW1.x(2)-H52r.x(2))/(comW1toH52*comW1toH52*comW1toH52)*((1-SpecialFunctions.erfc(comW1toH52/sqrtHMsigmas))-Math.sqrt(2)*comW1toH52/sqrtPiHMsigmas*Math.exp(-comW1toH52*comW1toH52/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq1ZcompW5 += chargeM5*(comW1.x(2)-M5r.x(2))/(comW1toM5*comW1toM5*comW1toM5)*((1-SpecialFunctions.erfc(comW1toM5/(2*sigmaM)))-Math.sqrt(2)*comW1toM5/sqrtPiMMsigmas*Math.exp(-comW1toM5*comW1toM5/(4*sigmaM*sigmaM)));

        
        Eq1.setX(0,Eq1XcompW2+Eq1XcompW3+Eq1XcompW4+Eq1XcompW5);
        Eq1.setX(1,Eq1YcompW2+Eq1YcompW3+Eq1YcompW4+Eq1YcompW5);
        Eq1.setX(2,Eq1ZcompW2+Eq1ZcompW3+Eq1ZcompW4+Eq1ZcompW5);

                
        double Eq2XcompW1 = 0.0;
        double Eq2YcompW1 = 0.0;
        double Eq2ZcompW1 = 0.0;
        double Eq2XcompW3 = 0.0;
        double Eq2YcompW3 = 0.0;
        double Eq2ZcompW3 = 0.0;
        double Eq2XcompW4 = 0.0;
        double Eq2YcompW4 = 0.0;
        double Eq2ZcompW4 = 0.0;
        double Eq2XcompW5 = 0.0;
        double Eq2YcompW5 = 0.0;
        double Eq2ZcompW5 = 0.0;

        
        double comW2toH11 = Math.sqrt(comW2.Mv1Squared(H11r));
        double comW2toH12 = Math.sqrt(comW2.Mv1Squared(H12r));
        double comW2toM1 = Math.sqrt(comW2.Mv1Squared(M1r));
        double comW2toH31 = Math.sqrt(comW2.Mv1Squared(H31r));
        double comW2toH32 = Math.sqrt(comW2.Mv1Squared(H32r));
        double comW2toM3 = Math.sqrt(comW2.Mv1Squared(M3r));
        double comW2toH41 = Math.sqrt(comW2.Mv1Squared(H41r));
        double comW2toH42 = Math.sqrt(comW2.Mv1Squared(H42r));
        double comW2toM4 = Math.sqrt(comW2.Mv1Squared(M4r));

        double comW2toH51 = Math.sqrt(comW2.Mv1Squared(H51r));
        double comW2toH52 = Math.sqrt(comW2.Mv1Squared(H52r));
        double comW2toM5 = Math.sqrt(comW2.Mv1Squared(M5r));

        
        // Contributions to sum from water molecule#1
        Eq2XcompW1 += chargeH11*(comW2.x(0)-H11r.x(0))/(comW2toH11*comW2toH11*comW2toH11)*((1-SpecialFunctions.erfc(comW2toH11/sqrtHMsigmas))-Math.sqrt(2)*comW2toH11/sqrtPiHMsigmas*Math.exp(-comW2toH11*comW2toH11/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq2XcompW1 += chargeH12*(comW2.x(0)-H12r.x(0))/(comW2toH12*comW2toH12*comW2toH12)*((1-SpecialFunctions.erfc(comW2toH12/sqrtHMsigmas))-Math.sqrt(2)*comW2toH12/sqrtPiHMsigmas*Math.exp(-comW2toH12*comW2toH12/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq2XcompW1 += chargeM1*(comW2.x(0)-M1r.x(0))/(comW2toM1*comW2toM1*comW2toM1)*((1-SpecialFunctions.erfc(comW2toM1/(2*sigmaM)))-Math.sqrt(2)*comW2toM1/sqrtPiMMsigmas*Math.exp(-comW2toM1*comW2toM1/(4*sigmaM*sigmaM)));

        Eq2YcompW1 += chargeH11*(comW2.x(1)-H11r.x(1))/(comW2toH11*comW2toH11*comW2toH11)*((1-SpecialFunctions.erfc(comW2toH11/sqrtHMsigmas))-Math.sqrt(2)*comW2toH11/sqrtPiHMsigmas*Math.exp(-comW2toH11*comW2toH11/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq2YcompW1 += chargeH12*(comW2.x(1)-H12r.x(1))/(comW2toH12*comW2toH12*comW2toH12)*((1-SpecialFunctions.erfc(comW2toH12/sqrtHMsigmas))-Math.sqrt(2)*comW2toH12/sqrtPiHMsigmas*Math.exp(-comW2toH12*comW2toH12/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq2YcompW1 += chargeM1*(comW2.x(1)-M1r.x(1))/(comW2toM1*comW2toM1*comW2toM1)*((1-SpecialFunctions.erfc(comW2toM1/(2*sigmaM)))-Math.sqrt(2)*comW2toM1/sqrtPiMMsigmas*Math.exp(-comW2toM1*comW2toM1/(4*sigmaM*sigmaM)));

        Eq2ZcompW1 += chargeH11*(comW2.x(2)-H11r.x(2))/(comW2toH11*comW2toH11*comW2toH11)*((1-SpecialFunctions.erfc(comW2toH11/sqrtHMsigmas))-Math.sqrt(2)*comW2toH11/sqrtPiHMsigmas*Math.exp(-comW2toH11*comW2toH11/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq2ZcompW1 += chargeH12*(comW2.x(2)-H12r.x(2))/(comW2toH12*comW2toH12*comW2toH12)*((1-SpecialFunctions.erfc(comW2toH12/sqrtHMsigmas))-Math.sqrt(2)*comW2toH12/sqrtPiHMsigmas*Math.exp(-comW2toH12*comW2toH12/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq2ZcompW1 += chargeM1*(comW2.x(2)-M1r.x(2))/(comW2toM1*comW2toM1*comW2toM1)*((1-SpecialFunctions.erfc(comW2toM1/(2*sigmaM)))-Math.sqrt(2)*comW2toM1/sqrtPiMMsigmas*Math.exp(-comW2toM1*comW2toM1/(4*sigmaM*sigmaM)));


        // Contributions to sum from water molecule#3
        Eq2XcompW3 += chargeH31*(comW2.x(0)-H31r.x(0))/(comW2toH31*comW2toH31*comW2toH31)*((1-SpecialFunctions.erfc(comW2toH31/sqrtHMsigmas))-Math.sqrt(2)*comW2toH31/sqrtPiHMsigmas*Math.exp(-comW2toH31*comW2toH31/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq2XcompW3 += chargeH32*(comW2.x(0)-H32r.x(0))/(comW2toH32*comW2toH32*comW2toH32)*((1-SpecialFunctions.erfc(comW2toH32/sqrtHMsigmas))-Math.sqrt(2)*comW2toH32/sqrtPiHMsigmas*Math.exp(-comW2toH32*comW2toH32/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq2XcompW3 += chargeM3*(comW2.x(0)-M3r.x(0))/(comW2toM3*comW2toM3*comW2toM3)*((1-SpecialFunctions.erfc(comW2toM3/(2*sigmaM)))-Math.sqrt(2)*comW2toM3/sqrtPiMMsigmas*Math.exp(-comW2toM3*comW2toM3/(4*sigmaM*sigmaM)));

        Eq2YcompW3 += chargeH31*(comW2.x(1)-H31r.x(1))/(comW2toH31*comW2toH31*comW2toH31)*((1-SpecialFunctions.erfc(comW2toH31/sqrtHMsigmas))-Math.sqrt(2)*comW2toH31/sqrtPiHMsigmas*Math.exp(-comW2toH31*comW2toH31/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq2YcompW3 += chargeH32*(comW2.x(1)-H32r.x(1))/(comW2toH32*comW2toH32*comW2toH32)*((1-SpecialFunctions.erfc(comW2toH32/sqrtHMsigmas))-Math.sqrt(2)*comW2toH32/sqrtPiHMsigmas*Math.exp(-comW2toH32*comW2toH32/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq2YcompW3 += chargeM3*(comW2.x(1)-M3r.x(1))/(comW2toM3*comW2toM3*comW2toM3)*((1-SpecialFunctions.erfc(comW2toM3/(2*sigmaM)))-Math.sqrt(2)*comW2toM3/sqrtPiMMsigmas*Math.exp(-comW2toM3*comW2toM3/(4*sigmaM*sigmaM)));

        Eq2ZcompW3 += chargeH31*(comW2.x(2)-H31r.x(2))/(comW2toH31*comW2toH31*comW2toH31)*((1-SpecialFunctions.erfc(comW2toH31/sqrtHMsigmas))-Math.sqrt(2)*comW2toH31/sqrtPiHMsigmas*Math.exp(-comW2toH31*comW2toH31/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq2ZcompW3 += chargeH32*(comW2.x(2)-H32r.x(2))/(comW2toH32*comW2toH32*comW2toH32)*((1-SpecialFunctions.erfc(comW2toH32/sqrtHMsigmas))-Math.sqrt(2)*comW2toH32/sqrtPiHMsigmas*Math.exp(-comW2toH32*comW2toH32/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq2ZcompW3 += chargeM3*(comW2.x(2)-M3r.x(2))/(comW2toM3*comW2toM3*comW2toM3)*((1-SpecialFunctions.erfc(comW2toM3/(2*sigmaM)))-Math.sqrt(2)*comW2toM3/sqrtPiMMsigmas*Math.exp(-comW2toM3*comW2toM3/(4*sigmaM*sigmaM)));
        

        // Contributions to sum from water molecule#4
        Eq2XcompW4 += chargeH41*(comW2.x(0)-H41r.x(0))/(comW2toH41*comW2toH41*comW2toH41)*((1-SpecialFunctions.erfc(comW2toH41/sqrtHMsigmas))-Math.sqrt(2)*comW2toH41/sqrtPiHMsigmas*Math.exp(-comW2toH41*comW2toH41/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq2XcompW4 += chargeH42*(comW2.x(0)-H42r.x(0))/(comW2toH42*comW2toH42*comW2toH42)*((1-SpecialFunctions.erfc(comW2toH42/sqrtHMsigmas))-Math.sqrt(2)*comW2toH42/sqrtPiHMsigmas*Math.exp(-comW2toH42*comW2toH42/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq2XcompW4 += chargeM4*(comW2.x(0)-M4r.x(0))/(comW2toM4*comW2toM4*comW2toM4)*((1-SpecialFunctions.erfc(comW2toM4/(2*sigmaM)))-Math.sqrt(2)*comW2toM4/sqrtPiMMsigmas*Math.exp(-comW2toM4*comW2toM4/(4*sigmaM*sigmaM)));

        Eq2YcompW4 += chargeH41*(comW2.x(1)-H41r.x(1))/(comW2toH41*comW2toH41*comW2toH41)*((1-SpecialFunctions.erfc(comW2toH41/sqrtHMsigmas))-Math.sqrt(2)*comW2toH41/sqrtPiHMsigmas*Math.exp(-comW2toH41*comW2toH41/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq2YcompW4 += chargeH42*(comW2.x(1)-H42r.x(1))/(comW2toH42*comW2toH42*comW2toH42)*((1-SpecialFunctions.erfc(comW2toH42/sqrtHMsigmas))-Math.sqrt(2)*comW2toH42/sqrtPiHMsigmas*Math.exp(-comW2toH42*comW2toH42/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq2YcompW4 += chargeM4*(comW2.x(1)-M4r.x(1))/(comW2toM4*comW2toM4*comW2toM4)*((1-SpecialFunctions.erfc(comW2toM4/(2*sigmaM)))-Math.sqrt(2)*comW2toM4/sqrtPiMMsigmas*Math.exp(-comW2toM4*comW2toM4/(4*sigmaM*sigmaM)));

        Eq2ZcompW4 += chargeH41*(comW2.x(2)-H41r.x(2))/(comW2toH41*comW2toH41*comW2toH41)*((1-SpecialFunctions.erfc(comW2toH41/sqrtHMsigmas))-Math.sqrt(2)*comW2toH41/sqrtPiHMsigmas*Math.exp(-comW2toH41*comW2toH41/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq2ZcompW4 += chargeH42*(comW2.x(2)-H42r.x(2))/(comW2toH42*comW2toH42*comW2toH42)*((1-SpecialFunctions.erfc(comW2toH42/sqrtHMsigmas))-Math.sqrt(2)*comW2toH42/sqrtPiHMsigmas*Math.exp(-comW2toH42*comW2toH42/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq2ZcompW4 += chargeM4*(comW2.x(2)-M4r.x(2))/(comW2toM4*comW2toM4*comW2toM4)*((1-SpecialFunctions.erfc(comW2toM4/(2*sigmaM)))-Math.sqrt(2)*comW2toM4/sqrtPiMMsigmas*Math.exp(-comW2toM4*comW2toM4/(4*sigmaM*sigmaM)));

        
        // Contributions to sum from water molecule#5
        Eq2XcompW5 += chargeH51*(comW2.x(0)-H51r.x(0))/(comW2toH51*comW2toH51*comW2toH51)*((1-SpecialFunctions.erfc(comW2toH51/sqrtHMsigmas))-Math.sqrt(2)*comW2toH51/sqrtPiHMsigmas*Math.exp(-comW2toH51*comW2toH51/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq2XcompW5 += chargeH52*(comW2.x(0)-H52r.x(0))/(comW2toH52*comW2toH52*comW2toH52)*((1-SpecialFunctions.erfc(comW2toH52/sqrtHMsigmas))-Math.sqrt(2)*comW2toH52/sqrtPiHMsigmas*Math.exp(-comW2toH52*comW2toH52/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq2XcompW5 += chargeM5*(comW2.x(0)-M5r.x(0))/(comW2toM5*comW2toM5*comW2toM5)*((1-SpecialFunctions.erfc(comW2toM5/(2*sigmaM)))-Math.sqrt(2)*comW2toM5/sqrtPiMMsigmas*Math.exp(-comW2toM5*comW2toM5/(4*sigmaM*sigmaM)));

        Eq2YcompW5 += chargeH51*(comW2.x(1)-H51r.x(1))/(comW2toH51*comW2toH51*comW2toH51)*((1-SpecialFunctions.erfc(comW2toH51/sqrtHMsigmas))-Math.sqrt(2)*comW2toH51/sqrtPiHMsigmas*Math.exp(-comW2toH51*comW2toH51/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq2YcompW5 += chargeH52*(comW2.x(1)-H52r.x(1))/(comW2toH52*comW2toH52*comW2toH52)*((1-SpecialFunctions.erfc(comW2toH52/sqrtHMsigmas))-Math.sqrt(2)*comW2toH52/sqrtPiHMsigmas*Math.exp(-comW2toH52*comW2toH52/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq2YcompW5 += chargeM5*(comW2.x(1)-M5r.x(1))/(comW2toM5*comW2toM5*comW2toM5)*((1-SpecialFunctions.erfc(comW2toM5/(2*sigmaM)))-Math.sqrt(2)*comW2toM5/sqrtPiMMsigmas*Math.exp(-comW2toM5*comW2toM5/(4*sigmaM*sigmaM)));

        Eq2ZcompW5 += chargeH51*(comW2.x(2)-H51r.x(2))/(comW2toH51*comW2toH51*comW2toH51)*((1-SpecialFunctions.erfc(comW2toH51/sqrtHMsigmas))-Math.sqrt(2)*comW2toH51/sqrtPiHMsigmas*Math.exp(-comW2toH51*comW2toH51/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq2ZcompW5 += chargeH52*(comW2.x(2)-H52r.x(2))/(comW2toH52*comW2toH52*comW2toH52)*((1-SpecialFunctions.erfc(comW2toH52/sqrtHMsigmas))-Math.sqrt(2)*comW2toH52/sqrtPiHMsigmas*Math.exp(-comW2toH52*comW2toH52/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq2ZcompW5 += chargeM5*(comW2.x(2)-M5r.x(2))/(comW2toM5*comW2toM5*comW2toM5)*((1-SpecialFunctions.erfc(comW2toM5/(2*sigmaM)))-Math.sqrt(2)*comW2toM5/sqrtPiMMsigmas*Math.exp(-comW2toM5*comW2toM5/(4*sigmaM*sigmaM)));

        
        Eq2.setX(0,Eq2XcompW1+Eq2XcompW3+Eq2XcompW4+Eq2XcompW5);
        Eq2.setX(1,Eq2YcompW1+Eq2YcompW3+Eq2YcompW4+Eq2YcompW5);
        Eq2.setX(2,Eq2ZcompW1+Eq2ZcompW3+Eq2ZcompW4+Eq2ZcompW5);


        // Find Eq3
        double Eq3XcompW1	 = 0.0;
        double Eq3YcompW1 = 0.0;
        double Eq3ZcompW1 = 0.0;
        double Eq3XcompW2 = 0.0;
        double Eq3YcompW2 = 0.0;
        double Eq3ZcompW2 = 0.0;
        double Eq3XcompW4 = 0.0;
        double Eq3YcompW4 = 0.0;
        double Eq3ZcompW4 = 0.0;
        double Eq3XcompW5 = 0.0;
        double Eq3YcompW5 = 0.0;
        double Eq3ZcompW5 = 0.0;

        
        
        
        double comW3toH11 = Math.sqrt(comW3.Mv1Squared(H11r));
        double comW3toH12 = Math.sqrt(comW3.Mv1Squared(H12r));
        double comW3toM1 = Math.sqrt(comW3.Mv1Squared(M1r));

        double comW3toH21 = Math.sqrt(comW3.Mv1Squared(H21r));
        double comW3toH22 = Math.sqrt(comW3.Mv1Squared(H22r));
        double comW3toM2 = Math.sqrt(comW3.Mv1Squared(M2r));

        double comW3toH41 = Math.sqrt(comW3.Mv1Squared(H41r));
        double comW3toH42 = Math.sqrt(comW3.Mv1Squared(H42r));
        double comW3toM4 = Math.sqrt(comW3.Mv1Squared(M4r));

        double comW3toH51 = Math.sqrt(comW3.Mv1Squared(H51r));
        double comW3toH52 = Math.sqrt(comW3.Mv1Squared(H52r));
        double comW3toM5 = Math.sqrt(comW3.Mv1Squared(M5r));

        
        
        // Contributions to sum from water molecule#1       
        Eq3XcompW1 += chargeH11*(comW3.x(0)-H11r.x(0))/(comW3toH11*comW3toH11*comW3toH11)*((1-SpecialFunctions.erfc(comW3toH11/sqrtHMsigmas))-Math.sqrt(2)*comW3toH11/sqrtPiHMsigmas*Math.exp(-comW3toH11*comW3toH11/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq3XcompW1 += chargeH12*(comW3.x(0)-H12r.x(0))/(comW3toH12*comW3toH12*comW3toH12)*((1-SpecialFunctions.erfc(comW3toH12/sqrtHMsigmas))-Math.sqrt(2)*comW3toH12/sqrtPiHMsigmas*Math.exp(-comW3toH12*comW3toH12/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq3XcompW1 += chargeM1*(comW3.x(0)-M1r.x(0))/(comW3toM1*comW3toM1*comW3toM1)*((1-SpecialFunctions.erfc(comW3toM1/(2*sigmaM)))-Math.sqrt(2)*comW3toM1/sqrtPiMMsigmas*Math.exp(-comW3toM1*comW3toM1/(4*sigmaM*sigmaM)));

        Eq3YcompW1 += chargeH11*(comW3.x(1)-H11r.x(1))/(comW3toH11*comW3toH11*comW3toH11)*((1-SpecialFunctions.erfc(comW3toH11/sqrtHMsigmas))-Math.sqrt(2)*comW3toH11/sqrtPiHMsigmas*Math.exp(-comW3toH11*comW3toH11/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq3YcompW1 += chargeH12*(comW3.x(1)-H12r.x(1))/(comW3toH12*comW3toH12*comW3toH12)*((1-SpecialFunctions.erfc(comW3toH12/sqrtHMsigmas))-Math.sqrt(2)*comW3toH12/sqrtPiHMsigmas*Math.exp(-comW3toH12*comW3toH12/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq3YcompW1 += chargeM1*(comW3.x(1)-M1r.x(1))/(comW3toM1*comW3toM1*comW3toM1)*((1-SpecialFunctions.erfc(comW3toM1/(2*sigmaM)))-Math.sqrt(2)*comW3toM1/sqrtPiMMsigmas*Math.exp(-comW3toM1*comW3toM1/(4*sigmaM*sigmaM)));

        Eq3ZcompW1 += chargeH11*(comW3.x(2)-H11r.x(2))/(comW3toH11*comW3toH11*comW3toH11)*((1-SpecialFunctions.erfc(comW3toH11/sqrtHMsigmas))-Math.sqrt(2)*comW3toH11/sqrtPiHMsigmas*Math.exp(-comW3toH11*comW3toH11/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq3ZcompW1 += chargeH12*(comW3.x(2)-H12r.x(2))/(comW3toH12*comW3toH12*comW3toH12)*((1-SpecialFunctions.erfc(comW3toH12/sqrtHMsigmas))-Math.sqrt(2)*comW3toH12/sqrtPiHMsigmas*Math.exp(-comW3toH12*comW3toH12/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq3ZcompW1 += chargeM1*(comW3.x(2)-M1r.x(2))/(comW3toM1*comW3toM1*comW3toM1)*((1-SpecialFunctions.erfc(comW3toM1/(2*sigmaM)))-Math.sqrt(2)*comW3toM1/sqrtPiMMsigmas*Math.exp(-comW3toM1*comW3toM1/(4*sigmaM*sigmaM)));

        
        // Contributions to sum from water molecule#2
        Eq3XcompW2 += chargeH21*(comW3.x(0)-H21r.x(0))/(comW3toH21*comW3toH21*comW3toH21)*((1-SpecialFunctions.erfc(comW3toH21/sqrtHMsigmas))-Math.sqrt(2)*comW3toH21/sqrtPiHMsigmas*Math.exp(-comW3toH21*comW3toH21/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq3XcompW2 += chargeH22*(comW3.x(0)-H22r.x(0))/(comW3toH22*comW3toH22*comW3toH22)*((1-SpecialFunctions.erfc(comW3toH22/sqrtHMsigmas))-Math.sqrt(2)*comW3toH22/sqrtPiHMsigmas*Math.exp(-comW3toH22*comW3toH22/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq3XcompW2 += chargeM2*(comW3.x(0)-M2r.x(0))/(comW3toM2*comW3toM2*comW3toM2)*((1-SpecialFunctions.erfc(comW3toM2/(2*sigmaM)))-Math.sqrt(2)*comW3toM2/sqrtPiMMsigmas*Math.exp(-comW3toM2*comW3toM2/(4*sigmaM*sigmaM)));

        Eq3YcompW2 += chargeH21*(comW3.x(1)-H21r.x(1))/(comW3toH21*comW3toH21*comW3toH21)*((1-SpecialFunctions.erfc(comW3toH21/sqrtHMsigmas))-Math.sqrt(2)*comW3toH21/sqrtPiHMsigmas*Math.exp(-comW3toH21*comW3toH21/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq3YcompW2 += chargeH22*(comW3.x(1)-H22r.x(1))/(comW3toH22*comW3toH22*comW3toH22)*((1-SpecialFunctions.erfc(comW3toH22/sqrtHMsigmas))-Math.sqrt(2)*comW3toH22/sqrtPiHMsigmas*Math.exp(-comW3toH22*comW3toH22/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq3YcompW2 += chargeM2*(comW3.x(1)-M2r.x(1))/(comW3toM2*comW3toM2*comW3toM2)*((1-SpecialFunctions.erfc(comW3toM2/(2*sigmaM)))-Math.sqrt(2)*comW3toM2/sqrtPiMMsigmas*Math.exp(-comW3toM2*comW3toM2/(4*sigmaM*sigmaM)));

        Eq3ZcompW2 += chargeH21*(comW3.x(2)-H21r.x(2))/(comW3toH21*comW3toH21*comW3toH21)*((1-SpecialFunctions.erfc(comW3toH21/sqrtHMsigmas))-Math.sqrt(2)*comW3toH21/sqrtPiHMsigmas*Math.exp(-comW3toH21*comW3toH21/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq3ZcompW2 += chargeH22*(comW3.x(2)-H22r.x(2))/(comW3toH22*comW3toH22*comW3toH22)*((1-SpecialFunctions.erfc(comW3toH22/sqrtHMsigmas))-Math.sqrt(2)*comW3toH22/sqrtPiHMsigmas*Math.exp(-comW3toH22*comW3toH22/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq3ZcompW2 += chargeM2*(comW3.x(2)-M2r.x(2))/(comW3toM2*comW3toM2*comW3toM2)*((1-SpecialFunctions.erfc(comW3toM2/(2*sigmaM)))-Math.sqrt(2)*comW3toM2/sqrtPiMMsigmas*Math.exp(-comW3toM2*comW3toM2/(4*sigmaM*sigmaM)));
        

        // Contributions to sum from water molecule#4
        Eq3XcompW4 += chargeH41*(comW3.x(0)-H41r.x(0))/(comW3toH41*comW3toH41*comW3toH41)*((1-SpecialFunctions.erfc(comW3toH41/sqrtHMsigmas))-Math.sqrt(2)*comW3toH41/sqrtPiHMsigmas*Math.exp(-comW3toH41*comW3toH41/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq3XcompW4 += chargeH42*(comW3.x(0)-H42r.x(0))/(comW3toH42*comW3toH42*comW3toH42)*((1-SpecialFunctions.erfc(comW3toH42/sqrtHMsigmas))-Math.sqrt(2)*comW3toH42/sqrtPiHMsigmas*Math.exp(-comW3toH42*comW3toH42/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq3XcompW4 += chargeM4*(comW3.x(0)-M4r.x(0))/(comW3toM4*comW3toM4*comW3toM4)*((1-SpecialFunctions.erfc(comW3toM4/(2*sigmaM)))-Math.sqrt(2)*comW3toM4/sqrtPiMMsigmas*Math.exp(-comW3toM4*comW3toM4/(4*sigmaM*sigmaM)));

        Eq3YcompW4 += chargeH41*(comW3.x(1)-H41r.x(1))/(comW3toH41*comW3toH41*comW3toH41)*((1-SpecialFunctions.erfc(comW3toH41/sqrtHMsigmas))-Math.sqrt(2)*comW3toH41/sqrtPiHMsigmas*Math.exp(-comW3toH41*comW3toH41/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq3YcompW4 += chargeH42*(comW3.x(1)-H42r.x(1))/(comW3toH42*comW3toH42*comW3toH42)*((1-SpecialFunctions.erfc(comW3toH42/sqrtHMsigmas))-Math.sqrt(2)*comW3toH42/sqrtPiHMsigmas*Math.exp(-comW3toH42*comW3toH42/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq3YcompW4 += chargeM4*(comW3.x(1)-M4r.x(1))/(comW3toM4*comW3toM4*comW3toM4)*((1-SpecialFunctions.erfc(comW3toM4/(2*sigmaM)))-Math.sqrt(2)*comW3toM4/sqrtPiMMsigmas*Math.exp(-comW3toM4*comW3toM4/(4*sigmaM*sigmaM)));

        Eq3ZcompW4 += chargeH41*(comW3.x(2)-H41r.x(2))/(comW3toH41*comW3toH41*comW3toH41)*((1-SpecialFunctions.erfc(comW3toH41/sqrtHMsigmas))-Math.sqrt(2)*comW3toH41/sqrtPiHMsigmas*Math.exp(-comW3toH41*comW3toH41/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq3ZcompW4 += chargeH42*(comW3.x(2)-H42r.x(2))/(comW3toH42*comW3toH42*comW3toH42)*((1-SpecialFunctions.erfc(comW3toH42/sqrtHMsigmas))-Math.sqrt(2)*comW3toH42/sqrtPiHMsigmas*Math.exp(-comW3toH42*comW3toH42/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq3ZcompW4 += chargeM4*(comW3.x(2)-M4r.x(2))/(comW3toM4*comW3toM4*comW3toM4)*((1-SpecialFunctions.erfc(comW3toM4/(2*sigmaM)))-Math.sqrt(2)*comW3toM4/sqrtPiMMsigmas*Math.exp(-comW3toM4*comW3toM4/(4*sigmaM*sigmaM)));


        // Contributions to sum from water molecule#5
        Eq3XcompW5 += chargeH51*(comW3.x(0)-H51r.x(0))/(comW3toH51*comW3toH51*comW3toH51)*((1-SpecialFunctions.erfc(comW3toH51/sqrtHMsigmas))-Math.sqrt(2)*comW3toH51/sqrtPiHMsigmas*Math.exp(-comW3toH51*comW3toH51/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq3XcompW5 += chargeH52*(comW3.x(0)-H52r.x(0))/(comW3toH52*comW3toH52*comW3toH52)*((1-SpecialFunctions.erfc(comW3toH52/sqrtHMsigmas))-Math.sqrt(2)*comW3toH52/sqrtPiHMsigmas*Math.exp(-comW3toH52*comW3toH52/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq3XcompW5 += chargeM5*(comW3.x(0)-M5r.x(0))/(comW3toM5*comW3toM5*comW3toM5)*((1-SpecialFunctions.erfc(comW3toM5/(2*sigmaM)))-Math.sqrt(2)*comW3toM5/sqrtPiMMsigmas*Math.exp(-comW3toM5*comW3toM5/(4*sigmaM*sigmaM)));

        Eq3YcompW5 += chargeH51*(comW3.x(1)-H51r.x(1))/(comW3toH51*comW3toH51*comW3toH51)*((1-SpecialFunctions.erfc(comW3toH51/sqrtHMsigmas))-Math.sqrt(2)*comW3toH51/sqrtPiHMsigmas*Math.exp(-comW3toH51*comW3toH51/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq3YcompW5 += chargeH52*(comW3.x(1)-H52r.x(1))/(comW3toH52*comW3toH52*comW3toH52)*((1-SpecialFunctions.erfc(comW3toH52/sqrtHMsigmas))-Math.sqrt(2)*comW3toH52/sqrtPiHMsigmas*Math.exp(-comW3toH52*comW3toH52/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq3YcompW5 += chargeM5*(comW3.x(1)-M5r.x(1))/(comW3toM5*comW3toM5*comW3toM5)*((1-SpecialFunctions.erfc(comW3toM5/(2*sigmaM)))-Math.sqrt(2)*comW3toM5/sqrtPiMMsigmas*Math.exp(-comW3toM5*comW3toM5/(4*sigmaM*sigmaM)));

        Eq3ZcompW5 += chargeH51*(comW3.x(2)-H51r.x(2))/(comW3toH51*comW3toH51*comW3toH51)*((1-SpecialFunctions.erfc(comW3toH51/sqrtHMsigmas))-Math.sqrt(2)*comW3toH51/sqrtPiHMsigmas*Math.exp(-comW3toH51*comW3toH51/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq3ZcompW5 += chargeH52*(comW3.x(2)-H52r.x(2))/(comW3toH52*comW3toH52*comW3toH52)*((1-SpecialFunctions.erfc(comW3toH52/sqrtHMsigmas))-Math.sqrt(2)*comW3toH52/sqrtPiHMsigmas*Math.exp(-comW3toH52*comW3toH52/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq3ZcompW5 += chargeM5*(comW3.x(2)-M5r.x(2))/(comW3toM5*comW3toM5*comW3toM5)*((1-SpecialFunctions.erfc(comW3toM5/(2*sigmaM)))-Math.sqrt(2)*comW3toM5/sqrtPiMMsigmas*Math.exp(-comW3toM5*comW3toM5/(4*sigmaM*sigmaM)));

        
        Eq3.setX(0,Eq3XcompW1+Eq3XcompW2+Eq3XcompW4+Eq3XcompW5);
        Eq3.setX(1,Eq3YcompW1+Eq3YcompW2+Eq3YcompW4+Eq3YcompW5);
        Eq3.setX(2,Eq3ZcompW1+Eq3ZcompW2+Eq3ZcompW4+Eq3ZcompW5);

        
        
        
        
        // Find Eq4
        double Eq4XcompW1	 = 0.0;
        double Eq4YcompW1 = 0.0;
        double Eq4ZcompW1 = 0.0;
        double Eq4XcompW2 = 0.0;
        double Eq4YcompW2 = 0.0;
        double Eq4ZcompW2 = 0.0;
        double Eq4XcompW3 = 0.0;
        double Eq4YcompW3 = 0.0;
        double Eq4ZcompW3 = 0.0;
        double Eq4XcompW5 = 0.0;
        double Eq4YcompW5 = 0.0;
        double Eq4ZcompW5 = 0.0;

        
        
        
        double comW4toH11 = Math.sqrt(comW4.Mv1Squared(H11r));
        double comW4toH12 = Math.sqrt(comW4.Mv1Squared(H12r));
        double comW4toM1 = Math.sqrt(comW4.Mv1Squared(M1r));

        double comW4toH21 = Math.sqrt(comW4.Mv1Squared(H21r));
        double comW4toH22 = Math.sqrt(comW4.Mv1Squared(H22r));
        double comW4toM2 = Math.sqrt(comW4.Mv1Squared(M2r));

        double comW4toH31 = Math.sqrt(comW4.Mv1Squared(H31r));
        double comW4toH32 = Math.sqrt(comW4.Mv1Squared(H32r));
        double comW4toM3 = Math.sqrt(comW4.Mv1Squared(M3r));

        double comW4toH51 = Math.sqrt(comW4.Mv1Squared(H51r));
        double comW4toH52 = Math.sqrt(comW4.Mv1Squared(H52r));
        double comW4toM5 = Math.sqrt(comW4.Mv1Squared(M5r));

        
        // Contributions to sum from water molecule#1       
        Eq4XcompW1 += chargeH11*(comW4.x(0)-H11r.x(0))/(comW4toH11*comW4toH11*comW4toH11)*((1-SpecialFunctions.erfc(comW4toH11/sqrtHMsigmas))-Math.sqrt(2)*comW4toH11/sqrtPiHMsigmas*Math.exp(-comW4toH11*comW4toH11/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq4XcompW1 += chargeH12*(comW4.x(0)-H12r.x(0))/(comW4toH12*comW4toH12*comW4toH12)*((1-SpecialFunctions.erfc(comW4toH12/sqrtHMsigmas))-Math.sqrt(2)*comW4toH12/sqrtPiHMsigmas*Math.exp(-comW4toH12*comW4toH12/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq4XcompW1 += chargeM1*(comW4.x(0)-M1r.x(0))/(comW4toM1*comW4toM1*comW4toM1)*((1-SpecialFunctions.erfc(comW4toM1/(2*sigmaM)))-Math.sqrt(2)*comW4toM1/sqrtPiMMsigmas*Math.exp(-comW4toM1*comW4toM1/(4*sigmaM*sigmaM)));

        Eq4YcompW1 += chargeH11*(comW4.x(1)-H11r.x(1))/(comW4toH11*comW4toH11*comW4toH11)*((1-SpecialFunctions.erfc(comW4toH11/sqrtHMsigmas))-Math.sqrt(2)*comW4toH11/sqrtPiHMsigmas*Math.exp(-comW4toH11*comW4toH11/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq4YcompW1 += chargeH12*(comW4.x(1)-H12r.x(1))/(comW4toH12*comW4toH12*comW4toH12)*((1-SpecialFunctions.erfc(comW4toH12/sqrtHMsigmas))-Math.sqrt(2)*comW4toH12/sqrtPiHMsigmas*Math.exp(-comW4toH12*comW4toH12/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq4YcompW1 += chargeM1*(comW4.x(1)-M1r.x(1))/(comW4toM1*comW4toM1*comW4toM1)*((1-SpecialFunctions.erfc(comW4toM1/(2*sigmaM)))-Math.sqrt(2)*comW4toM1/sqrtPiMMsigmas*Math.exp(-comW4toM1*comW4toM1/(4*sigmaM*sigmaM)));

        Eq4ZcompW1 += chargeH11*(comW4.x(2)-H11r.x(2))/(comW4toH11*comW4toH11*comW4toH11)*((1-SpecialFunctions.erfc(comW4toH11/sqrtHMsigmas))-Math.sqrt(2)*comW4toH11/sqrtPiHMsigmas*Math.exp(-comW4toH11*comW4toH11/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq4ZcompW1 += chargeH12*(comW4.x(2)-H12r.x(2))/(comW4toH12*comW4toH12*comW4toH12)*((1-SpecialFunctions.erfc(comW4toH12/sqrtHMsigmas))-Math.sqrt(2)*comW4toH12/sqrtPiHMsigmas*Math.exp(-comW4toH12*comW4toH12/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq4ZcompW1 += chargeM1*(comW4.x(2)-M1r.x(2))/(comW4toM1*comW4toM1*comW4toM1)*((1-SpecialFunctions.erfc(comW4toM1/(2*sigmaM)))-Math.sqrt(2)*comW4toM1/sqrtPiMMsigmas*Math.exp(-comW4toM1*comW4toM1/(4*sigmaM*sigmaM)));

        
        // Contributions to sum from water molecule#2
        Eq4XcompW2 += chargeH21*(comW4.x(0)-H21r.x(0))/(comW4toH21*comW4toH21*comW4toH21)*((1-SpecialFunctions.erfc(comW4toH21/sqrtHMsigmas))-Math.sqrt(2)*comW4toH21/sqrtPiHMsigmas*Math.exp(-comW4toH21*comW4toH21/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq4XcompW2 += chargeH22*(comW4.x(0)-H22r.x(0))/(comW4toH22*comW4toH22*comW4toH22)*((1-SpecialFunctions.erfc(comW4toH22/sqrtHMsigmas))-Math.sqrt(2)*comW4toH22/sqrtPiHMsigmas*Math.exp(-comW4toH22*comW4toH22/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq4XcompW2 += chargeM2*(comW4.x(0)-M2r.x(0))/(comW4toM2*comW4toM2*comW4toM2)*((1-SpecialFunctions.erfc(comW4toM2/(2*sigmaM)))-Math.sqrt(2)*comW4toM2/sqrtPiMMsigmas*Math.exp(-comW4toM2*comW4toM2/(4*sigmaM*sigmaM)));

        Eq4YcompW2 += chargeH21*(comW4.x(1)-H21r.x(1))/(comW4toH21*comW4toH21*comW4toH21)*((1-SpecialFunctions.erfc(comW4toH21/sqrtHMsigmas))-Math.sqrt(2)*comW4toH21/sqrtPiHMsigmas*Math.exp(-comW4toH21*comW4toH21/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq4YcompW2 += chargeH22*(comW4.x(1)-H22r.x(1))/(comW4toH22*comW4toH22*comW4toH22)*((1-SpecialFunctions.erfc(comW4toH22/sqrtHMsigmas))-Math.sqrt(2)*comW4toH22/sqrtPiHMsigmas*Math.exp(-comW4toH22*comW4toH22/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq4YcompW2 += chargeM2*(comW4.x(1)-M2r.x(1))/(comW4toM2*comW4toM2*comW4toM2)*((1-SpecialFunctions.erfc(comW4toM2/(2*sigmaM)))-Math.sqrt(2)*comW4toM2/sqrtPiMMsigmas*Math.exp(-comW4toM2*comW4toM2/(4*sigmaM*sigmaM)));

        Eq4ZcompW2 += chargeH21*(comW4.x(2)-H21r.x(2))/(comW4toH21*comW4toH21*comW4toH21)*((1-SpecialFunctions.erfc(comW4toH21/sqrtHMsigmas))-Math.sqrt(2)*comW4toH21/sqrtPiHMsigmas*Math.exp(-comW4toH21*comW4toH21/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq4ZcompW2 += chargeH22*(comW4.x(2)-H22r.x(2))/(comW4toH22*comW4toH22*comW4toH22)*((1-SpecialFunctions.erfc(comW4toH22/sqrtHMsigmas))-Math.sqrt(2)*comW4toH22/sqrtPiHMsigmas*Math.exp(-comW4toH22*comW4toH22/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq4ZcompW2 += chargeM2*(comW4.x(2)-M2r.x(2))/(comW4toM2*comW4toM2*comW4toM2)*((1-SpecialFunctions.erfc(comW4toM2/(2*sigmaM)))-Math.sqrt(2)*comW4toM2/sqrtPiMMsigmas*Math.exp(-comW4toM2*comW4toM2/(4*sigmaM*sigmaM)));
        

        // Contributions to sum from water molecule#3
        Eq4XcompW3 += chargeH31*(comW4.x(0)-H31r.x(0))/(comW4toH31*comW4toH31*comW4toH31)*((1-SpecialFunctions.erfc(comW4toH31/sqrtHMsigmas))-Math.sqrt(2)*comW4toH31/sqrtPiHMsigmas*Math.exp(-comW4toH31*comW4toH31/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq4XcompW3 += chargeH32*(comW4.x(0)-H32r.x(0))/(comW4toH32*comW4toH32*comW4toH32)*((1-SpecialFunctions.erfc(comW4toH32/sqrtHMsigmas))-Math.sqrt(2)*comW4toH32/sqrtPiHMsigmas*Math.exp(-comW4toH32*comW4toH32/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq4XcompW3 += chargeM3*(comW4.x(0)-M3r.x(0))/(comW4toM3*comW4toM3*comW4toM3)*((1-SpecialFunctions.erfc(comW4toM3/(2*sigmaM)))-Math.sqrt(2)*comW4toM3/sqrtPiMMsigmas*Math.exp(-comW4toM3*comW4toM3/(4*sigmaM*sigmaM)));

        Eq4YcompW3 += chargeH31*(comW4.x(1)-H31r.x(1))/(comW4toH31*comW4toH31*comW4toH31)*((1-SpecialFunctions.erfc(comW4toH31/sqrtHMsigmas))-Math.sqrt(2)*comW4toH31/sqrtPiHMsigmas*Math.exp(-comW4toH31*comW4toH31/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq4YcompW3 += chargeH32*(comW4.x(1)-H32r.x(1))/(comW4toH32*comW4toH32*comW4toH32)*((1-SpecialFunctions.erfc(comW4toH32/sqrtHMsigmas))-Math.sqrt(2)*comW4toH32/sqrtPiHMsigmas*Math.exp(-comW4toH32*comW4toH32/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq4YcompW3 += chargeM3*(comW4.x(1)-M3r.x(1))/(comW4toM3*comW4toM3*comW4toM3)*((1-SpecialFunctions.erfc(comW4toM3/(2*sigmaM)))-Math.sqrt(2)*comW4toM3/sqrtPiMMsigmas*Math.exp(-comW4toM3*comW4toM3/(4*sigmaM*sigmaM)));

        Eq4ZcompW3 += chargeH31*(comW4.x(2)-H31r.x(2))/(comW4toH31*comW4toH31*comW4toH31)*((1-SpecialFunctions.erfc(comW4toH31/sqrtHMsigmas))-Math.sqrt(2)*comW4toH31/sqrtPiHMsigmas*Math.exp(-comW4toH31*comW4toH31/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq4ZcompW3 += chargeH32*(comW4.x(2)-H32r.x(2))/(comW4toH32*comW4toH32*comW4toH32)*((1-SpecialFunctions.erfc(comW4toH32/sqrtHMsigmas))-Math.sqrt(2)*comW4toH32/sqrtPiHMsigmas*Math.exp(-comW4toH32*comW4toH32/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq4ZcompW3 += chargeM3*(comW4.x(2)-M3r.x(2))/(comW4toM3*comW4toM3*comW4toM3)*((1-SpecialFunctions.erfc(comW4toM3/(2*sigmaM)))-Math.sqrt(2)*comW4toM3/sqrtPiMMsigmas*Math.exp(-comW4toM3*comW4toM3/(4*sigmaM*sigmaM)));


        // Contributions to sum from water molecule#5
        Eq4XcompW5 += chargeH51*(comW4.x(0)-H51r.x(0))/(comW4toH51*comW4toH51*comW4toH51)*((1-SpecialFunctions.erfc(comW4toH51/sqrtHMsigmas))-Math.sqrt(2)*comW4toH51/sqrtPiHMsigmas*Math.exp(-comW4toH51*comW4toH51/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq4XcompW5 += chargeH52*(comW4.x(0)-H52r.x(0))/(comW4toH52*comW4toH52*comW4toH52)*((1-SpecialFunctions.erfc(comW4toH52/sqrtHMsigmas))-Math.sqrt(2)*comW4toH52/sqrtPiHMsigmas*Math.exp(-comW4toH52*comW4toH52/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq4XcompW5 += chargeM5*(comW4.x(0)-M5r.x(0))/(comW4toM5*comW4toM5*comW4toM5)*((1-SpecialFunctions.erfc(comW4toM5/(2*sigmaM)))-Math.sqrt(2)*comW4toM5/sqrtPiMMsigmas*Math.exp(-comW4toM5*comW4toM5/(4*sigmaM*sigmaM)));

        Eq4YcompW5 += chargeH51*(comW4.x(1)-H51r.x(1))/(comW4toH51*comW4toH51*comW4toH51)*((1-SpecialFunctions.erfc(comW4toH51/sqrtHMsigmas))-Math.sqrt(2)*comW4toH51/sqrtPiHMsigmas*Math.exp(-comW4toH51*comW4toH51/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq4YcompW5 += chargeH52*(comW4.x(1)-H52r.x(1))/(comW4toH52*comW4toH52*comW4toH52)*((1-SpecialFunctions.erfc(comW4toH52/sqrtHMsigmas))-Math.sqrt(2)*comW4toH52/sqrtPiHMsigmas*Math.exp(-comW4toH52*comW4toH52/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq4YcompW5 += chargeM5*(comW4.x(1)-M5r.x(1))/(comW4toM5*comW4toM5*comW4toM5)*((1-SpecialFunctions.erfc(comW4toM5/(2*sigmaM)))-Math.sqrt(2)*comW4toM5/sqrtPiMMsigmas*Math.exp(-comW4toM5*comW4toM5/(4*sigmaM*sigmaM)));

        Eq4ZcompW5 += chargeH51*(comW4.x(2)-H51r.x(2))/(comW4toH51*comW4toH51*comW4toH51)*((1-SpecialFunctions.erfc(comW4toH51/sqrtHMsigmas))-Math.sqrt(2)*comW4toH51/sqrtPiHMsigmas*Math.exp(-comW4toH51*comW4toH51/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq4ZcompW5 += chargeH52*(comW4.x(2)-H52r.x(2))/(comW4toH52*comW4toH52*comW4toH52)*((1-SpecialFunctions.erfc(comW4toH52/sqrtHMsigmas))-Math.sqrt(2)*comW4toH52/sqrtPiHMsigmas*Math.exp(-comW4toH52*comW4toH52/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq4ZcompW5 += chargeM5*(comW4.x(2)-M5r.x(2))/(comW4toM5*comW4toM5*comW4toM5)*((1-SpecialFunctions.erfc(comW4toM5/(2*sigmaM)))-Math.sqrt(2)*comW4toM5/sqrtPiMMsigmas*Math.exp(-comW4toM5*comW4toM5/(4*sigmaM*sigmaM)));

        
        Eq4.setX(0,Eq4XcompW1+Eq4XcompW2+Eq4XcompW3+Eq4XcompW5);
        Eq4.setX(1,Eq4YcompW1+Eq4YcompW2+Eq4YcompW3+Eq4YcompW5);
        Eq4.setX(2,Eq4ZcompW1+Eq4ZcompW2+Eq4ZcompW3+Eq4ZcompW5);

        
        
        
        // Find Eq5
        double Eq5XcompW1	 = 0.0;
        double Eq5YcompW1 = 0.0;
        double Eq5ZcompW1 = 0.0;
        double Eq5XcompW2 = 0.0;
        double Eq5YcompW2 = 0.0;
        double Eq5ZcompW2 = 0.0;
        double Eq5XcompW3 = 0.0;
        double Eq5YcompW3 = 0.0;
        double Eq5ZcompW3 = 0.0;
        double Eq5XcompW4 = 0.0;
        double Eq5YcompW4 = 0.0;
        double Eq5ZcompW4 = 0.0;

        
        
        
        double comW5toH11 = Math.sqrt(comW5.Mv1Squared(H11r));
        double comW5toH12 = Math.sqrt(comW5.Mv1Squared(H12r));
        double comW5toM1 = Math.sqrt(comW5.Mv1Squared(M1r));

        double comW5toH21 = Math.sqrt(comW5.Mv1Squared(H21r));
        double comW5toH22 = Math.sqrt(comW5.Mv1Squared(H22r));
        double comW5toM2 = Math.sqrt(comW5.Mv1Squared(M2r));

        double comW5toH31 = Math.sqrt(comW5.Mv1Squared(H31r));
        double comW5toH32 = Math.sqrt(comW5.Mv1Squared(H32r));
        double comW5toM3 = Math.sqrt(comW5.Mv1Squared(M3r));

        double comW5toH41 = Math.sqrt(comW5.Mv1Squared(H41r));
        double comW5toH42 = Math.sqrt(comW5.Mv1Squared(H42r));
        double comW5toM4 = Math.sqrt(comW5.Mv1Squared(M4r));

        
        // Contributions to sum from water molecule#1       
        Eq5XcompW1 += chargeH11*(comW5.x(0)-H11r.x(0))/(comW5toH11*comW5toH11*comW5toH11)*((1-SpecialFunctions.erfc(comW5toH11/sqrtHMsigmas))-Math.sqrt(2)*comW5toH11/sqrtPiHMsigmas*Math.exp(-comW5toH11*comW5toH11/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq5XcompW1 += chargeH12*(comW5.x(0)-H12r.x(0))/(comW5toH12*comW5toH12*comW5toH12)*((1-SpecialFunctions.erfc(comW5toH12/sqrtHMsigmas))-Math.sqrt(2)*comW5toH12/sqrtPiHMsigmas*Math.exp(-comW5toH12*comW5toH12/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq5XcompW1 += chargeM1*(comW5.x(0)-M1r.x(0))/(comW5toM1*comW5toM1*comW5toM1)*((1-SpecialFunctions.erfc(comW5toM1/(2*sigmaM)))-Math.sqrt(2)*comW5toM1/sqrtPiMMsigmas*Math.exp(-comW5toM1*comW5toM1/(4*sigmaM*sigmaM)));

        Eq5YcompW1 += chargeH11*(comW5.x(1)-H11r.x(1))/(comW5toH11*comW5toH11*comW5toH11)*((1-SpecialFunctions.erfc(comW5toH11/sqrtHMsigmas))-Math.sqrt(2)*comW5toH11/sqrtPiHMsigmas*Math.exp(-comW5toH11*comW5toH11/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq5YcompW1 += chargeH12*(comW5.x(1)-H12r.x(1))/(comW5toH12*comW5toH12*comW5toH12)*((1-SpecialFunctions.erfc(comW5toH12/sqrtHMsigmas))-Math.sqrt(2)*comW5toH12/sqrtPiHMsigmas*Math.exp(-comW5toH12*comW5toH12/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq5YcompW1 += chargeM1*(comW5.x(1)-M1r.x(1))/(comW5toM1*comW5toM1*comW5toM1)*((1-SpecialFunctions.erfc(comW5toM1/(2*sigmaM)))-Math.sqrt(2)*comW5toM1/sqrtPiMMsigmas*Math.exp(-comW5toM1*comW5toM1/(4*sigmaM*sigmaM)));

        Eq5ZcompW1 += chargeH11*(comW5.x(2)-H11r.x(2))/(comW5toH11*comW5toH11*comW5toH11)*((1-SpecialFunctions.erfc(comW5toH11/sqrtHMsigmas))-Math.sqrt(2)*comW5toH11/sqrtPiHMsigmas*Math.exp(-comW5toH11*comW5toH11/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq5ZcompW1 += chargeH12*(comW5.x(2)-H12r.x(2))/(comW5toH12*comW5toH12*comW5toH12)*((1-SpecialFunctions.erfc(comW5toH12/sqrtHMsigmas))-Math.sqrt(2)*comW5toH12/sqrtPiHMsigmas*Math.exp(-comW5toH12*comW5toH12/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq5ZcompW1 += chargeM1*(comW5.x(2)-M1r.x(2))/(comW5toM1*comW5toM1*comW5toM1)*((1-SpecialFunctions.erfc(comW5toM1/(2*sigmaM)))-Math.sqrt(2)*comW5toM1/sqrtPiMMsigmas*Math.exp(-comW5toM1*comW5toM1/(4*sigmaM*sigmaM)));

        
        // Contributions to sum from water molecule#2
        Eq5XcompW2 += chargeH21*(comW5.x(0)-H21r.x(0))/(comW5toH21*comW5toH21*comW5toH21)*((1-SpecialFunctions.erfc(comW5toH21/sqrtHMsigmas))-Math.sqrt(2)*comW5toH21/sqrtPiHMsigmas*Math.exp(-comW5toH21*comW5toH21/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq5XcompW2 += chargeH22*(comW5.x(0)-H22r.x(0))/(comW5toH22*comW5toH22*comW5toH22)*((1-SpecialFunctions.erfc(comW5toH22/sqrtHMsigmas))-Math.sqrt(2)*comW5toH22/sqrtPiHMsigmas*Math.exp(-comW5toH22*comW5toH22/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq5XcompW2 += chargeM2*(comW5.x(0)-M2r.x(0))/(comW5toM2*comW5toM2*comW5toM2)*((1-SpecialFunctions.erfc(comW5toM2/(2*sigmaM)))-Math.sqrt(2)*comW5toM2/sqrtPiMMsigmas*Math.exp(-comW5toM2*comW5toM2/(4*sigmaM*sigmaM)));

        Eq5YcompW2 += chargeH21*(comW5.x(1)-H21r.x(1))/(comW5toH21*comW5toH21*comW5toH21)*((1-SpecialFunctions.erfc(comW5toH21/sqrtHMsigmas))-Math.sqrt(2)*comW5toH21/sqrtPiHMsigmas*Math.exp(-comW5toH21*comW5toH21/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq5YcompW2 += chargeH22*(comW5.x(1)-H22r.x(1))/(comW5toH22*comW5toH22*comW5toH22)*((1-SpecialFunctions.erfc(comW5toH22/sqrtHMsigmas))-Math.sqrt(2)*comW5toH22/sqrtPiHMsigmas*Math.exp(-comW5toH22*comW5toH22/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq5YcompW2 += chargeM2*(comW5.x(1)-M2r.x(1))/(comW5toM2*comW5toM2*comW5toM2)*((1-SpecialFunctions.erfc(comW5toM2/(2*sigmaM)))-Math.sqrt(2)*comW5toM2/sqrtPiMMsigmas*Math.exp(-comW5toM2*comW5toM2/(4*sigmaM*sigmaM)));

        Eq5ZcompW2 += chargeH21*(comW5.x(2)-H21r.x(2))/(comW5toH21*comW5toH21*comW5toH21)*((1-SpecialFunctions.erfc(comW5toH21/sqrtHMsigmas))-Math.sqrt(2)*comW5toH21/sqrtPiHMsigmas*Math.exp(-comW5toH21*comW5toH21/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq5ZcompW2 += chargeH22*(comW5.x(2)-H22r.x(2))/(comW5toH22*comW5toH22*comW5toH22)*((1-SpecialFunctions.erfc(comW5toH22/sqrtHMsigmas))-Math.sqrt(2)*comW5toH22/sqrtPiHMsigmas*Math.exp(-comW5toH22*comW5toH22/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq5ZcompW2 += chargeM2*(comW5.x(2)-M2r.x(2))/(comW5toM2*comW5toM2*comW5toM2)*((1-SpecialFunctions.erfc(comW5toM2/(2*sigmaM)))-Math.sqrt(2)*comW5toM2/sqrtPiMMsigmas*Math.exp(-comW5toM2*comW5toM2/(4*sigmaM*sigmaM)));
        

        // Contributions to sum from water molecule#3
        Eq5XcompW3 += chargeH31*(comW5.x(0)-H31r.x(0))/(comW5toH31*comW5toH31*comW5toH31)*((1-SpecialFunctions.erfc(comW5toH31/sqrtHMsigmas))-Math.sqrt(2)*comW5toH31/sqrtPiHMsigmas*Math.exp(-comW5toH31*comW5toH31/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq5XcompW3 += chargeH32*(comW5.x(0)-H32r.x(0))/(comW5toH32*comW5toH32*comW5toH32)*((1-SpecialFunctions.erfc(comW5toH32/sqrtHMsigmas))-Math.sqrt(2)*comW5toH32/sqrtPiHMsigmas*Math.exp(-comW5toH32*comW5toH32/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq5XcompW3 += chargeM3*(comW5.x(0)-M3r.x(0))/(comW5toM3*comW5toM3*comW5toM3)*((1-SpecialFunctions.erfc(comW5toM3/(2*sigmaM)))-Math.sqrt(2)*comW5toM3/sqrtPiMMsigmas*Math.exp(-comW5toM3*comW5toM3/(4*sigmaM*sigmaM)));

        Eq5YcompW3 += chargeH31*(comW5.x(1)-H31r.x(1))/(comW5toH31*comW5toH31*comW5toH31)*((1-SpecialFunctions.erfc(comW5toH31/sqrtHMsigmas))-Math.sqrt(2)*comW5toH31/sqrtPiHMsigmas*Math.exp(-comW5toH31*comW5toH31/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq5YcompW3 += chargeH32*(comW5.x(1)-H32r.x(1))/(comW5toH32*comW5toH32*comW5toH32)*((1-SpecialFunctions.erfc(comW5toH32/sqrtHMsigmas))-Math.sqrt(2)*comW5toH32/sqrtPiHMsigmas*Math.exp(-comW5toH32*comW5toH32/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq5YcompW3 += chargeM3*(comW5.x(1)-M3r.x(1))/(comW5toM3*comW5toM3*comW5toM3)*((1-SpecialFunctions.erfc(comW5toM3/(2*sigmaM)))-Math.sqrt(2)*comW5toM3/sqrtPiMMsigmas*Math.exp(-comW5toM3*comW5toM3/(4*sigmaM*sigmaM)));

        Eq5ZcompW3 += chargeH31*(comW5.x(2)-H31r.x(2))/(comW5toH31*comW5toH31*comW5toH31)*((1-SpecialFunctions.erfc(comW5toH31/sqrtHMsigmas))-Math.sqrt(2)*comW5toH31/sqrtPiHMsigmas*Math.exp(-comW5toH31*comW5toH31/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq5ZcompW3 += chargeH32*(comW5.x(2)-H32r.x(2))/(comW5toH32*comW5toH32*comW5toH32)*((1-SpecialFunctions.erfc(comW5toH32/sqrtHMsigmas))-Math.sqrt(2)*comW5toH32/sqrtPiHMsigmas*Math.exp(-comW5toH32*comW5toH32/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq5ZcompW3 += chargeM3*(comW5.x(2)-M3r.x(2))/(comW5toM3*comW5toM3*comW5toM3)*((1-SpecialFunctions.erfc(comW5toM3/(2*sigmaM)))-Math.sqrt(2)*comW5toM3/sqrtPiMMsigmas*Math.exp(-comW5toM3*comW5toM3/(4*sigmaM*sigmaM)));


        // Contributions to sum from water molecule#4
        Eq5XcompW4 += chargeH41*(comW5.x(0)-H41r.x(0))/(comW5toH41*comW5toH41*comW5toH41)*((1-SpecialFunctions.erfc(comW5toH41/sqrtHMsigmas))-Math.sqrt(2)*comW5toH41/sqrtPiHMsigmas*Math.exp(-comW5toH41*comW5toH41/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq5XcompW4 += chargeH42*(comW5.x(0)-H42r.x(0))/(comW5toH42*comW5toH42*comW5toH42)*((1-SpecialFunctions.erfc(comW5toH42/sqrtHMsigmas))-Math.sqrt(2)*comW5toH42/sqrtPiHMsigmas*Math.exp(-comW5toH42*comW5toH42/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq5XcompW4 += chargeM4*(comW5.x(0)-M4r.x(0))/(comW5toM4*comW5toM4*comW5toM4)*((1-SpecialFunctions.erfc(comW5toM4/(2*sigmaM)))-Math.sqrt(2)*comW5toM4/sqrtPiMMsigmas*Math.exp(-comW5toM4*comW5toM4/(4*sigmaM*sigmaM)));

        Eq5YcompW4 += chargeH41*(comW5.x(1)-H41r.x(1))/(comW5toH41*comW5toH41*comW5toH41)*((1-SpecialFunctions.erfc(comW5toH41/sqrtHMsigmas))-Math.sqrt(2)*comW5toH41/sqrtPiHMsigmas*Math.exp(-comW5toH41*comW5toH41/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq5YcompW4 += chargeH42*(comW5.x(1)-H42r.x(1))/(comW5toH42*comW5toH42*comW5toH42)*((1-SpecialFunctions.erfc(comW5toH42/sqrtHMsigmas))-Math.sqrt(2)*comW5toH42/sqrtPiHMsigmas*Math.exp(-comW5toH42*comW5toH42/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq5YcompW4 += chargeM4*(comW5.x(1)-M4r.x(1))/(comW5toM4*comW5toM4*comW5toM4)*((1-SpecialFunctions.erfc(comW5toM4/(2*sigmaM)))-Math.sqrt(2)*comW5toM4/sqrtPiMMsigmas*Math.exp(-comW5toM4*comW5toM4/(4*sigmaM*sigmaM)));

        Eq5ZcompW4 += chargeH41*(comW5.x(2)-H41r.x(2))/(comW5toH41*comW5toH41*comW5toH41)*((1-SpecialFunctions.erfc(comW5toH41/sqrtHMsigmas))-Math.sqrt(2)*comW5toH41/sqrtPiHMsigmas*Math.exp(-comW5toH41*comW5toH41/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq5ZcompW4 += chargeH42*(comW5.x(2)-H42r.x(2))/(comW5toH42*comW5toH42*comW5toH42)*((1-SpecialFunctions.erfc(comW5toH42/sqrtHMsigmas))-Math.sqrt(2)*comW5toH42/sqrtPiHMsigmas*Math.exp(-comW5toH42*comW5toH42/(2*(sigmaM*sigmaM+sigmaH*sigmaH))));
        Eq5ZcompW4 += chargeM4*(comW5.x(2)-M4r.x(2))/(comW5toM4*comW5toM4*comW5toM4)*((1-SpecialFunctions.erfc(comW5toM4/(2*sigmaM)))-Math.sqrt(2)*comW5toM4/sqrtPiMMsigmas*Math.exp(-comW5toM4*comW5toM4/(4*sigmaM*sigmaM)));

        
        Eq5.setX(0,Eq5XcompW1+Eq5XcompW2+Eq5XcompW3+Eq5XcompW4);
        Eq5.setX(1,Eq5YcompW1+Eq5YcompW2+Eq5YcompW3+Eq5YcompW4);
        Eq5.setX(2,Eq5ZcompW1+Eq5ZcompW2+Eq5ZcompW3+Eq5ZcompW4);

        
         * Finding the tensor used to relate the induced dipole moment Pi with the induced electric field Epi.
         * kmb, 8/9/06
         
        
        double r12 = Math.sqrt(comW1.Mv1Squared(comW2));
        double r13 = Math.sqrt(comW1.Mv1Squared(comW3));
        double r14 = Math.sqrt(comW1.Mv1Squared(comW4));
        double r15 = Math.sqrt(comW1.Mv1Squared(comW5));
        double r23 = Math.sqrt(comW2.Mv1Squared(comW3));
        double r24 = Math.sqrt(comW2.Mv1Squared(comW4));
        double r25 = Math.sqrt(comW2.Mv1Squared(comW5));
        double r34 = Math.sqrt(comW3.Mv1Squared(comW4));
        double r35 = Math.sqrt(comW3.Mv1Squared(comW5));
        double r45 = Math.sqrt(comW4.Mv1Squared(comW5));
        
        Vector3D r12Vector = new Vector3D();
        Vector3D r13Vector = new Vector3D();
        Vector3D r14Vector = new Vector3D();
        Vector3D r15Vector = new Vector3D();
        Vector3D r23Vector = new Vector3D();
        Vector3D r24Vector = new Vector3D();
        Vector3D r25Vector = new Vector3D();
        Vector3D r34Vector = new Vector3D();
        Vector3D r35Vector = new Vector3D();
        Vector3D r45Vector = new Vector3D();
        
        
        r12Vector.Ev1Mv2(comW1,comW2);  // is this the correct direction? kmb, 8/7/06 / Direction doesn't matter; kmb, 8/10/06
        r13Vector.Ev1Mv2(comW1,comW3);  // is this the correct direction? kmb, 8/7/06 / Direction doesn't matter; kmb, 8/10/06
        r14Vector.Ev1Mv2(comW1,comW4);  // is this the correct direction? kmb, 8/7/06 / Direction doesn't matter; kmb, 8/10/06
        r15Vector.Ev1Mv2(comW1,comW5);  // is this the correct direction? kmb, 8/7/06 / Direction doesn't matter; kmb, 8/10/06
        r23Vector.Ev1Mv2(comW2,comW3);  // is this the correct direction? kmb, 8/7/06 / Direction doesn't matter; kmb, 8/10/06
        r24Vector.Ev1Mv2(comW2,comW4);  // is this the correct direction? kmb, 8/7/06 / Direction doesn't matter; kmb, 8/10/06
        r25Vector.Ev1Mv2(comW2,comW5);  // is this the correct direction? kmb, 8/7/06 / Direction doesn't matter; kmb, 8/10/06
        r34Vector.Ev1Mv2(comW3,comW4);  // is this the correct direction? kmb, 8/7/06 / Direction doesn't matter; kmb, 8/10/06
        r35Vector.Ev1Mv2(comW3,comW5);  // is this the correct direction? kmb, 8/7/06 / Direction doesn't matter; kmb, 8/10/06
        r45Vector.Ev1Mv2(comW4,comW5);  // is this the correct direction? kmb, 8/7/06 / Direction doesn't matter; kmb, 8/10/06
        
        
        double f12 = (1-SpecialFunctions.erfc(r12/(2*sigmaM)))-(r12/(sigmaM*Math.sqrt(Math.PI)) + (r12*r12*r12)/(6*Math.sqrt(Math.PI)*sigmaM*sigmaM*sigmaM))*Math.exp(-r12*r12/(4*sigmaM*sigmaM));
        double f13 = (1-SpecialFunctions.erfc(r13/(2*sigmaM)))-(r13/(sigmaM*Math.sqrt(Math.PI)) + (r13*r13*r13)/(6*Math.sqrt(Math.PI)*sigmaM*sigmaM*sigmaM))*Math.exp(-r13*r13/(4*sigmaM*sigmaM));
        double f14 = (1-SpecialFunctions.erfc(r14/(2*sigmaM)))-(r14/(sigmaM*Math.sqrt(Math.PI)) + (r14*r14*r14)/(6*Math.sqrt(Math.PI)*sigmaM*sigmaM*sigmaM))*Math.exp(-r14*r14/(4*sigmaM*sigmaM));
        double f15 = (1-SpecialFunctions.erfc(r15/(2*sigmaM)))-(r15/(sigmaM*Math.sqrt(Math.PI)) + (r15*r15*r15)/(6*Math.sqrt(Math.PI)*sigmaM*sigmaM*sigmaM))*Math.exp(-r15*r15/(4*sigmaM*sigmaM));
        double f23 = (1-SpecialFunctions.erfc(r23/(2*sigmaM)))-(r23/(sigmaM*Math.sqrt(Math.PI)) + (r23*r23*r23)/(6*Math.sqrt(Math.PI)*sigmaM*sigmaM*sigmaM))*Math.exp(-r23*r23/(4*sigmaM*sigmaM));
        double f24 = (1-SpecialFunctions.erfc(r24/(2*sigmaM)))-(r24/(sigmaM*Math.sqrt(Math.PI)) + (r24*r24*r24)/(6*Math.sqrt(Math.PI)*sigmaM*sigmaM*sigmaM))*Math.exp(-r24*r24/(4*sigmaM*sigmaM));
        double f25 = (1-SpecialFunctions.erfc(r25/(2*sigmaM)))-(r25/(sigmaM*Math.sqrt(Math.PI)) + (r25*r25*r25)/(6*Math.sqrt(Math.PI)*sigmaM*sigmaM*sigmaM))*Math.exp(-r25*r25/(4*sigmaM*sigmaM));
        double f34 = (1-SpecialFunctions.erfc(r34/(2*sigmaM)))-(r34/(sigmaM*Math.sqrt(Math.PI)) + (r34*r34*r34)/(6*Math.sqrt(Math.PI)*sigmaM*sigmaM*sigmaM))*Math.exp(-r34*r34/(4*sigmaM*sigmaM));
        double f35 = (1-SpecialFunctions.erfc(r35/(2*sigmaM)))-(r35/(sigmaM*Math.sqrt(Math.PI)) + (r35*r35*r35)/(6*Math.sqrt(Math.PI)*sigmaM*sigmaM*sigmaM))*Math.exp(-r35*r35/(4*sigmaM*sigmaM));
        double f45 = (1-SpecialFunctions.erfc(r45/(2*sigmaM)))-(r45/(sigmaM*Math.sqrt(Math.PI)) + (r45*r45*r45)/(6*Math.sqrt(Math.PI)*sigmaM*sigmaM*sigmaM))*Math.exp(-r45*r45/(4*sigmaM*sigmaM));
        
        double g12 = (1-SpecialFunctions.erfc(r12/(2*sigmaM)))-(r12/(sigmaM*Math.sqrt(Math.PI)))*Math.exp(-r12*r12/(4*sigmaM*sigmaM));
        double g13 = (1-SpecialFunctions.erfc(r13/(2*sigmaM)))-(r13/(sigmaM*Math.sqrt(Math.PI)))*Math.exp(-r13*r13/(4*sigmaM*sigmaM));
        double g14 = (1-SpecialFunctions.erfc(r14/(2*sigmaM)))-(r14/(sigmaM*Math.sqrt(Math.PI)))*Math.exp(-r14*r14/(4*sigmaM*sigmaM));
        double g15 = (1-SpecialFunctions.erfc(r15/(2*sigmaM)))-(r15/(sigmaM*Math.sqrt(Math.PI)))*Math.exp(-r15*r15/(4*sigmaM*sigmaM));
        double g23 = (1-SpecialFunctions.erfc(r23/(2*sigmaM)))-(r23/(sigmaM*Math.sqrt(Math.PI)))*Math.exp(-r23*r23/(4*sigmaM*sigmaM));
        double g24 = (1-SpecialFunctions.erfc(r24/(2*sigmaM)))-(r24/(sigmaM*Math.sqrt(Math.PI)))*Math.exp(-r24*r24/(4*sigmaM*sigmaM));
        double g25 = (1-SpecialFunctions.erfc(r25/(2*sigmaM)))-(r25/(sigmaM*Math.sqrt(Math.PI)))*Math.exp(-r25*r25/(4*sigmaM*sigmaM));
        double g34 = (1-SpecialFunctions.erfc(r34/(2*sigmaM)))-(r34/(sigmaM*Math.sqrt(Math.PI)))*Math.exp(-r34*r34/(4*sigmaM*sigmaM));
        double g35 = (1-SpecialFunctions.erfc(r35/(2*sigmaM)))-(r35/(sigmaM*Math.sqrt(Math.PI)))*Math.exp(-r35*r35/(4*sigmaM*sigmaM));
        double g45 = (1-SpecialFunctions.erfc(r45/(2*sigmaM)))-(r45/(sigmaM*Math.sqrt(Math.PI)))*Math.exp(-r45*r45/(4*sigmaM*sigmaM));
        
        // Filling the unit matrix I
        
        double[][] I = new double[3][3];
                
        int i = 0;
        int j = 0;
        
        while (i < 3) {
        		while (j < 3) {
        			I[i][j] = 1;
            		j = j + 1;
        		}
        		i = i + 1;
        }
        
        Tensor3D I = new Tensor3D();
        
        I.E(1);
        
//        double[][] T12 = new double[3][3];

        Tensor3D T12 = new Tensor3D();
        
        T12.E(r12Vector,r12Vector);
        T12.TE(3*f12/(r12*r12));
        
        I.TE(g12);
        
        T12.ME(I);
        T12.TE(1/(r12*r12*r12));
        
        // T12 = T21, so I can get by for now in the case of B2!

        
        I.E(1);
        
//        double[][] T12 = new double[3][3];

        Tensor3D T13 = new Tensor3D();
        
        T13.E(r13Vector,r13Vector);
        T13.TE(3*f13/(r13*r13));
        
        I.TE(g13);
        
        T13.ME(I);
        T13.TE(1/(r13*r13*r13));
        

        I.E(1);
        
//        double[][] T12 = new double[3][3];

        Tensor3D T14 = new Tensor3D();
        
        T14.E(r14Vector,r14Vector);
        T14.TE(3*f14/(r14*r14));
        
        I.TE(g14);
        
        T14.ME(I);
        T14.TE(1/(r14*r14*r14));

        
        
        
        I.E(1);
        

        
        Tensor3D T15 = new Tensor3D();
        
        T15.E(r15Vector,r15Vector);
        T15.TE(3*f15/(r15*r15));
        
        I.TE(g15);
        
        T15.ME(I);
        T15.TE(1/(r15*r15*r15));

        
        
        
        I.E(1);

        //        double[][] T12 = new double[3][3];

        Tensor3D T23 = new Tensor3D();
        
        T23.E(r23Vector,r23Vector);
        T23.TE(3*f23/(r23*r23));
        
        I.TE(g23);
        
        T23.ME(I);
        T23.TE(1/(r23*r23*r23));


        I.E(1);
        
//        double[][] T12 = new double[3][3];

        Tensor3D T24 = new Tensor3D();
        
        T24.E(r24Vector,r24Vector);
        T24.TE(3*f24/(r24*r24));
        
        I.TE(g24);
        
        T24.ME(I);
        T24.TE(1/(r24*r24*r24));

        
        I.E(1);
        

        
        Tensor3D T25 = new Tensor3D();
        
        T25.E(r25Vector,r25Vector);
        T25.TE(3*f25/(r25*r25));
        
        I.TE(g25);
        
        T25.ME(I);
        T25.TE(1/(r25*r25*r25));

        
        I.E(1);
        //        double[][] T12 = new double[3][3];

        Tensor3D T34 = new Tensor3D();
        
        T34.E(r34Vector,r34Vector);
        T34.TE(3*f34/(r34*r34));
        
        I.TE(g34);
        
        T34.ME(I);
        T34.TE(1/(r34*r34*r34));

        I.E(1);

        
        Tensor3D T35 = new Tensor3D();
        
        T35.E(r35Vector,r35Vector);
        T35.TE(3*f35/(r35*r35));
        
        I.TE(g35);
        
        T35.ME(I);
        T35.TE(1/(r35*r35*r35));

        I.E(1);
        
        
        Tensor3D T45 = new Tensor3D();
        
        T45.E(r45Vector,r45Vector);
        T45.TE(3*f45/(r45*r45));
        
        I.TE(g45);
        
        T45.ME(I);
        T45.TE(1/(r45*r45*r45));

        I.E(1);
        
        // Now distribute the elements of the tensor into 3 separate "row" vectors
        // so I can do dot products with etomica math methods
        // kmb, 8/7/06
        
        Vector3D T12row1 = new Vector3D();
        Vector3D T12row2 = new Vector3D();
        Vector3D T12row3 = new Vector3D();
        
        T12row1.setX(0,T12.component(0,0));
        T12row1.setX(1,T12.component(0,1));
        T12row1.setX(2,T12.component(0,2));
        T12row2.setX(0,T12.component(1,0));
        T12row2.setX(1,T12.component(1,1));
        T12row2.setX(2,T12.component(1,2));
        T12row3.setX(0,T12.component(2,0));
        T12row3.setX(1,T12.component(2,1));
        T12row3.setX(2,T12.component(2,2));
        

        Vector3D T13row1 = new Vector3D();
        Vector3D T13row2 = new Vector3D();
        Vector3D T13row3 = new Vector3D();
        
        T13row1.setX(0,T13.component(0,0));
        T13row1.setX(1,T13.component(0,1));
        T13row1.setX(2,T13.component(0,2));
        T13row2.setX(0,T13.component(1,0));
        T13row2.setX(1,T13.component(1,1));
        T13row2.setX(2,T13.component(1,2));
        T13row3.setX(0,T13.component(2,0));
        T13row3.setX(1,T13.component(2,1));
        T13row3.setX(2,T13.component(2,2));


        Vector3D T14row1 = new Vector3D();
        Vector3D T14row2 = new Vector3D();
        Vector3D T14row3 = new Vector3D();
        
        T14row1.setX(0,T14.component(0,0));
        T14row1.setX(1,T14.component(0,1));
        T14row1.setX(2,T14.component(0,2));
        T14row2.setX(0,T14.component(1,0));
        T14row2.setX(1,T14.component(1,1));
        T14row2.setX(2,T14.component(1,2));
        T14row3.setX(0,T14.component(2,0));
        T14row3.setX(1,T14.component(2,1));
        T14row3.setX(2,T14.component(2,2));
        
        Vector3D T15row1 = new Vector3D();
        Vector3D T15row2 = new Vector3D();
        Vector3D T15row3 = new Vector3D();
        
        T15row1.setX(0,T15.component(0,0));
        T15row1.setX(1,T15.component(0,1));
        T15row1.setX(2,T15.component(0,2));
        T15row2.setX(0,T15.component(1,0));
        T15row2.setX(1,T15.component(1,1));
        T15row2.setX(2,T15.component(1,2));
        T15row3.setX(0,T15.component(2,0));
        T15row3.setX(1,T15.component(2,1));
        T15row3.setX(2,T15.component(2,2));

        
        
        Vector3D T23row1 = new Vector3D();
        Vector3D T23row2 = new Vector3D();
        Vector3D T23row3 = new Vector3D();
        
        T23row1.setX(0,T23.component(0,0));
        T23row1.setX(1,T23.component(0,1));
        T23row1.setX(2,T23.component(0,2));
        T23row2.setX(0,T23.component(1,0));
        T23row2.setX(1,T23.component(1,1));
        T23row2.setX(2,T23.component(1,2));
        T23row3.setX(0,T23.component(2,0));
        T23row3.setX(1,T23.component(2,1));
        T23row3.setX(2,T23.component(2,2));


        Vector3D T24row1 = new Vector3D();
        Vector3D T24row2 = new Vector3D();
        Vector3D T24row3 = new Vector3D();
        
        T24row1.setX(0,T24.component(0,0));
        T24row1.setX(1,T24.component(0,1));
        T24row1.setX(2,T24.component(0,2));
        T24row2.setX(0,T24.component(1,0));
        T24row2.setX(1,T24.component(1,1));
        T24row2.setX(2,T24.component(1,2));
        T24row3.setX(0,T24.component(2,0));
        T24row3.setX(1,T24.component(2,1));
        T24row3.setX(2,T24.component(2,2));

        
        Vector3D T25row1 = new Vector3D();
        Vector3D T25row2 = new Vector3D();
        Vector3D T25row3 = new Vector3D();
        
        T25row1.setX(0,T25.component(0,0));
        T25row1.setX(1,T25.component(0,1));
        T25row1.setX(2,T25.component(0,2));
        T25row2.setX(0,T25.component(1,0));
        T25row2.setX(1,T25.component(1,1));
        T25row2.setX(2,T25.component(1,2));
        T25row3.setX(0,T25.component(2,0));
        T25row3.setX(1,T25.component(2,1));
        T25row3.setX(2,T25.component(2,2));
        
        Vector3D T34row1 = new Vector3D();
        Vector3D T34row2 = new Vector3D();
        Vector3D T34row3 = new Vector3D();
        
        T34row1.setX(0,T34.component(0,0));
        T34row1.setX(1,T34.component(0,1));
        T34row1.setX(2,T34.component(0,2));
        T34row2.setX(0,T34.component(1,0));
        T34row2.setX(1,T34.component(1,1));
        T34row2.setX(2,T34.component(1,2));
        T34row3.setX(0,T34.component(2,0));
        T34row3.setX(1,T34.component(2,1));
        T34row3.setX(2,T34.component(2,2));
        
        
        Vector3D T35row1 = new Vector3D();
        Vector3D T35row2 = new Vector3D();
        Vector3D T35row3 = new Vector3D();
        
        T35row1.setX(0,T35.component(0,0));
        T35row1.setX(1,T35.component(0,1));
        T35row1.setX(2,T35.component(0,2));
        T35row2.setX(0,T35.component(1,0));
        T35row2.setX(1,T35.component(1,1));
        T35row2.setX(2,T35.component(1,2));
        T35row3.setX(0,T35.component(2,0));
        T35row3.setX(1,T35.component(2,1));
        T35row3.setX(2,T35.component(2,2));

        Vector3D T45row1 = new Vector3D();
        Vector3D T45row2 = new Vector3D();
        Vector3D T45row3 = new Vector3D();
        
        T45row1.setX(0,T45.component(0,0));
        T45row1.setX(1,T45.component(0,1));
        T45row1.setX(2,T45.component(0,2));
        T45row2.setX(0,T45.component(1,0));
        T45row2.setX(1,T45.component(1,1));
        T45row2.setX(2,T45.component(1,2));
        T45row3.setX(0,T45.component(2,0));
        T45row3.setX(1,T45.component(2,1));
        T45row3.setX(2,T45.component(2,2));
        
        // Set the induced dipole moments equal to 10% of the permanent dipole value
        // kmb, 8/7/06
        P1.setX(0,14.3952507082236);
        P1.setX(1,14.3952507082236);
        P1.setX(2,14.3952507082236);
        P2.setX(0,14.3952507082236);
        P2.setX(1,14.3952507082236);
        P2.setX(2,14.3952507082236);
        P3.setX(0,14.3952507082236);
        P3.setX(1,14.3952507082236);
        P3.setX(2,14.3952507082236);
        P4.setX(0,14.3952507082236);
        P4.setX(1,14.3952507082236);
        P4.setX(2,14.3952507082236);
        P5.setX(0,14.3952507082236);
        P5.setX(1,14.3952507082236);
        P5.setX(2,14.3952507082236);

        
        
        P1old.E(P1);
        P2old.E(P2);
        P3old.E(P3);
        P4old.E(P4);
        P5old.E(P4);

        double deltaP1 = 1.0;
        double deltaP2 = 1.0;
        double deltaP3 = 1.0;
        double deltaP4 = 1.0;
        double deltaP5 = 1.0;
        
        
        while (noSCFforP1 || noSCFforP2 || noSCFforP3 || noSCFforP4 || noSCFforP5) {

    			// First calculate Ep1, based upon guess for P2 and P3 and P4 and P5
        	
        		Ep1.setX(0,T12row1.dot(P2)+T13row1.dot(P3)+T14row1.dot(P4)+T15row1.dot(P5));
        		Ep1.setX(1,T12row2.dot(P2)+T13row2.dot(P3)+T14row2.dot(P4)+T15row2.dot(P5));
        		Ep1.setX(2,T12row3.dot(P2)+T13row3.dot(P3)+T14row3.dot(P4)+T15row3.dot(P5));
        		
        		// Now calculate new P1 from the value of Ep1
        		
        		double alphaPol = 1.444;
        		
        		P1.setX(0,alphaPol*(Eq1.x(0) + Ep1.x(0)));
        		P1.setX(1,alphaPol*(Eq1.x(1) + Ep1.x(1)));
        		P1.setX(2,alphaPol*(Eq1.x(2) + Ep1.x(2)));

        		// Next calculate Ep2
        		
        		Ep2.setX(0,T12row1.dot(P1)+T23row1.dot(P3)+T24row1.dot(P4)+T25row1.dot(P5));
        		Ep2.setX(1,T12row2.dot(P1)+T23row2.dot(P3)+T24row2.dot(P4)+T25row2.dot(P5));
        		Ep2.setX(2,T12row3.dot(P1)+T23row3.dot(P3)+T24row3.dot(P4)+T25row3.dot(P5));
        		
        		// Now calculate new P2
        		
        		P2.setX(0,alphaPol*(Eq2.x(0) + Ep2.x(0)));
        		P2.setX(1,alphaPol*(Eq2.x(1) + Ep2.x(1)));
        		P2.setX(2,alphaPol*(Eq2.x(2) + Ep2.x(2)));
        		
        		// Next calculate Ep3
        		
        		Ep3.setX(0,T13row1.dot(P1)+T23row1.dot(P2)+T34row1.dot(P4)+T35row1.dot(P5));
        		Ep3.setX(1,T13row2.dot(P1)+T23row2.dot(P2)+T34row2.dot(P4)+T35row2.dot(P5));
        		Ep3.setX(2,T13row3.dot(P1)+T23row3.dot(P2)+T34row3.dot(P4)+T35row3.dot(P5));
        		
        		// Now calculate new P3
        		
        		P3.setX(0,alphaPol*(Eq3.x(0) + Ep3.x(0)));
        		P3.setX(1,alphaPol*(Eq3.x(1) + Ep3.x(1)));
        		P3.setX(2,alphaPol*(Eq3.x(2) + Ep3.x(2)));

        		
        		// Next calculate Ep4
        		
        		Ep4.setX(0,T14row1.dot(P1)+T24row1.dot(P2)+T34row1.dot(P3)+T45row1.dot(P5));
        		Ep4.setX(1,T14row2.dot(P1)+T24row2.dot(P2)+T34row2.dot(P3)+T45row2.dot(P5));
        		Ep4.setX(2,T14row3.dot(P1)+T24row3.dot(P2)+T34row3.dot(P3)+T45row3.dot(P5));
        		
        		// Now calculate new P4
        		
        		P4.setX(0,alphaPol*(Eq4.x(0) + Ep4.x(0)));
        		P4.setX(1,alphaPol*(Eq4.x(1) + Ep4.x(1)));
        		P4.setX(2,alphaPol*(Eq4.x(2) + Ep4.x(2)));


        		// Next calculate Ep5
        		
        		Ep5.setX(0,T15row1.dot(P1)+T25row1.dot(P2)+T35row1.dot(P3)+T45row1.dot(P4));
        		Ep5.setX(1,T15row2.dot(P1)+T25row2.dot(P2)+T35row2.dot(P3)+T45row2.dot(P4));
        		Ep5.setX(2,T15row3.dot(P1)+T25row3.dot(P2)+T35row3.dot(P3)+T45row3.dot(P4));
        		
        		// Now calculate new P5
        		
        		P5.setX(0,alphaPol*(Eq5.x(0) + Ep5.x(0)));
        		P5.setX(1,alphaPol*(Eq5.x(1) + Ep5.x(1)));
        		P5.setX(2,alphaPol*(Eq5.x(2) + Ep5.x(2)));

        		// Evaluate the criteria
        		
	        	deltaP1 = Math.sqrt(P1.Mv1Squared(P1old));
	        	deltaP2 = Math.sqrt(P2.Mv1Squared(P2old));
	        	deltaP3 = Math.sqrt(P3.Mv1Squared(P3old));
	        	deltaP4 = Math.sqrt(P4.Mv1Squared(P4old));
	        	deltaP5 = Math.sqrt(P5.Mv1Squared(P5old));
	        	
	        	counterSCFloop = counterSCFloop + 1;
	    
	        	
	        	if (deltaP1<1e-15) {
	        		noSCFforP1 = false;
	        	}
	        	else {
	        		P1old.E(P1);
	        	}
	        	if (deltaP2<1e-15) {
	        		noSCFforP2 = false;
	        	}
	        	else {
	        		P2old.E(P2);
	        	}
	        	if (deltaP3<1e-15) {
	        		noSCFforP3 = false;
	        	}
	        	else {
	        		P3old.E(P3);
	        	}
	        	if (deltaP4<1e-15) {
	        		noSCFforP4 = false;
	        	}
	        	else {
	        		P4old.E(P4);
	        	}
	        	if (deltaP5<1e-15) {
	        		noSCFforP5 = false;
	        	}
	        	else {
	        		P5old.E(P5);
	        	}

	        	
	        	if (counterSCFloop >= 1000) {
	        		counterSCFloopOK = false;
	        		noSCFforP1 = false;
	        		noSCFforP2 = false;
	        		noSCFforP3 = false;
	        		noSCFforP4 = false;
	        		noSCFforP5 = false;
	        		loopFailures = loopFailures + 1;
//	        		System.out.println("counterSCFloop = " + counterSCFloop + ", exiting SCF loop due to likely large repulsion");
    		        chargeH11 = Electron.UNIT.toSim(0.6113);
    		        chargeH12 = Electron.UNIT.toSim(0.6113);
    		        chargeM1 = Electron.UNIT.toSim(-1.2226);
    		        chargeH21 = Electron.UNIT.toSim(0.6113);
    		        chargeH22 = Electron.UNIT.toSim(0.6113);
    		        chargeM2 = Electron.UNIT.toSim(-1.2226);
    		        chargeH31 = Electron.UNIT.toSim(0.6113);
    		        chargeH32 = Electron.UNIT.toSim(0.6113);
    		        chargeM3 = Electron.UNIT.toSim(-1.2226);
    		        chargeH41 = Electron.UNIT.toSim(0.6113);
    		        chargeH42 = Electron.UNIT.toSim(0.6113);
    		        chargeM4 = Electron.UNIT.toSim(-1.2226);
    		        chargeH51 = Electron.UNIT.toSim(0.6113);
    		        chargeH52 = Electron.UNIT.toSim(0.6113);
    		        chargeM5 = Electron.UNIT.toSim(-1.2226);
	        	}
	     
	        // REPEAT LOOP HERE.  RE-EVALUATE CHARGES ON W1 AND THEN REPEAT.
        }
        
        
        
         * Here is where I need to add the polarization term to the energy sum.
         * kmb 5/4/06
                 
                
                double alpha = 1.444;
                
                double uW1, uW2, uW3, uW4, uW5, UpolAtkins;
                        
                UpolAtkins = -0.5*(P1.dot(Eq1)+P2.dot(Eq2)+P3.dot(Eq3)+P4.dot(Eq4)+P5.dot(Eq5));
                       
                
                if (rO1O2 <= 2.1 || rO1O3 <= 2.1 || rO1O4 <= 2.1 || rO1O5 <= 2.1 || rO2O3 <= 2.1 || rO2O4 <= 2.1 || rO2O5 <= 2.1 || rO3O4 <= 2.1 || rO3O5 <= 2.1 || rO4O5 <= 2.1) UpolAtkins = 0; // value from Cummings
                
                sumSCF += UpolAtkins;
                
                
                uW1 = 1.855 + Debye.UNIT.fromSim(Math.sqrt(P1.squared()));
                uW2 = 1.855 + Debye.UNIT.fromSim(Math.sqrt(P2.squared()));
                uW3 = 1.855 + Debye.UNIT.fromSim(Math.sqrt(P3.squared()));
                uW4 = 1.855 + Debye.UNIT.fromSim(Math.sqrt(P4.squared()));
                uW5 = 1.855 + Debye.UNIT.fromSim(Math.sqrt(P5.squared()));
                
            //System.out.println("Actually going to return an energy not equal to infinity" + sum);
            //if (rOO >= 500) System.out.println("O-O distance = " + rOO + ", sum = " + sum);

                //System.out.println("uW1 = " + uW1 + "uW2 = " + uW2 + "uW3 = " + uW3 + "uW4 = " + uW4 + "uW5 = " + uW5);
                //System.out.println("tetramer energy = " + sumSCF);
                
        return sumSCF;
        
//        return scfSums;
    }//end of energy
*/    
    
    
    public IVector[] gradient(AtomSet pair){
        throw new MethodNotImplementedException();
    }
    public double hyperVirial(AtomSet pair){
        throw new MethodNotImplementedException();
    }
    public double integral(double rC){
        throw new MethodNotImplementedException();
    }
    public double virial(AtomSet pair){
        throw new MethodNotImplementedException();
    }

    public double getSigma() {return sigma;}

    private final void setSigma(double s) {
        sigma = s;
        sigma2 = s*s;
    }

    public final double getRange() {
        return Double.POSITIVE_INFINITY;
    }
    
    public double getEpsilon() {return epsilon;}
    
    private final void setEpsilon(double eps) {
        epsilon = eps;
        epsilon4 = 4*epsilon;
    }
    
/*    // Original setCharges() method; meant for non-polarizable, pair-wise additive potentials
    private final void setCharges() {
        chargeOO = chargeO * chargeO;
        chargeOH = chargeO * chargeH;
        chargeHH = chargeH * chargeH;
        chargeMM = chargeM * chargeM;
        chargeOM = chargeO * chargeM;
        chargeHM = chargeH * chargeM;
    }
*/
    
    

    public void setBox(Box box) {
    }

/*    public void setPhase(Phase phase) {
        cPair.setNearestImageTransformer(phase.boundary());
    }
*/
    
    public final double getPolarizationEnergy() {
        return UpolAtkins;
    }
    
    public double sigma , sigma2, sumO2LJ;
    public double epsilon, epsilon4, gamma;
//    public int counterSCFloop;
//    public boolean counterSCFloopOK;
    private double chargeH11;
	private double chargeH12;//= Electron.UNIT.toSim(0.52);
    private double chargeM1; //= Electron.UNIT.toSim(-1.04);
    private double chargeH21;
    private double chargeH22;//= Electron.UNIT.toSim(0.52);
    private double chargeM2; //= Electron.UNIT.toSim(-1.04);
    private double chargeH31;
    private double chargeH32;
    private double chargeM3;
    private double chargeH41;
    private double chargeH42;
    private double chargeM4;
    private double chargeH51;
    private double chargeH52;
    private double chargeM5;
    private IVector Eq1, Eq2, Eq3, Eq4, Eq5;
    private IVector Ep1, Ep2, Ep3, Ep4, Ep5;
    private IVector P1, P2, P3, P4, P5;
    private IVector P1old, P2old, P3old, P4old, P5old;
    private IVector comW1, comW2, comW3, comW4, comW5;
    private final IVector r12Vector;
    private final IVector T12P1, T12P2;
    private final IVector T12row1, T12row2, T12row3;
    private final IVector B1, T12Eq2;
    private final Tensor Tunit, T12, A1, T12T12, inverseA1;
    private final double core; // = 4.41; //4.41 = 2.1^2; value according to Cummings
    private final double sigmaM;
    private final double sigmaH;
    private final double sqrtHMsigmas;
    private final double massH;
    private final double massO;
    private final double totalMass;
    private final double sqrtPiHMsigmas;
    private final double sqrtPiMMsigmas;
    public double UpolAtkins;
	/* (non-Javadoc)
	 * @see etomica.potential.PotentialSoft#gradient(etomica.atom.AtomSet, etomica.space.Tensor)
	 */
	public IVector[] gradient(AtomSet atoms, Tensor pressureTensor) {
		// TODO Auto-generated method stub
		return null;
	}



    
//    private double deltaP1, deltaP2, deltaP3;

}