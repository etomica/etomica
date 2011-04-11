package etomica.potential;

import etomica.api.IAtom;
import etomica.api.IAtomList;
import etomica.api.IBoundary;
import etomica.api.IBox;
import etomica.api.IVector;
import etomica.api.IVectorMutable;
import etomica.atom.Atom;
import etomica.atom.AtomArrayList;
import etomica.space.ISpace;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.space3d.Space3D;
import etomica.units.Angle;
import etomica.units.Dimension;
import etomica.units.Energy;

/**
 * Ab initio non-additive trimer potential for He developed by Cencek, Patkowski, and Szalewicz JCP 131 064105 2009. 
 * @author kate
 */
public class P3CPSNonAdditiveHe extends Potential implements PotentialSoft {

    public P3CPSNonAdditiveHe(ISpace space) {
        super(3, space);
        drAB = space.makeVector();
        drBC = space.makeVector();
        drAC = space.makeVector();
        setAngle(Math.PI);
        gradient = new IVectorMutable[3];
        gradient[0] = space.makeVector();
        gradient[1] = space.makeVector();
        gradient[2] = space.makeVector();
        
        
    }

    public void setBox(IBox box) {
        boundary = box.getBoundary();
    }

    public double energy(IAtomList atomSet) {
    	
    	setA();
        setAlpha();
        setBeta3();
        setZ3();
        setBeta4220();
        setBeta4211();
        setZ4220();
        setZ4211();
        
        IAtom atomA = atomSet.getAtom(0);
        IAtom atomB = atomSet.getAtom(1);
        IAtom atomC = atomSet.getAtom(2);
        
        atomA.getPosition().TE(1.0/AngstromPerBohrRadius);
        atomB.getPosition().TE(1.0/AngstromPerBohrRadius);
        atomC.getPosition().TE(1.0/AngstromPerBohrRadius);
        
        drAB.Ev1Mv2(atomA.getPosition(),atomB.getPosition());
        drAC.Ev1Mv2(atomA.getPosition(),atomC.getPosition());
        drBC.Ev1Mv2(atomB.getPosition(),atomC.getPosition());
        
        double RAB = Math.sqrt(drAB.squared());
        double RAC = Math.sqrt(drAC.squared());
        double RBC = Math.sqrt(drBC.squared());
        System.out.println (RAB + " " + RBC+ " " + RAC);
        
        double costhetaA =  drAB.dot(drAC)/(RAB*RAC);
        double costhetaB = -drAB.dot(drBC)/(RAB*RBC);
        double costhetaC =  drAC.dot(drBC)/(RAC*RBC);
        
        
        
        double theta1; double theta2; double theta3;
        
        if (costhetaA > 1) { theta1 = 0;}
        else if (costhetaA < -1) { theta1 = Math.PI;}
        else {theta1 = Math.acos(costhetaA);}
        
        if (costhetaB > 1) {theta2 = 0;}
        else if (costhetaB < -1) {theta2 = Math.PI;}
        else {theta2 = Math.acos(costhetaB);}
        
        if (costhetaC > 1) {theta3 = 0;}
        else if (costhetaC < -1) {theta3 = Math.PI;}
        else {theta3 = Math.acos(costhetaC);}
        
        System.out.println(theta1 + " " + theta2 + " " + theta3 + "  " + (theta1+theta2+theta3)/Math.PI);
        
        double Vexp = 0;

        for (int k3 = 0; k3<=4; k3++) {
        	for (int k2 = 0; k2<=k3; k2++) {
        		for (int k1 = 0; k1<=k2; k1++) {
        			
        			double P = legendreP(k1,costhetaA)*legendreP(k2,costhetaB)*legendreP(k3,costhetaC);
        			P =    P + legendreP(k1,costhetaA)*legendreP(k2,costhetaC)*legendreP(k3,costhetaB);
        			P =    P + legendreP(k1,costhetaB)*legendreP(k2,costhetaA)*legendreP(k3,costhetaC);
        			P =    P + legendreP(k1,costhetaB)*legendreP(k2,costhetaC)*legendreP(k3,costhetaA);
        			P =    P + legendreP(k1,costhetaC)*legendreP(k2,costhetaA)*legendreP(k3,costhetaB);
        			P =    P + legendreP(k1,costhetaC)*legendreP(k2,costhetaB)*legendreP(k3,costhetaA);
        			
        			Vexp = Vexp + A[k1][k2][k3]*Math.exp(-alpha[k1][k2][k3]*(RAB+RBC+RAC))*P;
        			
        			//System.out.println(k1 + " " + k2 + " " + k3 + "  " + alpha[k1][k2][k3]);
        		}
        	}
        }
        
        //V3Disp
        
        double W111 = 3.0*Math.pow(RAB,-3)*Math.pow(RBC,-3)*Math.pow(RAC,-3);
               W111 = W111*(1.0 + (3.0*Math.cos(theta1)*Math.cos(theta2)*Math.cos(theta3) ));
        double n111AB = 3; double n111BC = 3; double n111AC = 3;
      
        double W112 = 3.0/16.0*Math.pow(RAB,-3)*Math.pow(RBC,-4)*Math.pow(RAC,-4);
        double term112 = (9.0*Math.cos(theta3)-25.0*Math.cos(3.0*theta3));
               term112 = term112 + (6.0*Math.cos(theta1-theta2)*(3.0+5.0*Math.cos(2.0*theta3)));
               W112 = W112*term112;
        double n112AB = 3; double n112BC = 4; double n112AC = 4;
        
       
        double W122 = 15.0/64.0*Math.pow(RAB,-4)*Math.pow(RBC,-5)*Math.pow(RAC,-4);
        double term122 = 3*(Math.cos(theta1)+5.0*Math.cos(3*theta1));
        	   // As in Bukowski and Szalewicz 2001
        	   //term122 = term122 + 20.0*Math.cos(theta2-theta3)*(1.0-3.0*Math.cos(2*theta1));
               //term122 = term122 + 70.0*Math.cos(2.0*(theta2-theta3))*Math.cos(theta1);
               // As in Cencek's et al's code
        	   term122 = term122 +   20.0*Math.cos(theta3-theta2)*(1.0-3.0*Math.cos(2*theta1));
        	   term122 = term122 + 70.0*Math.cos(2.0*(theta3-theta2))*Math.cos(theta1);
        	   W122 = W122*term122;
        double n122AB = 4; double n122BC = 5; double n122AC = 4;	   
        
        double W222 = 15.0/128.0*Math.pow(RAB,-5)*Math.pow(RBC,-5)*Math.pow(RAC,-5);
        double term222 = -27.0 + 220.0*Math.cos(theta1)*Math.cos(theta2)*Math.cos(theta3);
        	   // As in Bukowski and Szalewicz 2001
               //term222 = term222 + 490.0*Math.cos(2.0*theta1)*Math.cos(theta2)*Math.cos(theta3);
               // As in Cencek's et al's code
        	   term222 = term222 +   490.0*Math.cos(2.0*theta1)*Math.cos(2.0*theta2)*Math.cos(2.0*theta3);
               term222 = term222 + 175.0*(Math.cos(2.0*(theta1-theta2))+Math.cos(2.0*(theta2-theta3))+Math.cos(2.0*(theta3-theta1)));
               W222 = W222*term222;
        double n222AB = 5; double n222BC = 5; double n222AC = 5;	
        
        double W113 = 5.0/32.0*Math.pow(RAB,-3)*Math.pow(RBC,-5)*Math.pow(RAC,-5);
        double term113 = 9.0 + 8.0*Math.cos(2.0*theta3) - 49.0*Math.cos(4.0*theta3);
               term113 = term113 + 6.0*Math.cos(theta1-theta2)*(9.0*Math.cos(theta3)+7.0*Math.cos(3.0*theta3));
               W113 = W113*term113;
        double n113AB = 3; double n113BC = 5; double n113AC = 5;
       
        double D111 = getD(beta3[1][1][1],RAB,n111AB)*getD(beta3[1][1][1],RBC,n111BC)*getD(beta3[1][1][1],RAC,n111AC);
        double D112 = getD(beta3[1][1][2],RAB,n112AB)*getD(beta3[1][1][2],RBC,n112BC)*getD(beta3[1][1][2],RAC,n112AC);
        double D122 = getD(beta3[1][2][2],RAB,n122AB)*getD(beta3[1][2][2],RBC,n122BC)*getD(beta3[1][2][2],RAC,n122AC);
        double D222 = getD(beta3[2][2][2],RAB,n222AB)*getD(beta3[2][2][2],RBC,n222BC)*getD(beta3[2][2][2],RAC,n222AC);
        double D113 = getD(beta3[1][1][3],RAB,n113AB)*getD(beta3[1][1][3],RBC,n113BC)*getD(beta3[1][1][3],RAC,n113AC);
       
   
        double V3disp =   D111*W111*Z3[1][1][1];
        V3disp = V3disp + 3.0*D112*W112*Z3[1][1][2];
        V3disp = V3disp + 3.0*D122*W122*Z3[1][2][2];
        V3disp = V3disp + D222*W222*Z3[2][2][2];
        V3disp = V3disp + 3.0*D113*W113*Z3[1][1][3];
		
        //V4Disp
        double W1111_211 = 36.0/128.0*Math.pow(RAB,-6)*Math.pow(RBC,-3)*Math.pow(RAC,-3);
        double term1111_211 = -3.0 + Math.cos(theta2+theta3)*Math.cos(theta2+theta3);
               term1111_211 = term1111_211 + Math.cos(theta1+theta3)*Math.cos(theta1+theta3);
               term1111_211 = term1111_211 + 5.0*Math.cos(theta1+theta2)*Math.cos(theta1+theta2);
               W1111_211 = W1111_211*term1111_211;
   	    double D1111_211 = getD(beta4211[1][1][1][1],RAB,6)*getD(beta4211[1][1][1][1],RBC,3)*getD(beta4211[1][1][1][1],RAC,3);     
   	    
   	    
     	// As in Bukowski and Szalewicz 2001
   	   /* double W1111_220 =             (1.0+Math.cos(theta1)*Math.cos(theta1))*Math.pow(RAB,-6)*Math.pow(RBC,-6);
   	           W1111_220 = W1111_220 + (1.0+Math.cos(theta2)*Math.cos(theta2))*Math.pow(RBC,-6)*Math.pow(RAC,-6);
   	           W1111_220 = W1111_220 + (1.0+Math.cos(theta3)*Math.cos(theta3))*Math.pow(RAB,-6)*Math.pow(RAC,-6);
               W1111_220 = W1111_220*9.0;
        double D1111_220 = getD(beta4220[1][1][1][1],RAB,6)*getD(beta4220[1][1][1][1],RBC,6)*getD(beta4220[1][1][1][1],RAC,6);
	    */
   	    
   	    // As Cencek et al code:
    	double W1111_220a = 9.0*(1.0+Math.cos(theta1)*Math.cos(theta1))*Math.pow(RAB,-6)*Math.pow(RBC,-6);
    	double W1111_220b = 9.0*(1.0+Math.cos(theta2)*Math.cos(theta2))*Math.pow(RBC,-6)*Math.pow(RAC,-6);
    	double W1111_220c =  9.0*(1.0+Math.cos(theta3)*Math.cos(theta3))*Math.pow(RAB,-6)*Math.pow(RAC,-6);
               
        double D1111_220a = getD(beta4220[1][1][1][1],RAB,6)*getD(beta4220[1][1][1][1],RBC,6)*getD(beta4220[1][1][1][1],RAC,0);
        double D1111_220b = getD(beta4220[1][1][1][1],RAB,0)*getD(beta4220[1][1][1][1],RBC,6)*getD(beta4220[1][1][1][1],RAC,6);
        double D1111_220c = getD(beta4220[1][1][1][1],RAB,6)*getD(beta4220[1][1][1][1],RBC,0)*getD(beta4220[1][1][1][1],RAC,6);
        
        double W1112_211 = 1.0/32.0*Math.pow(RAB,-7)*Math.pow(RBC,-3)*Math.pow(RAC,-4);
        double term1112_211 = -144.0*Math.cos(theta1) + 36.0*Math.cos(theta2 + theta3);
               term1112_211 = term1112_211 + 216.0*Math.cos(theta2-theta3) - 120.0*Math.cos(3.0*theta3);
               term1112_211 = term1112_211 - 720.0*Math.cos(theta1-2.0*theta3) - 72.0*Math.cos(theta1-2*theta2);
               W1112_211 = W1112_211*term1112_211;
   	    double D1112_211 = getD(beta4211[1][1][1][2],RAB,7)*getD(beta4211[1][1][1][2],RBC,3)*getD(beta4211[1][1][1][2],RAC,4); 
   	    
   	    
   	    double W1121_211 = 1.0/32.0*Math.pow(RAB,-6)*Math.pow(RBC,-4)*Math.pow(RAC,-4);
   	    double term1121_211 = -111.0*Math.cos(theta3) - 750.0*Math.cos(3.0*theta3);
               term1121_211 = term1121_211 + 180.0*Math.cos(theta1+theta3) + 108.0*Math.cos(theta1-theta2);
               term1121_211 = term1121_211 - 90.0*Math.cos(theta3-2.0*theta1) - 90.0*Math.cos(theta3-2*theta2);
                  W1121_211 = W1121_211*term1121_211;
	    double D1121_211 = getD(beta4211[1][1][2][1],RAB,6)*getD(beta4211[1][1][2][1],RBC,4)*getD(beta4211[1][1][2][1],RAC,4); 
	    
   	    
   	    double W2111_211 = -9.0/2.0*Math.pow(RAB,-8)*Math.pow(RBC,-3)*Math.pow(RAC,-3);
   	    double term2111_211 = Math.cos(2.0*theta1) + Math.cos(2.0*theta2) + 6.0*Math.cos(2.0*theta3);
                  W2111_211 = W2111_211*term2111_211;
	    double D2111_211 = getD(beta4211[2][1][1][1],RAB,8)*getD(beta4211[2][1][1][1],RBC,3)*getD(beta4211[2][1][1][1],RAC,3);     
               
	    double W1211_220 = -1.0/64.0*Math.pow(RAB,-7)*Math.pow(RBC,-7);
               W1211_220 = W1211_220 * (1485.0*Math.cos(theta3)+ 384.0*Math.cos(3.0*theta3));
        double D1211_220 = getD(beta4220[1][2][1][1],RAB,7)*getD(beta4220[1][2][1][1],RBC,7);
 
        double W2111_220 = 0.25*Math.pow(RAB,-8)*Math.pow(RBC,-6);
        	   W2111_220 = W2111_220 * (369.0 + 288.0*Math.cos(theta3)*Math.cos(theta3));
        double D2111_220 = getD(beta4220[2][1][1][1],RAB,8)*getD(beta4220[2][1][1][1],RBC,6)*getD(beta4220[2][1][1][1],RAC,0);

        // As in Bukowski and Szalewicz 2001
        /*
        double V4disp =   D1111_211*W1111_211*Z4211[1][1][1][1]; 
        V4disp = V4disp + D1111_220*W1111_220*Z4220[1][1][1][1];

        V4disp = V4disp + D1112_211*W1112_211*Z4211[1][1][1][2];
        
        V4disp = V4disp + D1121_211*W1121_211*Z4211[1][1][2][1];
        
        V4disp = V4disp + D2111_211*W2111_211*Z4211[2][1][1][1];
        V4disp = V4disp + D2111_220*W2111_220*Z4220[2][1][1][1];
        
        V4disp = V4disp + D1211_220*W1211_220*Z4220[1][2][1][1];
        */
        
        //As in Cencek's code
        double V4disp =   0;
        V4disp = V4disp + D1111_220a*W1111_220a*Z4220[1][1][1][1];
        V4disp = V4disp + D1111_220b*W1111_220b*Z4220[1][1][1][1];
        V4disp = V4disp + D1111_220c*W1111_220c*Z4220[1][1][1][1];

        V4disp = V4disp + 6.0*D1112_211*W1112_211*Z4211[1][1][1][2];
        
        V4disp = V4disp + 3.0*D1121_211*W1121_211*Z4211[1][1][2][1];
        
        V4disp = V4disp + 3.0*D2111_211*W2111_211*Z4211[2][1][1][1];
        V4disp = V4disp + 6.0*D2111_220*W2111_220*Z4220[2][1][1][1];
        
        V4disp = V4disp + 3.0*D1211_220*W1211_220*Z4220[1][2][1][1];
        
        //return (Vexp+V3disp+V4disp)*KPerHartree; //Kelvin
        return (Vexp+V3disp+V4disp)*KPerHartree; //Kelvin
        
    }
    

    
    public double getD(double beta, double RXY, double nXY) {
    	
    	double D = 1.0;	
    	double factorial = 1.0;
    	for (int n=1;n<=nXY;n++) {
    		
    		factorial = factorial*n;
    		
    		D = D + Math.pow(beta*RXY,n)/factorial;
    	}
    	
    	D = 1.0 - (Math.exp(-beta*RXY)*D);

    	return D;
    }
    
    
    
    public double legendreP (int i, double x) {
    	
    	if (i == 0) {
    		return 1.0;
    	} else if (i == 1) {
    		return x;
    	} else if (i == 2) {
    		//1/2(3x^2-1)	
    		return 0.5*(3.0*x*x-1.0);	
    	} else if (i == 3) {
    		//1/2(5x^3-3x)
    		return 0.5*(5.0*x*x*x-3.0*x);
    	} else if (i == 4) {
    		//1/8(35x^4-30x^2+3)	
    		return 1.0/8.0*(35.0*x*x*x*x-30.0*x*x+3.0);	
    	} else if (i == 5) {
    		// 1/8(63x^5-70x^3+15x)
    		return 1.0/8.0*(63.0*x*x*x*x*x - 70.0*x*x*x + 15.0*x);
    	} else {
    		throw new RuntimeException("Cannot do that order of Legendre polynomial");
    	}
    }
    
    public void setAlpha() {
    	alpha[0][0][0]=1.16406382984624;
        alpha[0][0][1]=0.593775469740520;
        alpha[0][0][2]=0.579991620360140;
        alpha[0][0][3]=2.16547015214329;
        alpha[0][0][4]=2.98314825911820;
        alpha[0][1][1]=1.41509085761783;
        alpha[0][1][2]=1.68952244399482;
        alpha[0][1][3]=4.01020983745692;
        alpha[0][1][4]=0.918500573699908;
        alpha[0][2][2]=1.24672879488806;
        alpha[0][2][3]=2.05681119414852;
        alpha[0][2][4]=1.31995140236607;
        alpha[0][3][3]=3.73423077881244;
        alpha[0][3][4]=2.03094608487497;
        alpha[0][4][4]=2.16349588270954;
        alpha[1][1][1]=1.14951445564366;
        alpha[1][1][2]=1.13500643353835;
        alpha[1][1][3]=0.706124946374114;
        alpha[1][1][4]=3.38198341068258;
        alpha[1][2][2]=0.586893129829570;
        alpha[1][2][3]=2.50706679697449;
        alpha[1][2][4]=0.947942951899726;
        alpha[1][3][3]=1.32340552994745;
        alpha[1][3][4]=0.400000000000000;
        alpha[1][4][4]=3.08785463812471;
        alpha[2][2][2]=0.958505100822341;
        alpha[2][2][3]=1.18369957176966;
        alpha[2][2][4]=2.88051393389548;
        alpha[2][3][3]=3.89009793612663;
        alpha[2][3][4]=1.75941163648810;
        alpha[2][4][4]=0.663027178763952;
        alpha[3][3][3]=3.85457885867751;
        alpha[3][3][4]=1.10245630997729;
        alpha[3][4][4]=0.928262289570450;
        alpha[4][4][4]=1.10072817300826;
       
    }
    
    public void setA() {
    	A[0][0][0]=7.33779142427628;
    	A[0][0][1]=0.207562291096732E-02;
    	A[0][0][2]=0.201257721465354E-02;
    	A[0][0][3]=-162280.921808245;
    	A[0][0][4]=140426930.572053;
    	A[0][1][1]=3676.77127295706;
    	A[0][1][2]=40721.4699374272;
    	A[0][1][3]=-145478302663.879;
    	A[0][1][4]=-2.38575533677150;
    	A[0][2][2]=24.9222823036941;
    	A[0][2][3]=-136123.231454042;
    	A[0][2][4]=-2036.85749653572;
    	A[0][3][3]=65196789563.7813;
    	A[0][3][4]=-246369.934348762;
    	A[0][4][4]=337609.745545195;
    	A[1][1][1]=-303.389041542620;
    	A[1][1][2]=419.847355533508;
    	A[1][1][3]=0.428281282658233E-01;
    	A[1][1][4]=-6376463380.36163;
    	A[1][2][2]=0.121898031059658E-01;
    	A[1][2][3]=25447929.5505372;
    	A[1][2][4]=1.81213617023749;
    	A[1][3][3]=-2092.82734569787;
    	A[1][3][4]=0.681066685456400E-05;
    	A[1][4][4]=-146635205.381890;
    	A[2][2][2]=6.60802549511257;
    	A[2][2][3]=-344.693694202977;
    	A[2][2][4]=143064734.239803;
    	A[2][3][3]=-1008432134063.46;
    	A[2][3][4]=-9312.13937942436;
    	A[2][4][4]=-0.386910250300881E-02;
    	A[3][3][3]=711864422051.715;
    	A[3][3][4]=18.0536606525932;
    	A[3][4][4]=0.922564385936391;
    	A[4][4][4]=-4.52285479839575;

    }

    public void setBeta3() {
    	beta3[1][1][1]=0.850816031004730;
    	beta3[1][1][2]=1.03935993289613;
    	beta3[1][1][3]=2.35163790098234;
    	beta3[1][2][2]=20.0000000000000;
    	beta3[2][2][2]=7.74979337816275;
    }
    
    public void setZ3(){
    	Z3[1][1][1]=0.49311;
    	Z3[1][1][2]=0.92372;
    	Z3[1][1][3]=4.1241;
    	Z3[1][2][2]=1.7377;
    	Z3[2][2][2]=3.2839;
    }
    
    
    public void setBeta4211() {	
/*    	beta4211[1][1][2][2]=2.22023197004267;
    	beta4211[1][2][2][1]=2.33977220590245;
    	beta4211[2][1][1][2]=1.96782469219456;*/
    	
    	beta4211[1][1][1][2]=2.22023197004267;
    	beta4211[1][1][2][1]=2.33977220590245;
    	beta4211[2][1][1][1]=1.96782469219456;
    }
    
    public void setZ4211(){
    	Z4211[1][1][1][2]=-370.838300778413;
    	Z4211[1][1][2][1]=673.766716043939;
    	Z4211[2][1][1][1]=-553.474291722504;
    }
    
	public void setBeta4220(){
		beta4220[1][1][1][1]=1.76277419240966;
		beta4220[1][2][1][1]=2.13546395662687;
		beta4220[2][1][1][1]=0.959706781068175;

	}
	
	public void setZ4220(){
		Z4220[1][1][1][1]=-15.2910806164061;
		Z4220[1][2][1][1]=158.205832955569;
		Z4220[2][1][1][1]=112.479143795999;
	}
    
    
    /**
     * Sets the nominal bond angle (in radians)
     */
    public void setAngle(double newAngle) {
        angle = newAngle;
    }
    
    /**
     * Returns the nominal bond angle (in radians)
     */
    public double getAngle() {
        return angle;
    }
    
    public Dimension getAngleDimension() {
        return Angle.DIMENSION;
    }

    /**
     * Sets the characteristic energy of the potential
     */
    public void setEpsilon(double newEpsilon) {
        epsilon = newEpsilon;
    }
    
    /**
     * Returns the characteristic energy of the potential
     */
    public double getEpsilon() {
        return epsilon;
    }
    
    public Dimension getEpsilonDimension() {
        return Energy.DIMENSION;
    }
    
    public double getRange() {
        return Double.POSITIVE_INFINITY;
    }

    public IVector[] gradient(IAtomList atoms) {
       throw new RuntimeException("Sorry, no gradient available yet");
    }

    public IVector[] gradient(IAtomList atoms, Tensor pressureTensor) {
        return gradient(atoms);
    }

    public double virial(IAtomList atoms) {
        return 0;
    }

    protected final IVectorMutable drAB, drAC, drBC;
    protected IBoundary boundary;
    protected double angle;
    protected double epsilon;
    private static final long serialVersionUID = 1L;
    protected final IVectorMutable[] gradient;
    public static boolean bigAngle;
    protected double[][][] alpha = new double [5][5][5];
    protected double[][][] A = new double [5][5][5];
    protected double[][][] beta3 = new double [3][3][4];
    protected double[][][] Z3 = new double [3][3][4];
    protected double[][][][] beta4220= new double [3][3][3][3];
    protected double[][][][] Z4220= new double [3][3][3][3];
    protected double[][][][] beta4211= new double [3][3][3][3];
    protected double[][][][] Z4211= new double [3][3][3][3];
    private static final double AngstromPerBohrRadius = 0.529177; // Rounding provided by Pryzbytek et al. 2010
    private static final double KPerHartree = 315774.65; // Rounding provided by Pryzbytek et al. 2010
    
    
    
    public static void main(String[] args) {
        Space space = Space3D.getInstance();

        P3CPSNonAdditiveHe potential = new P3CPSNonAdditiveHe(space);
      
        Atom atom0 = new Atom(space);
        Atom atom1 = new Atom(space);
        Atom atom2 = new Atom(space);
        
        AtomArrayList atoms = new AtomArrayList(3);
        atoms.add(atom0);
        atoms.add(atom1);
        atoms.add(atom2);
        
        // Equilateral triangle 
        double a = 7.0*AngstromPerBohrRadius;
        IVector r0 = (IVector)space.makeVector(new double[] {0,0,0});
        IVector r1 = (IVector)space.makeVector(new double[] {a,0,0});
        IVector r2 = (IVector)space.makeVector(new double[] {a/2.0,a/2.0*Math.sqrt(3),0});
        
        // 4th config
        a = 5.6*AngstromPerBohrRadius;
        r0 = (IVector)space.makeVector(new double[] {0,0,0});
        r1 = (IVector)space.makeVector(new double[] {a,0,0});
        r2 = (IVector)space.makeVector(new double[] {2*a,0,0});
       
        
        atom0.getPosition().E(r0);
        atom1.getPosition().E(r1);
        atom2.getPosition().E(r2);
        
        
            
        double U = potential.energy(atoms);

        System.out.println(U*1000); //millikelvin
            
            
        
    }
    

}


