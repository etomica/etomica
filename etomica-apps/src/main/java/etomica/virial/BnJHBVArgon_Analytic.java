/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;


/**
 * Computes virial coefficients for argon using the analytical expressions 
 * developed by Jaeger, Hellmann, Bich, and Vogel (2011), JCP, 135: 084308.
 * 
 * The authors computed B4-B7 using MSMC.  They employed only a pair potential and 
 * a nonadditive trimer potential, the later of which they developed in this work. 
 * The potentials are ab initio.
 * 
 * A first-order quantum correction is included on B2-B4.
 * 
 * The nonadditive trimer contribution to B7 is neglected.
 *  
 * @author Kate Shaul
 */
public class BnJHBVArgon_Analytic {

    public BnJHBVArgon_Analytic() {
       
    	setParams();
    	
    }

    public void setParams() {
    	c2[-9+9]=0.00000000000E+00;
    	c2[-8+9]=0.00000000000E+00;
    	c2[-7+9]=0.00000000000E+00;
    	c2[-6+9]=0.00000000000E+00;
    	c2[-5+9]=-6.43967919993E-05;
    	c2[-4+9]=1.23670554715E-03;
    	c2[-3+9]=-4.71320933474E-02;
    	c2[-2+9]=-2.02956133576E-01;
    	c2[-1+9]=-1.91084609343E+01;
    	c2[0+9]=3.36620466678E+01;
    	c2[1+9]=8.42765093583E-01;
    	c2[2+9]=-8.30397887940E-03;
    	c2[3+9]=0.00000000000E+00;
    	c2_sqrt[(int)(0.5+0.5)]=-7.06033409390E+00;
    	c2_sqrt[(int)(-0.5+0.5)]=1.22377682977E+01;
   
    	c3Add[-9+9]=0.00000000000E+00;   
    	c3Add[-8+9]=-2.80484912948E-05;  
    	c3Add[-7+9]=9.97664367737E-05;   
    	c3Add[-6+9]=3.96382358075E-03;   
    	c3Add[-5+9]=-1.94439533716E-01;  
    	c3Add[-4+9]=1.62369882602E+00;   
    	c3Add[-3+9]=-1.07450939747E+01;  
    	c3Add[-2+9]=1.32613577121E+02;  
    	c3Add[-1+9]=-9.33231361502E+02;  
    	c3Add[0+9]=-1.35067272455E+02;  
    	c3Add[1+9]=0.00000000000E+00;   
    	c3Add[2+9]=0.00000000000E+00;   
    	c3Add[3+9]=0.00000000000E+00;   
    	c3Add_sqrt[(int)(0.5+0.5)]=3.80063418193E+00;   
    	c3Add_sqrt[(int)(-0.5+0.5)]=1.59181702940E+03; 
   
    	c3[-9+9]=0.00000000000E+00;
    	c3[-8+9]=-2.59962699927E-05;
    	c3[-7+9]=9.71352136025E-05;
    	c3[-6+9]= 4.01081195695E-03;
    	c3[-5+9]=-1.88028806073E-01;
    	c3[-4+9]=1.64722353365E+00;
    	c3[-3+9]=-1.01735017150E+01;
    	c3[-2+9]=1.37319632897E+02;
    	c3[-1+9]=-9.06798523084E+02;
    	c3[0+9]=-1.34839966042E+02;
    	c3[1+9]=0.00000000000E+00;
    	c3[2+9]=0.00000000000E+00;
    	c3[3+9]=0.00000000000E+00;
    	c3_sqrt[(int)(0.5+0.5)]=3.17008319513E+00;
    	c3_sqrt[(int)(-0.5+0.5)]=1.60786304709E+03;
    
	    c4Add[-9+9]=0.00000000000E+00;
	    c4Add[-8+9]=-2.00127668660E-01;
	    c4Add[-7+9]=5.16098878170E+00;
	    c4Add[-6+9]=-6.35049729674E+01; 
	    c4Add[-5+9]=4.50965629405E+02;
	    c4Add[-4+9]=-1.90825035081E+03;
	    c4Add[-3+9]=5.62937470293E+03;
	    c4Add[-2+9]=-1.29955934123E+04;
	    c4Add[-1+9]=1.47211820405E+04;
	    c4Add[0+9]=-1.44634498664E+03;
	    c4Add[1+9]=0.00000000000E+00;
	    c4Add[2+9]=0.00000000000E+00;
	    c4Add[3+9]=0.00000000000E+00;
	    c4Add_sqrt[(int)(0.5+0.5)]=0.00000000000E+00;
	    c4Add_sqrt[(int)(-0.5+0.5)]=8.68494433126E+03;
    
	   c4[-9+9]=0.00000000000E+00; 
	   c4[-8+9]=-1.67135881019E-01; 
	   c4[-7+9]=4.43369449980E+00; 
	   c4[-6+9]=-5.57705182924E+01; 
	   c4[-5+9]=4.08432542355E+02; 
	   c4[-4+9]=-1.77692888041E+03; 
	   c4[-3+9]=5.34599537332E+03; 
	   c4[-2+9]=-1.32195737029E+04; 
	   c4[-1+9]=1.54637920020E+04; 
	   c4[0+9]=-1.61307216584E+03; 
	   c4[1+9]=0.00000000000E+00; 
	   c4[2+9]=0.00000000000E+00; 
	   c4[3+9]=0.00000000000E+00; 
	   c4_sqrt[(int)(0.5+0.5)]=0.00000000000E+00; 
	   c4_sqrt[(int)(-0.5+0.5)]=9.18220801482E+03; 
	   
	   c5Add[-9+9]=-1.10204929946E+01;  
	   c5Add[-8+9]=2.78496482705E+02;    
	   c5Add[-7+9]=-3.13675356346E+03; 
	   c5Add[-6+9]=2.07849524109E+04;  
	   c5Add[-5+9]=-8.74911585695E+04;  
	   c5Add[-4+9]=2.43050769313E+05;
	   c5Add[-3+9]=-4.83090592400E+05;
	   c5Add[-2+9]=6.92127652177E+05;
	   c5Add[-1+9]=-1.18633976542E+06;
	   c5Add[ 0+9]=-1.07771673399E+06;
	   c5Add[1+9]=-3.95912565103E+04; 
	   c5Add[2+9]=3.29557972763E+02;  
	   c5Add[3+9]=0.00000000000E+00;
	   c5Add_sqrt[(int)(0.5+0.5)]=3.11271507733E+05;
	   c5Add_sqrt[(int)(-0.5+0.5)]=1.84249116502E+06;

	   c5[-9+9]=-8.12931771976E+00;  
	   c5[-8+9]=2.10732135345E+02;  
	   c5[-7+9]=-2.41658003885E+03;  
	   c5[-6+9]=1.63429112656E+04;  
	   c5[-5+9]=-7.05138154213E+04;  
	   c5[-4+9]=1.99415483631E+05;  
	   c5[-3+9]=-4.04948314511E+05;  
	   c5[-2+9]=5.97093981469E+05;  
	   c5[-1+9]=-1.01810828653E+06;  
	   c5[0+9]=-1.01881422623E+06;  
	   c5[1+9]=-3.83451453437E+04;  
	   c5[2+9]=3.30053974723E+02;  
	   c5[3+9]=0.00000000000E+00;  
	   c5_sqrt[(int)(0.5+0.5)]=2.97931080782E+05;  
	   c5_sqrt[(int)(-0.5+0.5)]=1.69700346091E+06;  
	   
	   c6Add[-9+9]=-2.46669683494E+03;  
	   c6Add[-8+9]=6.39699153867E+04;   
	   c6Add[-7+9]=-7.09098062427E+05;   
	   c6Add[-6+9]=4.49042269338E+06; 
	   c6Add[-5+9]=-1.83523054410E+07; 
	   c6Add[-4+9]=5.16617504712E+07; 
	   c6Add[-3+9]=-1.06548568219E+08;    
	   c6Add[-2+9]=1.81814766661E+08;   
	   c6Add[-1+9]=-4.12093357833E+08;  
	   c6Add[0+9]=-3.56150463200E+08; 
	   c6Add[1+9]=-1.99599126696E+07;  
	   c6Add[2+9]=3.55746274518E+05;  
	   c6Add[3+9]=-5.51255214190E+03;  
	   c6Add_sqrt[(int)(0.5+0.5)]=1.23201832960E+08;  
	   c6Add_sqrt[(int)(-0.5+0.5)]=5.55240091541E+08;    

	  c6[-9+9]=-1.90203554800E+03;  
	  c6[-8+9]=5.32262589262E+04;  
	  c6[-7+9]=-6.36172289974E+05;  
	  c6[-6+9]=4.34610669722E+06;  
	  c6[-5+9]=-1.92314811535E+07;  
	  c6[-4+9]=5.89977956039E+07;  
	  c6[-3+9]=-1.33251075455E+08;  
	  c6[-2+9]=2.50157778349E+08;  
	  c6[-1+9]=-6.28197313299E+08;  
	  c6[0+9]=-5.94655339923E+08;  
	  c6[1+9]=-3.62460403720E+07;  
	  c6[2+9]=6.94728540417E+05;  
	  c6[3+9]=-1.14158498664E+04;  
	  c6_sqrt[(int)(0.5+0.5)]=2.14777131405E+08;  
	  c6_sqrt[(int)(-0.5+0.5)]=8.86342715726E+08;  



	 c7Add[-9+9]=-1.47083149394E+05;  
	 c7Add[-8+9]=4.15178399999E+06;  
	 c7Add[-7+9]=-4.86501552311E+07;  
	 c7Add[-6+9]=3.09077440740E+08;  
	 c7Add[-5+9]=-1.17768612888E+09;  
	 c7Add[-4+9]=2.84360158630E+09;  
	 c7Add[-3+9]=-4.52704994280E+09;  
	 c7Add[-2+9]=5.10551112913E+09;  
	 c7Add[-1+9]=-5.82791714953E+09;  
	 c7Add[0+9]=-1.40927512001E+09;  
	 c7Add[1+9]=0.00000000000E+00;  
	 c7Add[2+9]=0.00000000000E+00;  
	 c7Add[3+9]=0.00000000000E+00;  
	 c7Add_sqrt[(int)(0.5+0.5)]=1.55873889035E+08; 
	 c7Add_sqrt[(int)(-0.5+0.5)]=4.59052283453E+09;  




    }
    
    public static  double getBn(double T, double[] c, double[] c_sqrt) {
    	double Bn = 0;
    	double reducedT=T/1000;
    	for (int k=-9;k<=3;k++) {
    		Bn = Bn + c[k+9]*Math.pow(reducedT,k);
    	}
    	
    	Bn = Bn + c_sqrt[1]*Math.sqrt(reducedT) + c_sqrt[0]/Math.sqrt(reducedT);
    	
        return Bn;
    }

    public static void main(String[] args) {
        
        System.out.println("B2 with first quantum correction:");
        System.out.println();
        double T=100; //Kelvin
        for (int i=1;i<12;i++) {
        	
        	double B2=getBn(T, c2, c2_sqrt);
        	System.out.println(T+"  "+B2);
        	T=T+10;
        }
        
        System.out.println();
        System.out.println("B3Add with first quantum correction:");
        System.out.println();
        T=100;
        for (int i=1;i<12;i++) {
        	
        	double B3Add=getBn(T, c3Add, c3Add_sqrt);
        	System.out.println(T+"  "+B3Add);
        	T=T+10;
        }
      
        System.out.println();
        System.out.println("B3 with first quantum correction:");
        System.out.println();
        T=100;
        for (int i=1;i<12;i++) {
        	
        	double B3=getBn(T, c3, c3_sqrt);
        	System.out.println(T+"  "+B3);
        	T=T+10;
        }
        
        System.out.println();
        System.out.println("B4Add  with first quantum correction:");
        System.out.println();
        T=100;
        for (int i=1;i<12;i++) {
        	
        	double B4Add=getBn(T, c4Add, c4Add_sqrt);
        	System.out.println(T+"  "+B4Add);
        	T=T+10;
        }
      
        System.out.println();
        System.out.println("B4 with first quantum correction:");
        System.out.println();
        T=100;
        for (int i=1;i<12;i++) {

        	double B4=getBn(T, c4, c4_sqrt);
        	System.out.println(T+"  "+B4);
        	T=T+10;
        }
        
        System.out.println();
        System.out.println("B5:");
        System.out.println();
        T=100;
        for (int i=1;i<12;i++) {

        	double B5=getBn(T, c5, c5_sqrt);
        	System.out.println(T+"  "+B5);
        	T=T+10;
        }
        
        System.out.println();
        System.out.println("B6:");
        System.out.println();
        T=100;
        for (int i=1;i<12;i++) {

        	double B6=getBn(T, c6, c6_sqrt);
        	System.out.println(T+"  "+B6);
        	T=T+10;
        }
        
        System.out.println();
        System.out.println("B7Add:");
        System.out.println();
        T=100;
        for (int i=1;i<12;i++) {

        	double B7Add=getBn(T, c7Add, c7Add_sqrt);
        	System.out.println(T+"  "+B7Add);
        	T=T+10;
        }
       
        	
       
    }
    

    protected final static double[] c2 = new double [13];
    protected final static double[] c3 = new double [13];
    protected final static double[] c4 = new double [13];
    protected final static double[] c5 = new double [13];
    protected final static double[] c6 = new double [13];
    
    protected final static double[] c3Add = new double [13];
    protected final static double[] c4Add = new double [13];
    protected final static double[] c5Add = new double [13];
    protected final static double[] c6Add = new double [13];
    protected final static double[] c7Add = new double [13];

    protected final static double[] c2_sqrt = new double [2];
    protected final static double[] c3_sqrt = new double [2];
    protected final static double[] c4_sqrt = new double [2];
    protected final static double[] c5_sqrt = new double [2];
    protected final static double[] c6_sqrt = new double [2];
    
    protected final static double[] c3Add_sqrt = new double [2];
    protected final static double[] c4Add_sqrt = new double [2];
    protected final static double[] c5Add_sqrt = new double [2];
    protected final static double[] c6Add_sqrt = new double [2];
    protected final static double[] c7Add_sqrt = new double [2];

}


