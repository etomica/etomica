package etomica.normalmode;

/**
 * 
 * Soft-Sphere Solid Equation of State
 *  for n = 12, 9 and 6
 * 
 * Reference: Tan, Schultz and Kofke, J Chem. Phys.,133(2010)
 * 			   Appendix: Free Energy Formulas
 * 
 * @author Tai Boon Tan
 *
 */
public class SoftSphereSolidEOS {
	
	public SoftSphereSolidEOS(){
		
	}
		
	public double getUlat0(int id, int n){
		return ulat0[id][n/3-2];
	}
	
	public double getVib0(int id, int n){
		return vib0[id][n/3-2];
	}
	
	public double[] getCoeffAc(int id , int n){
		return c[id][n/3-2];
	}
	
	public double[] getAcQuantity(int id, int n, double temp, double rho){
		double bAc = 0.0;
		double dAcdrho = 0.0;
		double dbAcdbeta = 0.0;
		
		double[] c = getCoeffAc(id, n);
		double upsilon = getUpsilon(n, temp, rho);
		double ups = 1; 
		
		for(int i=1; i<=7; i++){
			ups *= upsilon;
			bAc += c[i]*ups;
			
			dAcdrho += temp*c[i]*ups*(-n*i/(3.0*rho));
			dbAcdbeta += c[i]*ups*(-i*temp);
		}
		
		double[] values = new double[]{bAc, dAcdrho, dbAcdbeta};
		
		return values;
	}
	
	public double[] getAllQuantity(int id, int n, double temp, double rho){
		double upsilon = getUpsilon(n, temp, rho);
		double ulat0 = getUlat0(id, n);
		double Vib0 = getVib0(id, n);
		double[] AcValues = getAcQuantity(id, n, temp, rho); 
		
		double bUlat = ulat0/upsilon; 			            // Equation (A2)
		double bAVib = Vib0 - 3.0/2.0*Math.log(upsilon) 
		                + Math.log(rho*sigma*sigma*sigma);  // Equation (A3)
		double bAc   = AcValues[0];			 				// Equation ( 9)
		double bA    = bUlat + bAVib + bAc; 				// Equation (A1) 
		
		double p = rho*rho*((n/3.0)*(temp/(upsilon*rho))*ulat0  	
				   + AcValues[1]
				   + (n*temp)/(2*rho) + 1/rho);
		
		double u = temp/upsilon*ulat0
				   + AcValues[1]
				   + (3.0/2.0)*temp;
		double[] quantity = new double[]{bA, p/(rho*temp), u/temp};
		return quantity;
	}
	
	public double getUpsilon(int n, double temp, double rho){
		return (temp/epsilon)*Math.pow((rho*sigma*sigma*sigma), -n/3.0);
	}
	
	public static void main (String[] args){
		double temp = 1.0;
		double rho = 2.20617;
		int n = 6;
		
		SoftSphereSolidEOS ssA = new SoftSphereSolidEOS();
		
		double[] qFCC = ssA.getAllQuantity(FCC, n, temp, rho);
		double[] qBCC = ssA.getAllQuantity(BCC, n, temp, rho);
		
		
		System.out.println("[FCC] betaA: " + qFCC[0] + "; Z: " + qFCC[1] + "; betaU: "+ qFCC[2]);
		System.out.println("[BCC] betaA: " + qBCC[0] + "; Z: " + qBCC[1] + "; betaU: "+ qBCC[2]);
	}
	
    public final static int FCC = 0;
    public final static int HCP = 1;
    public final static int BCC = 2;
    protected double epsilon = 1.0;
    protected double sigma = 1.0;
    
    protected final double[][] ulat0 = new double[][]{{3.613475382,2.208391122,1.516485025},{3.613718378,2.208528128,1.516536721},{3.630711328}};
    protected final double[][] vib0 = new double[][]{{2.69819834,3.48217524,3.88845245},{2.70269996,3.48590058,3.89158992},{2.56473110}};
    protected final double[][][] c = new double[][][]{
    		{{0.0, -0.024270, -0.385254, -1.989540, 24.366184, -208.856577, 860.698810,-1469.431853},
    		 {0.0,  0.249526, -0.245724, -0.500979,  4.258323,  -16.542027,  30.811960,  -23.208566},
    		 {0.0,  0.464148, -0.423821,  0.466322, -0.526094,    0.148125,   0.401893,   -0.384095}},
    		
    		{{0.0, -0.023687, -0.351963, -2.934478, 37.070278, -280.028273, 953.287096, -1140.820226},
    		 {0.0, 0.250099,  -0.282838,  0.007132,  0.922672,   -5.218832,  10.890988,    -7.369508},
    		 {0.0, 0.462065,  -0.413562,  0.429175, -0.461716,    0.111491,   0.368595,    -0.352244}},
    		
    		{{0.0, 0.283913, -2.421458, 14.558777, -77.105689, 229.711416, -391.161702, 335.089915}},
    			 
    };
    
}
