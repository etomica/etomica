package etomica.virial;

import etomica.api.IMoleculeList;

public class MayerSWComponent {
	
	private double sigma;
	private double sigma2;
	private double lamda;
	private double lamda2;
	public static final double Y_Value = 1.0;
	
	public MayerSWComponent(){
		this(1.0, 2.0);
	}
	
	public MayerSWComponent(double sigma, double lamda){
		setSigma(sigma);
		setLamda(lamda);
	}
	
	public double f1(IMoleculeList pair, double r2, double beta){
//		System.out.print("r2="+r2+"  ");
		if (r2 < sigma2){
			
//			System.out.println("1<r2, f1/y return 0;");
			return 0.0;
		}else if (r2 < sigma2*lamda2){
//			System.out.println("1<=r2<4, f1/y return Y_Value");
			return Y_Value;
		}else{
//			System.out.println("4<r2, f1/y return 0");
			return 0.0;
		}
	}
	
	public double e2(IMoleculeList pair, double r2, double beta){
//		System.out.print("r2="+r2+"  ");
		if (r2 < sigma2){
//			System.out.println("r2<sigma2, e2 return 0");
			return 0.0;
		}else{
//			System.out.println("r2>sigma2, e2 return 1");
			return 1.0;
		}
	}
	
	/**
     * Returns the SW diameter.
     * @return double
     */
    public double getSigma() {
        return sigma;
    }
    
    /**
     * Sets the SW diameter.
     * @param sigma The sigma to set
     */
	public void setSigma(double sigma) {
        this.sigma = sigma;
        sigma2 = sigma*sigma;
    }
	
	/**
     * Returns the SW lamda.
     * @return double
     */
    public double getLamda() {
        return lamda;
    }
	
    /**
     * Sets the SW lamda.
     * @param lamda The lamda to set
     */
	public void setLamda(double lamda) {
        this.lamda = lamda;
        lamda2 = lamda*lamda;
    }

}
