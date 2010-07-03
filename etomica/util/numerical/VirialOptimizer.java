package etomica.util.numerical;

public class VirialOptimizer {
	
	public VirialOptimizer(String filenameP, String filenameB){
		
		pSimData = ArrayReader1D.getFromFile(filenameP);
		bVirialFromFile = ArrayReader1D.getFromFile(filenameB);
		bVirial = new double[bVirialFromFile.length + 1];
		rho = new double[pSimData.length];
		pPade = new double[pSimData.length];
		
		for (int i=0; i<bVirialFromFile.length;i++){
			bVirial[i] = bVirialFromFile[i][0];
		}
		
		for (int i=0; i<pSimData.length;i++){
			rho[i] = pSimData[i][0];
		}
		
	}
	
	public void calcpPade(int K, int L, double bGuess){
		/*
		 * Pade Approximation
		 */
		double pPadeNum, pPadeDenom;
		
		/*
		 * Optimizing higher-order B-Coefficient 
		 */
		bVirial[bVirial.length-1] = bGuess;
		
		PadeApproximation pade = new PadeApproximation(bVirial, K, L);
		pade.solveCoefficients();
		double[] A = pade.getA();
		double[] C = pade.getC();
		
		/*
		 * Pressure Calculation from Pade
		 */
		
		for (int irho=0; irho< rho.length; irho++){
			pPadeNum = 0.0;
			pPadeDenom = 0.0;
						
			for (int i=0; i<A.length;i++){
				pPadeNum += A[i]*Math.pow(rho[irho],i);
			}
			
			for (int i=0; i<C.length;i++){
				pPadeDenom += C[i]*Math.pow(rho[irho],i);
			}
			
			pPade[irho] = pPadeNum/pPadeDenom;
			
		}
	}
	
	public void optimizeHigherbVirial(int K, int L, double bGuess){
		
		//double bGuess = 3.1;
		
		calcpPade(K, L, bGuess);
		
		double TSS = 0.0;
		
		for (int i=0; i<pPade.length; i++){
			diff = (pPade[i]-pSimData[i][1]);
		//	System.out.println(i +" " + diff);
			TSS += diff*diff;
			
		}
		System.out.print( TSS);
	}
	
	
	public static void main(String[] args){
		String filenameP = "/tmp/dataExtrapolated1.dat";
		String filenameB = "/tmp/Bn6";
		
		VirialOptimizer vOpt = new VirialOptimizer(filenameP, filenameB);
		
		for (double i = 3.0; i < 3.1; i+=0.001){
			System.out.print("\n" + i +" ");
			vOpt.optimizeHigherbVirial(5, 3, i);		
		}
	
	}
	
	protected double diff;
	protected double[] bVirial, rho, pPade;
	protected double[][] pSimData, bVirialFromFile;
}
