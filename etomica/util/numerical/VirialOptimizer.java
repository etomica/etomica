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
	
	public void optimizeHigherbVirial(int K, int L, double initbGuess){
		
		double interval = 0.0001;
		bGuess = initbGuess;
		calcpPade(K, L, initbGuess);
		
		double TSS = 0.0;
		
		for (int i=0; i<pPade.length; i++){
			diff = (pPade[i]-pSimData[i][1]);
			TSS += diff*diff;
			
		}
				
		while(TSS > tol){
			double TSSleft = 0.0;
			double TSSright = 0.0;
			double bGuessleft  = bGuess - interval;
			double bGuessright = bGuess + interval;
			
			calcpPade(K, L, bGuessleft);
			for (int i=0; i<pPade.length; i++){
				diff = (pPade[i]-pSimData[i][1]);
				TSSleft += diff*diff;	
			}
			
			calcpPade(K, L, bGuessright);
			for (int i=0; i<pPade.length; i++){
				diff = (pPade[i]-pSimData[i][1]);
				TSSright += diff*diff;	
			}

			if(Math.abs(TSSright) <= Math.abs(TSSleft)){
				TSS = TSSright;
				bGuess = bGuessright;
			} else {
				TSS = TSSleft;
				bGuess = bGuessleft;
			}
			/*
			 * extra check
			 * break out from the while loop if the difference is small
			 *  although TSS does not meet the tolerance condition
			 */
			if (Math.abs(TSSright-TSSleft) < 1e-9) break;
		}
		
		
		System.out.println(bGuess +" " + TSS);
	}
	
	
	public static void main(String[] args){
		String filenameP = "/tmp/dataExtrapolated1.dat";
		String filenameB = "/tmp/Bn6";
		
		VirialOptimizer vOpt = new VirialOptimizer(filenameP, filenameB);
		
		vOpt.optimizeHigherbVirial(5, 3, 2.0);		
	
	}
	
	protected double diff;
	protected double[] bVirial, rho, pPade;
	protected double[][] pSimData, bVirialFromFile;
	protected double tol = 1e-4;
	protected double bGuess;
	
}
