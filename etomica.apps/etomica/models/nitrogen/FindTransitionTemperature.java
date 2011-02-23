package etomica.models.nitrogen;

import etomica.units.Pascal;
import etomica.util.Constants;
import etomica.util.numerical.AkimaSpline;
import etomica.util.numerical.PolynomialFit;

/**
 * Class that determines the transition properties for alpha and beta nitrogen
 * 
 * @author taitan
 *
 */


public class FindTransitionTemperature {
	
	public FindTransitionTemperature(double[] density, double desiredPressure, double tolerance){
		akimaSpline = new AkimaSpline();
		this.density = density;
		this.desiredPressure = desiredPressure;
		this.tolerance = tolerance;
		
	}
	
	public double[] AandPAlpha(double rho, double temperature){
		double[] a = new double[density.length];
		
		for(int i=0; i<a.length; i++){
			a[i] = AAlphaFunction(i, temperature);
//			System.out.println("a[" + i+"]: " + a[i]);
//			System.out.println(density[i]+" " + a[i]);

		}
		
		/*
		 * Smoothen out the curve by fitting a[i] to 
		 *  2nd order polynomial 
		 */
		double[] r = PolynomialFit.doFit(2, density, a);
		for (int i=0; i<a.length; i++){
			a[i] = r[0] + r[1]*density[i] + r[2]*density[i]*density[i];
//			System.out.println(density[i]+" " + a[i]);
		}
//		System.exit(1);
		akimaSpline.setInputData(density, a);
		double[] aAlpha =  akimaSpline.doInterpolation(new double[]{rho});
		double[] pAlpha = akimaSpline.doInterpolationDy(new double[]{rho});
		double[] AandP = new double[]{aAlpha[0], pAlpha[0]*rho*rho}; 
		
		return  AandP;
	}
	
	public double[] AandPBeta(double rho, double temperature){
		double[] a = new double[density.length];

//		System.out.println("\n"+rho +" ***A**");
		for(int i=0; i<a.length; i++){
			a[i] = ABetaFunction(i, temperature);
//			System.out.println("a[" + i+"]: " + a[i]);
//			System.out.println(density[i]+" " + a[i]);
		}

		/*
		 * Smoothen out the curve by fitting a[i] to 
		 *  2nd order polynomial 
		 */
		double[] r = PolynomialFit.doFit(2, density, a);

		for (int i=0; i<a.length; i++){
			a[i] = r[0] + r[1]*density[i] + r[2]*density[i]*density[i];
//			System.out.println(density[i]+" " + a[i]);
		}
//		System.exit(1);
		akimaSpline.setInputData(density, a);
		double[] aBeta =  akimaSpline.doInterpolation(new double[]{rho});
		double[] pBeta = akimaSpline.doInterpolationDy(new double[]{rho});
		double[] AandP = new double[]{aBeta[0], pBeta[0]*rho*rho}; 
		
		return  AandP;
	}
	
	
	
	
	
	
	public double UAlpha(double TransRho, double TransT){
		double[] betaA = new double[41];
		double[] invTemp = new double[41];
		
		for(int i=0; i<betaA.length; i++){
			double temp = 50.0-i*((50.0-40.0)/(betaA.length-1));
			invTemp[i] = 1.0/temp;
			betaA[i] = AandPAlpha(TransRho, temp)[0]*invTemp[i];
//			System.out.println(invTemp[i] + " " + betaA[i]);
		}
		
		double invT = 1.0/TransT;
		akimaSpline.setInputData(invTemp, betaA);
		double[] u = akimaSpline.doInterpolationDy(new double[]{invT});
		
//		for(int i=0; i<u.length; i++){
//			System.out.println(invTemp[i]+" " + u[i]);
//		}
		
		System.out.println("alpha energy at " + invT  + " " + u[0]);
		return  u[0];
	}
	
	
	public double UBeta(double TransRho, double TransT){
		double[] betaA = new double[41];
		double[] invTemp = new double[41];
		
		for(int i=0; i<betaA.length; i++){
			double temp = 50.0-i*((50.0-40.0)/(betaA.length-1));
			invTemp[i] = 1.0/temp;
			betaA[i] = AandPBeta(TransRho, temp)[0]*invTemp[i];
//			System.out.println(invTemp[i] + " " + betaA[i]);
		}
		
		double invT = 1.0/TransT;
		akimaSpline.setInputData(invTemp, betaA);
		double[] u = akimaSpline.doInterpolationDy(invTemp);
			
		System.out.println("beta energy at " + invT  + " " + u[0]);
		return  u[0];
	}
	
	
	
	
	
	
	
	
	
	public double AAlphaFunction(int i, double temperature){
		double T=temperature;
		double A = 0.0;
		
		if(i==0){
			double[] coeff= new double[]{-6.261357685927708e+02,
					-1.891741383390460e+01,
					1.039551984062851e+00,
					-2.001837643940738e-02,
					1.884494783755010e-05,
					4.048174846846921e-06,
					-3.776690295594055e-08};
			A = coeff[0];
			for (int k=1; k<coeff.length; k++){
				A += coeff[k]*T;
				T *=temperature;
			}

			return A;
			
		}else if(i==1){
			double[] coeff= new double[]{-8.599214805531283e+02,
					9.559413759156875e+00,
					-3.282519910280455e-01,
					1.173516346710990e-02,
					-3.136023961965519e-04,
					4.711082305293204e-06,
					-2.929548653295880e-08};
			A = coeff[0];
			for (int k=1; k<coeff.length; k++){
				A += coeff[k]*T;
				T *=temperature;
			}
			return A;
			
		}else if(i==2){
			double[] coeff= new double[]{1.942246710669904e+04,
					-2.970259498775863e+03,
					1.818887903951307e+02,
					-5.925209043839509e+00,
					1.083993200856117e-01,
					-1.056147337013756e-03,
					4.281191203094416e-06};
			A = coeff[0];
			for (int k=1; k<coeff.length; k++){
				A += coeff[k]*T;
				T *=temperature;
			}
			
			return A;
			
		}else if(i==3){
			double[] coeff= new double[]{-6.807932793560021e+02,
					-2.663560776061603e+00,
					-4.283360715115391e-01,
					4.168387012890450e-02,
					-1.320899095444043e-03,
					1.863008993129741e-05,
					-1.003813962081565e-07};
			A = coeff[0];
			for (int k=1; k<coeff.length; k++){
				A += coeff[k]*T;
				T *=temperature;
			}
	
			return A;
			
		}else if(i==4){
			double[] coeff= new double[]{-1.776146471555584e+02,
					-6.816571472764382e+01,
					3.002069004397685e+00,
					-4.927793870259439e-02,
					-7.340571135126849e-05,
					1.089853416945077e-05,
					-8.843559532947079e-08};
			A = coeff[0];
			for (int k=1; k<coeff.length; k++){
				A += coeff[k]*T;
				T *=temperature;
			}
			
			return A;
			
		}else if(i==5){
			double[] coeff= new double[]{-9.895611525976120e+02,
					3.587237584069572e+01,
					-2.349237867194948e+00,
					8.967565407403655e-02,
					-1.923348808027034e-03,
					2.173590527671551e-05,
					-1.014492837196539e-07};
			A = coeff[0];
			for (int k=1; k<coeff.length; k++){
				A += coeff[k]*T;
				T *=temperature;
			}
			
			return A;
			
		}else if(i==6){
			double[] coeff= new double[]{-8.516756968452440e+02,
					5.891999677075453e+00,
					6.808373049917926e-02,
					-7.592728862855116e-03,
					1.876178869773401e-04,
					-1.987199801728271e-06,
					7.235827231753336e-09};
			A = coeff[0];
			for (int k=1; k<coeff.length; k++){
				A += coeff[k]*T;
				T *=temperature;
			}
			
			return A;
			
		}else if(i==7){
			double[] coeff= new double[]{-1.259040011890240e+03,
					7.190013705627723e+01,
					-4.362104167943449e+00,
					1.502019935062358e-01,
					-2.955884393927122e-03,
					3.121776374460380e-05,
					-1.380334039622383e-07};
			A = coeff[0];
			for (int k=1; k<coeff.length; k++){
				A += coeff[k]*T;
				T *=temperature;
			}
			
			return A;
			
		}else if(i==8){
			double[] coeff= new double[]{1.011789637188801e+03,
					-2.553008121143495e+02,
					1.528617508680717e+01,
					-4.791644955706074e-01,
					8.390063140837461e-03,
					-7.795739297036732e-05,
					3.001482850386152e-07};
			A = coeff[0];
			for (int k=1; k<coeff.length; k++){
				A += coeff[k]*T;
				T *=temperature;
			}
			
			return A;
			
		}else if(i==9){
			double[] coeff= new double[]{3.212288460846401e+01,
					-1.147088595778359e+02,
					6.897096159477881e+00,
					-2.124384454824485e-01,
					3.623795829618482e-03,
					-3.255478382709135e-05,
					1.199790676093946e-07};
			A = coeff[0];
			for (int k=1; k<coeff.length; k++){
				A += coeff[k]*T;
				T *=temperature;
			}
			
			return A;
			
		} else {
			throw new RuntimeException("");
		}
		
	}
	
	public double ABetaFunction(int i, double temperature){
		double T=temperature;
		double A = 0.0;
		
		if(i==0){
			double[] coeff= new double[]{
					7.618616472534738e+01,
					-1.162947174692604e+02,
					6.960496566079405e+00,
					-2.167097389314131e-01,
					3.754508984908238e-03,
					-3.448459739666921e-05,
					1.313590099604301e-07
			
				};
			A = coeff[0];
			for (int k=1; k<coeff.length; k++){
				A += coeff[k]*T;
				T *=temperature;
			}
			
			return A;
			
		}else if(i==1){
			double[] coeff= new double[]{
					-1.313397017073949e+03,
					7.840228820021632e+01,
					-4.374389304305200e+00,
					1.343057863167351e-01,
					-2.342581584091469e-03,
					2.183897379029391e-05,
					-8.483136898112068e-08
			};
			A = coeff[0];
			for (int k=1; k<coeff.length; k++){
				A += coeff[k]*T;
				T *=temperature;
			}
			
			return A;
			
		}else if(i==2){
			double[] coeff= new double[]{
					-1.127623864637319e+03,
					4.915010049813907e+01,
					-2.479201052122164e+00,
					6.969493356132701e-02,
					-1.115260289615615e-03,
					9.510900966651905e-06,
					-3.364017465026284e-08,
			};
			A = coeff[0];
			for (int k=1; k<coeff.length; k++){
				A += coeff[k]*T;
				T *=temperature;
			}
			
			return A;
			
		}else if(i==3){
			double[] coeff= new double[]{
					-2.482194272913129e+03,
					2.389140556940391e+02,
					-1.353272265527601e+01,
					4.124943590116977e-01,
					-7.083602825625894e-03,
					6.482207457822456e-05,
					-2.468018720952815e-07
			};
			A = coeff[0];
			for (int k=1; k<coeff.length; k++){
				A += coeff[k]*T;
				T *=temperature;
			}
			
			return A;
			
		}else if(i==4){
			double[] coeff= new double[]{
					-4.478004075135981e+01,
					-9.994428032689008e+01,
					6.048596737591539e+00,
					-1.893799297333427e-01,
					3.296457631788535e-03,
					-3.041742429955601e-05,
					1.164117940531552e-07
			};
			A = coeff[0];
			for (int k=1; k<coeff.length; k++){
				A += coeff[k]*T;
				T *=temperature;
			}
			
			return A;
			
		}else if(i==5){
			double[] coeff= new double[]{
					-4.641592983468033e+02,
					-3.853057077959862e+01,
					2.316202821216171e+00,
					-6.871185605401343e-02,
					1.109716262497792e-03,
					-9.362063613731818e-06,
					3.227357959446576e-08
			};
			A = coeff[0];
			for (int k=1; k<coeff.length; k++){
				A += coeff[k]*T;
				T *=temperature;
			}
			
			return A;
			
		}else if(i==6){
			double[] coeff= new double[]{
					-1.558913266582600e+03,
					1.113184285749814e+02,
					-6.188446888088212e+00,
					1.876597158215231e-01,
					-3.218200723483583e-03,
					2.942999922235004e-05,
					-1.119534617990490e-07
			};
			A = coeff[0];
			for (int k=1; k<coeff.length; k++){
				A += coeff[k]*T;
				T *=temperature;
			}
			
			return A;
			
		}else if(i==7){
			double[] coeff= new double[]{
					-5.839871208431333e+02,
					-2.812239006614151e+01,
					2.110055954006458e+00,
					-7.508344531239587e-02,
					1.450611693132853e-03,
					-1.471479568477133e-05,
					6.155436270339370e-08
			};
			A = coeff[0];
			for (int k=1; k<coeff.length; k++){
				A += coeff[k]*T;
				T *=temperature;
			}
			
			return A;
			
		}else if(i==8){
			double[] coeff= new double[]{
					6.725618360046647e+01,
					-1.168236142838188e+02,
					7.120678001507573e+00,
					-2.250575518454532e-01,
					3.959990385671562e-03,
					-3.696732556318715e-05,
					1.432489468536373e-07
			};
			A = coeff[0];
			for (int k=1; k<coeff.length; k++){
				A += coeff[k]*T;
				T *=temperature;
			}
			
			return A;
			
		}else if(i==9){
			double[] coeff= new double[]{
					-7.817483814763191e+02,
					2.105113995984065e+00,
					2.111096409893683e-01,
					-1.161056479009565e-02,
					2.632940520678518e-04,
					-2.931618541082707e-06,
					1.308385249353614e-08
			};
			A = coeff[0];
			for (int k=1; k<coeff.length; k++){
				A += coeff[k]*T;
				T *=temperature;
			}
			
			return A;
			
		} else {
			throw new RuntimeException("");
		}
		
	}
		
	public double findDensityAlpha(double temperature){
	
	    double [] p = new double[density.length];

        for(int i=0; i<p.length; i++){
        	p[i] = AandPAlpha(density[i], temperature)[1];
//        	System.out.println("p["+i+"]: " + p[i]);
        }

        akimaSpline.setInputData(p, density);
        double[] rho =  akimaSpline.doInterpolation(new double[]{desiredPressure});
        double del = 0.0003;
        
        while (del > 1e-10) {
          p = new double[5];
          double[] newRho = new double[]{rho[0]-2*del, rho[0]-del, rho[0], rho[0]+del, rho[0]+2*del};

          for(int i=0; i<p.length; i++){
        	  p[i] = AandPAlpha(newRho[i], temperature)[1];
//        	  System.out.println("p["+i+"]: " + p[i]);
          }

          akimaSpline.setInputData(p, newRho);
          rho =  akimaSpline.doInterpolation(new double[]{desiredPressure});
//          System.out.println("rho: " + rho[0]);
//          System.exit(1);
          del *= 0.1;
        }
//        System.out.println("density: " + rho[0]);
		return rho[0];
	}
	
	public double findDensityBeta(double temperature){
		
	    double [] p = new double[density.length];

        for(int i=0; i<p.length; i++){
        	p[i] = AandPBeta(density[i], temperature)[1];
//        	System.out.println("p["+i+"]: " + p[i]);
//        	 System.out.println(density[i] +" " + p[i]);
        }
//        System.exit(1);
        akimaSpline.setInputData(p, density);
        double[] rho =  akimaSpline.doInterpolation(new double[]{desiredPressure});
        double del = 0.0003;

        while (del > 1e-10) {
          p = new double[5];
          double[] newRho = new double[]{rho[0]-2*del, rho[0]-del, rho[0], rho[0]+del, rho[0]+2*del};
       
          for(int i=0; i<p.length; i++){
        	  p[i] = AandPBeta(newRho[i], temperature)[1];
//        	  System.out.println("p["+i+"]: " + p[i]);
//          	  System.out.println(density[i] +" " + p[i]);
          	  
          }
//          System.exit(1);
          akimaSpline.setInputData(p, newRho);
          rho =  akimaSpline.doInterpolation(new double[]{desiredPressure});
          del *= 0.1;
        }
//        System.out.println(temperature + " density: " + rho[0]);
		return rho[0];
	}
	
	public double findTransTemperature(double T0){
		
		//Alpha-phase Gibbs free energy
		// working with per particle quantity) 
		// GAlpha = G/N or newAlphaA = A/N
		double T1;
		
		
//		for (int temp=49; temp<50; temp++){	
//			double AlphaRho0 = findDensityAlpha(temp);
//			double[] AandPAlpha0 = AandPAlpha(AlphaRho0, temp);
//			double AlphaG0 = AandPAlpha0[0] + AandPAlpha0[1]/AlphaRho0; 
//			
////			double BetaRho0 = findDensityBeta(temp);
////			double[] AandPBeta0 = AandPBeta(BetaRho0, temp);
////			double BetaG0 = AandPBeta0[0] + AandPBeta0[1]/BetaRho0; 
//		
////			double deltaG0 = BetaG0 - AlphaG0;
//			System.out.println(temp + " "+AandPAlpha0[0] + " " + AandPAlpha0[1] + " "+ AlphaRho0+ " "+AlphaG0);
////			System.out.println(AandPBeta0[0] + " " + AandPBeta0[1] + " "+ BetaRho0+ " "+BetaG0);
////			System.out.println(temp + " deltaG0: " + deltaG0);
//		}
//		System.exit(1);
		
		
		double AlphaRho0 = findDensityAlpha(T0);
		double[] AandPAlpha0 = AandPAlpha(AlphaRho0, T0);
		double AlphaG0 = AandPAlpha0[0] + AandPAlpha0[1]/AlphaRho0; 
		
		double BetaRho0 = findDensityBeta(T0);
		double[] AandPBeta0 = AandPBeta(BetaRho0, T0);
		double BetaG0 = AandPBeta0[0] + AandPBeta0[1]/BetaRho0; 
	
		double deltaG0 = BetaG0 - AlphaG0;
//		System.out.println(AandPAlpha0[0] + " " + AandPAlpha0[1] + " "+ AlphaRho0+ " "+AlphaG0);
//		System.out.println(AandPBeta0[0] + " " + AandPBeta0[1] + " "+ BetaRho0+ " "+BetaG0);
//		System.out.println("deltaG0: " + deltaG0);
//		System.exit(1);
		
		if(deltaG0 >=0.0){
			T1 = 47;
			
		} else {
			T1 = 36;
		}
//		System.out.println("T1: " + T1);
		double AlphaRho1 = findDensityAlpha(T1);
		double[] AandPAlpha1 = AandPAlpha(AlphaRho1, T1);
		double AlphaG1 = AandPAlpha1[0] + AandPAlpha1[1]/AlphaRho1; 
		
		double BetaRho1 = findDensityBeta(T1);
		double[] AandPBeta1 = AandPBeta(BetaRho1, T1);
		double BetaG1 = AandPBeta1[0] + AandPBeta1[1]/BetaRho1; 
		
		
//		System.exit(1);
		double deltaG1 = BetaG1 - AlphaG1;
//		System.out.println(AandPAlpha1[0] + " " + AandPAlpha1[1] + " "+ AlphaRho1+ " "+AlphaG1);
//		System.out.println(AandPBeta1[0] + " " + AandPBeta1[1] + " "+ BetaRho1+ " "+BetaG1);
//		System.out.println("deltaG1: " + deltaG1);
		
		double T=Double.NaN;
		double AlphaRhoNew = Double.NaN;
		double BetaRhoNew = Double.NaN;
		double deltaGNew = 1000;
		int i =0;
		while (Math.abs(deltaGNew) > tolerance){
			
			// Interpolation to determine T at deltaG = 0.0
			T = T0 + ((T1-T0)/(deltaG1-deltaG0))*(0.0 - deltaG0);
		
			AlphaRhoNew = findDensityAlpha(T);
			double[] AandPAlphaNew = AandPAlpha(AlphaRhoNew, T);
			double AlphaGNew = AandPAlphaNew[0] + AandPAlphaNew[1]/AlphaRhoNew; 
			
			BetaRhoNew = findDensityBeta(T);
			double[] AandPBetaNew = AandPBeta(BetaRhoNew, T);
			double BetaGNew = AandPBetaNew[0] + AandPBetaNew[1]/BetaRhoNew; 
			
			System.out.println("********  " + i + " ***********" );
			System.out.println("Temperature: " + T);
			System.out.println(AandPAlphaNew[0] + " " + AandPAlphaNew[1] + " "+ AlphaRhoNew+ " "+AlphaGNew);
			System.out.println(AandPBetaNew[0] + " " + AandPBetaNew[1] + " "+ BetaRhoNew+ " "+BetaGNew);
//			System.exit(1);
			
//			System.out.println(AandPAlphaNew[0]+" "+AandPAlphaNew[1]);
//			System.out.println(AandPBetaNew[0]+" "+AandPBetaNew[1]);
			deltaGNew = BetaGNew - AlphaGNew;

//			System.exit(1);
			if(deltaGNew > 0.0){
				deltaG0 = deltaGNew;
				T0 = T;
			} else {
				deltaG1 = deltaGNew;
				T1 = T;
			}	
			++i;
		}
		
		this.rhoAlphaTrans = AlphaRhoNew;
		this.rhoBetaTrans = BetaRhoNew;
		
//		System.out.println("rho_alpha: "+AlphaRhoNew + " ;rho_Beta: " + BetaRhoNew);
		return T;
	}
	
	public static void main (String[] args){
		
		double tolerance = 1e-8;
		double pGPa = 0.01;
		double desiredPressure = Pascal.UNIT.toSim(pGPa*1e9);
		double [] density = new double[]{0.0220, 0.0221, 0.0223, 0.0224, 0.0225, 0.0226, 0.0227, 0.0228, 0.0229, 0.0230};
		
		FindTransitionTemperature findTransTemp = new FindTransitionTemperature(density, desiredPressure, tolerance);
		System.out.println(pGPa + " GPa; Transition Temperature found(K): " + findTransTemp.findTransTemperature(33) 
				+ " rhoAlpha: " + findTransTemp.rhoAlphaTrans + " rhoBeta: " + findTransTemp.rhoBetaTrans);
		System.exit(1);
		
		double TransTemp = 45.1448;
		double TransRhoAlpha = 0.022387;
		double TransRhoBeta = 0.022137;
		
//		double TransTemp = 41.744886537646046;
//		double TransRhoAlpha = 0.02168671838160813;
//		double TransRhoBeta = 0.02149379604805435;
		
		 double uAlpha = findTransTemp.UAlpha(TransRhoAlpha, TransTemp);
		 double uBeta = findTransTemp.UBeta(TransRhoBeta, TransTemp);
		 
		 double[] aAlpha= findTransTemp.AandPAlpha(TransRhoAlpha, TransTemp);
		 double[] aBeta= findTransTemp.AandPBeta(TransRhoBeta, TransTemp);
		 
		 double vmAlpha = Constants.AVOGADRO*1e-24/TransRhoAlpha;
		 double vmBeta = Constants.AVOGADRO*1e-24/TransRhoBeta;

		 System.out.println("vmAlpha: "+vmAlpha + " ; vmBeta: " + vmBeta + " ; " + (vmBeta-vmAlpha));
		 double TS = (uBeta - uAlpha) - (aBeta[0]-aAlpha[0]);
		 System.out.println("TS: " + TS);
		 System.out.println("TS (cal/mol): " + TS*2.390057243);
		 System.out.println("S (cal/mol.K): " + TS*2.390057243/TransTemp);
		 
	}
	
	protected AkimaSpline akimaSpline;
	protected double[] density; 
	protected double desiredPressure;
	protected double tolerance;
	protected double rhoAlphaTrans, rhoBetaTrans;
}
