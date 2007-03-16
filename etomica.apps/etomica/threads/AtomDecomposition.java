package etomica.threads;

/**
 * 
 * @author msellers
 *
 */


public class AtomDecomposition {

		
		private static class IntegratorThread extends Thread {
			
			double deltat = 0.001;
			int box;
			int threadID;
			int Tatoms;
			int atoms;
			double[][] Tcoord;
			double[][] Tvelocity;
			double[][] Tforce;
			
			public IntegratorThread(int t, int n, int z, int b, double[][] c, double[][] v, double[][] f){
				this.Tcoord = c;
				this.Tvelocity = v;
				this.Tforce = f;
				this.threadID = t;
				this.Tatoms = n;
				this.atoms = z;
				this.box = b;
			}
			
			public void run() {
				
				double[] r = new double[3];
				
				int startAtom = threadID*Tatoms;
				int stopAtom = (threadID+1)*Tatoms;
				
				//Compute positions and half update velocities
				for (int i=startAtom; i<stopAtom; i++){
					for (int j=0; j<3; j++){
						Tcoord[i][j] += deltat*Tvelocity[i][j] + deltat*deltat*Tforce[i][j] * 0.5;
						Tvelocity[i][j] = Tvelocity[i][j] + deltat*Tforce[i][j] * 0.5;
					}
				}
				
				//Zero old forces, compute new
				for (int i=startAtom; i<stopAtom; i++){
					for (int j=0; j<atoms; j++){
						
						if (j-i==0){ 
							continue;
						}
						
						double r2 = 0.0;
						
						for(int l=0; l<3; l++){
							Tforce[i][l] = 0.0;
							
							r[l] = Tcoord[i][l] - Tcoord[j][l];
							r[l] = r[l] - box*Math.round(r[l]/box);
							r2 += r[l]*r[l];
						}
						
							//Force loop, w/ cuttoff
							if (r2 < 6.25){
								double r2i = 1.0/r2;
								double r6i = r2i*r2i*r2i;
								double ff = 48*r2i*r6i*(r6i-0.5);
							
								for(int l=0; l<3; l++) {
									Tforce[i][l] += ff * r[l];
								}
							}
					
					}
				}
				
				for (int i=startAtom; i<stopAtom; i++){
					for(int l=0; l<3; l++){
						
						Tvelocity[i][l] += deltat*Tforce[i][l]*0.5;
					}
				}
				
			}
			
			
			public double[][] returnCoord(){
				return Tcoord;
			}
			
			public double[][] returnVelocity(){
				return Tvelocity;
			}
			
			public double[][] returnForce(){
				return Tforce;
			}
			
		}
		
		public static void main(String [] args) {
			
//			User variables
			int numThreads = 4;
			int boxSide = 20;
			int numAtoms = (boxSide/2) * (boxSide/2) * (boxSide/2);
			
			int atomsPerThread = numAtoms / numThreads;
			
			IntegratorThread[] threads = new IntegratorThread[numThreads]; 
			
			// Create position, velocity, and force matrix
			double[][] coord = new double[numAtoms][3];
			double[][] velocity = new double[numAtoms][3];
			double[][] force = new double[numAtoms][3];
			double[] velsum = new double[3];
			
			int n = 0;
			for (int k=0; k<(boxSide/2); k++) {
				for (int l=0; l<(boxSide/2); l++) {
					for (int m=0; m<(boxSide/2); m++) {
						coord[n][0] = 2*k;
						coord[n][1] = 2*l;
						coord[n][2] = 2*m;							
					
						n++;
					}
				}	
			}
			
			velsum[0] = 0.0;
			velsum[1] = 0.0;
			velsum[2] = 0.0;
			
			for (int i=0; i<numAtoms; i++) {
				
				velocity[i][0] = 2*Math.random() - 1;
				velocity[i][1] = 2*Math.random() - 1;
				velocity[i][2] = 2*Math.random() - 1;
				
				velsum[0] += velocity[i][0];
				velsum[1] += velocity[i][1];
				velsum[2] += velocity[i][2];
				
				force[i][0] = 0.0;
				force[i][1] = 0.0;
				force[i][2] = 0.0;
			}
			
			for (int i=0; i<numAtoms; i++) {
				velocity[i][0] -= velsum[0];
				velocity[i][1] -= velsum[1];
				velocity[i][2] -= velsum[2];
			}
			
			//Main Simulation, dishing atoms to threads.
			int nsteps = 10000;
			int istep = 0;
			double[][] coordPass = new double[atomsPerThread][3];
			double[][] velocityPass = new double[atomsPerThread][3];
			double[][] forcePass = new double[atomsPerThread][3];
			
			for (istep=0; istep<nsteps; istep++) {
								
				// Give each thread a slice of the coordinate matrix to work with
				for (int i=0; i<numThreads; i++) {				
						
					//Fill arrays to pass
					coordPass = coord;
					velocityPass = velocity;
					forcePass = force;
					
					//Create thread, pass array, start thread
					threads[i] = new IntegratorThread(i, atomsPerThread, numAtoms, boxSide, coordPass, velocityPass, forcePass);
					threads[i].start();
				}
				
				// Each thread integrates 1 timestep, and reports array back
				try {
					for (int i=0; i<numThreads; i++) {				
						
						threads[i].join();
						
						coordPass = threads[i].returnCoord();
						velocityPass = threads[i].returnVelocity();
						forcePass = threads[i].returnForce();
						
						int start = i*atomsPerThread;
						int stop = (i+1)*atomsPerThread;
						
						for (int j=start; j<stop; j++){
							coord[j] = coordPass[j];
							velocity[j] = velocityPass[j];
							force[j] = forcePass[j];
						}
						
					}
					
				}
				catch (InterruptedException e) {
					// fall through
				}
				
			System.out.println(istep + ".  Atom 0's x-position is " + coord[0][0]);	
			System.out.println(istep + ".  Atom 932's y-position is " + coord[932][1]);
			}
			
			
		}
}
