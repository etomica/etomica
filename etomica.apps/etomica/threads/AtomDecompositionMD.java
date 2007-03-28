package etomica.threads;

public class AtomDecompositionMD {
	
	//	 User variables
	int numThreads = 1;
	int boxSide = 20;
	double deltat = 0.001;
	int numAtoms = (boxSide/2) * (boxSide/2) * (boxSide/2);
	
	int atomsPerThread = numAtoms / numThreads;
	int ready1 = 0;
	int ready2 = 0;
	double[][] coord = new double[numAtoms][3];
	double[][] velocity = new double[numAtoms][3];
	double[][] force = new double[numAtoms][3];
	double[] velsum = new double[3];
	
	ThreadsMD[] threads = new ThreadsMD[numThreads]; 
	
	public void initializeSim(){
		
		//	 Create position, velocity, and force matrix
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
	}
	
	public void createThreads(){
		// Create threads and pointers
		for (int i=0; i<numThreads; i++) {				
			threads[i] = new ThreadsMD(i, atomsPerThread, numAtoms, boxSide, deltat, coord, velocity, force, this);
			threads[i].start();
		}
	}
	
	public void runMD(){

		int nsteps = 10000;
		int istep = 0;
		
		initializeSim();
		createThreads();
		
		for (istep=0; istep<nsteps; istep++) {
			// VERLET INTEGRATION
			synchronized(this){
				ready1++;
				notifyAll();
			}
			
			for(int i=0; i<numThreads; i++){
				synchronized(threads[i]){
					try{
						if (threads[i].stepcounter1<ready1){
							threads[i].wait();	
						}
					}
					catch(InterruptedException e){
					}
			
				}
			}
			
			// THREAD WRITE
			synchronized(this){
				ready2++;
				notifyAll();
			}
			
			for(int i=0; i<numThreads; i++){
				synchronized(threads[i]){
					try{
						if (threads[i].stepcounter2<ready2){
							threads[i].wait();	
						}
					}
					catch(InterruptedException e){
					}
			
				}
			}
		
		
		}
		
		synchronized(this){
			for(int i=0; i<numThreads; i++){
				threads[i].run = false;
			}
			notifyAll();
		}
		
	
		
	}
	
	public static void main(String [] args) {
		
		double t1 = System.currentTimeMillis();
	
		AtomDecompositionMD AD = new AtomDecompositionMD();
		
		AD.runMD();
		
		double t2 = System.currentTimeMillis();
		t2 =  (t2-t1)/1000;
		
		System.out.println(t2 +" seconds.");
	}

	
}
