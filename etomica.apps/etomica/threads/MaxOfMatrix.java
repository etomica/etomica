package etomica.threads;

/**
 * 		
 * 		The MAIN method gets a large matrix, and creates 4 threads to operate on
 * portions of that matrix.  The threads are started, and the MAIN method now tries
 * to 'join' the threads.  Once the threads are joined, another way of being sure
 * they have finished, the maximum value from the threads array is returned.  Each
 * return is compared against the global maximum.  When all the threads have returned
 * their value, the global max is printed.
 *  
 * @author msellers
 * 
 */

public class MaxOfMatrix {
	
	public static int[][] getBigMatrix(){
		
		int[][] matrix = new int[4][20];
		int filler = 37;
		
		for (int i=0; i<4; i++){
			for (int j=0; j<20; j++){
				
				matrix[i][j] = filler;
				filler++;
			}
		}
				
		return matrix;
	}
	
	private static class WorkerThread extends Thread {
		int max = Integer.MIN_VALUE;
		int[] threadArray;
		
		// A method to take the array from MAIN and let it be used in WorderThread
		public WorkerThread(int[] threadArray) {
			this.threadArray = threadArray;
		}
		
		// Find the maximum value in the thread's own array
		public void run() {
			for (int i = 0; i < threadArray.length; i++)
				max = Math.max(max, threadArray[i]);
		}
		
		public int getMax() {
			return max;
		}
	}
	
	public static void main(String[] args) {
		WorkerThread[] threads = new WorkerThread[4];
		int[][] bigMatrix = getBigMatrix();
		int max = Integer.MIN_VALUE;
		
		//Give each thread a row of the matrix to work on
		for (int i=0; i<4; i++) {
			threads[i] = new WorkerThread(bigMatrix[i]);
			threads[i].start();
		}
		
		//Wait for each thread to finish
		try {
			for (int i=0; i < 4; i++) {
				threads[i].join();
				System.out.println("Thread maximum was: " + threads[i].getMax());
				max = Math.max(max, threads[i].getMax());
			}
		}
		catch (InterruptedException e) {
			//fall through
		}
		
		System.out.println("Global maximum value was " + max);
	}

}
