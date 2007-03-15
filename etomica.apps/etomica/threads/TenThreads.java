package etomica.threads;

/**
 * JAVA THREADS EXAMPLE 3
 * 
 * TenThreads -- split up a matrix among ten threads, threads find local max, then main finds global max.
 * 		
 * 		This example shows a program that creates ten threads, each of which do some work.  It waits
 * for them all to finish, then gathers the results.  Notice how we can create multiple threads with
 * a FOR loop.
 * 
 * @author msellers
 * 
 * Example taken from "Introduction to Java threads", IBM developerWorks. http://ibm.com/developerWorks
 */

public class TenThreads {
		// Pay no attention to the value behind the method.
		public static int[][] getBigHairyMatrix() {
			return null;
		}
	
		private static class WorkerThread extends Thread {
			int max = Integer.MIN_VALUE;
			int[] ourArray;
			
			public WorkerThread(int[] ourArray) {
				this.ourArray = ourArray;
			}
			
			// Find the max value in our particular piece of the array
			public void run() {
				for (int i = 0; i<ourArray.length; i++)
					max = Math.max(max, ourArray[i]);
			}
			
			public int getMax() {
				return max;
			}
		}
		
		public static void main (String[] args) {
			WorkerThread[] threads = new WorkerThread[10];
			int[][] bigMatrix = getBigHairyMatrix();
			int max = Integer.MIN_VALUE;
			
		
			
			
			
			// Give each thread a slice of the matrix to work with
			for (int i=0; i<10; i++) {
				threads[i] = new WorkerThread(bigMatrix[i]);
				threads[i].start();
			}
			
			// Wait for each thread to finish
			try {
				for (int i=0; i<10; i++) {
					threads[i].join();
					max = Math.max(max, threads[i].getMax());
				}
			}
			catch (InterruptedException e) {
				// fall through
			}
			
			System.out.println("Maximum value was " + max);
		}
}
