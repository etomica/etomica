package etomica.threads;

/**
 * JAVA THREADS EXAMPLE 2
 * 
 * SyncExample -- Two threads are writing to the same variables.  Have no fear, the are SYNCHRONIZED.
 * 
 * 	Using synchronized blocks allows you to perform a group of related updates as a set without
 * worrying about other threads interrupting or seeing the intermediate results of a computation.  The 
 * following example code with either print "1 0" or "0 1".  In the absence of synchronization, it 
 * could also print "1 1" (or even "0 0", believe it or not).
 * 
 * @author msellers
 * 
 * Example taken from "Introduction to Java threads", IBM developerWorks. http://ibm.com/developerWorks
 */

public class SyncExample {
	private static Object lockObject = new Object();
	
	private static class Thread1 extends Thread {
		public void run() {
			synchronized (lockObject) {
				x = y = 0;
				System.out.println(x);
			}
		}
	}
	
	private static class Thread2 extends Thread {
		public void run() {
			synchronized (lockObject) {
				x = y = 1;
				System.out.println(y);
			}
		}
	}
	
	public static void main(String[] args) {
		new Thread1().start();
		new Thread2().start();
	}

}
