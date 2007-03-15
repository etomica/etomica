package etomica.threads;

/**
 * JAVA THREADS EXAMPLE 1
 * 
 * CalculatePrimes -- calculate as many primes as we can in ten seconds
 * 		
 * 		This example uses two threads, one for timing and one to do actual work.  The main thread
 * calculates prime number using a very straightforward algorithm.  Before it starts, it creates and
 * starts a timer thread, which will sleep for ten seconds, and then set a flag that the main thread
 * will check.  After ten second, the main thread wil stop.  Note that the shared flag is declared
 * volatile.
 * @author msellers
 * 
 * Example taken from "Introduction to Java threads", IBM developerWorks. http://ibm.com/developerWorks
 */

public class CalculatePrimes extends Thread {
	
	public static final int MAX_PRIMES = 1000000;
	public static final int TEN_SECONDS = 10000;
	
	public volatile boolean finished = false;
	
	public void run() {
		//array to store prime numbers
		int[] primes = new int[MAX_PRIMES];
		int count = 0;
		
		for (int i=2; count<MAX_PRIMES; i++) {
			// Check to see if the timer has expired
			if (finished) {
				break;
			}
			
			boolean prime = true;
			for (int j=0; j<count; j++) {
				if (i % primes[j] == 0) {
					prime = false;
					break;
				}
			}
			
			if (prime) {
				primes[count++] = i;
				System.out.println("Found prime: " + i);
			}
		}
	}
	
	//main thread spawns and starts calculator thread. sleeps for 10sec (timer). sets finished = true, to break and stop calculator thread.
	public static void main(String[] args) {
		CalculatePrimes calculator = new CalculatePrimes();
		calculator.start();
		try{
			Thread.sleep(TEN_SECONDS);
		}
		catch (InterruptedException e) {
			//fall through
		}
		calculator.finished = true;
	}
}
