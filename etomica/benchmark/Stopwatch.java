package etomica.benchmark;

/**
 * A class to help benchmark code.
 * It simulates a real stop watch.
 * 
 * Taken from "Java Platform Performance", by S. Wilson and J. Kesselman.
 */
 
public class Stopwatch {

    private long startTime = -1;
    private long stopTime = -1;
    private boolean running = false;
    
    /**
     * Starts the watch by recording the current time.
     * If running already, performs no action.
     */
    public Stopwatch start() {
        if(running) return this;
        startTime = System.currentTimeMillis();
        running = true;
        return this;
    }
    
    /**
     * Stops the watch by recording the current time.
     * Performs no action if not running.
     */
    public Stopwatch stop() {
        if(!running) return this;
        stopTime = System.currentTimeMillis();
        running = false;
        return this;
    }
    /**
     * Returns the elapsed time (in milliseconds) between start and stop method calls,
     * of between start and current time, if stop was not yet called.
     * If the watch has never been started, then return zero.
     */
     public long getElapsedTime() {
        if(startTime == -1) {
            return 0;
        }
        if(running) {
            return System.currentTimeMillis() - startTime;
        } else {
            return stopTime - startTime;
        }
     }
     
     /**
      * Sets the watch to its initial state.
      */
     public Stopwatch reset() {
        startTime = -1;
        stopTime = -1;
        running = false;
        return this;
     }
     
     /**
      * Demonstrates how this class is used.
      */
     public static void main(String[] args) {
        
        Stopwatch timer = new Stopwatch().start();
        
        long total = 0;
        for(int i=0; i < 10000000; i++) {
            total += i;
        }
        timer.stop();
        System.out.println("Uploop elapsed time (milliseconds): "+timer.getElapsedTime());

        timer.reset();
        timer.start();
        total = 0;
        for(int i=10000000-1; i>=0; i--) {
            total += i;
        }
        timer.stop();
        System.out.println("Dnloop elapsed time (milliseconds): "+timer.getElapsedTime());
        
     }
}