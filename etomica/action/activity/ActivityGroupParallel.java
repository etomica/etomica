/*
 * History
 * Created on Oct 28, 2004 by kofke
 */
package etomica.action.activity;

import java.util.LinkedList;

import etomica.Action;
import etomica.Activity;
import etomica.utility.Arrays;

/**
 * Organizer of simulation actions to be executed in parallel,
 * each on its own thread.
 */
public class ActivityGroupParallel extends Activity {

	/**
	 *  
	 */
	public ActivityGroupParallel() {
		this(new Action[0]);
	}

	public ActivityGroupParallel(Action[] actions) {
		super();
		setActions(actions);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see etomica.Activity#run()
	 */
	public void run() {
		if (numActions == 0)
			return;
		synchronized (this) {
			firstThread = new MyThread(actions[0]);
			lastThread = firstThread;
			firstThread.start();
			for (int i = 1; i < numActions; i++) {
				lastThread.nextThread = new MyThread(actions[i]);
				lastThread = lastThread.nextThread;
				lastThread.start();
			}
		}
		MyThread thread;
		synchronized (this) {
			thread = firstThread;
		}
		while (true) {
			synchronized (this) {
				if (thread == null) {
					lastThread = null;
					for (thread = firstThread; thread != null; thread = thread.nextThread) {
						if (!thread.completedNormally) {
							throw new etomica.exception.AbnormalCompletionException();
						}
					}
					return;
				}
			}
			try {
				thread.join();
			} catch (InterruptedException ex) {
			}
			synchronized (this) {
				thread = thread.nextThread;
			}
		}
	}

	public synchronized void halt() {
		for (int i=0; i<numActions; i++) {
			if (actions[i] instanceof Activity) ((Activity)actions[i]).halt();
		}
	}

	public synchronized void unPause() {
		for (int i=0; i<numActions; i++) {
			if (actions[i] instanceof Activity) ((Activity)actions[i]).unPause();
		}
	}

	public synchronized void pause() {
		for (int i=0; i<numActions; i++) {
			if (actions[i] instanceof Activity) ((Activity)actions[i]).pause();
		}
	}

	public synchronized Action[] getActions() {
		return actions;
	}

	public synchronized Action[] getCompletedActions() {
		if (firstThread == null) return new Action[0];
		LinkedList completedActions = new LinkedList();
		for (MyThread thread = firstThread; thread != null; thread = thread.nextThread) {
			if (!thread.isAlive()) {
				completedActions.add(thread.action);
			}
		}
		return (Action[])completedActions.toArray();
	}
	
	public synchronized void setActions(Action[] actions) {
		if (firstThread != null)
			return;
		this.actions = actions;
		numActions = actions.length;
	}

	public synchronized void addAction(Action newAction) {
		if (lastThread != null) {
			lastThread.nextThread = new MyThread(newAction);
			lastThread = lastThread.nextThread;
			lastThread.start();
		}
		else if (firstThread != null) {
			// If firstThread is not null, it means actionPerformed was called
			// and if lastThread is null it means actionPerformed finished.
			return;
		}
		actions = (Action[]) Arrays.addObject(actions, newAction);
		numActions++;
	}

	public synchronized boolean removeAction(Action action) {
		if (firstThread != null)
			return false;
		actions = (Action[]) Arrays.removeObject(actions, action);
		int newNumActions = actions.length;
		if (newNumActions == numActions)
			return false;
		numActions = newNumActions;
		return true;
	}

	protected int numActions;

	protected Action[] actions;

	protected MyThread lastThread, firstThread;

	private class MyThread extends Thread {
		public MyThread(Action action) {
			this.action = action;
		}

		public void run() {
			completedNormally = true;
			try {
				action.actionPerformed();
			} catch (Exception e) {
				e.printStackTrace();
				completedNormally = false;
			}
		}

		private Action action;

		private MyThread nextThread;

		private boolean completedNormally;
	}
}