/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.action.activity;

import java.util.LinkedList;

import etomica.action.Activity;
import etomica.action.IAction;
import etomica.util.Arrays;

/**
 * Organizer of simulation actions to be executed in parallel,
 * each on its own thread.
 */
public class ActivityGroupParallel extends Activity implements ActivityGroup {

    /**
	 *  
	 */
	public ActivityGroupParallel() {
		this(new IAction[0]);
	}

	public ActivityGroupParallel(IAction[] actions) {
		super();
		setActions(actions);
	}
    
	/*
	 * (non-Javadoc)
	 * 
	 * @see etomica.Activity#run()
	 */
	protected void run() {
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

	public synchronized IAction[] getAllActions() {
		return actions.clone();
	}

	public synchronized IAction[] getCompletedActions() {
		if (firstThread == null) return new IAction[0];
		LinkedList completedActions = new LinkedList();
		for (MyThread thread = firstThread; thread != null; thread = thread.nextThread) {
			if (!thread.isAlive()) {
				completedActions.add(thread.action);
			}
		}
		return (IAction[])completedActions.toArray();
	}
	
	public synchronized IAction[] getPendingActions() {
		if (firstThread == null) return actions.clone();
		return new IAction[0];
	}
	
	public synchronized IAction[] getCurrentActions() {
		if (firstThread == null) return new IAction[0];
		LinkedList currentActions = new LinkedList();
		for (MyThread thread = firstThread; thread != null; thread = thread.nextThread) {
			if (thread.isAlive()) {
				currentActions.add(thread.action);
			}
		}
		return (IAction[])currentActions.toArray();
	}
	
	public synchronized void setActions(IAction[] actions) {
		if (firstThread != null)
			return;
		this.actions = actions;
		numActions = actions.length;
	}

	public synchronized void addAction(IAction newAction) {
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
		actions = (IAction[]) Arrays.addObject(actions, newAction);
		numActions++;
	}

	public synchronized boolean removeAction(IAction action) {
		if (firstThread != null)
			return false;
		actions = (IAction[]) Arrays.removeObject(actions, action);
		int newNumActions = actions.length;
		if (newNumActions == numActions)
			return false;
		numActions = newNumActions;
		return true;
	}

    private static final long serialVersionUID = 1L;
	protected int numActions;
	protected IAction[] actions;
	protected MyThread lastThread, firstThread;

	private static class MyThread extends Thread {
		public MyThread(IAction action) {
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

		protected IAction action;
		protected MyThread nextThread;
		protected boolean completedNormally;
	}
}
