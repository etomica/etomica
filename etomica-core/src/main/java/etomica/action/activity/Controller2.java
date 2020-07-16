package etomica.action.activity;

import etomica.action.IAction;
import etomica.util.Debug;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.locks.LockSupport;

public class Controller2 {
    private final ExecutorService executor = Executors.newSingleThreadExecutor();
    private ActivityTask currentTask;
    private final ArrayDeque<ActivityTask> pendingActivities = new ArrayDeque<>(1);

    public Controller2() {
    }

    public void setActivity(Activity2 activity, long maxSteps, double sleepPeriodMillis) {
        ActivityTask task = new ActivityTask(activity, new LinkedBlockingQueue<>());
        task.maxSteps = maxSteps;
        this.currentTask = task;
        this.currentTask.sleepPeriodMillis = sleepPeriodMillis;
    }

    public void start() {
        executor.submit(this.currentTask);
    }

    public void pause() {
        this.currentTask.actionQueue.add(this.currentTask.pauseAction);
    }


    public void submitAction(IAction action) {
        this.currentTask.submitAction(action);
    }

    public double getSleepPeriod() {
        return this.currentTask.sleepPeriodMillis;
    }

    public void setSleepPeriod(double sleepPeriodMillis) {
        this.currentTask.sleepPeriodMillis = sleepPeriodMillis;
    }

    public long getMaxSteps() {
        return this.currentTask.maxSteps;
    }

    public void setMaxSteps(long maxSteps) {
        this.currentTask.maxSteps = maxSteps;
    }

    public void addAction(IAction action) {
        if (action instanceof Activity2) {
//            this.setActivity();
        } else {
            ActivityTask task = new ActivityTask(new Activity2() {
                @Override
                public void preAction() {

                }

                @Override
                public void postAction() {

                }

                @Override
                public void actionPerformed() {
                    action.actionPerformed();
                }
            }, new LinkedBlockingQueue<>());
        }

    }


    private static class ActivityTask implements Runnable {
        final Activity2 activity;
        long maxSteps;
        double sleepPeriodMillis;
        final LinkedBlockingQueue<IAction> actionQueue;

        private double sleepCarryover = 0.0;
        private boolean pauseFlag;
        private final List<IAction> actionsToRun = new ArrayList<>();

        final IAction pauseAction = () -> pauseFlag = true;
        final IAction resumeAction = () -> pauseFlag = false;
        final IAction toggleAction = () -> pauseFlag = !pauseFlag;

        ActivityTask(Activity2 activity, LinkedBlockingQueue<IAction> actionQueue) {
            this.activity = activity;
            this.pauseFlag = false;
            this.actionQueue = actionQueue;
        }

        public void submitAction(IAction action) {
            this.actionQueue.add(action);
        }

        @Override
        public void run() {
            activity.preAction();

            for (long currentStep = 0; currentStep < maxSteps; ) {
                if (pauseFlag) {
                    // If the activity is paused, just block waiting for Action events in the queue but don't invoke
                    // the main activity. Since the only way to unpause the activity is to send an "unpause" action,
                    // this effectively suspends the thread until it is resumed while still handling action events.
                    try {
                        IAction action = actionQueue.take();
                        action.actionPerformed();
                    } catch (InterruptedException e) {
                        throw new RuntimeException(e);
                    }

                } else {
                    // When the activity is not paused, run all the queued actions but don't wait if there are none.
                    actionQueue.drainTo(this.actionsToRun);
                    this.actionsToRun.forEach(IAction::actionPerformed);
                    this.actionsToRun.clear();
                }

                // An action event may have paused the activity so check again.
                if (!pauseFlag) {
                    if (Debug.ON) {
                        if (currentStep == Debug.START) { Debug.DEBUG_NOW = true; }
                        if (currentStep == Debug.STOP) { break; }
                        if (Debug.DEBUG_NOW && Debug.LEVEL > 0) {
                            System.out.println("*** integrator step " + currentStep);
                        }
                        Debug.stepCount = currentStep;
                    }

                    activity.actionPerformed();
                    currentStep++;

                    this.doSleepPeriod();
                }
            }

            activity.postAction();
        }

        private void doSleepPeriod() {
            if (this.sleepPeriodMillis > 0) {
                double nowSleep = sleepCarryover + this.sleepPeriodMillis;
                int sleepWholeMillis = (int) nowSleep;
                sleepCarryover = nowSleep - sleepWholeMillis;
                if (sleepWholeMillis > 0) {
                    try {
                        Thread.sleep(sleepWholeMillis);
                    } catch (InterruptedException e) {
                        throw new RuntimeException(e);
                    }
                }
            }
        }
    }
}
