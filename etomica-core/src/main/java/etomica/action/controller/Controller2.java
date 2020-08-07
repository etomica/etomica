package etomica.action.controller;

import etomica.action.IAction;
import etomica.action.activity.Activity2;
import etomica.action.activity.ActivityIntegrate2;
import etomica.integrator.Integrator;
import etomica.util.Debug;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.*;
import java.util.concurrent.atomic.AtomicBoolean;

public class Controller2 {
    private final ExecutorService executor = Executors.newSingleThreadExecutor();
    private ActivityTask currentTask;
    private final ArrayDeque<ActivityHandle> pendingActivities = new ArrayDeque<>(1);
    private final LinkedBlockingQueue<IAction> actionQueue = new LinkedBlockingQueue<>();
    private boolean started = false;
    private final Object lock = new Object();

    private final AtomicBoolean pauseFlag = new AtomicBoolean(true);
    private static final IAction WAKE_UP = () -> {};
    private Integrator integrator;


    public Controller2() {
    }

    /**
     * Adds an Activity to be executed by the Controller on a separate thread. Activities will be run in the
     * order they are added, with the next one starting after the current one completes. The returned
     * `CompletableFuture` will be completed when the activity finishes.
     * @param activity the repeated action.
     * @param maxSteps
     * @param sleepPeriodMillis
     */
    public ActivityHandle addActivity(Activity2 activity, long maxSteps, double sleepPeriodMillis) {
        this.tryToGetIntegrator(activity);
        ActivityTask task = new ActivityTask(activity, this.actionQueue, this.pauseFlag);
        if (this.currentTask == null) {
            this.currentTask = task;
        }
        task.maxSteps = maxSteps;
        task.sleepPeriodMillis = sleepPeriodMillis;
        CompletableFuture<Void> future = new CompletableFuture<>();
        ActivityHandle activityHandle = new ActivityHandle(future, task);
        synchronized (this.lock) {
            this.pendingActivities.add(activityHandle);
            this.processActivities();
        }
        return activityHandle;
    }

    public ActivityHandle addActivity(Activity2 activity) {
        return this.addActivity(activity, Long.MAX_VALUE, 0.0);
    }

    /**
     * Runs an activity to completion on the current thread.
     * @param activity
     * @param maxSteps
     */
    public void runActivityBlocking(Activity2 activity, long maxSteps) {
        ActivityTask task = new ActivityTask(activity, this.actionQueue, this.pauseFlag);
        task.maxSteps = maxSteps;
        task.run();
    }

    private void tryToGetIntegrator(Activity2 activity) {
        if (activity instanceof ActivityIntegrate2) {
            this.integrator = ((ActivityIntegrate2) activity).getIntegrator();
        }
    }

    public Integrator getIntegrator() {
        return integrator;
    }

    public void start() {
        synchronized (this.lock) {
            this.started = true;
            this.processActivities();
        }
    }

    private void processActivities() {
        if (this.started) {
            for (ActivityHandle activity : pendingActivities) {
                this.executor.submit(() -> {
                    this.currentTask = activity.task;
                    try {
                        activity.task.run();
                        activity.future.complete(null);
                    } catch (Throwable e) {
                        e.printStackTrace();
                        activity.future.completeExceptionally(e);
                    }
                });
            }
            this.pendingActivities.clear();
        }

    }

    public CompletableFuture<Void> pause() {
        this.pauseFlag.set(true);
        return submitActionInterrupt(WAKE_UP);
    }

    public CompletableFuture<Void> unpause() {
        this.pauseFlag.set(false);
        return submitActionInterrupt(WAKE_UP);
    }

    public CompletableFuture<Void> toggle() {
        this.pauseFlag.set(!this.pauseFlag.get());
        return submitActionInterrupt(WAKE_UP);
    }

    public CompletableFuture<Void> restartCurrentActivity() {
        return submitActionInterrupt(() -> {
            this.currentTask.activity.restart();
        });
    }

    public void completeActivities() {
        this.start();
        this.unpause();
        this.executor.shutdown();
        try {
            this.executor.awaitTermination(Long.MAX_VALUE, TimeUnit.SECONDS);
        } catch (InterruptedException e) {
            throw new RuntimeException(e);
        }
    }

    public boolean isPaused() {
        return this.pauseFlag.get();
    }

    public CompletableFuture<Void> submitActionInterrupt(IAction action) {
        CompletableFuture<Void> future = new CompletableFuture<>();
        if (this.currentTask == null) {
            this.executor.submit(() -> {
                action.actionPerformed();
                future.complete(null);
            });
            return future;
        }
        this.actionQueue.add(() -> {
            if (future.isCancelled()) {
                System.out.println("canceled");
                return;
            }
            try {
                action.actionPerformed();
                future.complete(null);
            } catch (Throwable ex) {
                future.completeExceptionally(ex);
            }
            System.out.println(action);
        });
        return future;
    }

    public CompletableFuture<Void> addActionSequential(IAction action) {
        Activity2 activity = new Activity2() {
            @Override
            public void preAction() {

            }

            @Override
            public void postAction() {

            }

            @Override
            public void restart() {

            }

            @Override
            public void actionPerformed() {
                action.actionPerformed();
            }
        };
        ActivityHandle handle = this.addActivity(activity, 1, 0);
        return handle.future;
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

    private static class ActivityTask implements Runnable {
        final Activity2 activity;
        long maxSteps;
        double sleepPeriodMillis;
        final LinkedBlockingQueue<IAction> actionQueue;

        private double sleepCarryover = 0.0;
        private final List<IAction> actionsToRun = new ArrayList<>();
        private final AtomicBoolean pauseFlag;

        ActivityTask(Activity2 activity, LinkedBlockingQueue<IAction> actionQueue, AtomicBoolean pauseFlag) {
            this.activity = activity;
            this.pauseFlag = pauseFlag;
            this.sleepPeriodMillis = 0.0;
            this.maxSteps = Long.MAX_VALUE;
            this.actionQueue = actionQueue;
        }

        @Override
        public void run() {
            activity.preAction();

            for (long currentStep = 0; currentStep < maxSteps; ) {
                if (pauseFlag.get()) {
                    // If the activity is paused, just block waiting for Action events in the queue but don't invoke
                    // the main activity. Since the only way to unpause the activity is to send an "unpause" action,
                    // this effectively suspends the thread until it is resumed while still handling action events.
                    try {
                        while (!actionQueue.isEmpty()) {
                            IAction action = actionQueue.take();
                            action.actionPerformed();
                        }
                    } catch (InterruptedException e) {
                        throw new RuntimeException(e);
                    }

                } else {
                    // When the activity is not paused, run all the queued actions but don't wait if there are none.
                    while (!actionQueue.isEmpty()) {
                        actionQueue.drainTo(this.actionsToRun);
                        this.actionsToRun.forEach(IAction::actionPerformed);
                        this.actionsToRun.clear();
                    }
                }

                // An action event may have paused the activity so check again.
                if (!pauseFlag.get()) {
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

    public static class ActivityHandle {
        public final CompletableFuture<Void> future;
        private final ActivityTask task;

        private ActivityHandle(CompletableFuture<Void> future, ActivityTask task) {
            this.future = future;
            this.task = task;
        }

        public ActivityHandle setMaxSteps(long maxSteps) {
            this.task.maxSteps = maxSteps;
            return this;
        }

        public long getMaxSteps() {
            return this.task.maxSteps;
        }

        public ActivityHandle setSleepPeriod(double sleepPeriodMillis) {
            this.task.sleepPeriodMillis = sleepPeriodMillis;
            return this;
        }

        public double getSleepPeriod() {
            return this.task.sleepPeriodMillis;
        }

        public void blockUntilComplete() {
            try {
                this.future.get();
            } catch (InterruptedException | ExecutionException e) {
                throw new RuntimeException(e);
            }
        }
    }

//    public static class ControllerFuture<T> extends CompletableFuture<T> {
//        @Override
//        public <U> CompletableFuture<U> newIncompleteFuture() {
//            return super.newIncompleteFuture();
//        }
//    }
}
