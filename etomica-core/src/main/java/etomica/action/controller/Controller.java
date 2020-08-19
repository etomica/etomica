package etomica.action.controller;

import etomica.action.IAction;
import etomica.action.activity.ActivityIntegrate;
import etomica.integrator.Integrator;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.*;
import java.util.concurrent.atomic.AtomicBoolean;

public class Controller {
    private final ExecutorService executor = Executors.newSingleThreadExecutor();
    private final ArrayDeque<ActivityHandle<?>> pendingActivities = new ArrayDeque<>(1);
    private boolean started = false;
    private final Object lock = new Object();

    private static final IAction WAKE_UP = () -> {};
    private Integrator integrator;

    private final ControllerHandle handle = new ControllerHandle();
    private ActivityHandle<?> currentTask;


    public Controller() {
    }

    /**
     * Adds an Activity to be executed by the Controller on a separate thread. Activities will be run in the
     * order they are added, with the next one starting after the current one completes. The returned
     * `CompletableFuture` will be completed when the activity finishes.
     * @param activity the repeated action.
     * @param maxSteps
     */

    public ActivityHandle<ActivityIntegrate> addActivity(ActivityIntegrate activity, long maxSteps) {
        activity.setMaxSteps(maxSteps);
        return addActivity(activity);
    }

    public ActivityHandle<ActivityIntegrate> addActivity(ActivityIntegrate activity, long maxSteps, double sleepPeriodMillis) {
        activity.setMaxSteps(maxSteps);
        this.handle.sleepPeriodMillis = sleepPeriodMillis;
        return addActivity(activity);
    }

    public <A extends Activity> ActivityHandle<A> addActivity(A activity) {
        tryToGetIntegrator(activity);
        CompletableFuture<Void> future = new CompletableFuture<>();
        ActivityHandle<A> handle = new ActivityHandle<>(activity, future);
        synchronized (this.lock) {
            this.pendingActivities.add(handle);
            this.processActivities();
        }
        return handle;
    }

    /**
     * Runs an activity to completion on the current thread.
     * @param activity
     */
    public void runActivityBlocking(Activity activity) {
        this.handle.pauseFlag.set(false);
        activity.runActivity(this.handle);
    }

    private void tryToGetIntegrator(Activity activity) {
        if (activity instanceof ActivityIntegrate) {
            this.integrator = ((ActivityIntegrate) activity).getIntegrator();
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
            for (ActivityHandle<?> handle : this.pendingActivities) {
                this.executor.submit(() -> {
                    this.currentTask = handle;
                    try {
                        handle.activity.runActivity(this.handle);
                        handle.future.complete(null);
                    } catch (Throwable e) {
                        e.printStackTrace();
                        handle.future.completeExceptionally(e);
                    }
                });
            }
            this.pendingActivities.clear();
        }

    }

    public CompletableFuture<Void> pause() {
        this.handle.pauseFlag.set(true);
        return submitActionInterrupt(WAKE_UP);
    }

    public CompletableFuture<Void> unpause() {
        this.handle.pauseFlag.set(false);
        return submitActionInterrupt(WAKE_UP);
    }

    public CompletableFuture<Void> toggle() {
        this.handle.pauseFlag.set(!this.handle.pauseFlag.get());
        return submitActionInterrupt(WAKE_UP);
    }

    public void halt() {
        this.start();
        this.unpause();
        this.executor.shutdown();
        try {
            this.executor.awaitTermination(Long.MAX_VALUE, TimeUnit.SECONDS);
        } catch (InterruptedException e) {
            throw new RuntimeException(e);
        }
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
        return this.handle.pauseFlag.get();
    }

    public boolean isRunningActivityStep() {
        return this.handle.isRunningActivityStep.get();
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
        this.handle.actionQueue.add(() -> {
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
        Activity activity = new Activity() {
            @Override
            public void runActivity(ControllerHandle handle) {
                handle.yield(action);
            }
        };
        ActivityHandle<?> handle = this.addActivity(activity);
        return handle.future;
    }

    public double getSleepPeriod() {
        return this.handle.sleepPeriodMillis;
    }

    public void setSleepPeriod(double sleepPeriodMillis) {
        this.handle.sleepPeriodMillis = sleepPeriodMillis;
    }

    public static class ActivityHandle<A extends Activity> {
        public final A activity;
        public final CompletableFuture<Void> future;


        public ActivityHandle(A activity, CompletableFuture<Void> future) {
            this.activity = activity;
            this.future = future;
        }
    }

    public static class ControllerHandle {
        private final LinkedBlockingQueue<IAction> actionQueue = new LinkedBlockingQueue<>();
        private final AtomicBoolean pauseFlag = new AtomicBoolean(true);
        private final AtomicBoolean isRunningActivityStep = new AtomicBoolean(false);

        private double sleepPeriodMillis = 0.0;
        private double sleepCarryover = 0.0;
        private final List<IAction> actionsToRun = new ArrayList<>();

        public final void yield(IAction action) {
            while (true) {
                if (pauseFlag.get()) {
                    // If the activity is paused, just block waiting for Action events in the queue but don't invoke
                    // the main activity. Since the only way to unpause the activity is to send an "unpause" action,
                    // this effectively suspends the thread until it is resumed while still handling action events.
                    try {
                        while (!actionQueue.isEmpty()) {
                            IAction a = actionQueue.take();
                            a.actionPerformed();
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
                    this.isRunningActivityStep.set(true);
                    action.actionPerformed();
                    this.isRunningActivityStep.set(false);

                    this.doSleepPeriod();
                    break;
                }
            }
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

        public LinkedBlockingQueue<IAction> getActionQueue() {
            return actionQueue;
        }

        public AtomicBoolean getPauseFlag() {
            return pauseFlag;
        }

        public AtomicBoolean getIsRunningActivityStep() {
            return isRunningActivityStep;
        }
    }
}
