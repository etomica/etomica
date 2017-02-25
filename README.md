![etomica](http://www.etomica.org/images/etomicanew.jpg)

# Instructions for Building

## Using Eclipse

### Option 1: Do everything in Eclipse

This requires Eclipse to be set up with Git and Gradle, which should be the
case if you are using "Eclipse IDE for Java Developers". (TODO: explain how to install
plugins if they didn't come with eclipse) 

* Click on "Checkout projects from Git" on the welcome screen, or go to
`File > Import`, then `Git > Projects from Git` and click Next

* Select "Clone URI" and click Next

* Enter `https://github.com/etomica/etomica` in the `URI` field and click Next without changing anything else.

* Make sure all branches are selected and click Next.

* Use the default destination as long as you haven't downloaded etomica to there already, and click Next.

* Wait for it to download, then select **_Import as general project_** and click Next, then Finish.

* Find the etomica folder in the Package Explorer, right click and go to `Configure > Add Gradle nature`. Wait while
Eclipse compiles and configures the project.

* You should now see each of the etomica subprojects listed in the Package or Project explorer. You can run all the Java classes from
Eclipse normally.

### Option 2: Use the command line and Eclipse

* From the command line, run `git clone https://github.com/etomica/etomica`

* `cd etomica` and `./gradlew eclipse`, which will download the gradle wrapper, compile the project, and generate the
    files for Eclipse.

* Open eclipse and go to `File > Open Projects from File System`, then select the etomica directory you cloned using the
    "Directory..." button or by typing in the path.
    
* Make sure "Search for nested projects" and "Detect and configure project natures" checkboxes are selected, and that
    all the entries in the Folder list are selected.
    * TODO: explain errors caused by conflicting project names
    
* Click Finish, the projects should now be ready to use.
    

## Using Gradle and the command line without an IDE

* `git clone https://github.com/etomica/etomica && cd etomica`

* `./gradlew build` to download the gradle wrapper, download the project dependencies, and build the project.

* `./gradlew tasks` to see which tasks are available.

* TODO: add more interesting tasks.

## Using Intellij IDEA

* Clone the project on the command line, then go to `File > Open...` and select it.
    * Alternatively, go to `File > New > Project from Version Control > git` and use `https:/github.com/etomica/etomica`.
        If you use this option, you will need to click "Import Gradle project" in the notification popup before proceeding.

* At this point you should be in the "Import Project from Gradle" dialog. Make sure 
"Create separate module per source set" and "Use default gradle wrapper" are selected (which should be the defaults).

* Make sure all of the modules are checked in the next dialog to appear and click Ok.
    

    

