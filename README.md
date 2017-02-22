![etomica](http://www.etomica.org/images/etomicanew.jpg)

# Instructions for Building

## Using Eclipse

### Option 1: Do everything in Eclipse

This requires Eclipse to be set up with Git and Gradle, which should be the
case if you are using "Eclipse IDE for Java Developers". (TODO: explain how to install
plugins if they didn't come with eclipse) 

1. Click on "Checkout projects from Git" on the welcome screen, or go to
`File > Import`, then `Git > Projects from Git` and click Next

2. Select "Clone URI" and click Next

3. Enter `https://github.com/etomica/etomica` in the `URI` field and click Next without changing anything else.

4. Make sure all branches are selected and click Next.

5. Use the default destination as long as you haven't downloaded etomica to there already, and click Next.

6. Wait for it to download, then select **_Import as general project_** and click Next, then Finish.

7. Find the etomica folder in the Package Explorer, right click and go to `Configure > Add Gradle nature`. Wait while
Eclipse compiles and configures the project.

8. You should now see each of the etomica subprojects listed in the Package or Project explorer. You can run all the Java classes from
Eclipse normally.

### Option 2: Use the command line and Eclipse

TODO

## Using Gradle and the command line without an IDE

TODO

## Using Intellij IDEA

TODO
