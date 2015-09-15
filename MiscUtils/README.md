# Description

Classes and methods to convenietly handle bioinformatics data. See individual classes for usage.

# Installation

You have created a java project where you want to use the classes in **MiscUtils**. The steps below assume you are using Eclipse.

## Option 1: Using compiled jar

Download the **MiscUtils** jar from [realeses](https://github.com/dariober/Java-cafe/releases). Configure you project's build path to include this jar 

Right-click your Project `->` Configure Build Path `->` select Libraries tab `->` Add external jars `->` Select `MiscUtils.x.x.x.jar`

## Option 2: From source

* Check out this repository *This retrieves the whole lot, not just MiscUtils but that's fine*
```
svn co https://github.com/dariober/Java-cafe/
```

* In the **src** directory of your Java project create a package called **miscUtils**.
* Using the GUI Finder or Explorer navigate to the **src** directory inside donwloaded **MiscUtils**. It might be somewehere like `~/.../Java-cafe/trunk/MiscUtils/src`
* Drag and drop some or all the content of `src` into the package **miscUtils**. The main dir in `MiscUtils/src` can (should) be skipped.
* Configure the project build path to add possible libraries required by MiscUtils. 
For example you might need [htsjdk-x.xxx.jar](https://github.com/broadinstitute/picard/releases/). 
So download the required jar files. Then right-click on the project `->` Configure Build Path `->` Libraries tab `->` Add external jars `->` Select `htsjdk-x.xxx.jar`

Now you should be ready to go.


