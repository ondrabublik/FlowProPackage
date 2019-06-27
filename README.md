# FlowProPackage
This is a root projects for [FlowPro](https://github.com/ondrabublik/FlowPro) and its modules. FlowProPackage includes [FlowPro](https://github.com/ondrabublik/FlowPro) itself, an interface [FlowProAPI](https://github.com/ondrabublik/FlowProAPI) according to which a user can create his or her modules and finally some basic modules which are implemented according to [FlowProAPI](https://github.com/ondrabublik/FlowProAPI).

FlowProPackage is not necessary for the use of 
[FlowPro](https://github.com/ondrabublik/FlowPro). 
[FlowPro](https://github.com/ondrabublik/FlowPro) can be used on its own, however FlowProPackage makes it easier to build [FlowPro](https://github.com/ondrabublik/FlowPro) along with its modules all at once.

## Built With
* [Gradle](https://gradle.org/)

## Prerequisites
* [Git](https://git-scm.com/)  
* Java 8  
* [Gradle](https://gradle.org/)

## Clone, build and run
```
git clone --recursive https://github.com/ondrabublik/FlowProPackage.git FlowProPackage
cd FlowProPackage
gradle build
cd FlowPro
java -jar FlowPro.jar master 0
```
