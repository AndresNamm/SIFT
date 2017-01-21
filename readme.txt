Youtube -  https://www.youtube.com/watch?v=123ZXYUJvu8

Overall parameters:
JDK 1.8
Opencv 2.4.13 and commons-math3-3.6.1


// This code is rather for me to test and understand the algorithm. I havent considered readibility too much in mind. Sorry :)



% Eclipse --------------------------------------------------------

Compile:

1. Import FeatureExtractionApplication project to eclipse

2. Right click FeatureExtractionApplication -> Properties -> Java Build Path -> Libraries -> Add External JAR ->Add commons-math3-3.6.1.jar and opencv-2413.jar

3. Go to Run Configuration -> Arguments -> add following to VM arguments:

-Djava.library.path=<path of opencv_java2413.dll>

Example for x64 dll, filename is not required: 
-Djava.library.path=C:/opencv/x64

4. Compile and run the program


Usage:

1. Load Image 1 and 2 (For image stitching, 1 has to be left image and 2 to be right image). Or Simply choose image from dropdown list

2. Click 1. Extract Feature 

3. Click 2. Describe Feature

4. Go back to Main tab and Click 3 Panorama

PS: All three buttons have to be clicked in sequence 

% INTELLIJJ --------------------------------------------------------

README for INTELLIJJ users

1) Select File->New->Project From existing sources
2) From the opened meny select the source folder our project named FeatureExtractionApplicationFinal .. 
3) Create project from source
4) Okay, okay ->  until it offers to open a new window with the project.
5) Add libraries Opencv 2.4.13 and commons-math3-3.6.1  - TUTORIAL-https://medium.com/@aadimator/how-to-set-up-opencv-in-intellij-idea-6eb103c1d45c#.jwuzxnckx



** The location of the images can be a source of the problem.
** If program is in D: drive and OpenCv ON another,  there might occur some issues. 
** Programmi asukoht viitses ei tohi olla täpitöhti. 
