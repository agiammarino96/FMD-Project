# FMD-Project
Project part of the course Functional Mechanical Design, Politecnico di Milano, academic year 2019/2020. Read FinalReport.pdf for the final report of the project. This project aims in the design of a multi-stage cam mechanism for an automatic machine. During this project two main tools have been used: MatLab programming and MSC Adams (software for multibody dynamics simulation). 

The steps for the design are summarized as follows:
1. Design of the motion laws involved in the mechanism using MatLab.
2. Design of the Cam profiles taking as input the motion laws of the previous step using MatLab.
3. Creation of the 3D model of the mechanism on MSC Adams, importing the profiles of the two cams from MatLab.
4. Kinematic and Dynamic simulation on MSC Adams for testing the design.
5. Post-processing of the results obtained from simulation.

The MatLab scripts of this project have been used for steps 1 and 2 (folder: src/To do before model Adams/) and for step 5 (folder: src/To do after model Adams/).

For obtaining the two cam profiles to be imported in MSC Adams, open folder "To do before model Adams":
1. Run DesignFirstCam.m --> generates internalCam.txt (internal profile first cam), externalCam.txt (external profile first cam), pitchProfileFirstCam.txt (pitch profile first cam), data_x.mat (motion law implemented by first cam) and curvature_radii_firstCam.mat (curvature radii internal and external profile first cam, needed later for computing the contact pressure).
2. Run motionlaw_Theta.m --> generates data_theta.mat (motion law implemented by second cam) using the function implemented in motionlaw.m.
3. Run DesignSecondCam.m --> generates SecondCam.txt (profile second cam), pitchProfileSecondCam.txt (pitch profile second cam) and curvature_radii_secondCam.mat (curvature radii second cam, needed later for computing the contact pressure).

For running the postprocessing open folder "To do after model Adams":
Run PostProcessing.m, the folder DataForPostprocessing contains data coming from the simulation performed in MSC Adams and the files curvature_radii_firstCam.mat and curvature_radii_secondCam.mat previously generated. 

The description of the relevant plots obtained by running the scripts can be found in FinalReport.pdf.
