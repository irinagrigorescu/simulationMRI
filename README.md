An educational tool for simulating MRI
======================================

This software is developed by Irina Grigorescu 
---------------------------------------------

Email: <irina.grigorescu.15@ucl.ac.uk>

This work is based on the book called "Magnetic Resonance Imaging: Physical Principles and Sequence Design" by Robert W. Brown, Yu-Chung Cheng and E. Mark Haacke.

The `project/learning/` folder contains different things I would like to learn.

The `report/` folder has Latex files with solved problems and notes.

### TODO
- [ ] Fix the phase stuff because phase is related to B1 field, and modelling it as a rotation about the z axis does not make sense physically.
- [ ] In the laboratory frame, make sure simulation ends with mu vector lying along the corresponding axis as it does in the rotating frame
- [ ] Remove Update buttons and update fields after enter https://practicalmatlab.wordpress.com/2012/07/09/textboxes/
- [ ] Add plot with Mz and Mxy in both rotating and laboratory frames
- [ ] Investigate off-resonance B1 

### Structure of the project
```
.
├── LICENSE                      # License file
├── README.md                    # This readme file
├── project/                     # MATLAB project files
│   ├── data/                    # Data files (*.mat)
│   ├── figures/                 # Saved figures
│   │   ├── plots/               
│   │   └── videos/
│   ├── src/                     # Functions for each chapter
│   │   ├── ch1/
│   │   ├── ch11/
│   │   ├── ch18/
│   │   ├── ch28/
│   │   └── ch9/
│   ├── helpers/                 # Helper functions
│   ├── learning/                # Scripts for when I want to learn about something
│   │   ├── guicreation/         # Learnt about how to create a GUI
│   │   │   ├── algs/
│   │   │   │   └── ch3/
│   │   │   └── gui/
│   │   │       └── ch3/
│   │   ├── imagingA2DObject/    # Simulating the entire process of imaging
│   │   └── spiral/              # Simulating spiral
│   └── test/                    # Tests
└── report                       # Latex report
    ├── Appendix.tex
    ├── MRIbook.pdf
    ├── MRIbook.tex              # Main report file
    ├── MRIbook.toc
    ├── chapters/                # Files for each chapter
    │   ├── Ch1.tex
    │   ├── Ch11.tex
    │   ├── Ch11Exercises.tex
    │   ├── Ch18.tex
    │   ├── Ch18Exercises.tex
    │   ├── Ch1Exercises.tex
    │   ├── Ch2.tex
    │   ├── Ch2Exercises.tex
    │   ├── Ch3.tex
    │   ├── Ch3Exercises.tex
    │   ├── Ch4.tex
    │   └── Ch4Exercises.tex
    ├── compile.sh               # Script to compile from Latex to Pdf file
    ├── img/                     # Images for the pdf 
    ├── mcode.sty
    └── references.bib           # References
 ```


