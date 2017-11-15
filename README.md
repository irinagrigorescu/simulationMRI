An educational tool for simulating MRI
======================================

This software is developed by Irina Grigorescu 
---------------------------------------------

Email: <irina.grigorescu.15@ucl.ac.uk>

This work is based on the famous _Green Bible_: Magnetic Resonance Imaging: Physical Principles and Sequence Design by Robert W. Brown, Yu-Chung Cheng and E. Mark Haacke

Additionally, the learning/ folder may contain different things I like to try.
For example, I am currently experimenting with the extended phase graph concept, so you will find a multi-pulse experiment that can simulate echo formation.

### TODO
- [ ] Fix the phase stuff because phase is related to B1 field, and modelling it as a rotation about the z axis does not make sense physically.
- [ ] In the laboratory frame, make sure simulation ends with mu vector lying along the corresponding axis as it does in the rotating frame
- [ ] Remove Update buttons and update fields after enter https://practicalmatlab.wordpress.com/2012/07/09/textboxes/
- [ ] Add plot with Mz and Mxy in both rotating and laboratory frames
- [ ] Investigate off-resonance B1 

```
.
├── LICENSE
├── README.md
├── project
│   ├── data/
│   ├── figures/
│   │   ├── plots/
│   │   └── videos/
│   ├── functions
│   │   ├── ch1/
│   │   ├── ch11/
│   │   ├── ch18/
│   │   ├── ch28/
│   │   └── ch9/
│   ├── helpers/
│   ├── learning
│   │   ├── guicreation/
│   │   │   ├── algs/
│   │   │   │   └── ch3/
│   │   │   └── gui/
│   │   │       └── ch3/
│   │   ├── imagingA2DObject/
│   │   └── spiral/
│   └── test/
└── report
    ├── Appendix.tex
    ├── MRIbook.pdf
    ├── MRIbook.tex
    ├── MRIbook.toc
    ├── chapters/
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
    ├── compile.sh
    ├── img/
    ├── mcode.sty
    └── references.bib
 ```


