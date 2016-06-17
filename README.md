# Face Perception SSVEP

This is an experiment from the [Face Categorisation Lab](http://face-categorization-lab.webnode.com/). The experiment is designed to detect a signature of face processing using steady-state visual evoked potentials. It does so by flashing random images rhythmically. Periodically, the image flashed will be a face, and at the frequency at which faces are shown, a clear response can be detected.

The images are not included in this repository, as there are copyright issues associated with sharing them freely. You can email the Face Categorisation Lab to get them (or get in touch with me).

Simply run face_cat.m to present images (make sure the images are in the right subdirectories).

Stimulation parameters are documented at the top of the script.

Relies on the psychtoolbox (for presentation) and fieldtrip (for analysis).
