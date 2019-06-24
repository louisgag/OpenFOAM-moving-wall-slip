# OpenFOAM-moving-wall-slip
This boundary conditon for OpenFOAM allows to have a moving wall which imposes the normal velocity while letting the tangent velocities slip

It is currently compatible with OpenFOAM-6.x from [OpenFOAM.org](http://www.openfoam.org)
and constitutes an updated version of the code posted on the [forum](https://www.cfd-online.com/Forums/openfoam-solving/105274-free-slip-moving-wall-bc.html).

To compile this new AMI motion function, simply put the *finiteVolume* folder into your user directory and run *wmake* from inside the folder.

This boundary condition can be especially useful when running an AMI case or an [embedded-AMI](http://github.com/louisgag/OpenFOAM-embedded-AMI) case.
