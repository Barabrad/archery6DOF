# archery6DOF
My first 6-DOF simulation: an arrow fired at a target

This simulation accounts for gravity and drag, and it allows the user to input a wind speed and direction.
Quaternions are used to keep track of the arrow's orientation, and to simplify the body-fixed moments calculations, which introduces small errors due to the irrational numbers used in the quaternions.

The input parameters are documented in archery6FOD.m. The output is the difference in coordinates from the target bulls-eye and the arrow tip upon either crossing the target's plane or hitting the ground (whichever comes first).

The following plots are generated:
(See pictures of the first and last in [my portfolio](https://barabrad.github.io/portfolio/personal/archery6DOF/))
- A 3D path of the arrow tip and center of mass (COM), along with the x, y, and z components of the two curves
- A 3D "path" of the arrow COM's velocity, along with the x, y, and z components of the curve, which helps with troubleshooting
- A 3D "path" of the arrow's body-fixed angular velocity, along with the x, y, and z components of the curve, which helps with troubleshooting
- A 3D plot showing the arrow's final state relative to the target's bulls-eye

The following assumptions were made to simplify the simulation:
- Constant air density, wind speed, and gravity over the arrow's flight
- The archer and target are in an inertial frame
- The archer won't impart any rotation on the arrow when firing, such that the initial velocity vector matches the arrow's pointing vector
- Fins are right triangles with negligible thickness with their bases at at the end of the arrow shaft
- The nock's COM is located at the end of the shaft (since there is an inserted part)
- The nock and tip are close enough to axisymmetric
- The nock's length is negligible for purposes of calculating drag
- Moments from drag on a component can be modeled by having the drag force act on the component's COM
- The following drag values will be used from [https://www.engineeringtoolbox.com/drag-coefficient-d_627.html](https://www.engineeringtoolbox.com/drag-coefficient-d_627.html):
  - Streamlined body (arrow front): 0.04
  - Wires and cables (arrow shaft side): 1.0-1.3 -> 1.15
  - Solid hemisphere flow normal to curved side (side of tip): 0.47
  - Solid hemisphere flow normal to flat side (arrow back): 1.17
  - Squared flat plate at 90 deg (fins normal to wind): 1.17